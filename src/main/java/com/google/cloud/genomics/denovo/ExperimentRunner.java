/*
 *Copyright 2014 Google Inc. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
 * in compliance with the License. You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software distributed under the License
 * is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
 * or implied. See the License for the specific language governing permissions and limitations under
 * the License.
 */
package com.google.cloud.genomics.denovo;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.Callset;
import com.google.api.services.genomics.model.ContigBound;
import com.google.api.services.genomics.model.Dataset;
import com.google.api.services.genomics.model.Read;
import com.google.api.services.genomics.model.Variant;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;


/**
 * Holds all the experiments in the project
 *
 *  The two stages :
 *
 *  stage1 : Filter to get candidate denovo mutation sites using variantStore variants stage2 :
 * Bayesian Inference engine to select mutation sites from candidate sites
 */
public class ExperimentRunner {

  static public final long PROJECT_ID = 1085016379660L;
  static public final int TOT_CHROMOSOMES = 24;
  static public final long MAX_VARIANT_RESULTS = 10000L;
  static public final long DEFAULT_START_POS = 1L;
  static public final Float GQX_THRESH = Float.valueOf((float) 30.0);
  static public final Float QD_THRESH = Float.valueOf((float) 2.0);
  static public final Float MQ_THRESH = Float.valueOf((float) 20.0);
  static public Genomics genomics;
  public String candidatesFile;
  private CommandLine cmdLine;
  

  public ExperimentRunner(Genomics _genomics, CommandLine _cmdLine) {
    genomics = _genomics;
    cmdLine = _cmdLine;

    // Check command line for candidates file
    checkAndAddCandidatesFile();

  }

  /**
   * Check to see that candidatesFile is defined for experiments
   */
  private void checkAndAddCandidatesFile() {
    if (cmdLine.stageId == "stage1" || cmdLine.stageId == "stage2") {
      if (cmdLine.candidatesFile == null) {
        cmdLine.getUsage();
        throw new RuntimeException("Candidates File required");
      }
    }
    candidatesFile = cmdLine.candidatesFile;
  }

  /*
   * Stage 1 : Get a list of the candidates from VCF file
   */
  public void stage1() throws IOException {
    System.out.println("---------Starting Stage 1 VCF Filter -----------");

    // Define Experiment Specific Constant Values
    final String TRIO_DATASET_ID = "2315870033780478914";
    
    Map<TrioIndividual,String> individualCallsetNameMap = new HashMap<>();
    individualCallsetNameMap.put(TrioIndividual.DAD,"NA12877");
    individualCallsetNameMap.put(TrioIndividual.MOM,"NA12878");
    individualCallsetNameMap.put(TrioIndividual.CHILD,"NA12879");
    individualCallsetNameMap = Collections.unmodifiableMap(individualCallsetNameMap);

    final File outdir = new File(System.getProperty("user.home"), ".denovo_experiments");
    DenovoUtil.helperCreateDirectory(outdir);
    final File logFile = new File(outdir, "exp1.log");
    final File callFile = new File(outdir, candidatesFile);
    final File variantFile = new File(outdir, "exp1.vars");

    // Open File Outout handles
    try (PrintWriter logWriter = new PrintWriter(logFile);
        PrintWriter callWriter = new PrintWriter(callFile);
        PrintWriter variantWriter = new PrintWriter(variantFile);) {

      @SuppressWarnings("unused")
      List<Dataset> allDatasetsInProject;
      List<Variant> variants;
      List<ContigBound> contigBounds;
      List<Callset> callsets;
      List<VariantContigStream> variantContigStreams = new ArrayList<VariantContigStream>();
      Map<TrioIndividual, String> dictRelationCallsetId = new HashMap<>();

      /* Get a list of all the datasets */
      allDatasetsInProject = DenovoUtil.getAllDatasets();

      /* Get all the callsets for the trio dataset */
      callsets = DenovoUtil.getCallsets(TRIO_DATASET_ID);

      // Create a family person type to callset id map
      for (Callset callset : callsets) {
        String callsetName = callset.getName(); 
        for(TrioIndividual individual : TrioIndividual.values() ) {
          if(callsetName.equals(individualCallsetNameMap.get(individual))) {
            dictRelationCallsetId.put(individual,callset.getId());
            break;
          }
        }
      }

      // Check that the mapping has the correct keys
      // One time Sanity check ; could be replaced with testing
      if (!new HashSet<TrioIndividual>(dictRelationCallsetId.keySet()).equals(
          new HashSet<TrioIndividual>(Arrays.asList(TrioIndividual.values())))) {
        throw new RuntimeException("Callsets not found");
      }

      /* Get a list of all the Variants per contig */
      contigBounds = DenovoUtil.getVariantsSummary(TRIO_DATASET_ID);

      long startTime = System.currentTimeMillis();
      long prevTime = startTime;

      /* Iterate through each contig and do variant filtering for each contig */
      for (ContigBound currentContig : contigBounds) {
        System.out.println();
        System.out.println("Currently processing contig : " + currentContig.getContig());

        VariantContigStream variantContigStream =
            new VariantContigStream(currentContig, TRIO_DATASET_ID);
        variantContigStreams.add(variantContigStream);

        DenovoCaller denovoCaller = new DenovoCaller(dictRelationCallsetId);

        long denovoCount = 0;
        long variantCount = 0;

        while (variantContigStream.hasMore()) {
          variants = variantContigStream.getVariants();
          for (Variant variant : variants) {

            variantCount++;

            DenovoResult denovoCallResult = denovoCaller.callDenovoVariantIteration1(variant);
            if (denovoCallResult != null) {

              denovoCount++;
              // callWriter.println("denovo candidate at " + currentContig.getContig() +
              // ":position "
              // + variant.getPosition() + " ; Details " + denovoCallResult.details);
              callWriter.println(currentContig.getContig() + "," + variant.getPosition());
              callWriter.flush();
            }
          }
        }
        System.out.println();
        System.out.println("Denovo calls/variant Calls" + Long.valueOf(denovoCount) + "/"
            + Long.valueOf(variantCount));

        long contigTime = System.currentTimeMillis();
        System.out.println();
        System.out.println(
            "Time elapsed at Chromosome " + (contigTime - prevTime) + " milliseconds");
        prevTime = contigTime;

      }

      long endTime = System.currentTimeMillis();
      System.out.println();
      System.out.println("Total time taken " + (endTime - startTime) + " milliseconds");
    }
  }


  /*
   * Experiment 2 : Reads in candidate calls from stage 1 output file and then refines the
   * candidates
   */
  public void stage2() throws IOException {

    System.out.println("---- Starting Stage2 Bayesian Caller -----");

    // Constant Values Needed for stage 2 experiments
    Map<TrioIndividual, String> datasetIdMap = new HashMap<>();
    datasetIdMap.put(TrioIndividual.DAD, "4140720988704892492");
    datasetIdMap.put(TrioIndividual.MOM, "2778297328698497799");
    datasetIdMap.put(TrioIndividual.CHILD, "6141326619449450766");
    datasetIdMap = Collections.unmodifiableMap(datasetIdMap);

    Map<TrioIndividual, String> callsetIdMap = new HashMap<>();
    callsetIdMap.put(TrioIndividual.DAD, "NA12877");
    callsetIdMap.put(TrioIndividual.MOM, "NA12878");
    callsetIdMap.put(TrioIndividual.CHILD, "NA12879");
    callsetIdMap = Collections.unmodifiableMap(callsetIdMap);


    final File outdir = new File(System.getProperty("user.home"), ".denovo_experiments");
    DenovoUtil.helperCreateDirectory(outdir);
    final File exp1CallsFile = new File(outdir, candidatesFile);

    /* Find the readset Ids associated with the datasets */
    Map<TrioIndividual, String> readsetIdMap = DenovoUtil.createReadsetIdMap(datasetIdMap, callsetIdMap);

    System.out.println();
    System.out.println("Readset Ids Found");
    System.out.println(readsetIdMap.toString());

    try (BufferedReader callCandidateReader = new BufferedReader(new FileReader(exp1CallsFile))) {
      for (String line; (line = callCandidateReader.readLine()) != null;) {
        String[] splitLine = line.split(",");
        if (splitLine.length != 2) {
          throw new RuntimeException("Could not parse line : " + line);
        }
        String chromosome = splitLine[0];
        Long candidatePosition = Long.valueOf(splitLine[1]);

        System.out.println("Processing " + chromosome + ":" + String.valueOf(candidatePosition));


        /* Get reads for the current position */
        Map<TrioIndividual, List<Read>> readMap = new HashMap<>();
        for (TrioIndividual trioIndividual : TrioIndividual.values()) {
          List<Read> reads = DenovoUtil
              .getReads(readsetIdMap.get(trioIndividual), chromosome,
              candidatePosition, candidatePosition);
          readMap.put(trioIndividual, reads);
        }

        /*
         * Extract the relevant bases for the currrent position
         */
        Map<TrioIndividual, ReadSummary> readSummaryMap = new HashMap<>();
        for (TrioIndividual trioIndividual : TrioIndividual.values()) {
          readSummaryMap.put(trioIndividual,
              new ReadSummary(readMap.get(trioIndividual), candidatePosition));
        }

        /*
         * Call the bayes inference algorithm to generate likelihood
         */
        boolean isDenovo = BayesInfer.infer(readSummaryMap,cmdLine);
        
        if (isDenovo) {
          System.out.println("######### Denovo detected ########");
        }


        /* Set threshold and write candidate if passed */
        // TODO
      }

    } catch (IOException e) {
      e.printStackTrace();
    }
  }
}
