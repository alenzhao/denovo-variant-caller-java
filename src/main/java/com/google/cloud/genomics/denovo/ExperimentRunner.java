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

import static com.google.cloud.genomics.denovo.DenovoUtil.individualCallsetNameMap;
import static com.google.cloud.genomics.denovo.DenovoUtil.readsetIdMap;

import com.google.api.services.genomics.model.Callset;
import com.google.api.services.genomics.model.ContigBound;
import com.google.api.services.genomics.model.Dataset;
import com.google.api.services.genomics.model.Read;
import com.google.api.services.genomics.model.Variant;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;
import com.google.common.base.Optional;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
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


  public String candidatesFile;
  private CommandLine cmdLine;

  public ExperimentRunner(CommandLine _cmdLine) {
    cmdLine = _cmdLine;

    // Check command line for candidates file
    checkAndAddCandidatesFile();
  }

  public void execute() throws IOException {
    if (cmdLine.stageId.equals("stage1")) {
      stage1();
    } else if (cmdLine.stageId.equals("stage2")) {
      stage2();
    } else {
      throw new IllegalArgumentException("Unknown stage : " + cmdLine.stageId);
    }
  }
  
  /**
   * Check to see that candidatesFile is defined for experiments
   */
  private void checkAndAddCandidatesFile() {
    if (cmdLine.stageId == "stage1" || cmdLine.stageId == "stage2") {
      if (cmdLine.candidatesFile == null) {
        cmdLine.getUsage();
        throw new IllegalArgumentException("Candidates File required");
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
      callsets = DenovoUtil.getCallsets(DenovoUtil.TRIO_DATASET_ID);

      // Create a family person type to callset id map
      for (Callset callset : callsets) {
        String callsetName = callset.getName();
        for (TrioIndividual individual : TrioIndividual.values()) {
          if (callsetName.equals(individualCallsetNameMap.get(individual))) {
            dictRelationCallsetId.put(individual, callset.getId());
            break;
          }
        }
      }

      // Check that the mapping has the correct keys
      // One time Sanity check ; could be replaced with testing
      if (!new HashSet<TrioIndividual>(dictRelationCallsetId.keySet()).equals(
          new HashSet<TrioIndividual>(Arrays.asList(TrioIndividual.values())))) {
        throw new IllegalStateException("Callsets not found");
      }

      /* Get a list of all the Variants per contig */
      contigBounds = DenovoUtil.getVariantsSummary(DenovoUtil.TRIO_DATASET_ID);

      long startTime = System.currentTimeMillis();
      long prevTime = startTime;

      /* Iterate through each contig and do variant filtering for each contig */
      for (ContigBound currentContig : contigBounds) {
        System.out.println();
        System.out.println("Currently processing contig : " + currentContig.getContig());

        VariantContigStream variantContigStream =
            new VariantContigStream(currentContig, DenovoUtil.TRIO_DATASET_ID);
        variantContigStreams.add(variantContigStream);

        DenovoCaller denovoCaller = new DenovoCaller(dictRelationCallsetId);

        long denovoCount = 0;
        long variantCount = 0;

        while (variantContigStream.hasMore()) {
          variants = variantContigStream.getVariants();
          for (Variant variant : variants) {

            variantCount++;

            Optional<String> denovoCallResultOptional =
                denovoCaller.callDenovoFromVarstore(variant);
            if (denovoCallResultOptional.isPresent()) {
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
   * Stage 2 : Reads in candidate calls from stage 1 output file and then refines the candidates
   */
  public void stage2() {

    System.out.println("---- Starting Stage2 Bayesian Caller -----");

    final File outdir = new File(System.getProperty("user.home"), ".denovo_experiments");
    DenovoUtil.helperCreateDirectory(outdir);
    final File stage1CallsFile = new File(outdir, candidatesFile);

    /* Find the readset Ids associated with the datasets */

    System.out.println();
    System.out.println("Readset Ids Found");
    System.out.println(readsetIdMap.toString());

    // Create the BayesNet inference object
    BayesInfer bayesInferrer = new BayesInfer(cmdLine.sequenceErrorRate, cmdLine.denovoMutationRate);
    
    /* Check if each candidate call is truly denovo by Bayesian denovo calling */
    
    try (BufferedReader callCandidateReader = new BufferedReader(new FileReader(stage1CallsFile))) {
      for (String line; (line = callCandidateReader.readLine()) != null;) {
        String[] splitLine = line.split(",");
        if (splitLine.length != 2) {
          throw new ParseException("Could not parse line : " + line, 0);
        }
        String chromosome = splitLine[0];
        Long candidatePosition = Long.valueOf(splitLine[1]);

        System.out.println("Processing " + chromosome + ":" + String.valueOf(candidatePosition));


        /* Get reads for the current position */
        Map<TrioIndividual, List<Read>> readMap = new HashMap<>();
        for (TrioIndividual trioIndividual : TrioIndividual.values()) {
          List<Read> reads = DenovoUtil.getReads(readsetIdMap.get(trioIndividual), chromosome,
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
        boolean isDenovo = bayesInferrer.infer(readSummaryMap);

        if (isDenovo) {
          System.out.println("######### Denovo detected ########");
        }

        /* Set threshold and write candidate if passed */
        // TODO
      }

    } catch (IOException | ParseException e) {
      e.printStackTrace();
    }
  }
}
