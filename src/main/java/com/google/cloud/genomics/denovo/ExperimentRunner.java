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

import static com.google.cloud.genomics.denovo.DenovoUtil.callsetIdMap;
import static com.google.cloud.genomics.denovo.DenovoUtil.datasetIdMap;
import static com.google.cloud.genomics.denovo.DenovoUtil.debugLevel;
import static com.google.cloud.genomics.denovo.DenovoUtil.individualCallsetNameMap;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.Callset;
import com.google.api.services.genomics.model.ContigBound;
import com.google.api.services.genomics.model.Read;
import com.google.api.services.genomics.model.Variant;
import com.google.cloud.genomics.denovo.DenovoUtil.InferenceMethod;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;
import com.google.common.base.Function;
import com.google.common.base.Optional;
import com.google.common.base.Predicate;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.ImmutableList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.ParseException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

/**
 * Holds all the experiments in the project
 *
 *  The two stages :
 *
 *  stage1 : Filter to get candidate denovo mutation sites using variantStore variants 
 *  stage2 : Bayesian Inference engine to select mutation sites from candidate sites
 */
public class ExperimentRunner {

  private CommandLine cmdLine;
  private Genomics genomics;
  private Map<TrioIndividual, String> individualCallsetIdMap = new HashMap<>();
  private final int numThreads;
  private Map<TrioIndividual, String> readsetIdMap = new HashMap<>();
  private final BayesInfer bayesInferrer;
  private final List<String> chromosomes;
  private final List<String> allChromosomes;
  private final InferenceMethod inferMethod;
  
  public ExperimentRunner(CommandLine cmdLine, Genomics genomics) throws IOException {
    this.cmdLine = cmdLine;
    this.genomics = genomics;
    numThreads = cmdLine.numThreads;
    
    // Set the overall System wide debug Level 
    DenovoUtil.debugLevel = cmdLine.debugLevel;
    DenovoUtil.LRT_THRESHOLD = cmdLine.lrtThreshold;
    
    try {
      readsetIdMap = DenovoUtil.createReadsetIdMap(datasetIdMap, callsetIdMap, genomics);
    } catch (IOException e) {
      // TODO(smoitra): Auto-generated catch block
      e.printStackTrace();
    }
    // Create the BayesNet inference object
    bayesInferrer = new BayesInfer(cmdLine.sequenceErrorRate, cmdLine.denovoMutationRate);
    allChromosomes = fetchAllChromosomes(genomics);
    chromosomes = verifyAndSetChromsomes(cmdLine.chromosomes);
    inferMethod = InferenceMethod.selectMethodfromString(cmdLine.inferMethod);
  }

  /*
   * Verify that the chromosomes supplied from the commandline are legal
   */
  private List<String> verifyAndSetChromsomes(List<String> chromosomes) {

    if (chromosomes == null || 
        chromosomes.size() == 1 && chromosomes.get(0).equals("all") ) {
      return allChromosomes;
    }

    for (String chr : chromosomes) {
      if (!allChromosomes.contains(chr)) {
        throw new IllegalArgumentException("Unknown chromosome " + chr);
      }
    }
    return chromosomes;
  }

  private ImmutableList<String> fetchAllChromosomes(Genomics genomics) throws IOException {
    // verify the commandline chromosomes
    List<ContigBound> contigBounds =
        DenovoUtil.getVariantsSummary(DenovoUtil.TRIO_DATASET_ID, genomics);

    ImmutableList<String> allChromosomes =
        FluentIterable.from(contigBounds).transform(new Function<ContigBound, String>() {
          @Override
          public String apply(ContigBound cb) {
            return cb.getContig();
          }
        }).toList();
    return allChromosomes;
  }

  public void execute() throws IOException, ParseException {
    if (cmdLine.stageId.equals("stage1")) {
      stage1();
    } else if (cmdLine.stageId.equals("stage2")) {
      stage2();
    } else {
      throw new IllegalArgumentException("Unknown stage : " + cmdLine.stageId);
    }
  }
  
  /*
   * Stage 1 : Get a list of the candidates from VCF file
   */
  public void stage1() throws IOException {
    if (debugLevel >= 0) {
      System.out.println("---------Starting Stage 1 VCF Filter -----------");  
    }
    
    // Define Experiment Specific Constant Values

    final File outdir = new File(System.getProperty("user.home"), ".denovo_experiments");
    DenovoUtil.helperCreateDirectory(outdir);
    final File outputFile = new File(outdir, cmdLine.outputFileName);

    // Open File Outout handles
    try (PrintWriter callWriter = new PrintWriter(outputFile);) {

      /* Get a list of all the datasets */
      DenovoUtil.getAllDatasets(genomics);

      /* Get all the callsets for the trio dataset */
      createCallsetIdMap(DenovoUtil.getCallsets(DenovoUtil.TRIO_DATASET_ID, genomics));

      /* Get a list of all the contigs */
      List<ContigBound> contigBounds = FluentIterable
          .from(DenovoUtil.getVariantsSummary(DenovoUtil.TRIO_DATASET_ID, genomics))
          .filter(new Predicate<ContigBound>() {

            @Override
            public boolean apply(ContigBound cb) {
              return chromosomes.contains(cb.getContig());
            }
          }).toList();

      ExecutorService executor = Executors.newFixedThreadPool(numThreads);
      /* Iterate through each contig and do variant filtering for each contig */
      for (ContigBound contig : contigBounds) {
        Runnable worker = new SimpleDenovoRunnable(callWriter, contig);
        executor.execute(worker);
      }

      executor.shutdown();
      while (!executor.isTerminated()) {
      }
      if (debugLevel >= 0) {
        System.out.println("All contigs processed");
      }
    }
  }

  private class SimpleDenovoRunnable implements Runnable {

    private final PrintWriter writer;
    private final ContigBound contig;

    public SimpleDenovoRunnable(PrintWriter writer, ContigBound contig) {
      this.writer = writer;
      this.contig = contig;
    }

    @Override
    public void run() {
      try {
        callSimpleDenovo(writer, contig);
      } catch (IOException e) {
        // TODO(smoitra): Auto-generated catch block
        e.printStackTrace();
      }
    }
  }

  /*
   * Writes calls to print stream
   */
  public synchronized void writeCalls(PrintWriter writer, String calls) {
    writer.print(calls);
    writer.flush();
  }

  /*
   * Creates a TrioType to callset ID map
   */
  private void createCallsetIdMap(List<Callset> callsets) {
    // Create a family person type to callset id map
    for (Callset callset : callsets) {
      String callsetName = callset.getName();
      for (TrioIndividual individual : TrioIndividual.values()) {
        if (callsetName.equals(individualCallsetNameMap.get(individual))) {
          individualCallsetIdMap.put(individual, callset.getId());
          break;
        }
      }
    }
  }

  /*
   * Check all the variants in the contig using a simple denovo caller
   */
  private void callSimpleDenovo(PrintWriter callWriter, ContigBound currentContig)
      throws IOException {
    // Create the caller object
    DenovoCaller denovoCaller = new DenovoCaller(individualCallsetIdMap);

    VariantContigStream variantContigStream = new VariantContigStream(genomics, currentContig,
        DenovoUtil.TRIO_DATASET_ID, DenovoUtil.MAX_VARIANT_RESULTS);

    while (variantContigStream.hasMore()) {
      List<Variant> variants = variantContigStream.getVariants();
      StringBuilder builder = new StringBuilder();

      for (Variant variant : variants) {

        Optional<String> denovoCallResultOptional = denovoCaller.callDenovoFromVarstore(variant);
        if (denovoCallResultOptional.isPresent()) {
          builder.append(
              String.format("%s,%d%n", currentContig.getContig(), variant.getPosition()));
        }
      }

      writeCalls(callWriter, builder.toString());
    }
  }

  /*
   * Stage 2 : Reads in candidate calls from stage 1 output file and then refines the candidates
   */
  public void stage2() throws FileNotFoundException, IOException, ParseException {

    if (debugLevel >= 0) {
      System.out.println("---- Starting Stage2 Bayesian Caller -----");  
    }

    final File outdir = new File(System.getProperty("user.home"), ".denovo_experiments");
    DenovoUtil.helperCreateDirectory(outdir);
    
    final File stage1CallsFile = new File(outdir, cmdLine.inputFileName);
    final File stage2CallsFile = new File(outdir, cmdLine.outputFileName);
    
    ExecutorService executor = new ThreadPoolExecutor(numThreads, // core thread pool size
        numThreads, // maximum thread pool size
        1, // time to wait before resizing pool
        TimeUnit.MINUTES, new ArrayBlockingQueue<Runnable>(numThreads, true), 
        new ThreadPoolExecutor.CallerRunsPolicy());

    /* Check if each candidate call is truly denovo by Bayesian denovo calling */
    try (BufferedReader callCandidateReader = new BufferedReader(new FileReader(stage1CallsFile));
        PrintWriter callWriter = new PrintWriter(stage2CallsFile);) {
      for (String line; (line = callCandidateReader.readLine()) != null;) {
        CallHolder callHolder = parseLine(line);

        /* Skip variant if chromosome does not match */
        if(!chromosomes.contains(callHolder.chromosome)) {
          continue;
        }
        
        if (numThreads >= 2) {
          executor.submit(new BayesDenovoRunnable(callHolder, callWriter));  
        } else {
          runBayesDenovoInference(callHolder, callWriter);
        }
      }
    } 
    // shutdown threadpool and wait
    executor.shutdown();
    while (!executor.isTerminated()) {
    }
  }

  private CallHolder parseLine(String line) throws ParseException {
    String[] splitLine = line.split(",");
    if (splitLine.length != 2) {
      throw new ParseException("Could not parse line : " + line, 0);
    }
    String chromosome = splitLine[0];
    Long candidatePosition = Long.valueOf(splitLine[1]);
    CallHolder callHolder = new CallHolder(chromosome, candidatePosition);
    return callHolder;
  }

  public Map<TrioIndividual, ReadSummary> getReadSummaryMap(Long candidatePosition,
      Map<TrioIndividual, List<Read>> readMap) {
    Map<TrioIndividual, ReadSummary> readSummaryMap = new HashMap<>();
    for (TrioIndividual trioIndividual : TrioIndividual.values()) {
      readSummaryMap.put(trioIndividual,
          new ReadSummary(readMap.get(trioIndividual), candidatePosition));
    }
    return readSummaryMap;
  }

  public Map<TrioIndividual, List<Read>> getReadMap(String chromosome, Long candidatePosition)
      throws IOException {
    /* Get reads for the current position */
    Map<TrioIndividual, List<Read>> readMap = new HashMap<>();
    for (TrioIndividual trioIndividual : TrioIndividual.values()) {
      List<Read> reads = DenovoUtil.getReads(readsetIdMap.get(trioIndividual), chromosome,
          candidatePosition, candidatePosition, genomics);
      readMap.put(trioIndividual, reads);
    }
    return readMap;
  }

  public Genomics getGenomics() {
    return genomics;
  }

  public void setGenomics(Genomics genomics) {
    this.genomics = genomics;
  }

  /*
   * Runs the Bayesian Denovo inference on a trio call 
   */
  private void runBayesDenovoInference(CallHolder callHolder, PrintWriter writer) {
    // get reads for chromosome and position
    Map<TrioIndividual, List<Read>> readMap = new HashMap<>();
    try {
      readMap = getReadMap(callHolder.chromosome, callHolder.position);
    } catch (IOException e) {
      // TODO(smoitra): Auto-generated catch block
      e.printStackTrace();
    }
    // Extract the relevant bases for the currrent position
    Map<TrioIndividual, ReadSummary> readSummaryMap =
        getReadSummaryMap(callHolder.position, readMap);

    // Call the bayes inference algorithm to generate likelihood
    BayesInfer.InferenceResult result = bayesInferrer.infer(readSummaryMap, inferMethod);

    if (debugLevel >= 1) {
      synchronized(this) {
        System.out.printf("%s,%d,%s", 
            callHolder.chromosome, callHolder.position, result.getDetails());
      }
    }
    
    if (result.isDenovo()) {
      writeCalls(writer, 
          String.format("%s,%d,%s", 
              callHolder.chromosome, callHolder.position, result.getDetails()));
    }
  }
  
  private class CallHolder {
    public final String chromosome;
    public final Long position;

    public CallHolder(String chromosome, Long position) {
      this.chromosome = chromosome;
      this.position = position;
    }
  }

  /*
   * Runs the Bayesian Denovo Calling code
   */
  private class BayesDenovoRunnable implements Runnable {

    private final CallHolder callHolder;
    private final PrintWriter writer;
    
    public BayesDenovoRunnable(CallHolder callHolder, PrintWriter writer) {
      this.callHolder = callHolder;
      this.writer = writer;
    }

    @Override
    public void run() {
      runBayesDenovoInference(callHolder, writer);
    }
  }
}
