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

import static com.google.cloud.genomics.denovo.DenovoUtil.debugLevel;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.CHILD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.DAD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.MOM;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.Call;
import com.google.api.services.genomics.model.Callset;
import com.google.api.services.genomics.model.ContigBound;
import com.google.api.services.genomics.model.Read;
import com.google.api.services.genomics.model.Readset;
import com.google.api.services.genomics.model.Variant;
import com.google.cloud.genomics.denovo.DenovoUtil.InferenceMethod;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;
import com.google.cloud.genomics.denovo.VariantsBuffer.PositionCall;
import com.google.common.base.Function;
import com.google.common.base.Optional;
import com.google.common.base.Predicate;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

import org.javatuples.Pair;

import java.io.BufferedReader;
import java.io.File;
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
 *  stage1 : Filter to get candidate denovo mutation sites using variantStore variants stage2 :
 * Bayesian Inference engine to select mutation sites from candidate sites
 */
public class ExperimentRunner {

  /**
   * Placeholder for running denovo variant calling
   */
  private CommandLine cmdLine;
  private Genomics genomics;
  private final int numThreads;
  private final Map<TrioIndividual, String> personToCallsetIdMap;
  private final Map<TrioIndividual, String> personToReadsetIdMap;
  private final Map<TrioIndividual, String> personToCallsetNameMap;
  private final Map<String, TrioIndividual> callsetIdToPersonMap;

  private final BayesInfer bayesInferrer;
  private final List<String> chromosomes;
  private final List<String> allChromosomes;
  private final InferenceMethod inferMethod;
  private final String datasetId;
  private final Long startPosition;
  private final Long endPosition;
  
  public ExperimentRunner(CommandLine cmdLine, Genomics genomics) throws IOException {
    this.cmdLine = cmdLine;
    this.genomics = genomics;
    numThreads = cmdLine.numThreads;

    // Set the overall System wide debug Level
    DenovoUtil.debugLevel = cmdLine.debugLevel;
    DenovoUtil.LRT_THRESHOLD = cmdLine.lrtThreshold;

    datasetId = cmdLine.datasetId;

    personToCallsetNameMap = createCallsetNameMap(cmdLine);
    personToReadsetIdMap = createReadsetIdMap(datasetId, personToCallsetNameMap, genomics);
    personToCallsetIdMap = createCallsetIdMap(DenovoUtil.getCallsets(datasetId, genomics));
    callsetIdToPersonMap = DenovoUtil.getReversedMap(personToCallsetIdMap);
    
    startPosition = cmdLine.startPosition;
    endPosition = cmdLine.endPosition;

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

    if (chromosomes == null || chromosomes.size() == 1 && chromosomes.get(0).equals("all")) {
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
    List<ContigBound> contigBounds = DenovoUtil.getVariantsSummary(datasetId, genomics);

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

      /* Get a list of all the contigs */
      List<ContigBound> contigBounds = FluentIterable
          .from(DenovoUtil.getVariantsSummary(datasetId, genomics))
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
    private int numTries = 0;

    public SimpleDenovoRunnable(PrintWriter writer, ContigBound contig) {
      this.writer = writer;
      this.contig = contig;
    }

    @Override
    public void run() {
      try {
        numTries++;
        callSimpleDenovo(writer, contig);
      } catch (IOException e) {
        e.printStackTrace();
        if (numTries < DenovoUtil.MAX_API_RETRIES) {
          try {
            Thread.sleep(DenovoUtil.API_WAIT_MILLISEC);
          } catch (InterruptedException e1) {
            e1.printStackTrace();
          }
          System.err.printf("Attempt #%d : contig %s%n", numTries + 1, contig);
          run();
        } else {
          System.err.printf("Failed to run contig : %s%n", contig);
        }
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
  private Map<TrioIndividual, String> createCallsetIdMap(List<Callset> callsets) {

    Map<TrioIndividual, String> callsetIdMap = new HashMap<>();
    // Create a family person type to callset id map
    for (Callset callset : callsets) {
      String callsetName = callset.getName();
      for (TrioIndividual individual : TrioIndividual.values()) {
        if (callsetName.equals(personToCallsetNameMap.get(individual))) {
          callsetIdMap.put(individual, callset.getId());
          break;
        }
      }
    }
    return callsetIdMap;
  }

  /**
   * @param datasetId
   * @param callsetNameMap
   * @return Map<String, String>
   * @throws IOException
   */
  public Map<TrioIndividual, String> createReadsetIdMap(String datasetId,
      Map<TrioIndividual, String> callsetNameMap, Genomics genomics) throws IOException {
    Map<TrioIndividual, String> readsetIdMap = new HashMap<>();

    List<Readset> readsets = DenovoUtil.getReadsets(datasetId, genomics);

    for (TrioIndividual individual : TrioIndividual.values()) {
      for (Readset readset : readsets) {
        String sampleName = readset.getName();
        String readsetId = readset.getId();

        if (sampleName.equals(callsetNameMap.get(individual))) {
          readsetIdMap.put(individual, readsetId);
        }
      }
    }
    // Check that the readsetIdMap is sane
    if (readsetIdMap.size() != 3) {
      throw new IllegalStateException("Borked readsetIdMap" + readsetIdMap);
    }
    return readsetIdMap;
  }

  public Map<TrioIndividual, String> createCallsetNameMap(CommandLine cmdLine) {
    Map<TrioIndividual, String> callsetNameMap = new HashMap<>();
    callsetNameMap.put(DAD, cmdLine.dadCallsetName);
    callsetNameMap.put(MOM, cmdLine.momCallsetName);
    callsetNameMap.put(CHILD, cmdLine.childCallsetName);
    return callsetNameMap;
  }

  /*
   * Check all the variants in the contig using a simple denovo caller
   */
  private void callSimpleDenovo(PrintWriter callWriter, ContigBound currentContig)
      throws IOException {
    // Create the caller object

    VariantsBuffer vbuffer = new VariantsBuffer();
    
    VariantContigStream variantContigStream = new VariantContigStream(genomics, 
        currentContig.getContig(),
        datasetId, 
        DenovoUtil.MAX_VARIANT_RESULTS, 
        startPosition, 
        endPosition == null ? currentContig.getUpperBound() : endPosition,
        Lists.newArrayList(personToCallsetIdMap.values()));

    while (variantContigStream.hasMore()) {
      StringBuilder builder = new StringBuilder();

      // Get a fresh batch of variants and filter those without calls
      Iterable<Variant> variants = FluentIterable
        .from(variantContigStream.getVariants())
        .filter(new Predicate<Variant>() {
          @Override
          public boolean apply(Variant variant) {
            return variant.getCalls() != null ? true : false;
          }});
      
      for (Variant variant : variants) {
        
        for (Call call : variant.getCalls()) {
          vbuffer.checkAndAdd(callsetIdToPersonMap.get(call.getCallsetId()), 
              Pair.with(variant, call));
          
          if(debugLevel > 1) {
            System.out.println(vbuffer);  
          }
        }
        
        // Try to see if buffer can be processed
        while(vbuffer.canProcess()) {
          if(debugLevel > 1) {
            System.out.println("canProcess1 : "+vbuffer.toString());  
          }

          Optional<PositionCall> nextCall = Optional.fromNullable(vbuffer.retrieveNextCall());
          if (nextCall.isPresent()) {
            System.out.println(nextCall);
            if (nextCall.get().isDenovo()) {
              builder.append(
                  String.format("%s,%d,%s%n", currentContig.getContig(), nextCall.get().position,
                      nextCall.get()));
            }
          }
          
          if(debugLevel > 1) {
            System.out.println("canProcess2 : "+vbuffer.toString());  
          }
          vbuffer.pop(CHILD);
          
          if(debugLevel > 1) {
            System.out.println("canProcess3 : "+vbuffer.toString());  
          }
        }
      }
      writeCalls(callWriter, builder.toString());
    }
  }

  /*
   * Stage 2 : Reads in candidate calls from stage 1 output file and then refines the candidates
   */
  public void stage2() throws ParseException, IOException {

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
        TimeUnit.MINUTES, new ArrayBlockingQueue<Runnable>(numThreads, true), new ThreadPoolExecutor.CallerRunsPolicy());

    /* Check if each candidate call is truly denovo by Bayesian denovo calling */
    try (BufferedReader callCandidateReader = new BufferedReader(new FileReader(stage1CallsFile));
        PrintWriter callWriter = new PrintWriter(stage2CallsFile);) {
      for (String line; (line = callCandidateReader.readLine()) != null;) {
        CallHolder callHolder = parseLine(line);

        /* Skip variant if chromosome does not match */
        if (!chromosomes.contains(callHolder.chromosome)) {
          continue;
        }

        if (numThreads >= 2) {
          executor.submit(new BayesDenovoRunnable(callHolder, callWriter));
        } else {
          new RunBayesDenovoWithRetries().run(callHolder, callWriter);
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
      List<Read> reads = DenovoUtil.getReads(personToReadsetIdMap.get(trioIndividual), chromosome,
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

  public Map<TrioIndividual, String> getPersonToCallsetIdMap() {
    return personToCallsetIdMap;
  }

  public Map<String, TrioIndividual> getCallsetIdToPersonMap() {
    return callsetIdToPersonMap;
  }
  
  /*
   * Runs the Bayesian Denovo inference on a trio call
   */
  private void runBayesDenovoInference(CallHolder callHolder, PrintWriter writer)
      throws IOException {
    // get reads for chromosome and position
    Map<TrioIndividual, List<Read>> readMap =
        getReadMap(callHolder.chromosome, callHolder.position);

    // Extract the relevant bases for the currrent position
    Map<TrioIndividual, ReadSummary> readSummaryMap =
        getReadSummaryMap(callHolder.position, readMap);

    // Call the bayes inference algorithm to generate likelihood
    BayesInfer.InferenceResult result = bayesInferrer.infer(readSummaryMap, inferMethod);

    if (debugLevel >= 1) {
      synchronized (this) {
        System.out.printf("%s,%d,%s", callHolder.chromosome, callHolder.position,
            result.getDetails());
      }
    }

    if (result.isDenovo()) {
      writeCalls(writer, String.format("%s,%d,%s", callHolder.chromosome, callHolder.position,
          result.getDetails()));
    }
  }

  private class CallHolder {
    public final String chromosome;
    public final Long position;

    public CallHolder(String chromosome, Long position) {
      this.chromosome = chromosome;
      this.position = position;
    }

    @Override
    public String toString() {
      return String.format("<%s,%s>", chromosome, position);
    }
  }

  /*
   * Wraps BayesDenovo Call so as to allow retries on IOException
   */
  private class RunBayesDenovoWithRetries {
    private int numTries = 0;

    public void run(CallHolder callHolder, PrintWriter writer) {
      try {
        numTries++;
        runBayesDenovoInference(callHolder, writer);
      } catch (IOException e) {
        e.printStackTrace();
        if (numTries < DenovoUtil.MAX_API_RETRIES) {
          try {
            Thread.sleep(DenovoUtil.API_WAIT_MILLISEC);
          } catch (InterruptedException e1) {
            e1.printStackTrace();
          }
          System.err.printf("Attempt #%d : call %s%n", numTries + 1, callHolder);
          run(callHolder, writer);
        } else {
          System.err.printf("Failed to run call : %s%n", callHolder);
        }
      }
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
      new RunBayesDenovoWithRetries().run(callHolder, writer);
    }
  }
}
