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

import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.CHILD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.DAD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.MOM;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.Call;
import com.google.api.services.genomics.model.Callset;
import com.google.api.services.genomics.model.ContigBound;
import com.google.api.services.genomics.model.Read;
import com.google.api.services.genomics.model.Readset;
import com.google.api.services.genomics.model.SearchCallsetsRequest;
import com.google.api.services.genomics.model.Variant;
import com.google.cloud.genomics.denovo.DenovoUtil.InferenceMethod;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;
import com.google.cloud.genomics.denovo.VariantsBuffer.PositionCall;
import com.google.cloud.genomics.utils.GenomicsFactory;
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
import java.security.GeneralSecurityException;
import java.text.ParseException;
import java.util.Collections;
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

  private final DenovoShared shared;
  private final BayesInfer bayesInferrer;
  
  
  public static ExperimentRunner initFromCommandLine(CommandLine cmdLine) 
      throws IOException, GeneralSecurityException {
    return new ExperimentRunner(cmdLine);
  }
  
  private ExperimentRunner(CommandLine cmdLine) throws IOException, GeneralSecurityException {

    Genomics genomics = GenomicsFactory.builder("genomics_denovo_caller").build()
        .fromClientSecretsFile(new File(cmdLine.clientSecretsFilename));

    Map<TrioIndividual, String> personToCallsetNameMap = createCallsetNameMap(cmdLine);
    Map<TrioIndividual, String> personToCallsetIdMap = createCallsetIdMap(
        getCallsets(cmdLine.datasetId, genomics), personToCallsetNameMap);
    List<String> allChromosomes = fetchAllChromosomes(cmdLine.datasetId, genomics);

    shared = new DenovoShared.Builder()
      .datasetId(cmdLine.datasetId)
      .numThreads(cmdLine.numThreads)
      .debugLevel(cmdLine.debugLevel)
      .lrtThreshold(cmdLine.lrtThreshold)
      .genomics(genomics)
      .personToCallsetNameMap(personToCallsetNameMap)
      .personToReadsetIdMap(createReadsetIdMap(cmdLine.datasetId, personToCallsetNameMap, genomics))
      .personToCallsetIdMap(personToCallsetIdMap)
      .callsetIdToPersonMap(DenovoUtil.getReversedMap(personToCallsetIdMap))
      .startPosition(cmdLine.startPosition)
      .endPosition(cmdLine.endPosition)
      .allChromosomes(allChromosomes)
      .chromosomes(verifyAndSetChromsomes(cmdLine.chromosomes, allChromosomes))
      .inferMethod(InferenceMethod.selectMethodfromString(cmdLine.inferMethod))     
      .stageId(cmdLine.stageId)
      .inputFileName(cmdLine.inputFileName)
      .outputFileName(cmdLine.outputFileName)
      .max_api_retries(cmdLine.maxApiRetries)
      .max_variant_results(cmdLine.maxVariantResults)
      .denovoMutationRate(cmdLine.denovoMutationRate)
      .sequenceErrorRate(cmdLine.sequenceErrorRate)
      .build();
    
      bayesInferrer = new BayesInfer(shared);
  }

  /**
   * @return List<Callset>
   * @throws IOException
   */
  public static List<Callset> getCallsets(String datasetId, Genomics genomics) throws IOException {
    List<Callset> callsets = genomics.callsets()
        .search(new SearchCallsetsRequest().setDatasetIds(Collections.singletonList(datasetId)))
        .execute()
        .getCallsets(); 
    return callsets;
  }

  
  /*
   * Verify that the chromosomes supplied from the commandline are legal
   */
  private List<String> verifyAndSetChromsomes(List<String> chromosomes, 
      List<String> allChromosomes) {

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

  private List<String> fetchAllChromosomes(String datasetId, Genomics genomics) throws IOException {
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
    if (shared.getStageId().equals("stage1")) {
      stage1();
    } else if (shared.getStageId().equals("stage2")) {
      stage2();
    } else {
      throw new IllegalArgumentException("Unknown stage : " + shared.getStageId());
    }
  }

  /*
   * Stage 1 : Get a list of the candidates from VCF file
   */
  public void stage1() throws IOException {
    if (shared.getDebugLevel() >= 0) {
      System.out.println("---------Starting Stage 1 VCF Filter -----------");
    }

    // Define Experiment Specific Constant Values

    final File outdir = new File(System.getProperty("user.home"), ".denovo_experiments");
    DenovoUtil.helperCreateDirectory(outdir);
    final File outputFile = new File(outdir, shared.getOutputFileName());

    // Open File Outout handles
    try (PrintWriter callWriter = new PrintWriter(outputFile);) {

      /* Get a list of all the contigs */
      List<ContigBound> contigBounds = FluentIterable
          .from(DenovoUtil.getVariantsSummary(shared.getDatasetId(), shared.genomics))
          .filter(new Predicate<ContigBound>() {

            @Override
            public boolean apply(ContigBound cb) {
              return shared.getChromosomes().contains(cb.getContig());
            }
          }).toList();

      ExecutorService executor = Executors.newFixedThreadPool(shared.getNumThreads());
      /* Iterate through each contig and do variant filtering for each contig */
      for (ContigBound contigBound : contigBounds) {
        Long startContigPos = shared.getStartPosition() == null ? 1L : shared.getStartPosition();
        Long endContigPos = shared.getEndPosition() == null 
            ? contigBound.getUpperBound() : shared.getEndPosition();
        Long strideLength = (endContigPos - startContigPos) / shared.getNumThreads();

        for (int threadIdx = 0 ; threadIdx < shared.getNumThreads() ; threadIdx++) {
          long start = startContigPos + threadIdx * strideLength;
          long end = threadIdx == shared.getNumThreads() - 1 ? endContigPos : start + strideLength - 1;
          Runnable worker = new SimpleDenovoRunnable(callWriter, contigBound.getContig(),
              start, end);
          executor.execute(worker);
        }
      }

      executor.shutdown();
      while (!executor.isTerminated()) {
      }
      if (shared.getDebugLevel() >= 1) {
        System.out.println("All contigs processed");
      }
    }
  }

  private class SimpleDenovoRunnable implements Runnable {

    private final PrintWriter writer;
    private final String contig;
    private final Long startPos;
    private final Long endPos;
    private int numTries = 0;

    public SimpleDenovoRunnable(PrintWriter writer, String contig, 
        Long startPosition, Long endPosition) {
      this.writer = writer;
      this.contig = contig;
      this.startPos = startPosition;
      this.endPos = endPosition;
    }

    @Override
    public void run() {
      try {
        numTries++;
        callSimpleDenovo(writer, contig, startPos, endPos);
      } catch (IOException e) {
        e.printStackTrace();
        if (numTries < shared.getMaxApiRetries()) {
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
  private Map<TrioIndividual, String> createCallsetIdMap(List<Callset> callsets, 
      Map<TrioIndividual, String> personToCallsetNameMap) {

    Map<TrioIndividual, String> callsetIdMap = new HashMap<>();
    // Create a family person type to callset id map
    for (Callset callset : callsets) {
      String callsetName = callset.getName();
      for (TrioIndividual person : TrioIndividual.values()) {
        if (callsetName.equals(personToCallsetNameMap.get(person))) {
          callsetIdMap.put(person, callset.getId());
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

    for (TrioIndividual person : TrioIndividual.values()) {
      for (Readset readset : readsets) {
        String sampleName = readset.getName();
        String readsetId = readset.getId();

        if (sampleName.equals(callsetNameMap.get(person))) {
          readsetIdMap.put(person, readsetId);
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
  private void callSimpleDenovo(PrintWriter callWriter, String contig, 
      Long startPosition, Long endPosition)
      throws IOException {
    // Create the caller object

    VariantsBuffer vbuffer = new VariantsBuffer();
    
    VariantContigStream variantContigStream = new VariantContigStream(contig,
        startPosition, 
        endPosition,
        Lists.newArrayList(shared.getPersonToCallsetIdMap().values()),
        shared);

    while (variantContigStream.hasMore()) {
      StringBuilder builder = new StringBuilder();

      // Get a fresh batch of variants and filter those without calls
      Optional<List<Variant>> variantsFromStream = 
          Optional.fromNullable(variantContigStream.getVariants());
      
      if (!variantsFromStream.isPresent()) {
        return; 
      }
      
      Iterable<Variant> variants = FluentIterable
        .from(variantsFromStream.get())
        .filter(new Predicate<Variant>() {
          @Override
          public boolean apply(Variant variant) {
            return variant.getCalls() != null ? true : false;
          }});
      
      for (Variant variant : variants) {
        
        for (Call call : variant.getCalls()) {
          vbuffer.checkAndAdd(shared.getCallsetIdToPersonMap().get(call.getCallsetId()), 
              Pair.with(variant, call));
          
        }
        if(shared.getDebugLevel() > 1) {
          System.out.println(vbuffer);  
        }
        
        // Try to see if buffer can be processed
        while(vbuffer.canProcess()) {
          if(shared.getDebugLevel() > 1) {
            System.out.println("canProcess1 : "+vbuffer.toString());  
          }

          Optional<PositionCall> nextCall = Optional.fromNullable(vbuffer.retrieveNextCall());
          if (nextCall.isPresent()) {
            if (shared.getDebugLevel() > 1) {
              System.out.println(nextCall);  
            }
            
            if (nextCall.get().isDenovo()) {
              builder.append(
                  String.format("%s,%d,%s%n", contig, nextCall.get().position,
                      nextCall.get()));
            }
          }
          
          if(shared.getDebugLevel() > 1) {
            System.out.println("canProcess2 : "+vbuffer.toString());  
          }
          vbuffer.pop(CHILD);
          
          if(shared.getDebugLevel() > 1) {
            System.out.println("canProcess3 : "+vbuffer.toString());  
          }
        }
      }
      writeCalls(callWriter, builder.toString());
    }
    
    // Flush remaining buffer
    StringBuilder builder = new StringBuilder();
    while (!vbuffer.isEmpty(CHILD)) {
      
      if(shared.getDebugLevel() > 1) {
        System.out.println("Flush1 : "+vbuffer.toString());  
      }
      
      Optional<PositionCall> nextCall = Optional.fromNullable(vbuffer.retrieveNextCall());
      if (nextCall.isPresent()) {
        System.out.println(nextCall);
        if (nextCall.get().isDenovo()) {
          builder.append(
              String.format("%s,%d,%s%n", contig, nextCall.get().position,
                  nextCall.get()));
        }
      }
      
      if(shared.getDebugLevel() > 1) {
        System.out.println("Flush2 : "+vbuffer.toString());  
      }
      vbuffer.pop(CHILD);
    }
    writeCalls(callWriter, builder.toString());
  }

  /*
   * Stage 2 : Reads in candidate calls from stage 1 output file and then refines the candidates
   */
  public void stage2() throws ParseException, IOException {

    if (shared.getDebugLevel() >= 0) {
      System.out.println("---- Starting Stage2 Bayesian Caller -----");
    }

    final File outdir = new File(System.getProperty("user.home"), ".denovo_experiments");
    DenovoUtil.helperCreateDirectory(outdir);

    final File stage1CallsFile = new File(outdir, shared.getInputFileName());
    final File stage2CallsFile = new File(outdir, shared.getOutputFileName());

    ExecutorService executor = new ThreadPoolExecutor(shared.getNumThreads(), // core thread pool size
        shared.getNumThreads(), // maximum thread pool size
        1, // time to wait before resizing pool
        TimeUnit.MINUTES, new ArrayBlockingQueue<Runnable>(shared.getNumThreads(), true),
        new ThreadPoolExecutor.CallerRunsPolicy());

    /* Check if each candidate call is truly denovo by Bayesian denovo calling */
    try (BufferedReader callCandidateReader = new BufferedReader(new FileReader(stage1CallsFile));
        PrintWriter callWriter = new PrintWriter(stage2CallsFile);) {
      for (String line; (line = callCandidateReader.readLine()) != null;) {
        CallHolder callHolder = parseLine(line);

        /* Skip variant if chromosome does not match */
        if (!shared.getChromosomes().contains(callHolder.chromosome)) {
          continue;
        }

        if (shared.getNumThreads() >= 2) {
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
    for (TrioIndividual person : TrioIndividual.values()) {
      readSummaryMap.put(person,
          new ReadSummary(readMap.get(person), candidatePosition));
    }
    return readSummaryMap;
  }

  public Map<TrioIndividual, List<Read>> getReadMap(String chromosome, Long candidatePosition)
      throws IOException {
    /* Get reads for the current position */
    Map<TrioIndividual, List<Read>> readMap = new HashMap<>();
    for (TrioIndividual person : TrioIndividual.values()) {
      List<Read> reads = DenovoUtil.getReads(shared.getPersonToReadsetIdMap().get(person), chromosome,
          candidatePosition, candidatePosition, shared.genomics);
      readMap.put(person, reads);
    }
    return readMap;
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
    BayesInfer.InferenceResult result = bayesInferrer.infer(readSummaryMap, 
        shared.getInferMethod());

    if (shared.getDebugLevel() >= 1) {
      synchronized (this) {
        System.out.printf("%s,%d,%s%n", callHolder.chromosome, callHolder.position,
            result.getDetails());
      }
    }

    if (result.isDenovo()) {
      writeCalls(writer, String.format("%s,%d,%s%n", callHolder.chromosome, callHolder.position,
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
        if (numTries < shared.getMaxApiRetries()) {
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
