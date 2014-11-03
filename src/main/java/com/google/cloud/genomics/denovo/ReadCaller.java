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
import com.google.api.services.genomics.model.Read;
import com.google.api.services.genomics.model.SearchReadsRequest;
import com.google.cloud.genomics.denovo.DenovoUtil.Chromosome;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioMember;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.ParseException;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

/**
 * Makes Denovo calls by examining reads at candidate position and applying Bayesian Inference
 */
public class ReadCaller extends DenovoCaller {

  private final DenovoShared shared;
  private final BayesInfer bayesInferrer;
  
  /**
   * @param shared shared project state
   */
  public ReadCaller(DenovoShared shared) {
    this.shared = shared;
    bayesInferrer = new BayesInfer(shared);
  }
  /* (non-Javadoc)
   * @see com.google.cloud.genomics.denovo.DenovoCaller#execute()
   */
  @Override
  public void execute() throws ParseException, IOException, IOException {
    shared.getLogger().info("---- Starting Bayesian Read Caller ----");

    final File inputFile = DenovoUtil.getNormalizedFile(shared.getInputFileName());
    final File outputFile = DenovoUtil.getNormalizedFile(shared.getOutputFileName());

    shared.getLogger().fine(String.format("Input File : %s", inputFile.getAbsolutePath()));
    shared.getLogger().fine(String.format("Output File : %s", outputFile.getAbsolutePath()));
    
    ExecutorService executor = new ThreadPoolExecutor(shared.getNumThreads(), // core thread pool size
        shared.getNumThreads(), // maximum thread pool size
        1, // time to wait before resizing pool
        TimeUnit.MINUTES, new ArrayBlockingQueue<Runnable>(shared.getNumThreads(), true),
        new ThreadPoolExecutor.CallerRunsPolicy());

    /* Check if each candidate call is truly denovo by Bayesian denovo calling */
    
    int lineCount = 0;
    try (BufferedReader inputReader = new BufferedReader(new FileReader(inputFile));
        PrintWriter outputWriter = new PrintWriter(outputFile);) {
      for (String line; (line = inputReader.readLine()) != null;) {
        CallHolder callHolder = parseLine(line);

        if(lineCount % DenovoUtil.READ_LOG_FREQ == 0 && lineCount > 0){
          shared.getLogger().info(String.format("%d Read candidates processed", lineCount));
        }
        lineCount++;
        
        /* Skip variant if chromosome does not match */
        if (!shared.getChromosomes().contains(
            Chromosome.valueOf(callHolder.chromosome.toUpperCase()))) {
          continue;
        }

        if (shared.getNumThreads() >= 2) {
          executor.submit(new BayesDenovoRunnable(callHolder, outputWriter));
        } else {
          new RunBayesDenovoWithRetries().run(callHolder, outputWriter);
        }
      }
      // shutdown threadpool and wait
      executor.shutdown();
      while (!executor.isTerminated()) {
      }
    }
    shared.getLogger().info("---- Read caller terminated ----");
  }
  
  /**
   * Retreives Reads and makes a call to inference engine
   * @param callHolder container for storing candidate calls
   * @param writer print stream 
   * @throws IOException
   */
  void runBayesDenovoInference(CallHolder callHolder, PrintWriter writer) throws IOException {
    // get reads for chromosome and position
    Map<TrioMember, List<Read>> readMap = getReadMap(callHolder.chromosome, callHolder.position);

    // Extract the relevant bases for the currrent position
    Map<TrioMember, ReadSummary> readSummaryMap = getReadSummaryMap(callHolder.position, readMap);

    // Call the bayes inference algorithm to generate likelihood
    BayesInfer.BayesCallResult result =
        bayesInferrer.infer(readSummaryMap, shared.getInferMethod());

    if (result.isDenovo()) {
      shared.getLogger().fine(String.format(
          "%s,%d,%s", callHolder.chromosome, callHolder.position, result.getDetails()));
      synchronized (this) {
        writeCalls(writer, String.format("%s,%d,%s%n", callHolder.chromosome, callHolder.position,
            result.getDetails()));
      }
    }
  }

  /**
   * Extract ReadMaps and genereates read summaries ; stores them in map
   * @param candidatePosition 
   * @param readMap A map from trio member to reads
   * @return a map from trio member to read summaries
   */
  Map<TrioMember, ReadSummary> getReadSummaryMap(Long candidatePosition,
      Map<TrioMember, List<Read>> readMap) {
    Map<TrioMember, ReadSummary> readSummaryMap = new TreeMap<>();
    for (TrioMember person : TrioMember.values()) {
      readSummaryMap.put(person,
          new ReadSummary(readMap.get(person), candidatePosition));
    }
    return readSummaryMap;
  }

  /**
   * Query API to retreive reads for each of the trio members
   * @param chromosome
   * @param candidatePosition
   * @return map from trio member to reads
   * @throws IOException
   */
  Map<TrioMember, List<Read>> getReadMap(String chromosome, Long candidatePosition)
      throws IOException {
    /* Get reads for the current position */
    Map<TrioMember, List<Read>> readMap = new HashMap<>();
    for (TrioMember person : TrioMember.values()) {
      List<Read> reads = getReads(shared.getPersonToReadGroupSetIdMap().get(person), chromosome,
          candidatePosition, candidatePosition + 1, shared.getGenomics());
      readMap.put(person, reads);
    }
    return readMap;
  }

  
  /** Make read search request
   * @param readGroupSetId
   * @param chromosomeName
   * @param startPos
   * @param endPos
   * @param genomics
   * @return list of reads at candidate position
   * @throws IOException
   */
  List<Read> getReads(String readGroupSetId, String chromosomeName, long startPos,
      long endPos, Genomics genomics) throws IOException {
    SearchReadsRequest request = new SearchReadsRequest()
        .setReadGroupSetIds(Collections.singletonList(readGroupSetId))
        .setReferenceName(chromosomeName)
        .setStart(startPos)
        .setEnd(endPos);

    return genomics.reads().search(request).execute().getAlignments();
  }

  /**
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

  /**
   * Wrapper for multithreaded calls
   */
  private class BayesDenovoRunnable implements Runnable {

    private final CallHolder callHolder;
    private final PrintWriter writer;

    BayesDenovoRunnable(CallHolder callHolder, PrintWriter writer) {
      this.callHolder = callHolder;
      this.writer = writer;
    }

    @Override
    public void run() {
      new RunBayesDenovoWithRetries().run(callHolder, writer);
    }
  }
}
