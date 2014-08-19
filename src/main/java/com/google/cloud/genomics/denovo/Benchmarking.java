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

import static com.google.cloud.genomics.denovo.DenovoUtil.callsetNameMap;
import static com.google.cloud.genomics.denovo.DenovoUtil.datasetIdMap;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.ContigBound;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;
import com.google.cloud.genomics.utils.GenomicsFactory;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.security.GeneralSecurityException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Runs benchmarking for read and variant search requests
 */
public class Benchmarking {

  final static String TRIO_DATASET_ID = "2315870033780478914";
  private int numRepeatExperiment;
  private int maxVarstoreRequests;
  private int maxReadstoreRequests;
  private int numThreads;
  private long maxVariantResults;
  private Genomics genomics;
  private Random random;
  final static long RANDOM_SEED = 42L;
  private final String benchmarkTarget;
  private final PrintStream logStream;
  private String startTimeStamp;
  private String endTimeStamp;
  private boolean contigousReads;
  
  public static class Builder {

    // Required parameters
    private final String benchmarkTarget;
    private final PrintStream logStream;

    // Optional parameters
    private int numThreads = 1;
    private int maxVarstoreRequests = 10;
    private int maxReadstoreRequests = 10;
    private int numRepeatExperiment = 2;
    private long maxVariantResults = DenovoUtil.MAX_VARIANT_RESULTS;
    private boolean contiguousReads = false;

    public Builder(String benchmarkTarget, PrintStream logStream) {
      this.benchmarkTarget = benchmarkTarget;
      this.logStream = logStream;
    }

    public Builder contiguousReads(boolean val) {
      contiguousReads = val;
      return this;
    }
    
    public Builder numThreads(int val) {
      numThreads = val;
      return this;
    }
    
    public Builder maxVarstoreRequests(int val) {
      maxVarstoreRequests = val;
      return this;
    }
    
    public Builder maxReadstoreRequests(int val) {
      maxReadstoreRequests = val;
      return this;
    }
    
    public Builder numRepeatExperiment(int val) {
      numRepeatExperiment = val;
      return this;
    }
    
    public Builder maxVariantResults(long val) {
      maxVariantResults = val;
      return this;
    }
    
    public Benchmarking build() {
      return new Benchmarking(this);
    }
  }

  @Override
  public String toString() {
    return String.format("benchmarkTarget : %s%n" +
        "numThreads : %d%n" + 
        "startTimeStamp : %s%n" +
        "endTimeStamp : %s%n" +
        "numRepeatExperiment : %d%n" +
        "maxReadstoreRequests : %d%n" +
        "contiguousReads : %b%n" +
        "maxVarstoreRequests : %d%n" +
        "maxVariantResults : %d"   
        , benchmarkTarget, numThreads, startTimeStamp, endTimeStamp, numRepeatExperiment, 
        maxReadstoreRequests, contigousReads, maxVarstoreRequests, maxVariantResults);
  }
  
  private Benchmarking(Builder builder) {
    benchmarkTarget = builder.benchmarkTarget;
    logStream = builder.logStream;
    maxVarstoreRequests = builder.maxVarstoreRequests;
    maxReadstoreRequests = builder.maxReadstoreRequests;
    maxVariantResults = builder.maxVariantResults;
    numRepeatExperiment = builder.numRepeatExperiment;
    numThreads = builder.numThreads;
    contigousReads = builder.contiguousReads;
    
    String homeDir = System.getProperty("user.home");
    String clientSecretsFilename = homeDir + "/Downloads/client_secrets.json";

    try {
      genomics = GenomicsFactory.builder("genomics_denovo_caller").build()
          .fromClientSecretsFile(new File(clientSecretsFilename));
    } catch (IOException | GeneralSecurityException e) {
      e.printStackTrace();
    }

    random = new Random(RANDOM_SEED);
  }

  public class BenchmarkCallable implements Callable<List<Long>> {

    @Override
    public List<Long> call() throws Exception {
      List<Long> executionTimes = new ArrayList<>();
      try {
        if (benchmarkTarget == "varstore") {
          executionTimes.addAll(timeVarstoreRequests());
        } else if (benchmarkTarget == "readstore") {
          executionTimes.addAll(timeReadstoreRequests());
        } else {
          throw new IllegalArgumentException("Unknown benchmark target" + benchmarkTarget);
        }
      } catch (IOException e) {
        e.printStackTrace();
      }

      return executionTimes;
    }
  }
  
  public void execute() {
    
    startTimeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").
        format(Calendar.getInstance().getTime());
    
    List<Long> executionTimes = null;
    if(numThreads == 1) {
      executionTimes = executeSingleThreaded(); 
    } else {
      executionTimes = executeMultiThreaded();
    }
    
    endTimeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").
        format(Calendar.getInstance().getTime());
    
    logStream.println(toString());
    printStats(executionTimes, logStream);
  }

  /**
   * Parallel process different threads 
   */
  private List<Long> executeMultiThreaded() {
    
    System.out.println("Executing on multiple threads");
    
    ExecutorService executor = Executors.newFixedThreadPool(numThreads);
    List<Future<List<Long>>> futures = new ArrayList<>();
    List<Long> executionTimes = new ArrayList<>();
    
    for (int i = 0; i < numRepeatExperiment; i++) {
      Callable<List<Long>> worker = new BenchmarkCallable();
      Future<List<Long>> submit = executor.submit(worker);
      futures.add(submit);
    }
    
    executor.shutdown();
    while (!executor.isTerminated()) {
    }
    
    for (Future<List<Long>> future : futures) {
      try {
        executionTimes.addAll(future.get());
      } catch (InterruptedException e) {
        e.printStackTrace();
      } catch (ExecutionException e) {
        e.printStackTrace();
      }
    }
    
    return executionTimes;
  }

  /**
   * Execute each experiment on a single thread
   */
  private List<Long> executeSingleThreaded() {
    
    System.out.println("Executing on a single thread");
    
    List<Long> executionTimes = new ArrayList<>();
    try {
      for (int i = 0; i < numRepeatExperiment; i++) {
        if (benchmarkTarget == "varstore") {
          executionTimes.addAll(timeVarstoreRequests());
        } else if (benchmarkTarget == "readstore") {
          executionTimes.addAll(timeReadstoreRequests());
        } else {
          throw new IllegalArgumentException("Unknown benchmark target" + benchmarkTarget);
        }
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
    return executionTimes;
  }

  private void printStats(List<Long> executionTimes, PrintStream out) {
    // Calculate summary stats
    DescriptiveStatistics stats = new DescriptiveStatistics();
    
    double totTime = 0.0;
    for (Long time : executionTimes) {
      stats.addValue(time);
      totTime += time;
    }

    out.printf("Mean/thread : %fms%n", totTime / executionTimes.size() / numThreads);
    out.printf("Min : %fms%n", stats.getMin());
    out.printf("Max : %fms%n", stats.getMax());
    out.printf("Mean : %fms%n", stats.getMean());
    out.printf("Median : %fms%n", stats.getPercentile(50));
    out.printf("Mean : %fms%n", stats.getMean());
    out.printf("Sdev : %fms%n", stats.getStandardDeviation());
  }

  /*
   * time the Varstore Requests
   */
  private List<Long> timeVarstoreRequests() throws IOException {

    // get all the contigbounds
    List<ContigBound> contigBounds = DenovoUtil.getVariantsSummary(TRIO_DATASET_ID, genomics);

    // extract a random contig
    ContigBound randomContig = contigBounds.get(random.nextInt(contigBounds.size()));

    VariantContigStream variantContigStream =
        new VariantContigStream(genomics, randomContig, TRIO_DATASET_ID, maxVariantResults);

    List<Long> timingResults = new ArrayList<>();

    for (int numRequests = 0; variantContigStream.hasMore() && numRequests < maxVarstoreRequests;
        numRequests++) {

      long startTime = Long.valueOf(System.currentTimeMillis());
      variantContigStream.getVariants();
      long endTime = Long.valueOf(System.currentTimeMillis());

      timingResults.add((endTime - startTime));
    }
    return timingResults;
  }

  /*
   * time the Varstore Requests
   */
  private List<Long> timeReadstoreRequests() throws IOException {

    // create the readsetIdMap
    Map<TrioIndividual, String> readsetIdMap =
        DenovoUtil.createReadsetIdMap(datasetIdMap, callsetNameMap, genomics);

    // get all the contigbounds
    List<ContigBound> contigBounds = DenovoUtil.getVariantsSummary(TRIO_DATASET_ID, genomics);

    List<Long> timingResults = new ArrayList<>();

    ContigBound contig = contigBounds.get(random.nextInt(contigBounds.size()));
    long candidatePosition = 1L;
    
    for (int numRequests = 0; numRequests < maxReadstoreRequests; numRequests++) {

      System.out.printf("Query #%d%n", numRequests);

      if(contigousReads) {
        candidatePosition += random.nextInt(20000);
      } else {
        // extract a random contig and position
        contig = contigBounds.get(random.nextInt(contigBounds.size()));
        candidatePosition = getRandomPosition(contig);
      }

      long startTime = Long.valueOf(System.currentTimeMillis());
      DenovoUtil.getReads(readsetIdMap.get(TrioIndividual.DAD), contig.getContig(),
          candidatePosition, candidatePosition, genomics);
      long endTime = Long.valueOf(System.currentTimeMillis());

      timingResults.add((endTime - startTime));
    }
    return timingResults;
  }

  private long getRandomPosition(ContigBound randomContig) {
    long nextLong = -1L;
    for (nextLong = random.nextLong(); nextLong < 0L; nextLong = random.nextLong()) {
    }
    long candidatePosition = nextLong % randomContig.getUpperBound() + 1;
    if (candidatePosition < 0L) {
      throw new IllegalStateException("candidatePos < 0");
    }
    return candidatePosition;
  }
}
