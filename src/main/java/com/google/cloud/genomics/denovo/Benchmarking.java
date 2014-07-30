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

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.ContigBound;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;
import com.google.cloud.genomics.utils.GenomicsFactory;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.security.GeneralSecurityException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Runs benchmarking for read and variant search requests
 */
public class Benchmarking {

  final static String TRIO_DATASET_ID = "2315870033780478914";
  private int maxVarstoreRequests;
  private int maxReadstoreRequests;
  private Genomics genomics;
  private Random random;
  final static long RANDOM_SEED = 42L;
  private final String benchmarkTarget;
  private final PrintStream logStream;

  public static class Builder {

    // Required parameters
    private final String benchmarkTarget;
    private final PrintStream logStream;

    // Optional parameters
    private int maxVarstoreRequests = 10;
    private int maxReadstoreRequests = 10;

    public Builder(String benchmarkTarget, PrintStream logStream) {
      this.benchmarkTarget = benchmarkTarget;
      this.logStream = logStream;
    }

    public Builder maxVarstoreRequests(int val) {
      maxVarstoreRequests = val;
      return this;
    }
    
    public Builder maxReadstoreRequests(int val) {
      maxReadstoreRequests = val;
      return this;
    }
    
    public Benchmarking build() {
      return new Benchmarking(this);
    }
  }

  @Override
  public String toString() {
    return String.format("benchmarkTarget : %s%n" +
        "maxVarstoreRequests : %d%n" + 
        "maxReadstoreRequests : %d" 
    , benchmarkTarget, maxVarstoreRequests, maxReadstoreRequests);
  }
  
  private Benchmarking(Builder builder) {
    benchmarkTarget = builder.benchmarkTarget;
    logStream = builder.logStream;
    maxVarstoreRequests = builder.maxVarstoreRequests;
    maxReadstoreRequests = builder.maxReadstoreRequests;
    
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

  public void execute() {
    List<Long> executionTimes = null;
    try {
      if (benchmarkTarget == "varstore") {
        executionTimes = timeVarstoreRequests();
      } else if (benchmarkTarget == "readstore") {
        executionTimes = timeReadstoreRequests();
      } else {
        throw new IllegalArgumentException("Unknown benchmark target" + benchmarkTarget);
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
    
    logStream.println(toString());
    
    printStats(executionTimes, logStream);
  }

  private void printStats(List<Long> executionTimes, PrintStream out) {
    // Calculate summary stats
    DescriptiveStatistics stats = new DescriptiveStatistics();
    for (Long time : executionTimes) {
      stats.addValue(time);
    }

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
        new VariantContigStream(genomics, randomContig, TRIO_DATASET_ID);

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
        DenovoUtil.createReadsetIdMap(datasetIdMap, callsetIdMap, genomics);

    // get all the contigbounds
    List<ContigBound> contigBounds = DenovoUtil.getVariantsSummary(TRIO_DATASET_ID, genomics);

    List<Long> timingResults = new ArrayList<>();

    for (int numRequests = 0; numRequests < maxReadstoreRequests; numRequests++) {

      System.out.printf("Query #%d%n", numRequests);
      // extract a random contig and position
      ContigBound randomContig = contigBounds.get(random.nextInt(contigBounds.size()));

      long candidatePosition = getRandomPosition(randomContig);

      long startTime = Long.valueOf(System.currentTimeMillis());
      DenovoUtil.getReads(readsetIdMap.get(TrioIndividual.DAD), randomContig.getContig(),
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
  
  public static void main(String[] args) throws FileNotFoundException {
    Benchmarking benchmarking ; 
    
    // Readstore benchmarking
    benchmarking = new Benchmarking.Builder("readstore", new PrintStream("/tmp/readstore"))
        .maxReadstoreRequests(20)
        .build();
    benchmarking.execute();
    
    // Varstore benchmarking
    benchmarking = new Benchmarking.Builder("varstore", new PrintStream("/tmp/varstore")).build();
    benchmarking.execute();
  }
}
