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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.Arrays;

/**
 * Performs benchmarking for readstore and variantstore search requests
 */
public class RunBenchmarks {

  private static File benchmarksDir; 

/*
 * Run benchmarking for maxResults field of variants/search
 */
  public static void benchmarkVarstoreMaxResults() throws FileNotFoundException {
    
    // varstore benchmarking
    for(long maxResults : Arrays.asList(1000L,5000L,10000L,20000L,30000L,40000L,50000L)) {

      String logFileString = benchmarksDir.getAbsolutePath() + "/varstore" +
          String.valueOf(maxResults/1000) + "K";
      Benchmarking benchmarking = new Benchmarking.Builder("varstore", new PrintStream(logFileString)).
          maxVariantResults(maxResults).
          maxVarstoreRequests(50).
          numRepeatExperiment(5).
          build();
      benchmarking.execute();
    }
  }
  
  /*
   * Run benchmarking for variants/search using threads
   */
  public static void benchmarkVarstoreThreaded() throws FileNotFoundException {

    int maxResults = 10000;
    for(int numThreads : Arrays.asList(2,5,10,20,50,100)) {
      String logFileString = benchmarksDir.getAbsolutePath() + "/varstore" +
          String.valueOf(maxResults/1000) + "K" + String.valueOf(numThreads);
      Benchmarking benchmarking = new Benchmarking.Builder("varstore", new PrintStream(logFileString)).
          maxVariantResults(maxResults).
          maxVarstoreRequests(10).
          numRepeatExperiment(numThreads * 5).
          numThreads(numThreads).
          build();
      benchmarking.execute();
    }
  }
  
  /*
   * Run benchmarking for reads/search
   */
  public static void benchmarkReadStore() throws FileNotFoundException {
    // Readstore benchmarking
    Benchmarking benchmarking = new Benchmarking.Builder("readstore", 
        new PrintStream(benchmarksDir.getAbsolutePath() + "/readstore")).
        maxReadstoreRequests(100).
        numRepeatExperiment(5).
        build();
    benchmarking.execute();
  }
  
  /*
   * Run benchmarking for reads/search
   */
  public static void benchmarkReadStoreContiguous() throws FileNotFoundException {
    // Readstore benchmarking
    Benchmarking benchmarking = new Benchmarking.Builder("readstore", 
        new PrintStream(benchmarksDir.getAbsolutePath() + "/readstoreContig")).
        maxReadstoreRequests(100).
        numRepeatExperiment(5).
        contiguousReads(true).
        build();
    benchmarking.execute();
  }
  /*
   * Run benchmarking for reads/search using multiple threads
   */
  public static void benchmarkReadstoreThreaded() throws FileNotFoundException {
    for(int numThreads : Arrays.asList(2,5,10,20,50,100)) {
        
      Benchmarking benchmarking = new Benchmarking.Builder("readstore", 
          new PrintStream(benchmarksDir.getAbsolutePath() + "/readstoreThreaded" + 
              String.valueOf(numThreads))).
          maxReadstoreRequests(100).
          numRepeatExperiment(numThreads * 5).
          numThreads(numThreads).
          build();
      benchmarking.execute();
    }
  }

  public static void benchmarkReadstoreContiguousThreaded() throws FileNotFoundException {
    for(int numThreads : Arrays.asList(2,5,10,20,50,100)) {
        
      Benchmarking benchmarking = new Benchmarking.Builder("readstore", 
          new PrintStream(benchmarksDir.getAbsolutePath() + "/readstoreContigThreaded" + 
              String.valueOf(numThreads))).
          maxReadstoreRequests(100).
          numRepeatExperiment(numThreads * 5).
          numThreads(numThreads).
          contiguousReads(true).
          build();
      benchmarking.execute();
    }
  }

  
  public static void main(String[] args) throws FileNotFoundException {
 
    benchmarksDir = new File(System.getProperty("user.home") + "/.benchmarks");
    DenovoUtil.helperCreateDirectory(benchmarksDir);

    benchmarkReadStoreContiguous();
    benchmarkReadstoreContiguousThreaded();
  }
}
