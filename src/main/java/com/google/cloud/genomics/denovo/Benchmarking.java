package com.google.cloud.genomics.denovo;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.ContigBound;
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
import com.google.cloud.genomics.utils.GenomicsFactory;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.security.GeneralSecurityException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Runs benchmarking for read and variant search requests
 */
public class Benchmarking {

  final static String TRIO_DATASET_ID = "2315870033780478914";
  final static Integer maxVarstoreRequests = 10;

  public static void main(String[] args) throws IOException, GeneralSecurityException {


    String homeDir = System.getProperty("user.home");
    String clientSecretsFilename = homeDir + "/Downloads/client_secrets.json";

    Genomics genomics = GenomicsFactory.builder("genomics_denovo_caller").build()
        .fromClientSecretsFile(new File(clientSecretsFilename));

    // Set genomics state so that all classes can use it
    DenovoUtil.setGenomics(genomics);

    List<Long> executionTimes = timeVarstoreRequests();

    printStats(executionTimes, System.out);

  }

  private static void printStats(List<Long> executionTimes, PrintStream out) {
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
  private static List<Long> timeVarstoreRequests() throws IOException {

    // get all the contigbounds
    List<ContigBound> contigBounds = DenovoUtil.getVariantsSummary(TRIO_DATASET_ID);

    // extract a random contig
    ContigBound randomContig = contigBounds.get(new Random().nextInt(contigBounds.size()));

    VariantContigStream variantContigStream =
        new VariantContigStream(randomContig, TRIO_DATASET_ID);

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
}
