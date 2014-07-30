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
import java.util.Map;
import java.util.Random;

/**
 * Runs benchmarking for read and variant search requests
 */
public class Benchmarking {

  final static String TRIO_DATASET_ID = "2315870033780478914";
  final static int MAX_VARSTORE_REQUESTS = 10;
  final static int MAX_READSTORE_REQUESTS = 10;
  private static Genomics genomics;
  private static Random random;
  final static long RANDOM_SEED = 42L;

  public static void main(String[] args) throws IOException, GeneralSecurityException {


    String homeDir = System.getProperty("user.home");
    String clientSecretsFilename = homeDir + "/Downloads/client_secrets.json";

    genomics = GenomicsFactory.builder("genomics_denovo_caller").build()
        .fromClientSecretsFile(new File(clientSecretsFilename));

    random = new Random(RANDOM_SEED);

    List<Long> executionTimes;
    executionTimes = timeVarstoreRequests();
    printStats(executionTimes, System.out);

    executionTimes = timeReadstoreRequests();
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
    List<ContigBound> contigBounds = DenovoUtil.getVariantsSummary(TRIO_DATASET_ID, genomics);

    // extract a random contig
    ContigBound randomContig = contigBounds.get(random.nextInt(contigBounds.size()));

    VariantContigStream variantContigStream =
        new VariantContigStream(genomics, randomContig, TRIO_DATASET_ID);

    List<Long> timingResults = new ArrayList<>();

    for (int numRequests = 0; variantContigStream.hasMore() && numRequests < MAX_VARSTORE_REQUESTS;
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
  private static List<Long> timeReadstoreRequests() throws IOException {

    // create the readsetIdMap
    Map<TrioIndividual, String> readsetIdMap =
        DenovoUtil.createReadsetIdMap(datasetIdMap, callsetIdMap, genomics);

    // get all the contigbounds
    List<ContigBound> contigBounds = DenovoUtil.getVariantsSummary(TRIO_DATASET_ID, genomics);

    List<Long> timingResults = new ArrayList<>();

    for (int numRequests = 0; numRequests < MAX_READSTORE_REQUESTS; numRequests++) {

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

  private static long getRandomPosition(ContigBound randomContig) {
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
