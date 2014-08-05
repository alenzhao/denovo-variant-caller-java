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

import static org.junit.Assert.*;

import com.google.api.services.genomics.Genomics;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;
import com.google.cloud.genomics.utils.GenomicsFactory;

import static com.google.cloud.genomics.denovo.DenovoUtil.Genotypes.*;

import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Map;

public class BayesInferTest {

  private static Genomics genomics;
  private static ExperimentRunner expRunner;
  private static BayesInfer bayesInferrer;

  @BeforeClass
  public static void setUp() throws Exception {

    String homeDir = System.getProperty("user.home");

    String argsString = "stage1 --candidates_file candidate.calls.tmp "
        + "--client_secrets_filename " + homeDir + "/Downloads/client_secrets.json "
        + "--seq_err_rate 1e-2 " + "--denovo_mut_rate 1e-8";
    String[] args = argsString.split(" ");

    CommandLine cmdLine = new CommandLine();
    cmdLine.setArgs(args);

    genomics = GenomicsFactory.builder("genomics_denovo_caller").build()
        .fromClientSecretsFile(new File(cmdLine.clientSecretsFilename));

    expRunner = new ExperimentRunner(cmdLine, genomics);
    
    bayesInferrer = new BayesInfer(cmdLine.sequenceErrorRate, cmdLine.denovoMutationRate);
  }

  @Test
  public void testGenomicsIsNotNull() {
    assertNotNull(genomics);
  }

  @Test
  public void testExpRunnerIsNotNull() {
    assertNotNull(expRunner);
  }
  
  @Test
  public void testTrioPos816785() throws IOException {
    Map<TrioIndividual, ReadSummary> readSummaryMap =
        expRunner.getReadSummaryMap(816785L, expRunner.getReadMap("chr1", 816785L));
    BayesInfer.InferResult result = bayesInferrer.infer(readSummaryMap);
    
    assertFalse(result.isDenovo());
    assertEquals("816785 => [CC,CC,CC]", Arrays.asList(CC,CC,CC), result.getMaxTrioGenoType());
  }
  
  @Test
  public void testTrioPos846600() throws IOException {
    Map<TrioIndividual, ReadSummary> readSummaryMap =
        expRunner.getReadSummaryMap(846600L, expRunner.getReadMap("chr1", 846600L));
    BayesInfer.InferResult result = bayesInferrer.infer(readSummaryMap);
    
    assertFalse(result.isDenovo());
    assertEquals("846600 => [CC,CC,CC]", Arrays.asList(CC,CC,CC), result.getMaxTrioGenoType());
  }

  @Test
  public void testTrioPos763769() throws IOException {
    Map<TrioIndividual, ReadSummary> readSummaryMap =
        expRunner.getReadSummaryMap(763769L, expRunner.getReadMap("chr1", 763769L));
    BayesInfer.InferResult result = bayesInferrer.infer(readSummaryMap);
    
    assertFalse(result.isDenovo());
    assertEquals("763769 => [AA,AA,AA]", Arrays.asList(AA,AA,AA), result.getMaxTrioGenoType());
  }

  @Test
  public void testTrioPos1298169() throws IOException {
    Map<TrioIndividual, ReadSummary> readSummaryMap =
        expRunner.getReadSummaryMap(1298169L, expRunner.getReadMap("chr1", 1298169L));
    BayesInfer.InferResult result = bayesInferrer.infer(readSummaryMap);

    assertFalse(result.isDenovo());
    assertEquals("1298169 => [TT,TT,TT]", Arrays.asList(TT, TT, TT), result.getMaxTrioGenoType());
  }
}
