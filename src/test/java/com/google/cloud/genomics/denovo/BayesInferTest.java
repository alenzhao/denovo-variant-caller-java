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

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;

import com.google.api.services.genomics.Genomics;
import com.google.cloud.genomics.utils.GenomicsFactory;

public class BayesInferTest {

  private static Genomics genomics;
  private static ExperimentRunner expRunner;

  @BeforeClass
  public static void setUp() throws Exception {

    String homeDir = System.getProperty("user.home");

    String argsString = "stage1 --candidates_file candidate.calls.tmp "
        + "--client_secrets_filename " + homeDir + "/Downloads/client_secrets.json "
        + "--require_all_scopes";
    String[] args = argsString.split(" ");

    CommandLine cmdLine = new CommandLine();
    cmdLine.setArgs(args);

    genomics = GenomicsFactory.builder("genomics_denovo_caller").build()
        .fromClientSecretsFile(new File(cmdLine.clientSecretsFilename));

    expRunner = new ExperimentRunner(genomics, cmdLine);
  }

  @AfterClass
  public static void tearDown() throws Exception {}

  @Test
  public void testGenomicsIsNotNull() {
    assertTrue(genomics != null);
  }

  @Test
  public void testExpRunnerIsNotNull() {
    assertTrue(expRunner != null);
  }


}
