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
import com.google.cloud.genomics.utils.GenomicsFactory;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

/**
 * Placeholder for running all Genomics Experiments.
 */
public class GenomicsExperiment {

  static CommandLine cmdLine;
  static ExperimentRunner expRunner;


  public static void main(String[] args) throws IOException {
    System.out.println("-------- Starting Genomics Experiment ---------");

    cmdLine = new CommandLine();

    try {
      // Parse the command line
      cmdLine.setArgs(args);
      Genomics genomics = GenomicsFactory.builder("genomics_denovo_caller").build()
          .fromClientSecretsFile(new File(cmdLine.clientSecretsFilename));

      expRunner = new ExperimentRunner(genomics, cmdLine);

      // Entry point for all Experiments
      executeExperiment(cmdLine.stageId);

    } catch (Exception e) {
      cmdLine.printHelp(e.getMessage() + "\n", System.err);
      e.printStackTrace();
      return;
    }
  }

  private static void executeExperiment(String stage_id) throws IllegalAccessException,
      IllegalArgumentException, InvocationTargetException, NoSuchMethodException {
    Method[] methods = ExperimentRunner.class.getDeclaredMethods();

    Method methodMatch = null;
    for (Method method : methods) {
      if (method.getName().equals(stage_id)) {
        methodMatch = method;
        break;
      }
    }
    if (methodMatch == null) {
      throw new NoSuchMethodException("No matching method found for Experiment : " + stage_id);
    } else {
      methodMatch.invoke(expRunner);
    }
  }

}
