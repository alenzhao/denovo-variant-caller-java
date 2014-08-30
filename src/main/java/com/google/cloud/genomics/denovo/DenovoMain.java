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
import java.io.IOException;
import java.security.GeneralSecurityException;
import java.text.ParseException;

/**
 * Placeholder for running all Genomics Experiments.
 */
public class DenovoMain {

  private static CommandLine cmdLine;

  public static void main(String[] args) throws IOException {
    
    cmdLine = new CommandLine();
    try {
      cmdLine.setArgs(args);
    } catch (Exception e) {
      cmdLine.printHelp(e.getMessage() + "\n", System.err);
      e.printStackTrace();
      System.exit(1);
    }
    
    
    
    try {
      DenovoRunner.initFromCommandLine(cmdLine).execute();
    } catch (GeneralSecurityException | ParseException e) {
      e.printStackTrace();
      return;
    }
    
  }
}
