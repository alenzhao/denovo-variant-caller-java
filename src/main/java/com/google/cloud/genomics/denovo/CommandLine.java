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

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.IOException;
import java.io.StringWriter;
import java.util.List;

/**
 * Command line options handler for GenomicsExperiment
 */
class CommandLine {

  CmdLineParser parser;

  @Argument(usage = "The stage of the calling pipeline ; usually stage1 or stage2",
      metaVar = "<stage id>", required = true)
  public String stageId = null;

  @Option(name = "--job_name", metaVar = "<job name>",
      usage = "Name of your job", required = true)
  public String jobName;

  @Option(name = "--output_file", metaVar = "<file>",
      usage = "File to write results")
  public String outputFileName = null;
  
  @Option(name = "--input_file", metaVar = "<file>",
      usage = "File to read from")
  public String inputFileName = null;
  
  @Option(name = "--client_secrets_filename", metaVar = "<client_secrets_filename>",
      usage = "Path to client_secrets.json")
  public String clientSecretsFilename = "client_secrets.json";

  @Option(name = "--seq_err_rate", metaVar = "<seq_err_rate>",
      usage = "Specify the sequence error rate (default 1e-2)")
  public double sequenceErrorRate = 1e-2;

  @Option(name = "--denovo_mut_rate", metaVar = "<denovo_mut_rate>",
      usage = "Specify the denovo mutation rate (default 1e-8)")
  public double denovoMutationRate = 1e-8;

  @Option(name = "--num_threads", metaVar = "<num_threads>",
      usage = "Specify the number of threads (default 1 ; 1 to 50 suggested)")
  public int numThreads = 1;

  @Option(name = "--debug_level", metaVar = "<debug_level>",
      usage = "specify the debug level (0 for no debug spew)")
  public int debugLevel = 0;

  @Option(name = "--chromosome", metaVar = "<chromosome>",
      usage = "specify the chromosomes to search (specify multiple times for multiple chromsomes)")
  public List<String> chromosomes;

  public CommandLine() {
    parser = new CmdLineParser(this);
  }

  public void setArgs(String[] args) throws CmdLineException {
    parser.parseArgument(args);
  }

  public void printHelp(String headline, Appendable out) throws IOException {
    out.append(headline).append("\n").append(getUsage());
  }

  public String getUsage() {
    StringWriter sw = new StringWriter();
    sw.append("Usage: GenomicsExperiment stage_id [flags...]\n");
    parser.printUsage(sw, null);
    return sw.toString();
  }

}
