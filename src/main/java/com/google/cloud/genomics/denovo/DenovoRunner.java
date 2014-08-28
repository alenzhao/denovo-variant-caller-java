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

import static com.google.cloud.genomics.denovo.DenovoUtil.Caller.READ;
import static com.google.cloud.genomics.denovo.DenovoUtil.Caller.VARIANT;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.CHILD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.DAD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.MOM;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.Callset;
import com.google.api.services.genomics.model.Readset;
import com.google.api.services.genomics.model.SearchCallsetsRequest;
import com.google.api.services.genomics.model.SearchReadsetsRequest;
import com.google.cloud.genomics.denovo.DenovoUtil.Chromosome;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;
import com.google.cloud.genomics.utils.GenomicsFactory;

import java.io.File;
import java.io.IOException;
import java.security.GeneralSecurityException;
import java.text.ParseException;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Holds all the experiments in the project
 *
 *  The two stages :
 *
 *  stage1 : Filter to get candidate denovo mutation sites using variantStore variants stage2 :
 * Bayesian Inference engine to select mutation sites from candidate sites
 */
public class DenovoRunner {

  private final DenovoShared shared;
  private CommandLine cmdLine;
  
  public static DenovoRunner initFromCommandLine(CommandLine cmdLine) 
      throws IOException, GeneralSecurityException {
    return new DenovoRunner(cmdLine);
  }
  
  private DenovoRunner(CommandLine cmdLine) throws IOException, GeneralSecurityException {

    Genomics genomics = GenomicsFactory.builder("genomics_denovo_caller").build()
        .fromClientSecretsFile(new File(cmdLine.clientSecretsFilename));

    Map<TrioIndividual, String> personToCallsetNameMap = createCallsetNameMap(cmdLine);
    Map<TrioIndividual, String> personToCallsetIdMap = createCallsetIdMap(
        getCallsets(cmdLine.datasetId, genomics), personToCallsetNameMap);
    this.cmdLine = cmdLine;
    
    // Get a list of all the contigBounds and the chromosomes
    
    shared = new DenovoShared.Builder()
      .datasetId(cmdLine.datasetId)
      .numThreads(cmdLine.numThreads)
      .debugLevel(cmdLine.debugLevel)
      .lrtThreshold(cmdLine.lrtThreshold)
      .genomics(genomics)
      .personToCallsetNameMap(personToCallsetNameMap)
      .personToReadsetIdMap(createReadsetIdMap(cmdLine.datasetId, personToCallsetNameMap, genomics))
      .personToCallsetIdMap(personToCallsetIdMap)
      .callsetIdToPersonMap(DenovoUtil.getReversedMap(personToCallsetIdMap))
      .startPosition(cmdLine.startPosition)
      .endPosition(cmdLine.endPosition)
      .chromosomes(cmdLine.chromosomes == null 
          ? Chromosome.ALL : EnumSet.copyOf(cmdLine.chromosomes))
      .inferMethod(cmdLine.inferMethod)
      .caller(cmdLine.caller)
      .inputFileName(cmdLine.inputFileName)
      .outputFileName(cmdLine.outputFileName)
      .max_api_retries(cmdLine.maxApiRetries)
      .max_variant_results(cmdLine.maxVariantResults)
      .denovoMutationRate(cmdLine.denovoMutationRate)
      .sequenceErrorRate(cmdLine.sequenceErrorRate)
      .build();
  }

  public void execute() throws IOException, ParseException, GeneralSecurityException {
    
    if (shared.getCaller() == VARIANT) {
      DenovoCallers.getVariantCaller(shared).execute();
    } else if (shared.getCaller() == READ && shared.getInputFileName() != null) {
      DenovoCallers.getReadCaller(shared).execute();
    } else if (shared.getCaller() == READ && shared.getInputFileName() == null) {
      String outFile = shared.getOutputFileName();
      String tempOutFile = outFile + ".tmp";
      cmdLine.outputFileName = tempOutFile;
      cmdLine.caller = VARIANT;
      DenovoRunner.initFromCommandLine(cmdLine).execute();
      cmdLine.inputFileName = tempOutFile;
      cmdLine.outputFileName = outFile;
      cmdLine.caller = READ;
      DenovoRunner.initFromCommandLine(cmdLine).execute();
    } else {
      throw new IllegalArgumentException("Unknown stage : " + shared.getCaller());
    }
  }

  /**
   * @return List<Callset>
   * @throws IOException
   */
  List<Callset> getCallsets(String datasetId, Genomics genomics) throws IOException {
    List<Callset> callsets = genomics.callsets()
        .search(new SearchCallsetsRequest().setDatasetIds(Collections.singletonList(datasetId)))
        .execute()
        .getCallsets(); 
    return callsets;
  }

  /*
   * Creates a TrioType to callset ID map
   */
  Map<TrioIndividual, String> createCallsetIdMap(List<Callset> callsets, 
      Map<TrioIndividual, String> personToCallsetNameMap) {

    Map<TrioIndividual, String> callsetIdMap = new HashMap<>();
    // Create a family person type to callset id map
    for (Callset callset : callsets) {
      String callsetName = callset.getName();
      for (TrioIndividual person : TrioIndividual.values()) {
        if (callsetName.equals(personToCallsetNameMap.get(person))) {
          callsetIdMap.put(person, callset.getId());
          break;
        }
      }
    }
    return callsetIdMap;
  }

  /**
   * @param datasetId
   * @param callsetNameMap
   * @return Map<String, String>
   * @throws IOException
   */
  Map<TrioIndividual, String> createReadsetIdMap(String datasetId,
      Map<TrioIndividual, String> callsetNameMap, Genomics genomics) throws IOException {
    Map<TrioIndividual, String> readsetIdMap = new HashMap<>();

    List<Readset> readsets = genomics.readsets()
        .search(
            new SearchReadsetsRequest().setDatasetIds(Collections.singletonList(datasetId)))
        .execute()
        .getReadsets();

    for (TrioIndividual person : TrioIndividual.values()) {
      for (Readset readset : readsets) {
        String sampleName = readset.getName();
        String readsetId = readset.getId();

        if (sampleName.equals(callsetNameMap.get(person))) {
          readsetIdMap.put(person, readsetId);
        }
      }
    }
    // Check that the readsetIdMap is sane
    if (readsetIdMap.size() != 3) {
      throw new IllegalStateException("Borked readsetIdMap" + readsetIdMap);
    }
    return readsetIdMap;
  }

  Map<TrioIndividual, String> createCallsetNameMap(CommandLine cmdLine) {
    Map<TrioIndividual, String> callsetNameMap = new HashMap<>();
    callsetNameMap.put(DAD, cmdLine.dadCallsetName);
    callsetNameMap.put(MOM, cmdLine.momCallsetName);
    callsetNameMap.put(CHILD, cmdLine.childCallsetName);
    return callsetNameMap;
  }

}
