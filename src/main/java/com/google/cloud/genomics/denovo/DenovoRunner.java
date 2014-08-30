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
import static com.google.cloud.genomics.denovo.DenovoUtil.Caller.FULL;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioMember.CHILD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioMember.DAD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioMember.MOM;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.Callset;
import com.google.api.services.genomics.model.Readset;
import com.google.api.services.genomics.model.SearchCallsetsRequest;
import com.google.api.services.genomics.model.SearchReadsetsRequest;
import com.google.cloud.genomics.denovo.DenovoUtil.Chromosome;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioMember;
import com.google.cloud.genomics.utils.GenomicsFactory;
import com.google.common.base.Function;
import com.google.common.collect.FluentIterable;

import java.io.File;
import java.io.IOException;
import java.security.GeneralSecurityException;
import java.text.ParseException;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Wrappper for caller objects
 * 
 * <ul>
 * <li>Initializes Shared state object from commandline options</li>
 * <li>Sets up calling pipeline and the caller objects</li>
 * </ul>
 */
public class DenovoRunner {

  private final DenovoShared shared;
  private CommandLine cmdLine;
  
  public static DenovoRunner initFromCommandLine(CommandLine cmdLine) 
      throws IOException, GeneralSecurityException {
    return new DenovoRunner(cmdLine);
  }
  
  /**
   * Initializes Engine params from command line args
   * 
   * @param cmdLine Commandline arguments
   * @throws IOException thrown when Maps cannot be created upon failure on querying APIs  
   * @throws GeneralSecurityException API security authorization failure
   */
  private DenovoRunner(CommandLine cmdLine) throws IOException, GeneralSecurityException {

    Genomics genomics = GenomicsFactory.builder("genomics_denovo_caller").build()
        .fromClientSecretsFile(new File(cmdLine.clientSecretsFilename));

    Map<TrioMember, String> personToCallsetNameMap = createCallsetNameMap(cmdLine);
    Map<TrioMember, String> personToCallsetIdMap = createCallsetIdMap(
        getCallsets(cmdLine.datasetId, genomics), personToCallsetNameMap);
    this.cmdLine = cmdLine;
    Set<Chromosome> chromosomes = cmdLine.chromosomes == null 
        ? Chromosome.ALL 
        : EnumSet.copyOf(FluentIterable
            .from(cmdLine.chromosomes)
            .transform(new Function<String, Chromosome>() {
                @Override
                public Chromosome apply(String input) {
                  return Chromosome.fromString(input);
                }
            }).toList());
    
    
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
      .chromosomes(chromosomes)
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

  /**
   * @throws IOException Various API related querying troubles
   * @throws ParseException Input file could not be properly parsed
   * @throws GeneralSecurityException API security authorization failure
   */
  public void execute() throws IOException, ParseException, GeneralSecurityException {
    
    if (shared.getCaller() == VARIANT) {
      DenovoCallers.getVariantCaller(shared).execute();
    } else if (shared.getCaller() == READ && shared.getInputFileName() != null) {
      DenovoCallers.getReadCaller(shared).execute();
    } else if (shared.getCaller() == READ && shared.getInputFileName() == null) {
      throw new IllegalArgumentException("Input calls file needed for read mode");
    } else if (shared.getCaller() == FULL) {
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
      throw new IllegalArgumentException("Unknown caller mode : " + shared.getCaller());
    }
  }

  /**
   * Get all the callsets in dataset
   * 
   * @param datasetId dataset under consideration
   * @param genomics genomics querying object 
   * @return list of all callsets
   * @throws IOException
   */
  List<Callset> getCallsets(String datasetId, Genomics genomics) throws IOException {
    List<Callset> callsets = genomics.callsets()
        .search(new SearchCallsetsRequest().setDatasetIds(Collections.singletonList(datasetId)))
        .execute()
        .getCallsets(); 
    return Collections.unmodifiableList(callsets);
  }

  /**
   * Map callset trio members to callset ids
   * 
   * @param callsets a list of all the callsets
   * @param personToCallsetNameMap a mapping from a tripo member to a callset name
   * @return a mapping from callset trio members to callset ids
   */
  Map<TrioMember, String> createCallsetIdMap(List<Callset> callsets, 
      Map<TrioMember, String> personToCallsetNameMap) {

    Map<TrioMember, String> callsetIdMap = new HashMap<>();
    // Create a family person type to callset id map
    for (Callset callset : callsets) {
      String callsetName = callset.getName();
      for (TrioMember person : TrioMember.values()) {
        if (callsetName.equals(personToCallsetNameMap.get(person))) {
          callsetIdMap.put(person, callset.getId());
          break;
        }
      }
    }
    return Collections.unmodifiableMap(callsetIdMap);
  }

  /**
   * Create a mapping from trio members to readset ids
   * 
   * @param datasetId The dataset under consideration
   * @param callsetNameMap A mapping from trio members to callset names
   * @param genomics The genomics querying object
   * @return A mapping from trio members to readset ids
   * @throws IOException
   */
  Map<TrioMember, String> createReadsetIdMap(String datasetId,
      Map<TrioMember, String> callsetNameMap, Genomics genomics) throws IOException {
    Map<TrioMember, String> readsetIdMap = new HashMap<>();

    List<Readset> readsets = genomics.readsets()
        .search(
            new SearchReadsetsRequest().setDatasetIds(Collections.singletonList(datasetId)))
        .execute()
        .getReadsets();

    for (TrioMember person : TrioMember.values()) {
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
    return Collections.unmodifiableMap(readsetIdMap);
  }

  /**
   * Extract callset names from cmdline and associate with trio members
   * 
   * @param cmdLine
   * @return a mapping from trio members to callset names
   */
  Map<TrioMember, String> createCallsetNameMap(CommandLine cmdLine) {
    Map<TrioMember, String> callsetNameMap = new HashMap<>();
    callsetNameMap.put(DAD, cmdLine.dadCallsetName);
    callsetNameMap.put(MOM, cmdLine.momCallsetName);
    callsetNameMap.put(CHILD, cmdLine.childCallsetName);
    return Collections.unmodifiableMap(callsetNameMap);
  }

}
