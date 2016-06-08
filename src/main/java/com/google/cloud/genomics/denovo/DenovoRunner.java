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

import static com.google.cloud.genomics.denovo.DenovoUtil.Caller.FULL;
import static com.google.cloud.genomics.denovo.DenovoUtil.Caller.READ;
import static com.google.cloud.genomics.denovo.DenovoUtil.Caller.VARIANT;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioMember.CHILD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioMember.DAD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioMember.MOM;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.CallSet;
import com.google.api.services.genomics.model.ReadGroupSet;
import com.google.api.services.genomics.model.SearchCallSetsRequest;
import com.google.api.services.genomics.model.SearchCallSetsResponse;
import com.google.api.services.genomics.model.SearchReadGroupSetsRequest;
import com.google.cloud.genomics.denovo.DenovoUtil.Chromosome;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioMember;
import com.google.cloud.genomics.utils.GenomicsFactory;
import com.google.common.base.Function;
import com.google.common.collect.FluentIterable;

import java.io.IOException;
import java.security.GeneralSecurityException;
import java.text.ParseException;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler;
import java.util.logging.Formatter;
import java.util.logging.LogManager;
import java.util.logging.LogRecord;
import java.util.logging.Logger;

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
        .fromApplicationDefaultCredential();

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
            }).toImmutableSet());

    // Initialize the shared object
    shared = new DenovoShared.Builder()
      .datasetId(cmdLine.datasetId)
      .numThreads(cmdLine.numThreads)
      .logger(setUpLogger(cmdLine))
      .lrtThreshold(cmdLine.lrtThreshold)
      .genomics(genomics)
      .personToCallsetNameMap(personToCallsetNameMap)
      .personToReadGroupSetIdMap(createReadGroupSetIdMap(cmdLine.datasetId, personToCallsetNameMap, genomics))
      .personToCallsetIdMap(personToCallsetIdMap)
      .callsetIdToPersonMap(DenovoUtil.getReversedMap(personToCallsetIdMap))
      .startPosition(cmdLine.startPosition)
      .endPosition(cmdLine.endPosition)
      .chromosomes(chromosomes)
      .inferMethod(cmdLine.inferMethod)
      .caller(cmdLine.caller)
      .inputFileName(cmdLine.inputFileName)
      .outputFileName(cmdLine.outputFileName)
      .max_variant_results(cmdLine.maxVariantResults)
      .denovoMutationRate(cmdLine.denovoMutationRate)
      .sequenceErrorRate(cmdLine.sequenceErrorRate)
      .build();
  }

  /** Set up the logger. Creates a new one each time.
   * @param cmdLine
   * @return the logger
   * @throws IOException
   */
  private Logger setUpLogger(CommandLine cmdLine) throws IOException {

    if (LogManager.getLogManager().getLogger(DenovoRunner.class.getName()) != null) {
      return LogManager.getLogManager().getLogger(DenovoRunner.class.getName());
    }

    Logger logger = Logger.getLogger(DenovoRunner.class.getName());
    ConsoleHandler consoleHandler = new ConsoleHandler();
    logger.setLevel(cmdLine.logLevel.getLevel());
    logger.setUseParentHandlers(false);
    Formatter conciseFormat = new Formatter(){
      @Override
      public String format(LogRecord record) {
        StringBuilder sb = new StringBuilder();
        sb.append(DenovoUtil.LogLevel.levelMap.get(record.getLevel()));
        sb.append(" : ");
        sb.append(record.getMessage());
        sb.append("\n");
        return sb.toString();
      }
    };
    consoleHandler.setFormatter(conciseFormat);
    consoleHandler.setLevel(cmdLine.logLevel.getLevel());
    logger.addHandler(consoleHandler);

    if (cmdLine.logFile != null) {
      FileHandler fileHandler = new FileHandler(
          DenovoUtil.getNormalizedFile(cmdLine.logFile).getAbsolutePath(), false);
      fileHandler.setFormatter(conciseFormat);
      fileHandler.setLevel(cmdLine.logLevel.getLevel());
      logger.addHandler(fileHandler);
    }
    return logger;
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
   * Get all the call sets in dataset
   *
   * @param datasetId dataset under consideration
   * @param genomics genomics querying object
   * @return list of all call sets
   * @throws IOException
   */
  List<CallSet> getCallsets(String datasetId, Genomics genomics) throws IOException {
    SearchCallSetsResponse response = genomics.callsets().search(new SearchCallSetsRequest()
        .setVariantSetIds(Collections.singletonList(datasetId))).execute();
    List<CallSet> callsets = response.getCallSets();
    return Collections.unmodifiableList(callsets);
  }

  /**
   * Map callset trio members to callset ids
   *
   * @param callsets a list of all the callsets
   * @param personToCallsetNameMap a mapping from a tripo member to a callset name
   * @return a mapping from callset trio members to callset ids
   */
  Map<TrioMember, String> createCallsetIdMap(List<CallSet> callsets,
      Map<TrioMember, String> personToCallsetNameMap) {

    Map<TrioMember, String> callsetIdMap = new HashMap<>();
    // Create a family person type to callset id map
    for (CallSet callset : callsets) {
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
   * Create a mapping from trio members to read group set ids
   *
   * @param datasetId The dataset under consideration
   * @param callsetNameMap A mapping from trio members to callset names
   * @param genomics The genomics querying object
   * @return A mapping from trio members to read group set ids
   * @throws IOException
   */
  Map<TrioMember, String> createReadGroupSetIdMap(String datasetId,
      Map<TrioMember, String> callsetNameMap, Genomics genomics) throws IOException {
    Map<TrioMember, String> readGroupSetIdMap = new HashMap<>();

    Genomics.Readgroupsets.Search search = genomics.readgroupsets()
        .search(new SearchReadGroupSetsRequest()
            .setDatasetIds(Collections.singletonList(datasetId)));

    List<ReadGroupSet> readGroupSets = search.execute().getReadGroupSets();

    for (TrioMember person : TrioMember.values()) {
      for (ReadGroupSet readGroupSet : readGroupSets) {
        String sampleName = readGroupSet.getName();
        String id = readGroupSet.getId();

        if (sampleName.equals(callsetNameMap.get(person))) {
          readGroupSetIdMap.put(person, id);
        }
      }
    }
    // Check that the readGroupSetIdMap is sane
    if (readGroupSetIdMap.size() != 3) {
      throw new IllegalStateException("Borked readGroupSetIdMap" + readGroupSetIdMap);
    }
    return Collections.unmodifiableMap(readGroupSetIdMap);
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
