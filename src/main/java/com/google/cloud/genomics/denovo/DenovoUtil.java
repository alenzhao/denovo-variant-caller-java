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

import java.io.File;
import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.Call;
import com.google.api.services.genomics.model.Callset;
import com.google.api.services.genomics.model.ContigBound;
import com.google.api.services.genomics.model.Dataset;
import com.google.api.services.genomics.model.GetVariantsSummaryResponse;
import com.google.api.services.genomics.model.ListDatasetsResponse;
import com.google.api.services.genomics.model.Read;
import com.google.api.services.genomics.model.Readset;
import com.google.api.services.genomics.model.SearchCallsetsRequest;
import com.google.api.services.genomics.model.SearchCallsetsResponse;
import com.google.api.services.genomics.model.SearchReadsRequest;
import com.google.api.services.genomics.model.SearchReadsResponse;
import com.google.api.services.genomics.model.SearchReadsetsRequest;
import com.google.api.services.genomics.model.SearchReadsetsResponse;
import com.google.api.services.genomics.model.SearchVariantsRequest;

/*
 * Utility functions shared by other classes in Denovo project
 */
public class DenovoUtil {

  public static final double EPS = 1e-12;

  public static enum TrioIndividual {
    DAD("DAD"), MOM("MOM"), CHILD("CHILD");

    private String name;

    private TrioIndividual(String _name) {
      name = _name;
    }

    public String getName() {
      return name;
    }

  }
  public static enum Genotypes {
    AA("AA"),
    AC("AC"),
    AT("AT"),
    AG("AG"),
    CC("CC"),
    CT("CT"),
    CG("CG"),
    TT("TT"),
    TG("TG"),
    GG("GG");

    private final String alleles;

    private Genotypes(String _alleles) {
      alleles = _alleles;
    }

    public final String getAlleles() {
      return alleles;
    }

  }

  public static SearchVariantsRequest createSearchVariantsRequest(SearchVariantsRequest oldRequest,
      ContigBound contig,
      long startPos,
      long endPos,
      String datasetId,
      String nextPageToken) {

    // Init searchRequest obj
    SearchVariantsRequest searchVariantsRequest;
    if (oldRequest == null) {
      searchVariantsRequest = new SearchVariantsRequest();
    } else {
      searchVariantsRequest = oldRequest;
    }

    searchVariantsRequest
        .setContig(contig.getContig())
        .setDatasetId(datasetId)
        .setStartPosition(startPos)
        .setEndPosition(endPos)
        .setMaxResults(BigInteger.valueOf(ExperimentRunner.MAX_VARIANT_RESULTS))
        .setPageToken(nextPageToken);

    return searchVariantsRequest;
  }

  private static SearchReadsRequest createSearchReadsRequest(String readsetId,
      String chromosomeName, long startPos, long endPos) {
    // Init searchRequest obj
    SearchReadsRequest searchReadsRequest = new SearchReadsRequest();
    List<String> readsetIdList = new ArrayList<String>();
    readsetIdList.add(readsetId);

    searchReadsRequest.setReadsetIds(readsetIdList).setSequenceName(chromosomeName)
        .setSequenceStart(BigInteger.valueOf(startPos)).setSequenceEnd(BigInteger.valueOf(endPos));

    return searchReadsRequest;
  }

  private static SearchReadsetsRequest createSearchReadsetsRequest(String datasetId) {

    // Init searchRequest obj
    SearchReadsetsRequest searchReadsetsRequest = new SearchReadsetsRequest();

    // pack the dataset Id into a list
    List<String> datasetIds = new ArrayList<String>();
    datasetIds.add(datasetId);

    searchReadsetsRequest.setDatasetIds(datasetIds);

    return searchReadsetsRequest;

  }

  private static SearchCallsetsRequest createSearchCallsetsRequest(String datasetId) {

    // Init searchRequest obj
    SearchCallsetsRequest searchCallsetsRequest = new SearchCallsetsRequest();

    // pack the dataset Id into a list
    List<String> datasetIds = new ArrayList<String>();
    datasetIds.add(datasetId);

    searchCallsetsRequest.setDatasetIds(datasetIds);

    return searchCallsetsRequest;

  }


  /**
   * @param datasetIdMap
   * @param callsetIdMap
   * @return Map<String, String>
   * @throws IOException
   */
  public static Map<TrioIndividual, String> createReadsetIdMap(
      Map<TrioIndividual, String> datasetIdMap, Map<TrioIndividual, String> callsetIdMap)
      throws IOException {
    Map<TrioIndividual, String> readsetIdMap = new HashMap<>();

    for (TrioIndividual trioIndividual : datasetIdMap.keySet()) {
      List<Readset> readsets = getReadsets(datasetIdMap.get(trioIndividual));

      for (Readset readset : readsets) {
        // System.out.println(
        // readset.getDatasetId() + ":" + readset.getName() + ":" + readset.getId());
        String sampleName = readset.getName();
        String readsetId = readset.getId();

        for (TrioIndividual individual : callsetIdMap.keySet()) {
          if (callsetIdMap.get(individual).equals(sampleName)) {
            readsetIdMap.put(individual, readsetId);
          }
        }
      }
    }
    // Check that the readsetIdMap is sane
    if (readsetIdMap.size() != 3) {
      throw new IllegalStateException("Borked readsetIdMap" + readsetIdMap.toString());
    }
    return readsetIdMap;
  }

  public static List<Integer> getGenotype(Call call) {
    List<Integer> genoType = new ArrayList<Integer>();

    String genoTypeString = call.getInfo().get("GT").get(0);
    String splitChar = null;
    if (genoTypeString.contains("/")) {
      splitChar = "/";
    } else if (genoTypeString.contains("|")) {
      splitChar = "|";
    } else {
      return null;
    }

    for (String allele : genoTypeString.split(splitChar)) {
      genoType.add(Integer.valueOf(allele));
    }
    return genoType;

  }

  public static List<Read> getReads(String readsetId, String chromosomeName, long startPos,
      long endPos) throws IOException {

    SearchReadsRequest searchReadsRequest =
        createSearchReadsRequest(readsetId, chromosomeName, startPos, endPos);


    SearchReadsResponse searchReadsExecuted =
        ExperimentRunner.genomics.reads().search(searchReadsRequest).execute();

    List<Read> reads = searchReadsExecuted.getReads();
    return reads;
  }



  private static List<Readset> getReadsets(String datasetId) throws IOException {
    SearchReadsetsRequest searchReadsetsRequest = createSearchReadsetsRequest(datasetId);
    Genomics.Readsets.Search search =
        ExperimentRunner.genomics.readsets().search(searchReadsetsRequest);
    SearchReadsetsResponse execute = search.execute();
    List<Readset> readsets = execute.getReadsets();

    return readsets;
  }


  /**
   * @return List<Callset>
   * @throws IOException
   */
  public static List<Callset> getCallsets(String datasetId) throws IOException {
    SearchCallsetsRequest searchCallsetsRequest = createSearchCallsetsRequest(datasetId);
    Genomics.Callsets.Search search =
        ExperimentRunner.genomics.callsets().search(searchCallsetsRequest);
    SearchCallsetsResponse execute = search.execute();
    List<Callset> callsets = execute.getCallsets();

    System.out.println();
    System.out.println("Callsets found for dataset: " + datasetId);
    System.out.print(execute.toPrettyString());
    return callsets;
  }

  /**
   * @return List<Dataset>
   * @throws IOException
   */
  public static List<Dataset> getAllDatasets() throws IOException {
    System.out.println();
    System.out.println("######## Datasets under Project ########");

    // Get a list of all the datasets associated with project id
    Genomics.Datasets.List datasetRequest =
        ExperimentRunner.genomics.datasets().list().setProjectId(ExperimentRunner.PROJECT_ID);
    datasetRequest.setDisableGZipContent(true);

    ListDatasetsResponse execute = datasetRequest.execute();
    List<Dataset> datasets = execute.getDatasets();

    System.out.println("Datasets : " + execute.toPrettyString());
    return datasets;
  }

  /**
   * @return List<ContigBound>
   * @throws IOException
   */
  public static List<ContigBound> getVariantsSummary(String datasetId) throws IOException {
    System.out.println();
    System.out.println("######## Variants in Dataset ########");
    System.out.println("Querying DatasetID : " + datasetId);

    Genomics.Variants.GetSummary variantsSummaryRequest =
        ExperimentRunner.genomics.variants().getSummary().setDatasetId(datasetId);
    variantsSummaryRequest.setDisableGZipContent(true);

    GetVariantsSummaryResponse execute = variantsSummaryRequest.execute();
    // System.out.println("Variants : "+execute.toPrettyString());


    List<ContigBound> contigBounds = execute.getContigBounds();
    BigInteger totBases = BigInteger.valueOf(0);
    for (ContigBound contigBound : contigBounds) {
      totBases = totBases.add(BigInteger.valueOf(contigBound.getUpperBound()));
    }
    System.out.println("Total Number of Bases : " + totBases.toString());
    return contigBounds;
  }


  public static void helperCreateDirectory(File theDir) {
    // if the directory does not exist, create it
    if (!theDir.exists()) {
      System.out.println("creating directory: " + theDir);
      theDir.mkdir();
    }
  }


}
