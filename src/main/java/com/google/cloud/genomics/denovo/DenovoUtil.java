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

import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.CHILD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.DAD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.MOM;

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
import com.google.common.base.Optional;

import java.io.File;
import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/*
 * Utility functions shared by other classes in Denovo project
 */
public class DenovoUtil {

  public static final double EPS = 1e-12;
  static public final long PROJECT_ID = 1085016379660L;
  static public final int TOT_CHROMOSOMES = 24;
  static public final long MAX_VARIANT_RESULTS = 10000L;
  static public final long DEFAULT_START_POS = 1L;
  static public final Float GQX_THRESH = Float.valueOf((float) 30.0);
  static public final Float QD_THRESH = Float.valueOf((float) 2.0);
  static public final Float MQ_THRESH = Float.valueOf((float) 20.0);
  static public final String TRIO_DATASET_ID = "2315870033780478914";

  static public Map<String, Float> qualityThresholdMap = new HashMap<>();
  public static Map<TrioIndividual, String> datasetIdMap = new HashMap<>();
  public static Map<TrioIndividual, String> callsetIdMap = new HashMap<>();  
  static public Map<TrioIndividual, String> individualCallsetNameMap = new HashMap<>();
  
  static public int debugLevel = 0; 
  
  static {
    // Constant Values Needed for stage 2 experiments
    datasetIdMap.put(DAD, "4140720988704892492");
    datasetIdMap.put(MOM, "2778297328698497799");
    datasetIdMap.put(CHILD, "6141326619449450766");
    datasetIdMap = Collections.unmodifiableMap(datasetIdMap);

    callsetIdMap.put(DAD, "NA12877");
    callsetIdMap.put(MOM, "NA12878");
    callsetIdMap.put(CHILD, "NA12879");
    callsetIdMap = Collections.unmodifiableMap(callsetIdMap);
    
    qualityThresholdMap.put("GQX", GQX_THRESH);
    qualityThresholdMap.put("QD", QD_THRESH);
    qualityThresholdMap.put("MQ", MQ_THRESH);
    qualityThresholdMap = Collections.unmodifiableMap(qualityThresholdMap);

    individualCallsetNameMap.put(DAD, "NA12877");
    individualCallsetNameMap.put(MOM, "NA12878");
    individualCallsetNameMap.put(CHILD, "NA12879");
    individualCallsetNameMap = Collections.unmodifiableMap(individualCallsetNameMap);
  }
  
  public enum TrioIndividual {
    DAD, MOM, CHILD;
  }

  public enum Haplotype {
    A, C, G, T;
    
    public static final EnumSet<Haplotype> allHaplotypes = EnumSet.allOf(Haplotype.class);
    
    public EnumSet<Haplotype> getMutants() {
      EnumSet<Haplotype> difference = allHaplotypes.clone();
      difference.remove(this);
      return difference;
    }
    
    public Haplotype getTransversion(Haplotype hap) {
      switch(hap) {
        case A : return G;
        case G : return A;
        case C : return T;
        case T : return C;
        default : throw new IllegalArgumentException("Unknown haplotype " + hap);
      }
    }
    
    public EnumSet<Haplotype> getTransition(Haplotype hap) {
      switch(hap) {
        case A : return EnumSet.of(C, T);
        case G : return EnumSet.of(C, T);
        case C : return EnumSet.of(A, G);
        case T : return EnumSet.of(A, G);
        default : throw new IllegalArgumentException("Unknown haplotype " + hap);
      }
    }
  }
  
  public enum Genotypes {
    AA(true),
    AC(false),
    AT(false),
    AG(false),
    CC(true),
    CT(false),
    CG(false),
    TT(true),
    TG(false),
    GG(true);

    private final boolean isDiploid;

    Genotypes(boolean isDiploid) {
      this.isDiploid = isDiploid;
    }

    public boolean isDiploid() {
      return isDiploid;
    }
  }

  public static SearchVariantsRequest createSearchVariantsRequest(SearchVariantsRequest oldRequest,
      ContigBound contig,
      long startPos,
      long endPos,
      String datasetId,
      String nextPageToken,
      long maxVariantResults) {

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
        .setMaxResults(BigInteger.valueOf(maxVariantResults))
        .setPageToken(nextPageToken);

    return searchVariantsRequest;
  }

  private static SearchReadsRequest createSearchReadsRequest(String readsetId,
      String chromosomeName, long startPos, long endPos) {

    return new SearchReadsRequest().setReadsetIds(Collections.singletonList(readsetId))
        .setSequenceName(chromosomeName).setSequenceStart(BigInteger.valueOf(startPos))
        .setSequenceEnd(BigInteger.valueOf(endPos));
  }

  private static SearchReadsetsRequest createSearchReadsetsRequest(String datasetId) {

    // Init searchRequest obj
    return new SearchReadsetsRequest().setDatasetIds(Collections.singletonList(datasetId));
  }

  private static SearchCallsetsRequest createSearchCallsetsRequest(String datasetId) {

    // Init searchRequest obj
    return new SearchCallsetsRequest().setDatasetIds(Collections.singletonList(datasetId));
  }

  /**
   * @param datasetIdMap
   * @param callsetIdMap
   * @return Map<String, String>
   * @throws IOException
   */
  public static Map<TrioIndividual, String> createReadsetIdMap(
      Map<TrioIndividual, String> datasetIdMap, Map<TrioIndividual, String> callsetIdMap, 
      Genomics genomics)
      throws IOException {
    Map<TrioIndividual, String> readsetIdMap = new HashMap<>();

    for (TrioIndividual trioIndividual : datasetIdMap.keySet()) {
      List<Readset> readsets = getReadsets(datasetIdMap.get(trioIndividual), genomics);

      for (Readset readset : readsets) {
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
      throw new IllegalStateException("Borked readsetIdMap" + readsetIdMap);
    }
    return readsetIdMap;
  }

  public static Optional<List<Integer>> getGenotypeFromInfoField(Call call) {
    String genoTypeString = call.getInfo().get("GT").get(0);
    String splitChar = null;
    if (genoTypeString.contains("/")) {
      splitChar = "/";
    } else if (genoTypeString.contains("|")) {
      splitChar = "\\|";
    } else {
      return Optional.absent();
    }

    List<Integer> genoType = new ArrayList<>();
    for (String allele : genoTypeString.split(splitChar)) {
      genoType.add(Integer.valueOf(allele));
    }
    return Optional.of(genoType);
  }

  public static List<Read> getReads(String readsetId, String chromosomeName, long startPos,
      long endPos, Genomics genomics) throws IOException {

    SearchReadsRequest searchReadsRequest =
        createSearchReadsRequest(readsetId, chromosomeName, startPos, endPos);

    SearchReadsResponse searchReadsExecuted =
        genomics.reads().search(searchReadsRequest).execute();

    List<Read> reads = searchReadsExecuted.getReads();
    return reads;
  }

  public static List<Readset> getReadsets(String datasetId, Genomics genomics) throws IOException {
    SearchReadsetsRequest searchReadsetsRequest = createSearchReadsetsRequest(datasetId);
    Genomics.Readsets.Search search =
        genomics.readsets().search(searchReadsetsRequest);
    SearchReadsetsResponse execute = search.execute();
    List<Readset> readsets = execute.getReadsets();

    return readsets;
  }

  /**
   * @return List<Callset>
   * @throws IOException
   */
  public static List<Callset> getCallsets(String datasetId, Genomics genomics) throws IOException {
    SearchCallsetsRequest searchCallsetsRequest = createSearchCallsetsRequest(datasetId);
    Genomics.Callsets.Search search =
        genomics.callsets().search(searchCallsetsRequest);
    SearchCallsetsResponse execute = search.execute();
    List<Callset> callsets = execute.getCallsets();

    return callsets;
  }

  /**
   * @return List<Dataset>
   * @throws IOException
   */
  public static List<Dataset> getAllDatasets(Genomics genomics) throws IOException {

    // Get a list of all the datasets associated with project id
    Genomics.Datasets.List datasetRequest =
        genomics.datasets().list().setProjectId(PROJECT_ID);
    datasetRequest.setDisableGZipContent(true);

    ListDatasetsResponse execute = datasetRequest.execute();
    List<Dataset> datasets = execute.getDatasets();

    return datasets;
  }

  /**
   * @return List<ContigBound>
   * @throws IOException
   */
  public static List<ContigBound> getVariantsSummary(String datasetId, Genomics genomics) 
      throws IOException {

    Genomics.Variants.GetSummary variantsSummaryRequest =
        genomics.variants().getSummary().setDatasetId(datasetId);
    variantsSummaryRequest.setDisableGZipContent(true);

    GetVariantsSummaryResponse execute = variantsSummaryRequest.execute();
    if (debugLevel >= 2) {
      System.out.println("Variants : " + execute.toPrettyString());
    }

    List<ContigBound> contigBounds = execute.getContigBounds();
    BigInteger totBases = BigInteger.valueOf(0);
    for (ContigBound contigBound : contigBounds) {
      totBases = totBases.add(BigInteger.valueOf(contigBound.getUpperBound()));
    }
    return contigBounds;
  }

  public static void helperCreateDirectory(File theDir) {
    // if the directory does not exist, create it
    if (!theDir.exists()) {
      System.err.println("creating directory: " + theDir);
      theDir.mkdir();
    }
  }

  /**
   * Check if the particular genotype is denovo i.e. present in kids but not in parents
   *
   * @param trioGenotypeList
   * @return isDenovo
   */
  public static boolean checkTrioGenoTypeIsDenovo(List<Genotypes> trioGenotypeList) {
    Genotypes genoTypeDad = trioGenotypeList.get(0);
    Genotypes genoTypeMom = trioGenotypeList.get(1);
    Genotypes genoTypeChild = trioGenotypeList.get(2);

    String childAlleles = genoTypeChild.name();
    String momAlleles = genoTypeMom.name();
    String dadAlleles = genoTypeDad.name();

    String c1 = childAlleles.substring(0, 1);
    String c2 = childAlleles.substring(1, 2);
    boolean predicate1 = momAlleles.contains(c1) & dadAlleles.contains(c2);
    boolean predicate2 = momAlleles.contains(c2) & dadAlleles.contains(c1);
    boolean predicate3 = !(predicate1 | predicate2);
    return predicate3;
  }
}
