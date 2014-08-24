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
import com.google.api.services.genomics.model.Call;
import com.google.api.services.genomics.model.Callset;
import com.google.api.services.genomics.model.ContigBound;
import com.google.api.services.genomics.model.GetVariantsSummaryResponse;
import com.google.api.services.genomics.model.Read;
import com.google.api.services.genomics.model.Readset;
import com.google.api.services.genomics.model.SearchCallsetsRequest;
import com.google.api.services.genomics.model.SearchCallsetsResponse;
import com.google.api.services.genomics.model.SearchReadsRequest;
import com.google.api.services.genomics.model.SearchReadsResponse;
import com.google.api.services.genomics.model.SearchReadsetsRequest;
import com.google.api.services.genomics.model.SearchReadsetsResponse;
import com.google.api.services.genomics.model.SearchVariantsRequest;
import com.google.api.services.genomics.model.Variant;
import com.google.cloud.genomics.denovo.DenovoUtil.Allele;
import com.google.common.base.Optional;

import org.javatuples.Triplet;

import java.io.File;
import java.io.IOException;
import java.math.BigInteger;
import java.util.Arrays;
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
  static public final int TOT_CHROMOSOMES = 24;
  static public final long MAX_VARIANT_RESULTS = 10000L;
  static public final Float GQX_THRESH = Float.valueOf((float) 30.0);
  static public final Float QD_THRESH = Float.valueOf((float) 2.0);
  static public final Float MQ_THRESH = Float.valueOf((float) 20.0);
  public static final int MAX_API_RETRIES = 5;
  public static final long API_WAIT_MILLISEC = 5000;

  public static double LRT_THRESHOLD = 1.0;
  static public Map<String, Float> qualityThresholdMap = new HashMap<>();
  static public Map<Triplet<Genotype, Genotype, Genotype>, Boolean> isDenovoMap = 
      new HashMap<>();
  
  static public int debugLevel = 0; 
  
  static {
    // Constant Values Needed for stage 2 experiments
    qualityThresholdMap.put("GQX", GQX_THRESH);
    qualityThresholdMap.put("QD", QD_THRESH);
    qualityThresholdMap.put("MQ", MQ_THRESH);
    qualityThresholdMap = Collections.unmodifiableMap(qualityThresholdMap);

    for (Genotype genotypeDad : Genotype.values()) {
      for (Genotype genotypeMom : Genotype.values()) {
        for (Genotype genotypeChild : Genotype.values()) {
          String childAlleles = genotypeChild.name();
          String momAlleles = genotypeMom.name();
          String dadAlleles = genotypeDad.name();

          String c1 = childAlleles.substring(0, 1);
          String c2 = childAlleles.substring(1, 2);
          boolean alleleInParents = momAlleles.contains(c1) & dadAlleles.contains(c2);
          boolean alleleInParentsMirror = momAlleles.contains(c2) & dadAlleles.contains(c1);
          boolean isDenovo = !(alleleInParents || alleleInParentsMirror);
          isDenovoMap.put(Triplet.with(genotypeDad, genotypeMom, genotypeChild), isDenovo);
        }
      }
    }
  }
  
  public enum TrioIndividual {
    DAD, MOM, CHILD;
  }

  public enum Allele {
    A, C, G, T;
    
    public static final EnumSet<Allele> allHaplotypes = EnumSet.allOf(Allele.class);
    
    public EnumSet<Allele> getMutants() {
      EnumSet<Allele> difference = allHaplotypes.clone();
      difference.remove(this);
      return difference;
    }
    
    public Allele getTransversion(Allele all) {
      switch(all) {
        case A : return G;
        case G : return A;
        case C : return T;
        case T : return C;
        default : throw new IllegalArgumentException("Unknown haplotype " + all);
      }
    }
    
    public EnumSet<Allele> getTransition(Allele all) {
      switch(all) {
        case A : return EnumSet.of(C, T);
        case G : return EnumSet.of(C, T);
        case C : return EnumSet.of(A, G);
        case T : return EnumSet.of(A, G);
        default : throw new IllegalArgumentException("Unknown haplotype " + all);
      }
    }
  }
  
  public enum Genotype {
    AA(Zygosity.HOMOZYGOUS, Allele.A),
    AC(Zygosity.HETEROZYGOUS, Allele.A, Allele.C),
    AG(Zygosity.HETEROZYGOUS, Allele.A, Allele.G),
    AT(Zygosity.HETEROZYGOUS, Allele.A, Allele.T),
    CC(Zygosity.HOMOZYGOUS, Allele.C),
    CG(Zygosity.HETEROZYGOUS, Allele.C, Allele.G),
    CT(Zygosity.HETEROZYGOUS, Allele.C, Allele.T),
    GG(Zygosity.HOMOZYGOUS, Allele.G),
    GT(Zygosity.HETEROZYGOUS, Allele.G, Allele.T),
    TT(Zygosity.HOMOZYGOUS, Allele.T);

    private enum Zygosity {HOMOZYGOUS, HETEROZYGOUS}

    private final Zygosity zygosity;
    private final EnumSet<Allele> alleles;
    
    Genotype(Zygosity zygosity, Allele a) {
      this.zygosity = zygosity;
      alleles = EnumSet.of(a);
    }

    Genotype(Zygosity zygosity, Allele a, Allele b) {
      this.zygosity= zygosity;
      alleles = EnumSet.of(a,b);
    }

    /*
     * Return a genotype from a pair of allele objects
     */
    public static Genotype valueOfPairAlleles(Allele a, Allele b) {
      Allele[] allelePair = new Allele[]{a, b};
      Arrays.sort(allelePair);
      return valueOf(allelePair[0].toString()+allelePair[1].toString());
    }

    public EnumSet<Allele> getAlleleSet() {
      return alleles;
    }
    
    public boolean isHomozygous() {
      return zygosity == Zygosity.HOMOZYGOUS;
    }

    /**
     * @param base the base to check for 
     * @return whether this contains the allele
     */
    public boolean containsAllele(Allele base) {
      return alleles.contains(base);
    }
  }

  public enum InferenceMethod {
    MAP, BAYES, LRT;
    
    public static InferenceMethod selectMethodfromString(String method) {
      for (InferenceMethod inferMethod : InferenceMethod.values()) {
        if(inferMethod.name().toLowerCase().equals(method.toLowerCase())) { return inferMethod; }  
      }
      throw new IllegalArgumentException("Unknown method " + method);
    }
  }
  
  /*
   * Reverse a dictionary
   */
  public static <K, V> Map<V, K> getReversedMap(Map<K, V> map) {
    Map<V, K> reversed = new HashMap<>();
    for (Map.Entry<K, V> entry : map.entrySet()) {
      reversed.put(entry.getValue(), entry.getKey());
    }
    return reversed;
  }
  
  public static SearchVariantsRequest createSearchVariantsRequest(SearchVariantsRequest oldRequest,
      String contig,
      long startPos,
      long endPos,
      String datasetId,
      String nextPageToken,
      long maxVariantResults,
      List<String> callsetIds) {

    // Init searchRequest obj
    SearchVariantsRequest searchVariantsRequest;
    if (oldRequest == null) {
      searchVariantsRequest = new SearchVariantsRequest();
    } else {
      searchVariantsRequest = oldRequest;
    }

    searchVariantsRequest
        .setContig(contig)
        .setDatasetId(datasetId)
        .setStartPosition(startPos)
        .setEndPosition(endPos)
        .setMaxResults(BigInteger.valueOf(maxVariantResults))
        .setPageToken(nextPageToken)
        .setCallsetIds(callsetIds);

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

  public static Optional<List<Integer>> getGenotype(Call call) {
    List<Integer> genoType = call.getGenotype();

    if (genoType.contains(-1)) {
      return Optional.absent();
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

  /*
   * overloads and forwards to function
   */
  public static boolean checkTrioGenoTypeIsDenovo(Genotype genotypeDad, Genotype genotypeMom, 
      Genotype genotypeChild) {
    return isDenovoMap.get(Triplet.with(genotypeDad, genotypeMom, genotypeChild));
  }
  
  /**
   * Check if the particular genotype is denovo i.e. present in kids but not in parents
   *
   * @param trioGenotypeList
   * @return isDenovo
   */
  public static boolean checkTrioGenoTypeIsDenovo(List<Genotype> trioGenotypeList) {
    Genotype genotypeDad = trioGenotypeList.get(0);
    Genotype genotypeMom = trioGenotypeList.get(1);
    Genotype genotypeChild = trioGenotypeList.get(2);
    return checkTrioGenoTypeIsDenovo(genotypeDad, genotypeMom, genotypeChild);
  }
  
  /*
   * Does the variant not contain any variations
   */
  public static boolean isConserved(Variant variant) {
    return variant.getInfo() != null && variant.getInfo().containsKey("BLOCKAVG_min30p3a");
  }
  
  public static Call getCallInVariant(Variant variant, TrioIndividual person, 
      Map<TrioIndividual, String> personToCallsetId) {
    for(Call call: variant.getCalls()) {
      if (personToCallsetId.get(person).equals(call.getCallsetId())) {
        return call;
      }
    }
    return null;
  }
}
