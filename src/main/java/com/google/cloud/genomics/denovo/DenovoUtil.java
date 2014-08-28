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

import org.javatuples.Triplet;

import java.io.File;
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
  static public Map<String, Float> qualityThresholdMap = new HashMap<>();
  static public Map<Triplet<Genotype, Genotype, Genotype>, Boolean> isDenovoMap = new HashMap<>();

  static {
    // Constant Values Needed for stage 2 experiments
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

  public enum Caller { VARIANT, READ }
  
  public enum TrioIndividual {
    DAD, MOM, CHILD;

    public static final EnumSet<TrioIndividual> PARENTS = EnumSet.of(DAD, MOM);
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
      switch (all) {
        case A:
          return G;
        case G:
          return A;
        case C:
          return T;
        case T:
          return C;
        default:
          throw new IllegalArgumentException("Unknown haplotype " + all);
      }
    }

    public EnumSet<Allele> getTransition(Allele all) {
      switch (all) {
        case A:
          return EnumSet.of(C, T);
        case G:
          return EnumSet.of(C, T);
        case C:
          return EnumSet.of(A, G);
        case T:
          return EnumSet.of(A, G);
        default:
          throw new IllegalArgumentException("Unknown haplotype " + all);
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

    private enum Zygosity {
      HOMOZYGOUS, HETEROZYGOUS
    }

    private final Zygosity zygosity;
    private final EnumSet<Allele> alleles;

    Genotype(Zygosity zygosity, Allele a) {
      this.zygosity = zygosity;
      alleles = EnumSet.of(a);
    }

    Genotype(Zygosity zygosity, Allele a, Allele b) {
      this.zygosity = zygosity;
      alleles = EnumSet.of(a, b);
    }

    /*
     * Return a genotype from a pair of allele objects
     */
    public static Genotype valueOfPairAlleles(Allele a, Allele b) {
      Allele[] allelePair = new Allele[] {a, b};
      Arrays.sort(allelePair);
      return valueOf(allelePair[0].toString() + allelePair[1].toString());
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
        if (inferMethod.name().toLowerCase().equals(method.toLowerCase())) {
          return inferMethod;
        }
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
}
