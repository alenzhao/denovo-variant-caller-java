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
  public static Map<Triplet<Genotype, Genotype, Genotype>, Boolean> isDenovoMap = new HashMap<>();

  static {
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
    isDenovoMap = Collections.unmodifiableMap(isDenovoMap);
  }

  /**
   * Type of Caller. Variant bases or read based (more expensive)
   */
  public enum Caller { VARIANT, READ, FULL }
  
  public enum Chromosome {
    CHR1,
    CHR2,
    CHR3,
    CHR4,
    CHR5,
    CHR6,
    CHR7,
    CHR8,
    CHR9,
    CHR10,
    CHR11,
    CHR12,
    CHR13,
    CHR14,
    CHR15,
    CHR16,
    CHR17,
    CHR18,
    CHR19,
    CHR20,
    CHR21,
    CHR22,
    CHRX,
    CHRY,
    CHRM;
    
    public static final EnumSet<Chromosome> ALL = EnumSet.allOf(Chromosome.class);
    
    public static Chromosome fromString(String str) {
      str = str.toUpperCase();
      try {
        return Chromosome.valueOf(str);
      } catch (IllegalArgumentException e) {
        for(int i = 1 ; i <= 22 ; i++) {
          if (str.equals(String.valueOf(i))) {
            return Chromosome.valueOf("CHR" + String.valueOf(i));
          }
        }
        if (Arrays.asList("X","Y","M").contains(str.toUpperCase())) {
          return Chromosome.valueOf("CHR" + str);
        }
        throw new IllegalArgumentException("Unknown Chromosome " + str);
      }
    }
  }

  /**
   * Mebers in the trio
   */
  public enum TrioMember {
    CHILD, MOM, DAD;

    public static final EnumSet<TrioMember> PARENTS = EnumSet.of(DAD, MOM);
  }

  public enum Allele {
    A, C, G, T;

    public static final EnumSet<Allele> allHaplotypes = EnumSet.allOf(Allele.class);

    /**
     * @return set of alleles that are different from this allele
     */
    public EnumSet<Allele> getMutants() {
      EnumSet<Allele> difference = allHaplotypes.clone();
      difference.remove(this);
      return difference;
    }

    /**
     * @param allele 
     * @return transverse allele
     */
    public Allele getTransversion(Allele allele) {
      switch (allele) {
        case A:
          return G;
        case G:
          return A;
        case C:
          return T;
        case T:
          return C;
        default:
          throw new IllegalArgumentException("Unknown haplotype " + allele);
      }
    }

    /**
     * @param allele
     * @return transition allele
     */
    public EnumSet<Allele> getTransition(Allele allele) {
      switch (allele) {
        case A:
          return EnumSet.of(C, T);
        case G:
          return EnumSet.of(C, T);
        case C:
          return EnumSet.of(A, G);
        case T:
          return EnumSet.of(A, G);
        default:
          throw new IllegalArgumentException("Unknown haplotype " + allele);
      }
    }
  }

  /**
   * All genotypes present by combining pairs of alleles
   */
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

    /**
     * @param a
     * @param b
     * @return a genotype from a pair of allele objects
     */
    public static Genotype valueOfPairAlleles(Allele a, Allele b) {
      Allele[] allelePair = new Allele[] {a, b};
      Arrays.sort(allelePair);
      return valueOf(allelePair[0].toString() + allelePair[1].toString());
    }

    /**
     * @return associated alleles with genotype
     */
    public EnumSet<Allele> getAlleleSet() {
      return alleles;
    }

    /**
     * @return zygosity of genotype
     */
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

  /**
   * Different Bayesian Inference Methods
   */
  public enum InferenceMethod {
    MAP {
      @Override
      boolean isDenovo(BayesInferenceResult result, DenovoShared shared) {
        return checkTrioGenoTypeIsDenovo(result.getMaxTrioGenotype());
        }
    }, 
    BAYES {
      @Override
      boolean isDenovo(BayesInferenceResult result, DenovoShared shared) {
        return result.getBayesDenovoProb() > 0.5;
      }
    }, 
    LRT {
      @Override
      boolean isDenovo(BayesInferenceResult result, DenovoShared shared) {
        return result.getLikelihoodRatio() > shared.getLrtThreshold();
      }
    };

    /**
     * Is the call denovo based on evidence from inference
     * @param result the result from bayes inference step
     * @param shared shared state in project
     * @return is the call denovo
     */
    abstract boolean isDenovo(BayesInferenceResult result, DenovoShared shared);
  }

  /**
   * Reverse a dictionary
   * @param map
   * @return revered map
   */
  static <K, V> Map<V, K> getReversedMap(Map<K, V> map) {
    Map<V, K> reversed = new HashMap<>();
    for (Map.Entry<K, V> entry : map.entrySet()) {
      reversed.put(entry.getValue(), entry.getKey());
    }
    return reversed;
  }

  static File getNormalizedFile(String fpath) {
    return new File(fpath).isAbsolute() 
        ? new File(fpath)
        : new File(System.getProperty("user.dir"), fpath);
  }

  /**
   * Check if the particular genotype is denovo i.e. present in kids but not in parents
   * @param genotypeDad
   * @param genotypeMom
   * @param genotypeChild
   * @return is it denovo
   */
  public static boolean checkTrioGenoTypeIsDenovo(Genotype genotypeDad, Genotype genotypeMom,
      Genotype genotypeChild) {
    return isDenovoMap.get(Triplet.with(genotypeDad, genotypeMom, genotypeChild));
  }

  /**
   * Overloads and forwards to {@link #checkTrioGenoTypeIsDenovo(Genotype, Genotype, Genotype) checkTrioTypeIsDenovo}
   * @param trioGenotypeList
   * @return is teh trio of genotypes denovo
   */
  static boolean checkTrioGenoTypeIsDenovo(List<Genotype> trioGenotypeList) {
    Genotype genotypeDad = trioGenotypeList.get(0);
    Genotype genotypeMom = trioGenotypeList.get(1);
    Genotype genotypeChild = trioGenotypeList.get(2);
    return checkTrioGenoTypeIsDenovo(genotypeDad, genotypeMom, genotypeChild);
  }
}
