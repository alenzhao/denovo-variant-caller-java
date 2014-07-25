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

import com.google.cloud.genomics.denovo.DenovoUtil.Genotypes;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;

import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.CHILD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.DAD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.MOM;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/*
 * DenovoBayesNet implements abstract BayesNet
 */
public class DenovoBayesNet implements BayesNet<TrioIndividual, Genotypes> {

  private double sequenceErrorRate;
  private double denovoMutationRate;
  private Map<TrioIndividual, Node<TrioIndividual, Genotypes>> nodeMap;

  public DenovoBayesNet(double sequenceErrorRate, double denovoMutationRate) {
    setNodeMap(new HashMap<TrioIndividual, Node<TrioIndividual, Genotypes>>());
    setSequenceErrorRate(sequenceErrorRate);
    setDenovoMutationRate(denovoMutationRate);
  }

  @Override
  public void addNode(Node<TrioIndividual, Genotypes> node) {
    getNodeMap().put(node.getId(), node);
  }

  /*
   * Prints the conditional probability to terminal
   */
  /**
   * @param conditionalProbabilityTable
   */
  public static void printConditionalProbabilityTable(PrintStream out,
      Map<List<Genotypes>, Double> conditionalProbabilityTable) {

    for (Genotypes dadGenotype : Genotypes.values()) {
      for (Genotypes momGenotype : Genotypes.values()) {
        for (Genotypes childGenotype : Genotypes.values()) {
          List<Genotypes> cptKey = Arrays.asList(dadGenotype, momGenotype, childGenotype);
          Double probVal = conditionalProbabilityTable.get(cptKey);
          out.printf("%s : %s%n", cptKey, probVal);
        }
      }
    }
  }

  /*
   * Creates the conditional probability table to be used in the bayes net One each for mom, dad and
   * child
   */
  /**
   * @param individual
   * @return conditionalProbabilityTable
   */
  Map<List<Genotypes>, Double> createConditionalProbabilityTable(TrioIndividual individual) {

    Map<List<Genotypes>, Double> conditionalProbabilityTable = new HashMap<>();
    if (individual == DAD || individual == MOM) {
      for (Genotypes genoType : Genotypes.values()) {
        conditionalProbabilityTable.put(Collections.singletonList(genoType),
            1.0 / Genotypes.values().length);
      }
    } else { // individual == TrioIndividuals.CHILD

      // Loops over parent Genotypes
      for (Genotypes genoTypeDad : Genotypes.values()) {
        for (Genotypes genoTypeMom : Genotypes.values()) {

          int validInheritanceCases = 0;
          // Initial pass
          for (Genotypes genoTypeChild : Genotypes.values()) {
            double prob;
            String dadAlleles = genoTypeDad.name();
            String momAlleles = genoTypeMom.name();
            String childAlleles = genoTypeChild.name();
            String c1 = childAlleles.substring(0, 1);
            String c2 = childAlleles.substring(1, 2);
            boolean predicate1 = momAlleles.contains(c1) & dadAlleles.contains(c2);
            boolean predicate2 = momAlleles.contains(c2) & dadAlleles.contains(c1);
            boolean predicate3 = predicate1 | predicate2;
            
            boolean isDenovo = DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(genoTypeDad,genoTypeMom,genoTypeChild));
            
            if (predicate3) {
              prob = 1.0;
              validInheritanceCases++;
            } else {
              prob = 0.0;
            }

            conditionalProbabilityTable.put(Arrays.asList(genoTypeDad, genoTypeMom, genoTypeChild),
                prob);
          }
          // Secondary Pass to normalize prob values
          for (Genotypes genoTypeChild : Genotypes.values()) {
            List<Genotypes> cptKey = Arrays.asList(genoTypeDad, genoTypeMom, genoTypeChild);

            boolean isNotInheritenceCase =
                -DenovoUtil.EPS <= conditionalProbabilityTable.get(cptKey)
                && conditionalProbabilityTable.get(cptKey) <= DenovoUtil.EPS;
            conditionalProbabilityTable.put(
                cptKey,
                isNotInheritenceCase
                    ? getDenovoMutationRate()
                    : 1.0 / validInheritanceCases
                        - getDenovoMutationRate()
                            * (Genotypes.values().length - validInheritanceCases) / validInheritanceCases);
          }

          // Sanity check - probabilities should add up to 1.0 (almost)
          double totProb = 0.0;
          for (Genotypes genoTypeChild : Genotypes.values()) {
            List<Genotypes> cptKey = Arrays.asList(genoTypeDad, genoTypeMom, genoTypeChild);
            totProb += conditionalProbabilityTable.get(cptKey);
          }
          if (Math.abs(totProb - 1.0) > DenovoUtil.EPS) {
            throw new IllegalStateException(
                "cpt probabilities not adding up : " + String.valueOf(totProb));
          }
        }
      }
    }
    return conditionalProbabilityTable;
  }

  /**
   * Get the log likelihood for a particular read base
   *
   * @param genoType
   * @param isHomozygous
   * @param base
   * @return logLikeliHood
   */
  public double getBaseLikelihood(Genotypes genoType, boolean isHomozygous, String base) {
    return isHomozygous
        ? genoType.name().contains(base)
            ? Math.log(1 - getSequenceErrorRate())
            : Math.log(getSequenceErrorRate()) - Math.log(3)
        : genoType.name().contains(base)
            ? Math.log(1 - 2 * getSequenceErrorRate() / 3) - Math.log(2)
            : Math.log(getSequenceErrorRate()) - Math.log(3);
  }

  /**
   * Get the log likelihood of reads for a particular individual in a trio for all their possible
   * Genotypes
   *
   * @param readSummaryMap
   * @return individualLogLikelihood
   */
  public Map<TrioIndividual, Map<Genotypes, Double>> getIndividualLogLikelihood(
      Map<TrioIndividual, ReadSummary> readSummaryMap) {
    Map<TrioIndividual, Map<Genotypes, Double>> individualLogLikelihood = new HashMap<>();
    for (TrioIndividual trioIndividual : TrioIndividual.values()) {

      ReadSummary readSummary = readSummaryMap.get(trioIndividual);
      Map<Genotypes, Double> genoTypeLogLikelihood = getGenoTypeLogLikelihood(readSummary);
      individualLogLikelihood.put(trioIndividual, genoTypeLogLikelihood);
    }
    return individualLogLikelihood;
  }

  /**
   * Get the log likelihood for all possible Genotypes for a set of reads
   *
   * @param readSummary
   * @return genotypeLogLikelihood
   */
  public Map<Genotypes, Double> getGenoTypeLogLikelihood(ReadSummary readSummary) {
    Map<Genotypes, Double> genotypeLogLikelihood = new HashMap<>();
    for (Genotypes genoType : Genotypes.values()) {
      Map<String, Integer> count = readSummary.getCount();
      boolean isHomozygous =
          genoType.name().substring(0, 1).equals(genoType.name().substring(1, 2));

      double readlogLikelihood = 0.0;
      for (String base : count.keySet()) {
        readlogLikelihood += getBaseLikelihood(genoType, isHomozygous, base);
      }
      genotypeLogLikelihood.put(genoType, readlogLikelihood);
    }
    return genotypeLogLikelihood;
  }

  /**
   * Infer the most likely genotype given likelihood for all the Genotypes of the trio
   *
   * @param individualLogLikelihood
   * @return maxgenoType
   */
  public List<Genotypes> getMaxGenoType(
      Map<TrioIndividual, Map<Genotypes, Double>> individualLogLikelihood) {
    double maxLogLikelihood = Double.NEGATIVE_INFINITY;
    List<Genotypes> maxGenoType = null;

    // Calculate overall bayes net log likelihood
    for (Genotypes genoTypeDad : Genotypes.values()) {
      for (Genotypes genoTypeMom : Genotypes.values()) {
        for (Genotypes genoTypeChild : Genotypes.values()) {
          double logLikelihood = 0.0;

          /* Get likelihood from the reads */
          logLikelihood += individualLogLikelihood.get(DAD).get(genoTypeDad);
          logLikelihood += individualLogLikelihood.get(MOM).get(genoTypeMom);
          logLikelihood += individualLogLikelihood.get(CHILD).get(genoTypeChild);

          /* Get likelihoods from the trio relationship */
          logLikelihood += getNodeMap().get(DAD).getConditionalProbabilityTable().get(
              Collections.singletonList(genoTypeDad));
          logLikelihood += getNodeMap().get(MOM).getConditionalProbabilityTable().get(
              Collections.singletonList(genoTypeMom));
          logLikelihood += getNodeMap().get(CHILD).getConditionalProbabilityTable().get(
              Arrays.asList(genoTypeDad, genoTypeMom, genoTypeChild));

          if (logLikelihood > maxLogLikelihood) {
            maxLogLikelihood = logLikelihood;
            maxGenoType = Arrays.asList(genoTypeDad, genoTypeMom, genoTypeChild);
          }
        }
      }
    }
    return maxGenoType;
  }

  /**
   * @return the sequenceErrorRate
   */
  public double getSequenceErrorRate() {
    return sequenceErrorRate;
  }

  /**
   * @param sequenceErrorRate the sequenceErrorRate to set
   */
  public void setSequenceErrorRate(double sequenceErrorRate) {
    this.sequenceErrorRate = sequenceErrorRate;
  }

  /**
   * @return the denovoMutationRate
   */
  public double getDenovoMutationRate() {
    return denovoMutationRate;
  }

  /**
   * @param denovoMutationRate the denovoMutationRate to set
   */
  public void setDenovoMutationRate(double denovoMutationRate) {
    this.denovoMutationRate = denovoMutationRate;
  }

  /**
   * @return the nodeMap
   */
  public Map<TrioIndividual, Node<TrioIndividual, Genotypes>> getNodeMap() {
    return nodeMap;
  }

  /**
   * @param hashMap the nodeMap to set
   */
  public void setNodeMap(Map<TrioIndividual, Node<TrioIndividual, Genotypes>> hashMap) {
    this.nodeMap = hashMap;
  }
}
