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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/*
 * DenovoBayesNet implements abstract BayesNet
 */
public class DenovoBayesNet implements BayesNet<TrioIndividual, Genotypes> {

  private final double sequenceErrorRate;
  private final double denovoMutationRate;
  private Map<TrioIndividual, Node<TrioIndividual, Genotypes>> nodeMap;

  public DenovoBayesNet(double sequenceErrorRate, double denovoMutationRate) {
    nodeMap = new HashMap<TrioIndividual, Node<TrioIndividual, Genotypes>>();
    this.sequenceErrorRate = sequenceErrorRate;
    this.denovoMutationRate = denovoMutationRate; 
    
    // Initialize the nodes and prob tables
    initializeTrioNodes();
  }

  /**
   * 
   */
  private void initializeTrioNodes() {
    addNode(new Node<>(DAD, null, createConditionalProbabilityTable(DAD)));
    addNode(new Node<>(MOM, null, createConditionalProbabilityTable(MOM)));
    List<Node<TrioIndividual, Genotypes>> childParents = new ArrayList<>();
    childParents.add(nodeMap.get(DAD));
    childParents.add(nodeMap.get(MOM));
    addNode(new Node<>(CHILD, childParents, createConditionalProbabilityTable(CHILD)));
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
  public Map<List<Genotypes>, Double> createConditionalProbabilityTable(TrioIndividual individual) {

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


          // Initial pass to count valid inheritance cases
          int validInheritanceCases = 0;
          for (Genotypes genoTypeChild : Genotypes.values()) {
            boolean isDenovo = DenovoUtil.checkTrioGenoTypeIsDenovo(
                Arrays.asList(genoTypeDad, genoTypeMom, genoTypeChild));

            if (!isDenovo) {
              validInheritanceCases++;
            }
          }

          // Secondary Pass to fill in probability values
          for (Genotypes genoTypeChild : Genotypes.values()) {
            List<Genotypes> cptKey = Arrays.asList(genoTypeDad, genoTypeMom, genoTypeChild);

            boolean isDenovo = DenovoUtil.checkTrioGenoTypeIsDenovo(cptKey);
            conditionalProbabilityTable.put(
                cptKey,
                isDenovo
                    ? getDenovoMutationRate()
                    : (1.0
                          - getDenovoMutationRate() 
                               * (Genotypes.values().length - validInheritanceCases)) 
                      / validInheritanceCases
                );
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
      for (Map.Entry<String, Integer> entry : count.entrySet()) {
        readlogLikelihood += entry.getValue() * 
            getBaseLikelihood(genoType, isHomozygous, entry.getKey());
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
          logLikelihood += getLogLikelihoodFromCPT(DAD, genoTypeDad);
          logLikelihood += getLogLikelihoodFromCPT(MOM, genoTypeMom);
          logLikelihood += getLogLikelihoodFromCPT(CHILD, genoTypeDad, genoTypeMom, genoTypeChild); 
              
          if (logLikelihood > maxLogLikelihood) {
            maxLogLikelihood = logLikelihood;
            maxGenoType = Arrays.asList(genoTypeDad, genoTypeMom, genoTypeChild);
          }
        }
      }
    }
    return maxGenoType;
  }

  private double getLogLikelihoodFromCPT(TrioIndividual individual, Genotypes... keyParts) {
    List<Genotypes> cptKey = Arrays.asList(keyParts);
    return Math.log(getNodeMap().
        get(individual).
        getConditionalProbabilityTable().
        get(cptKey));
  }

  /**
   * @return the sequenceErrorRate
   */
  public double getSequenceErrorRate() {
    return sequenceErrorRate;
  }

  /**
   * @return the denovoMutationRate
   */
  public double getDenovoMutationRate() {
    return denovoMutationRate;
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
