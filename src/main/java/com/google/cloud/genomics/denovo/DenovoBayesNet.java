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

import com.google.cloud.genomics.denovo.DenovoUtil.Allele;
import com.google.cloud.genomics.denovo.DenovoUtil.Genotype;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;

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
public class DenovoBayesNet implements BayesNet<TrioIndividual, Genotype> {

  public static class InferenceResult {
    public final List<Genotype> maxTrioGenotype;
    public final double bayesDenovoProb;
    public final double logOfLikelihoodRatio;
    
    public InferenceResult(List<Genotype> maxTrioGenotype,
        double bayesDenovoProb,
        double logOfLikelihoodRatio) {
      this.maxTrioGenotype = maxTrioGenotype;
      this.bayesDenovoProb = bayesDenovoProb;
      this.logOfLikelihoodRatio = logOfLikelihoodRatio;
    }
  }

  private final double sequenceErrorRate;
  private final double denovoMutationRate;
  private Map<TrioIndividual, Node<TrioIndividual, Genotype>> nodeMap;

  public DenovoBayesNet(double sequenceErrorRate, double denovoMutationRate) {
    nodeMap = new HashMap<TrioIndividual, Node<TrioIndividual, Genotype>>();
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
    List<Node<TrioIndividual, Genotype>> childParents = new ArrayList<>();
    childParents.add(nodeMap.get(DAD));
    childParents.add(nodeMap.get(MOM));
    addNode(new Node<>(CHILD, childParents, createConditionalProbabilityTable(CHILD)));
  }

  @Override
  public void addNode(Node<TrioIndividual, Genotype> node) {
    getNodeMap().put(node.getId(), node);
  }

  /*
   * Prints the conditional probability to terminal
   */
  /**
   * @param conditionalProbabilityTable
   */
  public static void printConditionalProbabilityTable(PrintStream out,
      Map<List<Genotype>, Double> conditionalProbabilityTable) {

    for (Genotype dadGenotype : Genotype.values()) {
      for (Genotype momGenotype : Genotype.values()) {
        for (Genotype childGenotype : Genotype.values()) {
          List<Genotype> cptKey = Arrays.asList(dadGenotype, momGenotype, childGenotype);
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
  public Map<List<Genotype>, Double> createConditionalProbabilityTable(TrioIndividual individual) {

    Map<List<Genotype>, Double> conditionalProbabilityTable = new HashMap<>();
    if (individual == DAD || individual == MOM) {
      for (Genotype genoType : Genotype.values()) {
        conditionalProbabilityTable.put(Collections.singletonList(genoType),
            (genoType.isHomozygous() ? 1.0 : 2.0) / (Allele.values().length * Allele.values().length));
      }
    } else { // individual == TrioIndividuals.CHILD

      // Loops over parent Genotypes
      for (Genotype genoTypeDad : Genotype.values()) {
        for (Genotype genoTypeMom : Genotype.values()) {

          // Get a map of mendelian genotypes and their frequencies of occurence
          Map<Genotype, Integer> mendelianAlleles = mendelianGenotypes(genoTypeDad, genoTypeMom);
          int numDenovoGenotypes = Genotype.values().length - mendelianAlleles.size();
          Integer mendelianCount = 0;
          for(Integer count : mendelianAlleles.values() ) { 
            mendelianCount += count;
          }
          
          // Initial pass to count valid inheritance cases
          for (Genotype genoTypeChild : Genotype.values()) {
            List<Genotype> cptKey = Arrays.asList(genoTypeDad, genoTypeMom, genoTypeChild);
            boolean isDenovo = DenovoUtil.checkTrioGenoTypeIsDenovo(cptKey);
          
            double value = isDenovo
                ? getDenovoMutationRate() / numDenovoGenotypes
                : (1.0 - getDenovoMutationRate()) / mendelianCount 
                    * mendelianAlleles.get(genoTypeChild);
            conditionalProbabilityTable.put(
                cptKey,
                value
                );
          }
        }
      }
    }
    return conditionalProbabilityTable;
  }

  /*
   * Returns a Map containing mendelian Genotypes and their counts in a mendelian inheritance 
   * scenario
   */
  private Map<Genotype, Integer> mendelianGenotypes(Genotype genoTypeDad, Genotype genoTypeMom) {
    Map<Genotype, Integer> mendelCases = new HashMap<>();
    for (int ii = 0 ; ii < genoTypeDad.name().length() ; ii++) {
      for (int jj= 0 ; jj < genoTypeMom.name().length() ; jj++) {
        char[] mendelianCharAlleles = new char[2];
        mendelianCharAlleles[0] = genoTypeDad.name().charAt(ii) ;
        mendelianCharAlleles[1] = genoTypeMom.name().charAt(jj) ;
        Arrays.sort(mendelianCharAlleles);
        Genotype mendelianAlleles = 
            Genotype.getGenoTypeFromString(new String(mendelianCharAlleles));
        mendelCases.put(mendelianAlleles, 
            (mendelCases.containsKey(mendelianAlleles) 
                ? mendelCases.get(mendelianAlleles)
                : 0) + 1 );
      }
    }
    return mendelCases;
    
  }

  /**
   * Get the log likelihood for a particular read base
   *
   * @param genoType
   * @param base
   * @return logLikeliHood
   */
  public double getBaseLogLikelihood(Genotype genoType, String base) {
    return genoType.isHomozygous()
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
  public Map<TrioIndividual, Map<Genotype, Double>> getIndividualLogLikelihood(
      Map<TrioIndividual, ReadSummary> readSummaryMap) {
    Map<TrioIndividual, Map<Genotype, Double>> individualLogLikelihood = new HashMap<>();
    for (TrioIndividual trioIndividual : TrioIndividual.values()) {

      ReadSummary readSummary = readSummaryMap.get(trioIndividual);
      Map<Genotype, Double> genoTypeLogLikelihood = getGenoTypeLogLikelihood(readSummary);
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
  public Map<Genotype, Double> getGenoTypeLogLikelihood(ReadSummary readSummary) {
    Map<Genotype, Double> genotypeLogLikelihood = new HashMap<>();
    for (Genotype genoType : Genotype.values()) {
      Map<String, Integer> count = readSummary.getCount();

      double readlogLikelihood = 0.0;
      for (Map.Entry<String, Integer> entry : count.entrySet()) {
        readlogLikelihood += entry.getValue() * 
            getBaseLogLikelihood(genoType, entry.getKey());
      }
      genotypeLogLikelihood.put(genoType, readlogLikelihood);
    }
    return genotypeLogLikelihood;
  }

  private double getLogLikelihoodFromCPT(TrioIndividual individual, Genotype... keyParts) {
    List<Genotype> cptKey = Arrays.asList(keyParts);
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
  public Map<TrioIndividual, Node<TrioIndividual, Genotype>> getNodeMap() {
    return nodeMap;
  }

  /**
   * @param hashMap the nodeMap to set
   */
  public void setNodeMap(Map<TrioIndividual, Node<TrioIndividual, Genotype>> hashMap) {
    this.nodeMap = hashMap;
  }

  public DenovoBayesNet.InferenceResult performInference(
      Map<TrioIndividual, ReadSummary> readSummaryMap) {

    // Calculate Likelihoods of the different reads
    Map<TrioIndividual, Map<Genotype, Double>> individualLogLikelihood =
        getIndividualLogLikelihood(readSummaryMap);
    
    double maxLogLikelihood = Double.NEGATIVE_INFINITY;
    double denovoLikelihood = 0.0;
    double mendelianLikelihood = 0.0;
    
    List<Genotype> maxGenoType = null;

    // Calculate overall bayes net log likelihood
    for (Genotype genoTypeDad : Genotype.values()) {
      for (Genotype genoTypeMom : Genotype.values()) {
        for (Genotype genoTypeChild : Genotype.values()) {
          double logLikelihood = 0.0;

          /* Get likelihood from the reads */
          logLikelihood += individualLogLikelihood.get(DAD).get(genoTypeDad);
          logLikelihood += individualLogLikelihood.get(MOM).get(genoTypeMom);
          logLikelihood += individualLogLikelihood.get(CHILD).get(genoTypeChild);

          /* Get likelihoods from the trio relationship */
          logLikelihood += getLogLikelihoodFromCPT(DAD, genoTypeDad);
          logLikelihood += getLogLikelihoodFromCPT(MOM, genoTypeMom);
          logLikelihood += getLogLikelihoodFromCPT(CHILD, genoTypeDad, genoTypeMom, genoTypeChild);

          if (DenovoUtil.checkTrioGenoTypeIsDenovo(genoTypeDad, genoTypeMom, genoTypeChild)) {
            denovoLikelihood += Math.exp(logLikelihood);
          } else {
            mendelianLikelihood += Math.exp(logLikelihood);
          }
          
          if (logLikelihood > maxLogLikelihood) {
            maxLogLikelihood = logLikelihood;
            maxGenoType = Arrays.asList(genoTypeDad, genoTypeMom, genoTypeChild);
          }
        }
      }
    }
    
    double bayesDenovoProb = denovoLikelihood  / (denovoLikelihood + mendelianLikelihood);
    
    // ln(likelihood null/likelihood alternate)
    double logOfLikelihoodRatio = mendelianLikelihood - denovoLikelihood;
    
    return new DenovoBayesNet.InferenceResult(maxGenoType, bayesDenovoProb, logOfLikelihoodRatio);
  }
}
