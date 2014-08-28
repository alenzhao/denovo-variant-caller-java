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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/*
 * DenovoBayesNet implements abstract BayesNet
 */
public class DenovoBayesNet {

    private final DenovoShared shared;
    private final Map<TrioIndividual, Node<TrioIndividual, Genotype>> nodeMap; 
    
  public DenovoBayesNet(DenovoShared shared) {
    this.shared = shared;
    nodeMap = new HashMap<TrioIndividual, Node<TrioIndividual, Genotype>>();

    // Initialize the conditional Probability table
    addNode(new Node<>(DAD, null, createConditionalProbabilityTable(DAD)));
    addNode(new Node<>(MOM, null, createConditionalProbabilityTable(MOM)));
    List<Node<TrioIndividual, Genotype>> childParents = new ArrayList<>();
    childParents.add(nodeMap.get(DAD));
    childParents.add(nodeMap.get(MOM));
    addNode(new Node<>(CHILD, childParents, createConditionalProbabilityTable(CHILD)));
  }

  void addNode(Node<TrioIndividual, Genotype> node) {
    getNodeMap().put(node.getId(), node);
  }

  /*
   * Creates the conditional probability table to be used in the bayes net One each for mom, dad and
   * child
   */
  Map<List<Genotype>, Double> createConditionalProbabilityTable(TrioIndividual individual) {

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
                ? shared.getDenovoMutationRate() / numDenovoGenotypes
                : (1.0 - shared.getDenovoMutationRate()) / mendelianCount 
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
            Genotype.valueOf(new String(mendelianCharAlleles));
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
   * @param genotype
   * @param base
   * @return logLikeliHood
   */
  double getBaseLogLikelihood(Genotype genotype, Allele base) {
    if (base == null || genotype == null) {
      throw new NullPointerException("Can't get base log likelihood of null members");
    }

    return genotype.isHomozygous()
        ? genotype.containsAllele(base)
            ? Math.log(1 - shared.getSequenceErrorRate())
            : Math.log(shared.getSequenceErrorRate()) - Math.log(3)
        : genotype.containsAllele(base)
            ? Math.log(1 - 2 * shared.getSequenceErrorRate() / 3) - Math.log(2)
            : Math.log(shared.getSequenceErrorRate()) - Math.log(3);
  }

  /**
   * Get the log likelihood of reads for a particular individual in a trio for all their possible
   * Genotypes
   *
   * @param readSummaryMap
   * @return individualLogLikelihood
   */
  Map<TrioIndividual, Map<Genotype, Double>> getIndividualLogLikelihood(
      Map<TrioIndividual, ReadSummary> readSummaryMap) {
    Map<TrioIndividual, Map<Genotype, Double>> individualLogLikelihood = new HashMap<>();
    for (TrioIndividual person : TrioIndividual.values()) {

      ReadSummary readSummary = readSummaryMap.get(person);
      Map<Genotype, Double> genoTypeLogLikelihood = getReadSummaryLogLikelihood(readSummary);
      individualLogLikelihood.put(person, genoTypeLogLikelihood);
    }
    return individualLogLikelihood;
  }

  /**
   * Get the log likelihood for all possible Genotypes for a set of reads
   *
   * @param readSummary
   * @return genotypeLogLikelihood
   */
  Map<Genotype, Double> getReadSummaryLogLikelihood(ReadSummary readSummary) {
    
    if (readSummary == null) {
      throw new NullPointerException("Did not expect ReadSummary to be null");
    }
    
    Map<Genotype, Double> genotypeLogLikelihood = new HashMap<>();
    for (Genotype genoType : Genotype.values()) {
      Map<Allele, Integer> count = readSummary.getCount();

      double readlogLikelihood = 0.0;
      for (Map.Entry<Allele, Integer> entry : count.entrySet()) {
        readlogLikelihood += entry.getValue() * 
            getBaseLogLikelihood(genoType, entry.getKey());
      }
      genotypeLogLikelihood.put(genoType, readlogLikelihood);
    }
    return genotypeLogLikelihood;
  }

  double getLogLikelihoodFromCPT(TrioIndividual person, Genotype... keyParts) {
    if (TrioIndividual.PARENTS.contains(person) && keyParts.length != 1) {
      throw new IllegalArgumentException("Expected one Genotype argument : got " + keyParts.length);
    }
    if (person == CHILD && keyParts.length != 3) {
      throw new IllegalArgumentException("Expected three Genotype argument : got " + keyParts.length);
    }
    
    List<Genotype> cptKey = Arrays.asList(keyParts);
    return Math.log(getNodeMap().
        get(person).
        getConditionalProbabilityTable().
        get(cptKey));
  }

  /**
   * @return the nodeMap
   */
  Map<TrioIndividual, Node<TrioIndividual, Genotype>> getNodeMap() {
    return nodeMap;
  }

  public BayesInferenceResult performInference(Map<TrioIndividual, ReadSummary> readSummaryMap) {

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
          
          double logLikelihood = 0;
          logLikelihood += getTrioGenotypeLogLikelihood(individualLogLikelihood, genoTypeDad,
              genoTypeMom, genoTypeChild);
          logLikelihood += getRelationshipLogLikelihood(genoTypeDad, genoTypeMom, genoTypeChild);
          
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
    double likelihoodRatio = denovoLikelihood / mendelianLikelihood;
    
    return new BayesInferenceResult(maxGenoType, maxLogLikelihood, 
        bayesDenovoProb, likelihoodRatio, Math.log(mendelianLikelihood), Math.log(denovoLikelihood));
  }

  /**
   * @param individualLogLikelihood
   * @param genoTypeDad
   * @param genoTypeMom
   * @param genoTypeChild
   * @return the log likelihood of the trio type
   */
  double getTrioGenotypeLogLikelihood(
      Map<TrioIndividual, Map<Genotype, Double>> individualLogLikelihood, 
      Genotype genoTypeDad, 
      Genotype genoTypeMom, 
      Genotype genoTypeChild) {
    double logLikelihood = 0.0;

    /* Get likelihood from the reads */
    logLikelihood += individualLogLikelihood.get(DAD).get(genoTypeDad);
    logLikelihood += individualLogLikelihood.get(MOM).get(genoTypeMom);
    logLikelihood += individualLogLikelihood.get(CHILD).get(genoTypeChild);
    return logLikelihood;
  }

  /**
   * @param genoTypeDad
   * @param genoTypeMom
   * @param genoTypeChild
   * @return likelihoods of the trio relationship
   */
  double getRelationshipLogLikelihood(Genotype genoTypeDad, 
      Genotype genoTypeMom, 
      Genotype genoTypeChild) {
    double logLikelihood = 0;
    logLikelihood += getLogLikelihoodFromCPT(DAD, genoTypeDad);
    logLikelihood += getLogLikelihoodFromCPT(MOM, genoTypeMom);
    logLikelihood += getLogLikelihoodFromCPT(CHILD, genoTypeDad, genoTypeMom, genoTypeChild);
    return logLikelihood;
  }
}
