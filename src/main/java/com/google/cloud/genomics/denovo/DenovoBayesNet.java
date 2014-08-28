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

import static com.google.cloud.genomics.denovo.DenovoUtil.TrioMember.CHILD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioMember.DAD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioMember.MOM;

import com.google.cloud.genomics.denovo.DenovoUtil.Allele;
import com.google.cloud.genomics.denovo.DenovoUtil.Genotype;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioMember;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Bayes Net with the following structure
 * 
 *    Dad       Mom
 *     |\       /|
 *     | \     / |
 *     |  \   /  |
 *     |  Child  |
 *     |    |    |
 *     Rd   Rc   Rm
 *     
 *  P(D,M,C,R) = P(Rd)P(Rd|D)P(Rm|M)P(Rc|C)P(C|D,M) 
 */ 
public class DenovoBayesNet {

    private final DenovoShared shared;
    private final Map<TrioMember, Node<TrioMember, Genotype>> nodeMap; 
    
  /**
   * @param shared Shared parameters for tool
   */
  public DenovoBayesNet(DenovoShared shared) {
    this.shared = shared;
    nodeMap = new HashMap<TrioMember, Node<TrioMember, Genotype>>();

    // Initialize the conditional Probability table
    addNode(new Node<>(DAD, null, createConditionalProbabilityTable(DAD)));
    addNode(new Node<>(MOM, null, createConditionalProbabilityTable(MOM)));
    List<Node<TrioMember, Genotype>> childParents = new ArrayList<>();
    childParents.add(nodeMap.get(DAD));
    childParents.add(nodeMap.get(MOM));
    addNode(new Node<>(CHILD, childParents, createConditionalProbabilityTable(CHILD)));
  }

  /**
   * add a new node to bayes net
   * @param node node to be added 
   */
  void addNode(Node<TrioMember, Genotype> node) {
    getNodeMap().put(node.getId(), node);
  }

  /**
   * Creates the conditional probability table to be used in the bayes net One each for mom, dad and
   * child
   * @param individual : Mom, dad or child
   * @return Conditional Probability Table
   */
  Map<List<Genotype>, Double> createConditionalProbabilityTable(TrioMember individual) {

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

  /**
   * Get the set of genotypes that are mendelian for a certain mom and dad genotype pair
   * @param genoTypeDad Dad genotype
   * @param genoTypeMom Mom genotype
   * @return a Map containing mendelian Genotypes and their counts by mendelian inheritance
   */
  Map<Genotype, Integer> mendelianGenotypes(Genotype genoTypeDad, Genotype genoTypeMom) {
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
   * @return logLikeliHood likelihood of genotype generating the base
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
   * @param readSummaryMap Summary stats for reads for all trio members
   * @return individualLogLikelihood likelihood of reads for a member in the trio
   */
  Map<TrioMember, Map<Genotype, Double>> getIndividualLogLikelihood(
      Map<TrioMember, ReadSummary> readSummaryMap) {
    Map<TrioMember, Map<Genotype, Double>> individualLogLikelihood = new HashMap<>();
    for (TrioMember person : TrioMember.values()) {

      ReadSummary readSummary = readSummaryMap.get(person);
      Map<Genotype, Double> genoTypeLogLikelihood = getReadSummaryLogLikelihood(readSummary);
      individualLogLikelihood.put(person, genoTypeLogLikelihood);
    }
    return individualLogLikelihood;
  }

  /**
   * Get the log likelihood for all possible Genotypes for a set of reads
   *
   * @param readSummary Summary stats for reads for all trio members
   * @return genotypeLogLikelihood likelihood of reads for a particular genotype
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

  /**
   * Extract likelihood from conditional probability table
   * 
   * @param person a member of the trio
   * @param keyParts varargs containing genotype keys
   * @return likelihood value
   */
  double getLogLikelihoodFromCPT(TrioMember person, Genotype... keyParts) {
    if (TrioMember.PARENTS.contains(person) && keyParts.length != 1) {
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

  Map<TrioMember, Node<TrioMember, Genotype>> getNodeMap() {
    return nodeMap;
  }

  /**
   * Performs Bayesian inference by adding iteraing over all possible genotypes in the trio and 
   * calculating the likelihood the of trip by adding the likelihood of the reads as well as the 
   * likelihood of the trio of a genotypes according to laws of inheritance 
   * 
   * @param readSummaryMap
   * @return A BayesInference result object storing likelihoods, genpotypes and other pertinent data
   */
  public BayesInferenceResult performInference(Map<TrioMember, ReadSummary> readSummaryMap) {

    // Calculate Likelihoods of the different reads
    Map<TrioMember, Map<Genotype, Double>> individualLogLikelihood =
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
   * Likelihood of the genotype trio according to probabilities of mendelian genetics
   * 
   * @param individualLogLikelihood likelihood of reads for a member in the trio
   * @param genoTypeDad
   * @param genoTypeMom
   * @param genoTypeChild
   * @return the log likelihood of the trio type
   */
  double getTrioGenotypeLogLikelihood(
      Map<TrioMember, Map<Genotype, Double>> individualLogLikelihood, 
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
   * likelihoods of the trio relationship
   * 
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
