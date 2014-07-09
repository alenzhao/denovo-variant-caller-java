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

import com.google.cloud.genomics.denovo.ReadSummary;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/*
 * Performs Bayesian Inference over reads Contains logic for creating the bayes net, calculating
 * likelihoods, performing Maximum A (MAP) inference and checking whether MAP candidate is indeed a
 * denovo variant
 */
public class BayesInfer {
  private static final double SEQ_ERR_RATE = 1e-8;
  private static final double DENOVO_MUT_RATE = 1e-2; // 1% error rate
  private static final String[] GENOTYPES =
      {"AA", "AC", "AT", "AG", "CC", "CT", "CG", "TT", "TG", "GG"};
  private static BayesNet bn;

  private BayesInfer() {}

  // static initialization for bayes net object
  private static void Init() {
    if (bn != null) {
      return;
    }
    bn = new BayesNet();
    bn.addNode(new Node("DAD", null, createConditionalProbabilityTable("DAD")));
    bn.addNode(new Node("MOM", null, createConditionalProbabilityTable("MOM")));
    List<Node> childParents = new ArrayList<Node>();
    childParents.add(bn.nodeMap.get("DAD"));
    childParents.add(bn.nodeMap.get("MOM"));
    bn.addNode(new Node("CHILD", childParents, createConditionalProbabilityTable("CHILD")));

  }

  /*
   * Creates the conditional probability table to be used in the bayes net One each for mom, dad and
   * child
   */
  /**
   * @param id
   * @return conditionalProbabilityTable
   */
  private static Map<String, Double> createConditionalProbabilityTable(String id) {

    Map<String, Double> conditionalProbabilityTable = new HashMap<String, Double>();
    switch (id) {
      case "DAD": // falls through
      case "MOM":
        for (String genoType : GENOTYPES) {
          conditionalProbabilityTable.put(genoType, 1.0 / GENOTYPES.length);
        }
        break;
      case "CHILD":
        // Loops over parent genotypes
        for (String genoTypeDad : GENOTYPES) {
          for (String genoTypeMom : GENOTYPES) {

            int validInheritanceCases = 0;
            // Initial pass
            for (String genoTypeChild : GENOTYPES) {
              double prob;
              String c1 = genoTypeChild.substring(0, 1);
              String c2 = genoTypeChild.substring(1, 2);
              boolean predicate1 = genoTypeMom.contains(c1) & genoTypeDad.contains(c2);
              boolean predicate2 = genoTypeMom.contains(c2) & genoTypeDad.contains(c1);
              boolean predicate3 = (predicate1 | predicate2);
              if (predicate3) {
                prob = 1.0;
                validInheritanceCases++;
              } else {
                prob = 0.0;
              }
              String cptKey = genoTypeDad + "|" + genoTypeMom + "|" + genoTypeChild;
              conditionalProbabilityTable.put(cptKey, prob);
            }
            // Secondary Pass to normalize prob values
            for (String genoTypeChild : GENOTYPES) {
              String cptKey = genoTypeDad + "|" + genoTypeMom + "|" + genoTypeChild;
              if (conditionalProbabilityTable.get(cptKey) == 0.0) {
                conditionalProbabilityTable.put(cptKey, DENOVO_MUT_RATE);
              } else {
                conditionalProbabilityTable.put(cptKey, 1.0 / validInheritanceCases
                    - DENOVO_MUT_RATE * (GENOTYPES.length - validInheritanceCases)
                    / (validInheritanceCases));
              }
            }
            // Sanity check - probabilities should add up to 1.0 (almost)
            double totProb = 0.0;
            for (String genoTypeChild : GENOTYPES) {
              String cptKey = genoTypeDad + "|" + genoTypeMom + "|" + genoTypeChild;
              totProb += conditionalProbabilityTable.get(cptKey);
            }
            if (Math.abs(totProb - 1.0) > 1e-12) {
              throw new RuntimeException(
                  "cpt probabilities not adding up : " + String.valueOf(totProb));
            }
          }
        }
        break;
      default:
        throw new RuntimeException("Unknown CPT id : " + id);

    }
    return conditionalProbabilityTable;
  }

  /*
   * Performs inference given a set of mom, dad and child reads to determine the most likely
   * genotype for the trio
   */
  public static boolean infer(Map<String, ReadSummary> readSummaryMap) {

    // Initialize the bayes net if not already done
    if (bn == null) {
      Init();
    }
    // Calculate Likelihoods of the different reads
    Map<String, Map<String, Double>> individualLogLikelihood =
        getIndividualLogLikelihood(readSummaryMap);

    // Get the trio genotype with the max likelihood
    String maxTrioGenoType = getMaxGenoType(individualLogLikelihood);

    // Check that the MAP genotype has indeed the highest likelihood
    boolean checkTrioGenoTypeIsDenovo = checkTrioGenoTypeIsDenovo(maxTrioGenoType);

    System.out.print(",");
    for (String p : readSummaryMap.keySet()) {
      System.out.print(p + ":" + readSummaryMap.get(p).getCount().toString() + ";");
    }
    System.out.print(",maxGenoType=" + maxTrioGenoType);
    System.out.print(",isDenovo=" + checkTrioGenoTypeIsDenovo);
    System.out.println();

    return checkTrioGenoTypeIsDenovo;
  }

  /**
   * Get the log likelihood of reads for a particular individual in a trio for all their possible
   * genotypes
   *
   * @param readSummaryMap
   * @return individualLogLikelihood
   */
  private static Map<String, Map<String, Double>> getIndividualLogLikelihood(
      Map<String, ReadSummary> readSummaryMap) {
    Map<String, Map<String, Double>> individualLogLikelihood = new HashMap<>();
    for (String trioIndividual : Arrays.asList("DAD", "MOM", "CHILD")) {

      ReadSummary readSummary = readSummaryMap.get(trioIndividual);
      Map<String, Double> genoTypeLogLikelihood = getGenoTypeLogLikelihood(readSummary);
      individualLogLikelihood.put(trioIndividual, genoTypeLogLikelihood);
    }
    return individualLogLikelihood;
  }

  /**
   * Get the log likelihood for all possible genotypes for a set of reads
   *
   * @param readSummary
   * @return genotypeLogLikelihood
   */
  private static Map<String, Double> getGenoTypeLogLikelihood(ReadSummary readSummary) {
    Map<String, Double> genotypeLogLikelihood = new HashMap<>();
    for (String genoType : GENOTYPES) {
      Map<String, Integer> count = readSummary.getCount();
      boolean isHomozygous = genoType.substring(0, 1).equals(genoType.substring(1, 2));

      double readlogLikelihood = 0.0;
      for (String base : count.keySet()) {
        double ll = getBaseLikelihood(genoType, isHomozygous, base);
        readlogLikelihood += ll;
      }
      genotypeLogLikelihood.put(genoType, readlogLikelihood);
    }
    return genotypeLogLikelihood;
  }

  /**
   * Check if the particular genotype is denovo i.e. present in kids but not in parents
   *
   * @param maxGenoType
   * @return isDenovo
   */
  private static boolean checkTrioGenoTypeIsDenovo(String maxGenoType) {
    String[] maxGenoTypeSplit = maxGenoType.split("\\|");
    String genoTypeDad = maxGenoTypeSplit[0];
    String genoTypeMom = maxGenoTypeSplit[1];
    String genoTypeChild = maxGenoTypeSplit[2];
    String c1 = genoTypeChild.substring(0, 1);
    String c2 = genoTypeChild.substring(1, 2);
    boolean predicate1 = genoTypeMom.contains(c1) & genoTypeDad.contains(c2);
    boolean predicate2 = genoTypeMom.contains(c2) & genoTypeDad.contains(c1);
    boolean predicate3 = !(predicate1 | predicate2);
    return predicate3;
  }

  /**
   * Infer the most likely genotype given likelihood for all the genotypes of the trio
   *
   * @param individualLogLikelihood
   * @return maxgenoType
   */
  private static String getMaxGenoType(Map<String, Map<String, Double>> individualLogLikelihood) {
    double maxLogLikelihood = Double.NEGATIVE_INFINITY;
    String maxGenoType = "AA|AA|AA";
    // Calculate overall bayes net log likelihood
    for (String genoTypeDad : GENOTYPES) {
      for (String genoTypeMom : GENOTYPES) {
        for (String genoTypeChild : GENOTYPES) {
          double ll = individualLogLikelihood.get("DAD").get(genoTypeDad)
              + individualLogLikelihood.get("MOM").get(genoTypeMom)
              + individualLogLikelihood.get("CHILD").get(genoTypeChild)
              + bn.nodeMap.get("DAD").conditionalProbabilityTable.get(genoTypeDad)
              + bn.nodeMap.get("MOM").conditionalProbabilityTable.get(genoTypeMom) + bn.nodeMap
                  .get("CHILD").conditionalProbabilityTable.get(
                  genoTypeDad + "|" + genoTypeMom + "|" + genoTypeChild);

          if (ll > maxLogLikelihood) {
            maxLogLikelihood = ll;
            maxGenoType = genoTypeDad + "|" + genoTypeMom + "|" + genoTypeChild;
          }
        }
      }
    }
    return maxGenoType;
  }

  /**
   * Get the log likelihood for a particular read base
   *
   * @param genoType
   * @param isHomozygous
   * @param base
   * @return logLikeliHood
   */
  private static double getBaseLikelihood(String genoType, boolean isHomozygous, String base) {
    double logLikeliHood = 0.0;
    if (isHomozygous) {
      if (genoType.contains(base)) {
        logLikeliHood = Math.log(1 - SEQ_ERR_RATE);
      } else {
        logLikeliHood = Math.log(SEQ_ERR_RATE) - Math.log(3);
      }
    } else {
      if (genoType.contains(base)) {
        logLikeliHood = Math.log(1 - 2 * SEQ_ERR_RATE / 3) - Math.log(2);
      } else {
        logLikeliHood = Math.log(SEQ_ERR_RATE) - Math.log(3);
      }
    }
    return logLikeliHood;
  }
}


/*
 * Bayes net data structure
 */
class BayesNet {
  public Map<String, Node> nodeMap;

  public BayesNet() {
    nodeMap = new HashMap<String, Node>();
  }

  public void addNode(Node node) {
    nodeMap.put(node.id, node);
  }
}


/*
 * Individual Node in the Bayes Net
 */
class Node {
  public String id;
  public List<Node> parents;
  public Map<String, Double> conditionalProbabilityTable;

  public Node(String id, List<Node> parents, Map<String, Double> cpt) {
    this.id = id;
    this.parents = parents;
    this.conditionalProbabilityTable = cpt;
  }
}
