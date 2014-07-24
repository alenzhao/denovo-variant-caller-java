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
import com.google.cloud.genomics.denovo.ReadSummary;
import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.collect.Iterables;

import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.CHILD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.DAD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.MOM;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;


/*
 * Performs Bayesian Inference over reads Contains logic for creating the bayes net, calculating
 * likelihoods, performing Maximum A (MAP) inference and checking whether MAP candidate is indeed a
 * denovo variant
 */
public class BayesInfer {
  private static DenovoBayesNet bn;
  private static boolean isInitialized = false;

  private BayesInfer() {}

  // static initialization for bayes net object
  public static void Init(CommandLine cmdLine) {

    // Check if the static class has already been initialized
    if (isInitialized) {
      return;
    }

    // Create a new Bayes net and fill in the params
    bn = new DenovoBayesNet(cmdLine.sequenceErrorRate, cmdLine.denovoMutationRate);
    bn.addNode(new Node<>(DAD, null, bn.createConditionalProbabilityTable(DAD)));
    bn.addNode(new Node<>(MOM, null, bn.createConditionalProbabilityTable(MOM)));
    List<Node<TrioIndividual, Genotypes>> childParents = new ArrayList<>();
    childParents.add(bn.nodeMap.get(DAD));
    childParents.add(bn.nodeMap.get(MOM));
    bn.addNode(new Node<>(CHILD, childParents, bn.createConditionalProbabilityTable(CHILD)));

    // Set the initialization flag
    isInitialized = true;
  }


  /*
   * Performs inference given a set of mom, dad and child reads to determine the most likely
   * genotype for the trio
   */
  public static boolean infer(Map<TrioIndividual, ReadSummary> readSummaryMap,
      CommandLine cmdLine) {

    // Initialize the bayes net if not already done
    if (!isInitialized) {
      Init(cmdLine);
    }
    // Calculate Likelihoods of the different reads
    Map<TrioIndividual, Map<Genotypes, Double>> individualLogLikelihood =
        bn.getIndividualLogLikelihood(readSummaryMap);

    // Get the trio genotype with the max likelihood
    List<Genotypes> maxTrioGenoType = bn.getMaxGenoType(individualLogLikelihood);

    // Check that the MAP genotype has indeed the highest likelihood
    boolean checkTrioGenoTypeIsDenovo = DenovoUtil.checkTrioGenoTypeIsDenovo(maxTrioGenoType);

    // Convert to Tree Map in order to order the keys
    TreeMap<TrioIndividual, ReadSummary> treeReadSummaryMap = new TreeMap<>();
    treeReadSummaryMap.putAll(readSummaryMap);

    String readCounts = Joiner.on(";").join(Iterables.transform(treeReadSummaryMap.entrySet(),
        new Function<Entry<TrioIndividual, ReadSummary>, String>() {
          @Override
          public String apply(Entry<TrioIndividual, ReadSummary> e) {
            return String.format("%s:%s", e.getKey().name(), e.getValue().getCount().toString());
          }
        }));

    System.out.format("readCounts=%s,maxGenoType=%s,isDenovo=%b%n", readCounts,
        maxTrioGenoType.toString(), checkTrioGenoTypeIsDenovo);

    return checkTrioGenoTypeIsDenovo;
  }
}
