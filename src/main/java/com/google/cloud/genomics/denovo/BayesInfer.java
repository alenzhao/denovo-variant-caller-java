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

import com.google.cloud.genomics.denovo.DenovoUtil.Genotypes;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;
import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.collect.Iterables;

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
  private DenovoBayesNet dbn;
    
  public BayesInfer(Double sequenceErrorRate, Double denovoMutationRate) {

    // Create a new Bayes net and fill in the params
    dbn = new DenovoBayesNet(sequenceErrorRate, denovoMutationRate);
    dbn.addNode(new Node<>(DAD, null, dbn.createConditionalProbabilityTable(DAD)));
    dbn.addNode(new Node<>(MOM, null, dbn.createConditionalProbabilityTable(MOM)));
    List<Node<TrioIndividual, Genotypes>> childParents = new ArrayList<>();
    childParents.add(dbn.getNodeMap().get(DAD));
    childParents.add(dbn.getNodeMap().get(MOM));
    dbn.addNode(new Node<>(CHILD, childParents, dbn.createConditionalProbabilityTable(CHILD)));
  }

  /*
   * Performs inference given a set of mom, dad and child reads to determine the most likely
   * genotype for the trio
   */
  public InferResult infer(Map<TrioIndividual, ReadSummary> readSummaryMap) {

    // Calculate Likelihoods of the different reads
    Map<TrioIndividual, Map<Genotypes, Double>> individualLogLikelihood =
        dbn.getIndividualLogLikelihood(readSummaryMap);

    // Get the trio genotype with the max likelihood
    List<Genotypes> maxTrioGenoType = dbn.getMaxGenoType(individualLogLikelihood);

    // Check that the MAP genotype has indeed the highest likelihood
    boolean checkTrioGenoTypeIsDenovo = DenovoUtil.checkTrioGenoTypeIsDenovo(maxTrioGenoType);

    // Convert to Tree Map in order to order the keys
    TreeMap<TrioIndividual, ReadSummary> treeReadSummaryMap = new TreeMap<>();
    treeReadSummaryMap.putAll(readSummaryMap);

    String readCounts = Joiner.on(";").join(Iterables.transform(treeReadSummaryMap.entrySet(),
        new Function<Entry<TrioIndividual, ReadSummary>, String>() {
          @Override
          public String apply(Entry<TrioIndividual, ReadSummary> e) {
            return String.format("%s:%s", e.getKey().name(), e.getValue().getCount());
          }
        }));

    InferResult result = new InferResult(checkTrioGenoTypeIsDenovo, 
        String.format("readCounts=%s,maxGenoType=%s,isDenovo=%b%n", readCounts,
            maxTrioGenoType, checkTrioGenoTypeIsDenovo));

    return result;
  }
  
  public static class InferResult {
    private final boolean isDenovo;
    private final String details;

    /**
     * @param checkTrioGenoTypeIsDenovo
     * @param format
     */
    public InferResult(boolean checkTrioGenoTypeIsDenovo, String format) {
      this.isDenovo = checkTrioGenoTypeIsDenovo;
      this.details = format;
    }

    public boolean isDenovo() {
      return isDenovo;
    }

    public String getDetails() {
      return details;
    }

  }
}
