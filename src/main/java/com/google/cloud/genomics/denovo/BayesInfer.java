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

import static com.google.cloud.genomics.denovo.DenovoUtil.InferenceMethod.BAYES;
import static com.google.cloud.genomics.denovo.DenovoUtil.InferenceMethod.LRT;
import static com.google.cloud.genomics.denovo.DenovoUtil.InferenceMethod.MAP;

import com.google.cloud.genomics.denovo.DenovoUtil.Genotype;
import com.google.cloud.genomics.denovo.DenovoUtil.InferenceMethod;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;
import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.collect.Iterables;

import org.javatuples.Pair;

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

    // Create a new Denovo BayesNet
    dbn = new DenovoBayesNet(sequenceErrorRate, denovoMutationRate);
  }

  /*
   * Performs inference given a set of mom, dad and child reads to determine the most likely
   * genotype for the trio
   */
  public InferenceResult infer(Map<TrioIndividual, ReadSummary> readSummaryMap,
      InferenceMethod inferMethod) {

    if (inferMethod == MAP) { return performMAPInference(readSummaryMap); }
    if (inferMethod == BAYES) { return performBayesInference(readSummaryMap); }
    if (inferMethod == LRT) { return performLRTInference(readSummaryMap); }
    
    throw new IllegalArgumentException("Unexpected method " + inferMethod);
  }

  private InferenceResult createInferenceResult(Pair<List<Genotype>, Boolean> pair, String readCounts) {
    InferenceResult result = new InferenceResult(pair.getValue1(), pair.getValue0(), 
        String.format("readCounts=%s,maxGenoType=%s,isDenovo=%b", readCounts,
            pair.getValue0(), pair.getValue1()));
    return result;
  }

  private String createReadCountString(Map<TrioIndividual, ReadSummary> readSummaryMap) {
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
    return readCounts;
  }

  private InferenceResult performMAPInference(Map<TrioIndividual, ReadSummary> readSummaryMap) {
    
    // Get the trio genotype with the max likelihood
    DenovoBayesNet.InferenceResult dbnResult = dbn.performInference(readSummaryMap);

    // Check that the MAP genotype has indeed the highest likelihood
    boolean isDenovo = DenovoUtil.checkTrioGenoTypeIsDenovo(dbnResult.maxTrioGenotype);

    return createInferenceResult(Pair.with(dbnResult.maxTrioGenotype, isDenovo), 
        createReadCountString(readSummaryMap));
  }
  
  private InferenceResult performBayesInference(Map<TrioIndividual, ReadSummary> readSummaryMap) {
    
    // Get the trio genotype with the max likelihood
    DenovoBayesNet.InferenceResult dbnResult = dbn.performInference(readSummaryMap);

    // Use the bayesian classifier rule
    boolean isDenovo = dbnResult.bayesDenovoProb > 0.5;

    return createInferenceResult(new Pair<>(dbnResult.maxTrioGenotype, isDenovo), 
        createReadCountString(readSummaryMap));
  }
  
  private InferenceResult performLRTInference(Map<TrioIndividual, ReadSummary> readSummaryMap) {
    
    // Get the trio genotype with the max likelihood
    DenovoBayesNet.InferenceResult dbnResult = dbn.performInference(readSummaryMap);

    // Use the likelihood ratio rule
    boolean isDenovo = dbnResult.likelihoodRatio > DenovoUtil.LRT_THRESHOLD;

    return createInferenceResult(new Pair<>(dbnResult.maxTrioGenotype, isDenovo), 
        createReadCountString(readSummaryMap));
  }

  
  public static class InferenceResult {
    private final boolean isDenovo;
    private final String details;
    private final List<Genotype> maxTrioGenoType;

    public InferenceResult(boolean isDenovo, List<Genotype> maxTrioGenoType, String format) {
      this.isDenovo = isDenovo;
      this.details = format;
      this.maxTrioGenoType = maxTrioGenoType;
    }

    @Override
    public String toString() {
      return "<" + details + ">";
    }
    
    public boolean isDenovo() {
      return isDenovo;
    }

    public String getDetails() {
      return details;
    }

    public List<Genotype> getMaxTrioGenoType() {
      return maxTrioGenoType;
    }
  }
}
