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

import com.google.cloud.genomics.denovo.DenovoUtil.Genotype;
import com.google.cloud.genomics.denovo.DenovoUtil.InferenceMethod;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioMember;
import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.collect.Iterables;

import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

/**
 * Performs Bayesian Inference over reads. Contains logic for creating the bayes net, calculating
 * likelihoods, performing Maximum A (MAP) inference and checking whether MAP candidate is indeed a
 * denovo variant
 */
class BayesInfer {
  private final DenovoBayesNet dbn;
  private final DenovoShared shared;
  
  /**
   * @param shared the shared parameters for the tool
   */
  BayesInfer(DenovoShared shared) {
    dbn = new DenovoBayesNet(shared);
    this.shared = shared;
  }
  
  /**
   * Delegates the inference task to the MAP, Bayes Rule or LRT handler methods
   * @param readSummaryMap the summary statistics of the reads at a position
   * @param inferMethod the chosen bayesian inference method
   * @return result of the inference procedure for that position
   */
  BayesCallResult infer(Map<TrioMember, ReadSummary> readSummaryMap,
      InferenceMethod inferMethod) {

    // Get the trio genotype with the max likelihood
    BayesInferenceResult result = dbn.performInference(readSummaryMap);
    boolean isDenovo = inferMethod.isDenovo(result, shared);

    String readCounts = Joiner.on(";").join(Iterables.transform(readSummaryMap.entrySet(),
        new Function<Entry<TrioMember, ReadSummary>, String>() {
          @Override
          public String apply(Entry<TrioMember, ReadSummary> e) {
            return String.format("%s:%s", e.getKey().name(), e.getValue().getCount());
          }
        }));

    return new BayesCallResult(isDenovo, result.getMaxTrioGenotype(), 
        String.format("readCounts=%s,maxGenoType=%s,isDenovo=%b", 
            readCounts,
            result.getMaxTrioGenotype(), isDenovo));
  }

  /**
   * This container holds the call result from the Bayesian Inference procedure  
   */
  public static class BayesCallResult {
    private final boolean denovo;
    private final String details;
    private final List<Genotype> maxTrioGenoType;

    /**
     * @param isDenovo if call is denovo
     * @param maxTrioGenoType the genotype of the trio with max likelihood
     * @param details additional details
     */
    public BayesCallResult(boolean isDenovo, List<Genotype> maxTrioGenoType, String details) {
      this.denovo = isDenovo;
      this.details = details;
      this.maxTrioGenoType = maxTrioGenoType;
    }

    @Override
    public String toString() {
      return "<" + details + ">";
    }
    
    public boolean isDenovo() {
      return denovo;
    }

    public String getDetails() {
      return details;
    }

    public List<Genotype> getMaxTrioGenoType() {
      return maxTrioGenoType;
    }
  }
}
