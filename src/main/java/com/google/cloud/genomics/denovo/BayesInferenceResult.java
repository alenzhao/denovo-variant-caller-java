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

import java.util.List;

/**
 * Container for results from the Bayesian Inference on the Bayes Net
 */
public class BayesInferenceResult {
  private final List<Genotype> maxTrioGenotype;
  private final double maxLikelihood;
  private final double bayesDenovoProb;
  private final double likelihoodRatio;
  private final double mendelianLogLikelihood;
  private final double denovoLogLikelihood;

  /**
   * @param maxTrioGenotype the genotype of the trio with max likelihood
   * @param maxLikelihood the max likelihood value
   * @param bayesDenovoProb the bayesian probability of the position being denovo
   * @param logOfLikelihoodRatio the log of the likelihood ratio of denovo over mendelian 
   * @param mendelianLogLikelihood the log likelihood of the mendelian case
   * @param denovoLogLikelihood the log likelihood of the denovo case
   */
  public BayesInferenceResult(List<Genotype> maxTrioGenotype,
      double maxLikelihood,
      double bayesDenovoProb,
      double logOfLikelihoodRatio,
      double mendelianLogLikelihood,
      double denovoLogLikelihood) {
    this.maxTrioGenotype = maxTrioGenotype;
    this.maxLikelihood = maxLikelihood;
    this.bayesDenovoProb = bayesDenovoProb;
    this.likelihoodRatio = logOfLikelihoodRatio;
    this.mendelianLogLikelihood = mendelianLogLikelihood;
    this.denovoLogLikelihood = denovoLogLikelihood;
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("<");
    sb.append("maxGT : " + getMaxTrioGenotype().toString());
    sb.append(", maxLL : " + String.valueOf(getMaxLikelihood()));
    sb.append(", mendelLL : " + String.valueOf(getMendelianLogLikelihood()));
    sb.append(", denovoLL : " + String.valueOf(getDenovoLogLikelihood()));
    sb.append(", bayesProb : " + String.valueOf(getBayesDenovoProb()));
    sb.append(", llRatio : " + String.valueOf(getLikelihoodRatio()));
    sb.append(">");
    return sb.toString();
  }

  /**
   * @return the maxTrioGenotype
   */
  public List<Genotype> getMaxTrioGenotype() {
    return maxTrioGenotype;
  }

  /**
   * @return the maxLikelihood
   */
  public double getMaxLikelihood() {
    return maxLikelihood;
  }

  /**
   * @return the bayesDenovoProb
   */
  public double getBayesDenovoProb() {
    return bayesDenovoProb;
  }

  /**
   * @return the likelihoodRatio
   */
  public double getLikelihoodRatio() {
    return likelihoodRatio;
  }

  /**
   * @return the mendelianLogLikelihood
   */
  public double getMendelianLogLikelihood() {
    return mendelianLogLikelihood;
  }

  /**
   * @return the denovoLogLikelihood
   */
  public double getDenovoLogLikelihood() {
    return denovoLogLikelihood;
  }
}
