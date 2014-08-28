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

import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.AG;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.CC;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.CT;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.GG;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.TT;
import static com.google.cloud.genomics.denovo.DenovoUtil.InferenceMethod.BAYES;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.cloud.genomics.denovo.DenovoUtil.TrioMember;

import org.junit.Test;

import java.util.Arrays;
import java.util.Map;

public class BayesInferBayesTest extends BayesInferTest {

  @Test
  public void testTrioChr1Pos1298169() {
    /*{MOM={T=30}, CHILD={A=6, T=39}, DAD={A=2, T=41}}*/
    Map<TrioMember, ReadSummary> readSummaryMap = createReadSummaryMapChr1Pos1298169();
    BayesInfer.BayesCallResult result = bayesInferrer.infer(readSummaryMap, BAYES);

    assertFalse(result.isDenovo());
    assertEquals("1298169 => [TT,TT,TT]", Arrays.asList(TT, TT, TT), result.getMaxTrioGenoType());
  }
  
  @Test
  public void testTrioPosChr170041751(){
    /*{DAD={T=2, C=58}, CHILD={T=8, C=28}, MOM={T=2, C=51}}*/
    Map<TrioMember, ReadSummary> readSummaryMap = createReadSummaryMapChr170041751();
    BayesInfer.BayesCallResult result = bayesInferrer.infer(readSummaryMap, BAYES);

    assertEquals("70041751 => [CC,CC,CC]", Arrays.asList(CC, CC, CT), result.getMaxTrioGenoType());
    assertTrue(result.isDenovo());
  }

  @Test
  public void testTrioChr1Pos149035163() {
    
    Map<TrioMember, ReadSummary> readSummaryMap = createReadSummaryMapChr1Pos149035163();
    BayesInfer.BayesCallResult result = bayesInferrer.infer(readSummaryMap, BAYES);

    assertEquals("149035163 => [CC,CC,CC]", Arrays.asList(CC, CC, CC), result.getMaxTrioGenoType());
    assertFalse(result.isDenovo());
  }
  
  /* Begin testing positive gold stanfdards */
  @Test
  public void testTrioPosChr1pos75884343() {
    Map<TrioMember, ReadSummary> readSummaryMap =
        createReadSummaryMapChr1Pos75884343();
    BayesInfer.BayesCallResult result = bayesInferrer.infer(readSummaryMap, BAYES);

    assertEquals("75884343 => [TT,TT,CT]", Arrays.asList(TT, TT, CT), result.getMaxTrioGenoType());
    assertTrue(result.isDenovo());
  }
  
  @Test
  public void testTrioPosChr1pos110583335() {
    Map<TrioMember, ReadSummary> readSummaryMap =
        createReadSummaryMapChr1Pos110583335();
    BayesInfer.BayesCallResult result = bayesInferrer.infer(readSummaryMap, BAYES);

    assertEquals("110583335 => [GG,GG,AG]", Arrays.asList(GG, GG, AG), result.getMaxTrioGenoType());
    assertTrue(result.isDenovo());
  }
}
