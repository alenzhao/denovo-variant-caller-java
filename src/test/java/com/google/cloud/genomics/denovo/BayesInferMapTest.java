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

import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.AA;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.AG;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.CC;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.CT;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.GG;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.TT;
import static com.google.cloud.genomics.denovo.DenovoUtil.InferenceMethod.MAP;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.cloud.genomics.denovo.DenovoUtil.TrioMember;

import org.junit.Test;

import java.util.Arrays;
import java.util.Map;

public class BayesInferMapTest extends BayesInferTest {

  @Test 
  public void testAllBasesSame() {
    ReadSummary summary = createSameReadSummary();
    Map<TrioMember, ReadSummary> summaryMap = createMapReadSummary(summary, summary, summary);
    
    BayesInfer.BayesCallResult result = bayesInferrer.infer(summaryMap, MAP);
    assertFalse(result.isDenovo());
    assertEquals(summaryMap.toString()+" => [AA,AA,AA]", 
        Arrays.asList(AA,AA,AA), result.getMaxTrioGenoType());
  }
  

  @Test
  public void testAlmostAllBasesSame() {
    ReadSummary summary = createAlmostSameReadSummary();
    Map<TrioMember, ReadSummary> summaryMap = createMapReadSummary(summary, summary, summary);
    
    BayesInfer.BayesCallResult result = bayesInferrer.infer(summaryMap, MAP);
    assertFalse(result.isDenovo());
    assertEquals(summaryMap.toString()+" => [AA,AA,AA]", 
        Arrays.asList(AA,AA,AA), result.getMaxTrioGenoType());
  }

  @Test
  /** 
   * A very interesting edge case 
   */
  public void testChrXPos154226820() {
    Map<TrioMember, ReadSummary> readSummaryMap = createReadSummaryMapChrXPos154226820();
    BayesInfer.BayesCallResult result = bayesInferrer.infer(readSummaryMap, MAP);
    assertFalse(result.isDenovo());
    assertEquals(readSummaryMap.toString()+" => [CT,TT,CT]", 
        Arrays.asList(CT,TT,CT), result.getMaxTrioGenoType());
  }

  @Test
  public void testTrioPos816785() {
    Map<TrioMember, ReadSummary> readSummaryMap = createReadSummaryMapChr1Pos816785();
    BayesInfer.BayesCallResult result = bayesInferrer.infer(readSummaryMap, MAP);
    
    assertFalse(result.isDenovo());
    assertEquals("816785 => [CC,CC,CC]", Arrays.asList(CC,CC,CC), result.getMaxTrioGenoType());
  }
  
  @Test
  public void testTrioPos846600(){
    Map<TrioMember, ReadSummary> readSummaryMap =
        createReadSummaryMapChr1Pos846600L();        
    BayesInfer.BayesCallResult result = bayesInferrer.infer(readSummaryMap, MAP);
    
    assertFalse(result.isDenovo());
    assertEquals("846600 => [CC,CC,CC]", Arrays.asList(CC,CC,CC), result.getMaxTrioGenoType());
  }

  @Test
  public void testTrioPos149035163(){
    Map<TrioMember, ReadSummary> readSummaryMap =
        createReadSummaryMapChr1Pos149035163L();
    BayesInfer.BayesCallResult result = bayesInferrer.infer(readSummaryMap, MAP);

    assertEquals("149035163 => [CC,CC,CC]", Arrays.asList(CC, CC, CC), result.getMaxTrioGenoType());
    assertFalse(result.isDenovo());
  }
  
  /* Begin testing positive gold stanfdards */
  @Test
  public void testTrioPosChr1pos75884343() {
    Map<TrioMember, ReadSummary> readSummaryMap =
        createReadSummaryMapChr1Pos75884343();
    BayesInfer.BayesCallResult result = bayesInferrer.infer(readSummaryMap, MAP);

    assertEquals("75884343 => [TT,TT,CT]", Arrays.asList(TT, TT, CT), result.getMaxTrioGenoType());
    assertTrue(result.isDenovo());
  }
  
  @Test
  public void testTrioPosChr1pos110583335() {
    Map<TrioMember, ReadSummary> readSummaryMap =
        createReadSummaryMapChr1Pos110583335();
    BayesInfer.BayesCallResult result = bayesInferrer.infer(readSummaryMap, MAP);

    assertEquals("110583335 => [GG,GG,AG]", Arrays.asList(GG, GG, AG), result.getMaxTrioGenoType());
    assertTrue(result.isDenovo());
  }
}
