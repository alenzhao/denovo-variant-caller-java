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
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.CT;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.TT;
import static com.google.cloud.genomics.denovo.DenovoUtil.InferenceMethod.MAP;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.CHILD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.DAD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.MOM;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.cloud.genomics.denovo.DenovoUtil.Allele;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;

import org.junit.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class BayesInferMapTest extends BayesInferTest {

  @Test 
  public void testAllBasesSame() {
    ReadSummary summary = createSameReadSummary();
    Map<TrioIndividual, ReadSummary> summaryMap = createMapReadSummary(summary, summary, summary);
    
    BayesInfer.InferenceResult result = bayesInferrer.infer(summaryMap, MAP);
    assertFalse(result.isDenovo());
    assertEquals(summaryMap.toString()+" => [AA,AA,AA]", 
        Arrays.asList(AA,AA,AA), result.getMaxTrioGenoType());
  }
  

  @Test
  public void testAlmostAllBasesSame() {
    ReadSummary summary = createAlmostSameReadSummary();
    Map<TrioIndividual, ReadSummary> summaryMap = createMapReadSummary(summary, summary, summary);
    
    BayesInfer.InferenceResult result = bayesInferrer.infer(summaryMap, MAP);
    assertFalse(result.isDenovo());
    assertEquals(summaryMap.toString()+" => [AA,AA,AA]", 
        Arrays.asList(AA,AA,AA), result.getMaxTrioGenoType());
  }

  @Test
  public void testChrXPos154226820() {
    Map<TrioIndividual, ReadSummary> readSummaryMap = new HashMap<>();
    for (TrioIndividual person : TrioIndividual.values()) {
      Map<Allele, Integer> baseCount = new HashMap<>();
      if (person == DAD ) baseCount.put(Allele.T,28);
      if (person == MOM ) baseCount.put(Allele.T,36);
      if (person == CHILD ) {
        baseCount.put(Allele.T,33);
        baseCount.put(Allele.C,15);
      }
      readSummaryMap.put(person, new ReadSummary().setCount(baseCount));
    }
    
    BayesInfer.InferenceResult result = bayesInferrer.infer(readSummaryMap, MAP);
    assertTrue(result.isDenovo());
    assertEquals(readSummaryMap.toString()+" => [TT,TT,CT]", 
        Arrays.asList(TT,TT,CT), result.getMaxTrioGenoType());
  }

  
//  @Test
//  public void testTrioPos816785MAP() throws IOException {
//    Map<TrioIndividual, ReadSummary> readSummaryMap =
//        expRunner.getReadSummaryMap(816785L, expRunner.getReadMap("chr1", 816785L));
//    BayesInfer.InferenceResult result = bayesInferrer.infer(readSummaryMap, MAP);
//    
//    assertFalse(result.isDenovo());
//    assertEquals("816785 => [CC,CC,CC]", Arrays.asList(CC,CC,CC), result.getMaxTrioGenoType());
//  }
//  
//  @Test
//  public void testTrioPos846600MAP() throws IOException {
//    Map<TrioIndividual, ReadSummary> readSummaryMap =
//        expRunner.getReadSummaryMap(846600L, expRunner.getReadMap("chr1", 846600L));
//    BayesInfer.InferenceResult result = bayesInferrer.infer(readSummaryMap, MAP);
//    
//    assertFalse(result.isDenovo());
//    assertEquals("846600 => [CC,CC,CC]", Arrays.asList(CC,CC,CC), result.getMaxTrioGenoType());
//  }
//
//  @Test
//  public void testTrioPos763769MAP() throws IOException {
//    Map<TrioIndividual, ReadSummary> readSummaryMap =
//        expRunner.getReadSummaryMap(763769L, expRunner.getReadMap("chr1", 763769L));
//    BayesInfer.InferenceResult result = bayesInferrer.infer(readSummaryMap, MAP);
//    
//    assertFalse(result.isDenovo());
//    assertEquals("763769 => [AA,AA,AA]", Arrays.asList(AA,AA,AA), result.getMaxTrioGenoType());
//  }
//
//  @Test
//  public void testTrioPos1298169MAP() throws IOException {
//    Map<TrioIndividual, ReadSummary> readSummaryMap =
//        expRunner.getReadSummaryMap(1298169L, expRunner.getReadMap("chr1", 1298169L));
//    BayesInfer.InferenceResult result = bayesInferrer.infer(readSummaryMap, MAP);
//
//    assertFalse(result.isDenovo());
//    assertEquals("1298169 => [TT,TT,TT]", Arrays.asList(TT, TT, TT), result.getMaxTrioGenoType());
//  }
//  
//  @Test
//  @Ignore("Known Borderline Failure")
//  /*chr1,70041751,readCounts=DAD:{T=2, C=58};MOM:{T=2, C=51};
//   * CHILD:{T=8, C=28},maxGenoType=[CC, CC, CT],isDenovo=true
//   */
//  public void testTrioPos70041751MAP() throws IOException {
//    Map<TrioIndividual, ReadSummary> readSummaryMap =
//        expRunner.getReadSummaryMap(70041751L, expRunner.getReadMap("chr1", 70041751L));
//    BayesInfer.InferenceResult result = bayesInferrer.infer(readSummaryMap, MAP);
//
//    assertEquals("70041751 => [CC,CC,CC]", Arrays.asList(CC, CC, CC), result.getMaxTrioGenoType());
//    assertFalse(result.isDenovo());
//  }
//
//  @Test
//  /*chr1,149035163,readCounts=DAD:{T=24, A=2, C=225, -=5};MOM:{T=22, G=3, A=6, C=223, -=2};
//   * CHILD:{T=34, G=1, A=2, C=218, -=1},maxGenoType=[CC, CC, CT],isDenovo=true
//   */
//  public void testTrioPos149035163MAP() throws IOException {
//    Map<TrioIndividual, ReadSummary> readSummaryMap =
//        expRunner.getReadSummaryMap(149035163L, expRunner.getReadMap("chr1", 149035163L));
//    BayesInfer.InferenceResult result = bayesInferrer.infer(readSummaryMap, MAP);
//
//    assertEquals("149035163 => [CC,CC,CC]", Arrays.asList(CC, CC, CC), result.getMaxTrioGenoType());
//    assertFalse(result.isDenovo());
//  }
}
