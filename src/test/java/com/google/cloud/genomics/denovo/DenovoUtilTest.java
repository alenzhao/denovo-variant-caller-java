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
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.AC;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.AG;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.AT;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.CC;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.CG;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.CT;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.GG;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.TG;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.TT;
import static org.junit.Assert.assertEquals;

import com.google.api.services.genomics.model.Call;

import org.junit.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Some tests for the DenovoUtil class
 */
public class DenovoUtilTest {

  @Test
  public void testGetGenotypeFrominfoFieldPhased() {
    List<Integer> genotype = Arrays.asList(0, 0);
    Call call =
        new Call().setInfo(Collections.singletonMap("GT", Collections.singletonList("0|0")));
    assertEquals(genotype, DenovoUtil.getGenotypeFromInfoField(call).get());
  }

  @Test
  public void testGetGenotypeFrominfoFieldUnPhased() {
    List<Integer> genotype = Arrays.asList(0, 1);
    Call call =
        new Call().setInfo(Collections.singletonMap("GT", Collections.singletonList("0/1")));
    assertEquals(genotype, DenovoUtil.getGenotypeFromInfoField(call).get());
  }
  
  @Test(expected = NumberFormatException.class)
  public void testGetGenotypeFrominfoFieldNotInt() {
    List<Integer> genotype = Arrays.asList(0, 1);
    Call call =
        new Call().setInfo(Collections.singletonMap("GT", Collections.singletonList("./1")));
    assertEquals(genotype, DenovoUtil.getGenotypeFromInfoField(call).get());
  }
  
  /**
   * Test method for {@link com.google.cloud.genomics.denovo.DenovoUtil#checkTrioGenoTypeIsDenovo(java.util.List)}.
   */
  @Test
  public void testCheckTrioGenoTypeIsDenovo() {
    
    // Denovo Inheritance cases
    assertEquals("TT|AA,AA", true,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AA, AA, TT)));
    assertEquals("AT|AA,AA", true,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AA, AA, AT)));
    assertEquals("TT|AA,AC", true,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AA, AC, TT)));
    assertEquals("CC|AA,AC", true,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AA, AC, CC)));
    assertEquals("AT|AA,AC", true,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AA, AC, AT)));
    assertEquals("AA|AA,CG", true,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AA, CG, AA)));
    assertEquals("GG|AC,TG", true,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AC, TG, GG)));
    assertEquals("AA|AC,TG", true,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AC, TG, AA)));
    assertEquals("AC|AC,TG", true,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AC, TG, AC)));
    assertEquals("TG|AC,TG", true,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AC, TG, TG)));

    // Normal Inheritance cases
    assertEquals("AT|AC,TG", false,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AC, TG, AT)));
    assertEquals("CT|AC,TG", false,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AC, TG, CT)));
    assertEquals("AG|AC,TG", false,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AC, TG, AG)));
    assertEquals("CG|AC,TG", false,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AC, TG, CG)));
    assertEquals("AA|AA,AC", false,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AA, AC, AA)));
    assertEquals("AC|AA,AC", false,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AA, AC, AC)));
    assertEquals("AC|AA,CG", true,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AA, CG, AA)));
    assertEquals("AG|AA,CG", true,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AA, CG, AA)));
    assertEquals("AA|AC,AC", false,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AC, AC, AA)));
    assertEquals("AC|AC,AC", false,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AC, AC, AC)));
    assertEquals("CC|AC,AC", false,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AC, AC, CC)));
  }

}
