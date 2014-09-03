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
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.GT;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotype.TT;
import static org.junit.Assert.assertEquals;

import com.google.cloud.genomics.denovo.DenovoUtil.Chromosome;

import org.junit.Test;

import java.util.Arrays;

/**
 * Some tests for the DenovoUtil class
 */
public class DenovoUtilTest extends DenovoTest {

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
    assertEquals("GG|AC,GT", true,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AC, GT, GG)));
    assertEquals("AA|AC,GT", true,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AC, GT, AA)));
    assertEquals("AC|AC,GT", true,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AC, GT, AC)));
    assertEquals("GT|AC,GT", true,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AC, GT, GT)));

    // Normal Inheritance cases
    assertEquals("AT|AC,GT", false,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AC, GT, AT)));
    assertEquals("CT|AC,GT", false,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AC, GT, CT)));
    assertEquals("AG|AC,GT", false,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AC, GT, AG)));
    assertEquals("CG|AC,GT", false,
        DenovoUtil.checkTrioGenoTypeIsDenovo(Arrays.asList(AC, GT, CG)));
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
  
  @Test(expected=IllegalArgumentException.class)
  public void testChromosome_chr100() {
     Chromosome.fromString("chr100");
  }

  @Test
  public void testChromosome_chr1() {
    assertEquals(Chromosome.CHR1, Chromosome.fromString("chr1"));
  }

  @Test
  public void testChromosome_CHR1() {
    assertEquals(Chromosome.CHR1, Chromosome.fromString("CHR1"));
  }
  
  @Test
  public void testChromosome_1() {
    assertEquals(Chromosome.CHR1, Chromosome.fromString("chr1"));
  }
  
  @Test
  public void testChromosome_X() {
    assertEquals(Chromosome.CHRX, Chromosome.fromString("chrX"));
  }

  @Test
  public void testChromosome_x() {
    assertEquals(Chromosome.CHRX, Chromosome.fromString("chrx"));
  }
}