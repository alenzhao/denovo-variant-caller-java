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
import static org.junit.Assert.assertNotNull;

import com.google.cloud.genomics.denovo.DenovoUtil.Allele;
import com.google.cloud.genomics.denovo.DenovoUtil.Genotype;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;

import org.junit.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

/**
 * Tests DenovoBayesNet
 */
public class DenovoBayesNetTest extends DenovoTest {

  private static final DenovoBayesNet dbn;
  private static final Map<List<Genotype>, Double> conditionalProbabilityTable;
  private static final double EPS = 1e-12;
  private static final double EPS_SMALL = 1e-20;

  static {
    dbn = new DenovoBayesNet(1e-2, 1e-8);

    conditionalProbabilityTable = new HashMap<>();
    int numGenotypes = Genotype.values().length;
    for (Genotype genotype : Genotype.values()) {
      conditionalProbabilityTable.put(Collections.singletonList(genotype),
          Double.valueOf(1.0 / numGenotypes));
    }

    // makes sure conditionalProbabilityTable is set up properly
    assertSumsToOne(conditionalProbabilityTable.values(), EPS);
  }
  
  @Test
  public void testDenovoBayesNet() {
    assertNotNull(dbn);
    assertEquals(1e-2, dbn.getSequenceErrorRate(), EPS);
    assertEquals(1e-8, dbn.getDenovoMutationRate(), EPS);
  }

  /**
   * Test method for {@link com.google.cloud.genomics.denovo.DenovoBayesNet#addNode(com.google.cloud.genomics.denovo.Node)}
   */
  @Test
  public void testAddNode_NodeOfTrioIndividualGenotypes() {
    Node<TrioIndividual, Genotype> dadNode =
        new Node<>(TrioIndividual.DAD, null, conditionalProbabilityTable);

    Node<TrioIndividual, Genotype> momNode =
        new Node<>(TrioIndividual.MOM, null, conditionalProbabilityTable);

    Node<TrioIndividual, Genotype> childNode = new Node<>(TrioIndividual.CHILD,
        Arrays.asList(dadNode, momNode), conditionalProbabilityTable);

    dbn.addNode(dadNode);
    dbn.addNode(momNode);
    dbn.addNode(childNode);

    assertEquals(dbn.getNodeMap().get(TrioIndividual.DAD), dadNode);
    assertEquals(dbn.getNodeMap().get(TrioIndividual.MOM), momNode);
    assertEquals(dbn.getNodeMap().get(TrioIndividual.CHILD), childNode);

    assertEquals(dbn.getNodeMap().get(TrioIndividual.CHILD).getParents(), Arrays.asList(dadNode, momNode));
    assertEquals(dbn.getNodeMap().get(TrioIndividual.DAD).getParents(), null);
    assertEquals(dbn.getNodeMap().get(TrioIndividual.MOM).getParents(), null);
  }

  /**
   * Test method for {@link com.google.cloud.genomics.denovo.DenovoBayesNet#createConditionalProbabilityTable(com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual)}
   * .
   */
  public void testParentCreateConditionalProbabilityTable(TrioIndividual individual) {
    Map<List<Genotype>, Double> cpt = dbn.createConditionalProbabilityTable(individual);

    // check keys
    assertEquals(
        new HashSet<>(
            Arrays.asList(
                Collections.singletonList(Genotype.AA),
                Collections.singletonList(Genotype.AC),
                Collections.singletonList(Genotype.AT),
                Collections.singletonList(Genotype.AG),
                Collections.singletonList(Genotype.CC),
                Collections.singletonList(Genotype.CT),
                Collections.singletonList(Genotype.CG),
                Collections.singletonList(Genotype.TT),
                Collections.singletonList(Genotype.GT),
                Collections.singletonList(Genotype.GG))),
        cpt.keySet());

    // check values
    assertEquals(0.0625, cpt.get(Collections.singletonList(Genotype.AA)), EPS);
    assertEquals(0.125, cpt.get(Collections.singletonList(Genotype.AC)), EPS);
    assertEquals(0.125, cpt.get(Collections.singletonList(Genotype.AT)), EPS);
    assertEquals(0.125, cpt.get(Collections.singletonList(Genotype.AG)), EPS);
    assertEquals(0.0625, cpt.get(Collections.singletonList(Genotype.CC)), EPS);
    assertEquals(0.125, cpt.get(Collections.singletonList(Genotype.CT)), EPS);
    assertEquals(0.125, cpt.get(Collections.singletonList(Genotype.CG)), EPS);
    assertEquals(0.0625, cpt.get(Collections.singletonList(Genotype.TT)), EPS);
    assertEquals(0.125, cpt.get(Collections.singletonList(Genotype.GT)), EPS);
    assertEquals(0.0625, cpt.get(Collections.singletonList(Genotype.GG)), EPS);
    
    // check total probability Values
    assertSumsToOne(cpt.values(), EPS);
  }

  @Test
  public void testDadCreateConditionalProbabilityTable() {
    testParentCreateConditionalProbabilityTable(TrioIndividual.DAD);
  }

  @Test
  public void testMomCreateConditionalProbabilityTable() {
    testParentCreateConditionalProbabilityTable(TrioIndividual.MOM);
  }

  @Test
  public void testChildCreateConditionalProbabilityTableValues() {
    Map<List<Genotype>, Double> cpt = dbn.createConditionalProbabilityTable(TrioIndividual.CHILD);

    // check some key values
    assertEquals("AA|AA,AA", 1.0, cpt.get(Arrays.asList(AA, AA, AA)), 1e-7);
    assertEquals("TT|AA,AA", 1e-9, cpt.get(Arrays.asList(AA, AA, TT)), 1e-9);
    assertEquals("AT|AA,AA", 1e-9, cpt.get(Arrays.asList(AA, AA, AT)), 1e-9);
    assertEquals("TT|AA,AC", 1e-9, cpt.get(Arrays.asList(AA, AC, TT)), 1e-9);
    assertEquals("CC|AA,AC", 1e-9, cpt.get(Arrays.asList(AA, AC, CC)), 1e-9);
    assertEquals("AT|AA,AC", 1e-9, cpt.get(Arrays.asList(AA, AC, AT)), 1e-9);
    assertEquals("AA|AA,AC", 0.5, cpt.get(Arrays.asList(AA, AC, AA)), 1e-7);
    assertEquals("AC|AA,AC", 0.5, cpt.get(Arrays.asList(AA, AC, AC)), 1e-7);
    assertEquals("GG|AC,GT", 1e-9, cpt.get(Arrays.asList(AC, GT, GG)), 1e-9);
    assertEquals("AA|AC,GT", 1e-9, cpt.get(Arrays.asList(AC, GT, AA)), 1e-9);
    assertEquals("AC|AC,GT", 1e-9, cpt.get(Arrays.asList(AC, GT, AC)), 1e-9);
    assertEquals("GT|AC,GT", 1e-9, cpt.get(Arrays.asList(AC, GT, GT)), 1e-9);
    assertEquals("AT|AC,GT", 0.25, cpt.get(Arrays.asList(AC, GT, AT)), 1e-7);
    assertEquals("CT|AC,GT", 0.25, cpt.get(Arrays.asList(AC, GT, CT)), 1e-7);
    assertEquals("AG|AC,GT", 0.25, cpt.get(Arrays.asList(AC, GT, AG)), 1e-7);
    assertEquals("CG|AC,GT", 0.25, cpt.get(Arrays.asList(AC, GT, CG)), 1e-7);
    assertEquals("AT|AA,GT", 0.5, cpt.get(Arrays.asList(AA, GT, AT)), 1e-7);
    assertEquals("AG|AA,GT", 0.5, cpt.get(Arrays.asList(AA, GT, AG)), 1e-7);
    assertEquals("GT|AA,GT", 1e-9, cpt.get(Arrays.asList(AA, GT, GT)), 1e-9);
    assertEquals("AC|AA,GT", 1e-9, cpt.get(Arrays.asList(AA, GT, AC)), 1e-9);
    
    // harder
    assertEquals("AA|AC,AC", 0.25, cpt.get(Arrays.asList(AC, AC, AA)), 1e-7);
    assertEquals("AC|AC,AC", 0.5, cpt.get(Arrays.asList(AC, AC, AC)), 1e-7);
    assertEquals("CC|AC,AC", 0.25, cpt.get(Arrays.asList(AC, AC, CC)), 1e-7);

  }

  @Test
  public void testChildCreateConditionalProbabilityTableTotalProbability() {
    Map<List<Genotype>, Double> cpt = dbn.createConditionalProbabilityTable(TrioIndividual.CHILD);

    // Sanity check - probabilities should add up to 1.0 (almost)
    for (Genotype genoTypeDad : Genotype.values()) {
      for (Genotype genoTypeMom : Genotype.values()) {
        double totProb = 0.0;
        for (Genotype genoTypeChild : Genotype.values()) {
          totProb += cpt.get(Arrays.asList(genoTypeDad, genoTypeMom, genoTypeChild));
        }
        assertEquals(1.0, totProb, EPS);
      }
    }
  }
  
  @Test
  public void testBaseLogLikelihood_HomozygousMatch() {
    for(Genotype genotype : Genotype.values()) {
      if (genotype.isHomozygous()) {
        for (Allele allele : genotype.getAlleleSet()) {
          assertEquals(Math.log(0.99), dbn.getBaseLogLikelihood(genotype, allele), EPS);
        }
      }
    }
  }

  @Test
  public void testBaseLogLikelihood_HomozygousNoMatch() {
    for(Genotype genotype : Genotype.values()) {
      if (genotype.isHomozygous()) {
        for (Allele allele : EnumSet.complementOf(genotype.getAlleleSet())) {
          assertEquals(Math.log(1e-2 / 3), dbn.getBaseLogLikelihood(genotype, allele), EPS);
        }
      }
    }
  }

  @Test
  public void testBaseLogLikelihood_HeterozygousMatch() {
    for(Genotype genotype : Genotype.values()) {
      if (!genotype.isHomozygous()) {
        for (Allele allele : genotype.getAlleleSet()) {
          assertEquals(Math.log(0.5 * (1 - 2 * 1e-2 / 3)), dbn.getBaseLogLikelihood(genotype, allele), EPS);
        }
      }
    }
  }

  @Test
  public void testBaseLogLikelihood_HeterozygousNoMatch() {
    for(Genotype genotype : Genotype.values()) {
      if (!genotype.isHomozygous()) {
        for (Allele allele : EnumSet.complementOf(genotype.getAlleleSet())) {
          assertEquals(Math.log(1e-2 / 3), dbn.getBaseLogLikelihood(genotype, allele), EPS);
        }
      }
    }
  }

  @Test
  public void testGetReadSummaryLogLikelihood_AllSame() {
    ReadSummary summary = BayesInferTest.createSameReadSummary();
    Map<Genotype, Double> llMap = 
        dbn.getReadSummaryLogLikelihood(summary);
    for (Genotype gt : Genotype.values()) {
      assertEquals(dbn.getBaseLogLikelihood(gt, Allele.A)*40, llMap.get(gt), EPS_SMALL);  
    }
  }
  
  @Test
  public void testGetReadSummaryLogLikelihood_AlmostSame() {
    ReadSummary summary = BayesInferTest.createAlmostSameReadSummary();
    Map<Genotype, Double> llMap = 
        dbn.getReadSummaryLogLikelihood(summary);
    for (Genotype gt : Genotype.values()) {
      double ll = 0;
      ll += dbn.getBaseLogLikelihood(gt, Allele.A) * 38;
      ll += dbn.getBaseLogLikelihood(gt, Allele.C) * 2;
      ll += dbn.getBaseLogLikelihood(gt, Allele.G) * 3;
      assertEquals(ll, llMap.get(gt), EPS);  
    }
  }
  
  @Test
  public void testGetReadSummaryLogLikelihood_AllGap() {
    Map<Genotype, Double> llMap = 
        dbn.getReadSummaryLogLikelihood(new ReadSummary(new HashMap<Allele, Integer>()));
    
    for(Genotype gt : Genotype.values()) {
      assertEquals(0.0, llMap.get(gt), EPS);  
    }
  }
  
  
  @Test(expected = NullPointerException.class)
  public void testBaseLogLikelihood_null() {
    dbn.getBaseLogLikelihood(AA, null);
  }
  
  @Test(expected = NullPointerException.class)
  public void testGetReadSummaryLogLikelihood_null() {
    dbn.getReadSummaryLogLikelihood(null);
  }
  
}
