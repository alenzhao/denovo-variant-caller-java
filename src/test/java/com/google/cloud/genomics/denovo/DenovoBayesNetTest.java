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
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.CHILD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.DAD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.MOM;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import com.google.cloud.genomics.denovo.DenovoBayesNet.InferenceResult;
import com.google.cloud.genomics.denovo.DenovoUtil.Allele;
import com.google.cloud.genomics.denovo.DenovoUtil.Genotype;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;

import org.junit.Before;
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

  private DenovoBayesNet dbn;
  private static final Map<List<Genotype>, Double> conditionalProbabilityTable;
  private static final double EPS_SMALL = 1e-12;
  private static final double EPS_SMALLER = 1e-20;
  private static final double EPS_MEDIUM = 1e-5;
  private static final double EPS_LARGE = 1;

  static {
    conditionalProbabilityTable = new HashMap<>();
    int numGenotypes = Genotype.values().length;
    for (Genotype genotype : Genotype.values()) {
      conditionalProbabilityTable.put(Collections.singletonList(genotype),
          Double.valueOf(1.0 / numGenotypes));
    }

    // makes sure conditionalProbabilityTable is set up properly
    assertSumsToOne(conditionalProbabilityTable.values(), EPS_SMALL);
  }

  @Before
  public void setUp() {
    dbn = new DenovoBayesNet(1e-2, 1e-8);
  }
  
  @Test
  public void testDenovoBayesNet() {
    assertNotNull(dbn);
    assertEquals(1e-2, dbn.getSequenceErrorRate(), EPS_SMALL);
    assertEquals(1e-8, dbn.getDenovoMutationRate(), EPS_SMALL);
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

    assertEquals(dbn.getNodeMap().get(TrioIndividual.CHILD).getParents(), 
        Arrays.asList(dadNode, momNode));
    assertEquals(dbn.getNodeMap().get(TrioIndividual.DAD).getParents(), null);
    assertEquals(dbn.getNodeMap().get(TrioIndividual.MOM).getParents(), null);
  }

  /**
   * Test method for {@link com.google.cloud.genomics.denovo.DenovoBayesNet#createConditionalProbabilityTable(com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual)}
   * .
   */
  public void testParentCreateConditionalProbabilityTable(TrioIndividual person) {
    Map<List<Genotype>, Double> cpt = dbn.createConditionalProbabilityTable(person);

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
    assertEquals(0.0625, cpt.get(Collections.singletonList(Genotype.AA)), EPS_SMALL);
    assertEquals(0.125, cpt.get(Collections.singletonList(Genotype.AC)), EPS_SMALL);
    assertEquals(0.125, cpt.get(Collections.singletonList(Genotype.AT)), EPS_SMALL);
    assertEquals(0.125, cpt.get(Collections.singletonList(Genotype.AG)), EPS_SMALL);
    assertEquals(0.0625, cpt.get(Collections.singletonList(Genotype.CC)), EPS_SMALL);
    assertEquals(0.125, cpt.get(Collections.singletonList(Genotype.CT)), EPS_SMALL);
    assertEquals(0.125, cpt.get(Collections.singletonList(Genotype.CG)), EPS_SMALL);
    assertEquals(0.0625, cpt.get(Collections.singletonList(Genotype.TT)), EPS_SMALL);
    assertEquals(0.125, cpt.get(Collections.singletonList(Genotype.GT)), EPS_SMALL);
    assertEquals(0.0625, cpt.get(Collections.singletonList(Genotype.GG)), EPS_SMALL);
    
    // check total probability Values
    assertSumsToOne(cpt.values(), EPS_SMALL);
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
        assertEquals(1.0, totProb, EPS_SMALL);
      }
    }
  }
  
  @Test
  public void testBaseLogLikelihood_HomozygousMatch() {
    for(Genotype genotype : Genotype.values()) {
      if (genotype.isHomozygous()) {
        for (Allele allele : genotype.getAlleleSet()) {
          assertEquals(Math.log(0.99), dbn.getBaseLogLikelihood(genotype, allele), EPS_SMALL);
        }
      }
    }
  }

  @Test
  public void testBaseLogLikelihood_HomozygousNoMatch() {
    for(Genotype genotype : Genotype.values()) {
      if (genotype.isHomozygous()) {
        for (Allele allele : EnumSet.complementOf(genotype.getAlleleSet())) {
          assertEquals(Math.log(1e-2 / 3), dbn.getBaseLogLikelihood(genotype, allele), EPS_SMALL);
        }
      }
    }
  }

  @Test
  public void testBaseLogLikelihood_HeterozygousMatch() {
    for(Genotype genotype : Genotype.values()) {
      if (!genotype.isHomozygous()) {
        for (Allele allele : genotype.getAlleleSet()) {
          assertEquals(Math.log(0.5 * (1 - 2 * 1e-2 / 3)), 
              dbn.getBaseLogLikelihood(genotype, allele), EPS_SMALL);
        }
      }
    }
  }

  @Test
  public void testBaseLogLikelihood_HeterozygousNoMatch() {
    for(Genotype genotype : Genotype.values()) {
      if (!genotype.isHomozygous()) {
        for (Allele allele : EnumSet.complementOf(genotype.getAlleleSet())) {
          assertEquals(Math.log(1e-2 / 3), dbn.getBaseLogLikelihood(genotype, allele), EPS_SMALL);
        }
      }
    }
  }

  @Test
  public void testGetReadSummaryLogLikelihood_AllSame() {
    ReadSummary summary = createSameReadSummary();
    Map<Genotype, Double> llMap = 
        dbn.getReadSummaryLogLikelihood(summary);
    for (Genotype gt : Genotype.values()) {
      assertEquals(dbn.getBaseLogLikelihood(gt, Allele.A)*40, llMap.get(gt), EPS_SMALLER);  
    }
  }
  
  @Test
  public void testGetReadSummaryLogLikelihood_AlmostSame() {
    ReadSummary summary = createAlmostSameReadSummary();
    Map<Genotype, Double> llMap = 
        dbn.getReadSummaryLogLikelihood(summary);
    for (Genotype gt : Genotype.values()) {
      double ll = 0;
      ll += dbn.getBaseLogLikelihood(gt, Allele.A) * 38;
      ll += dbn.getBaseLogLikelihood(gt, Allele.C) * 2;
      ll += dbn.getBaseLogLikelihood(gt, Allele.G) * 3;
      assertEquals(ll, llMap.get(gt), EPS_SMALL);  
    }
  }
  
  @Test
  public void testGetReadSummaryLogLikelihood_AllGap() {
    Map<Genotype, Double> llMap = 
        dbn.getReadSummaryLogLikelihood(new ReadSummary(new HashMap<Allele, Integer>()));
    
    for(Genotype gt : Genotype.values()) {
      assertEquals(0.0, llMap.get(gt), EPS_SMALL);  
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
  
  @Test
  public void testgetLogLikelihoodFromCPT_Parents() {
    // check values
    assertEquals(Math.log(0.0625),dbn.getLogLikelihoodFromCPT(DAD, AA), EPS_SMALL); 
    assertEquals(Math.log(0.125),dbn.getLogLikelihoodFromCPT(DAD, AC), EPS_SMALL);
    assertEquals(Math.log(0.0625),dbn.getLogLikelihoodFromCPT(MOM, CC), EPS_SMALL); 
    assertEquals(Math.log(0.125),dbn.getLogLikelihoodFromCPT(MOM, CG), EPS_SMALL);
  }

  @Test
  public void testgetLogLikelihoodFromCPT_Child() {
    // check values
    Map<List<Genotype>, Double> cpt = dbn.createConditionalProbabilityTable(CHILD);
    assertEquals("CC|AA,AC", 1e-9, cpt.get(Arrays.asList(AA, AC, CC)), EPS_LARGE);
    assertEquals("CC|AA,AC", Math.log(1e-9), 
        dbn.getLogLikelihoodFromCPT(CHILD, AA, AC, CC), EPS_LARGE);
    assertEquals("AT|AA,AC", Math.log(1e-9), 
        dbn.getLogLikelihoodFromCPT(CHILD, AA, AC, AT), EPS_LARGE);
    assertEquals("AA|AA,AC", Math.log(0.5), 
        dbn.getLogLikelihoodFromCPT(CHILD, AA, AC, AA), EPS_MEDIUM);
    assertEquals("AA|AC,AC", Math.log(0.25), 
        dbn.getLogLikelihoodFromCPT(CHILD, AC, AC, AA), EPS_MEDIUM);
    assertEquals("AC|AC,AC", Math.log(0.5), 
        dbn.getLogLikelihoodFromCPT(CHILD, AC, AC, AC), EPS_MEDIUM);
    assertEquals("CC|AC,AC", Math.log(0.25), 
        dbn.getLogLikelihoodFromCPT(CHILD, AC, AC, CC), EPS_MEDIUM);
  }
  
  @Test
  public void getTrioGenotypeLogLikelihood_AllSame() {
    ReadSummary summary = createSameReadSummary();
    Map<TrioIndividual, ReadSummary> summaryMap = 
        createMapReadSummary(summary, summary, summary);
    Map<TrioIndividual, Map<Genotype, Double>> individualLogLikelihood = 
        dbn.getIndividualLogLikelihood(summaryMap);
    Map<Genotype, Double> llMap = 
        dbn.getReadSummaryLogLikelihood(summary);

    assertEquals(llMap.get(AA) * 3, 
        dbn.getTrioGenotypeLogLikelihood(individualLogLikelihood, AA, AA, AA), EPS_MEDIUM);
    assertEquals(llMap.get(AA) * 2 + llMap.get(AT), 
        dbn.getTrioGenotypeLogLikelihood(individualLogLikelihood, AA, AT, AA), EPS_MEDIUM);
    assertEquals(llMap.get(AA) + llMap.get(AT) + llMap.get(CG), 
        dbn.getTrioGenotypeLogLikelihood(individualLogLikelihood, AA, AT, CG), EPS_MEDIUM);
    assertEquals(llMap.get(AA) + llMap.get(AT) + llMap.get(CG), 
        dbn.getTrioGenotypeLogLikelihood(individualLogLikelihood, CG, AT, AA), EPS_MEDIUM);
  }

  @Test
  public void getTrioGenotypeLogLikelihood_AlmostSame() {
    ReadSummary summary = createSameReadSummary();
    ReadSummary summary2 = createAlmostSameReadSummary();
    Map<TrioIndividual, ReadSummary> summaryMap = 
        createMapReadSummary(summary, summary2, summary);
    Map<TrioIndividual, Map<Genotype, Double>> individualLogLikelihood = 
        dbn.getIndividualLogLikelihood(summaryMap);
    Map<Genotype, Double> llMap = 
        dbn.getReadSummaryLogLikelihood(summary);
    Map<Genotype, Double> llMap2 = 
        dbn.getReadSummaryLogLikelihood(summary2);

    assertEquals(llMap.get(AA) * 2 + llMap2.get(AA), 
        dbn.getTrioGenotypeLogLikelihood(individualLogLikelihood, AA, AA, AA), EPS_MEDIUM);
    assertEquals(llMap.get(AA) * 2 + llMap2.get(AT), 
        dbn.getTrioGenotypeLogLikelihood(individualLogLikelihood, AA, AT, AA), EPS_MEDIUM);
    assertEquals(llMap.get(AA) + llMap2.get(AT) + llMap.get(CG), 
        dbn.getTrioGenotypeLogLikelihood(individualLogLikelihood, AA, AT, CG), EPS_MEDIUM);
    assertEquals(llMap.get(AT) + llMap2.get(AA) + llMap.get(CG), 
        dbn.getTrioGenotypeLogLikelihood(individualLogLikelihood, CG, AA, AT), EPS_MEDIUM);
  }

  @Test
  public void getTrioGenotypeLogLikelihood() {
    assertEquals(dbn.getLogLikelihoodFromCPT(DAD, AA) + 
        dbn.getLogLikelihoodFromCPT(MOM, AA) + 
        dbn.getLogLikelihoodFromCPT(CHILD, AA, AA, AA), 
        dbn.getRelationshipLogLikelihood(AA, AA, AA), EPS_MEDIUM);

    assertEquals(dbn.getLogLikelihoodFromCPT(DAD, AA) + 
        dbn.getLogLikelihoodFromCPT(MOM, AT) + 
        dbn.getLogLikelihoodFromCPT(CHILD, AA, AT, AA), 
        dbn.getRelationshipLogLikelihood(AA, AT, AA), EPS_MEDIUM);
    
    assertEquals(dbn.getLogLikelihoodFromCPT(DAD, AA) + 
        dbn.getLogLikelihoodFromCPT(MOM, AT) + 
        dbn.getLogLikelihoodFromCPT(CHILD, AA, AT, CG), 
        dbn.getRelationshipLogLikelihood(AA, AT, CG), EPS_MEDIUM);


    assertEquals(dbn.getLogLikelihoodFromCPT(DAD, CG) + 
        dbn.getLogLikelihoodFromCPT(MOM, AT) + 
        dbn.getLogLikelihoodFromCPT(CHILD, CG, AT, AA), 
        dbn.getRelationshipLogLikelihood(CG, AT, AA), EPS_MEDIUM);
  }
  
  @Test(expected = IllegalArgumentException.class)
  public void testgetLogLikelihoodFromCPT_IncorrectArgs1() {
    dbn.getLogLikelihoodFromCPT(CHILD, AA);
  }

  @Test(expected = IllegalArgumentException.class)
  public void testgetLogLikelihoodFromCPT_IncorrectArgs2() {
    dbn.getLogLikelihoodFromCPT(CHILD, AA, AC);
  }

  @Test(expected = IllegalArgumentException.class)
  public void testgetLogLikelihoodFromCPT_IncorrectArgs3() {
    dbn.getLogLikelihoodFromCPT(DAD, AA, AC, AG);
  }

  @Test(expected = IllegalArgumentException.class)
  public void testgetLogLikelihoodFromCPT_IncorrectArgs4() {
    dbn.getLogLikelihoodFromCPT(MOM, AA, AC, AG);
  }
  
  @Test
  public void testPeformInference_AllSame() {
    ReadSummary summary = createSameReadSummary();
    Map<Genotype, Double> llMap = 
        dbn.getReadSummaryLogLikelihood(summary);
    Map<TrioIndividual, ReadSummary> summaryMap = 
       createMapReadSummary(summary, summary, summary);
    InferenceResult result = dbn.performInference(summaryMap);
    
    double answer = llMap.get(AA) * 3 + dbn.getLogLikelihoodFromCPT(CHILD, AA,AA,AA) 
        + dbn.getLogLikelihoodFromCPT(DAD, AA) + dbn.getLogLikelihoodFromCPT(MOM, AA);
    
    assertEquals("ll(AA,AA,AA)", answer, result.maxLikelihood, EPS_SMALL);
  }

  @Test
  public void testPeformInference_AlmostSame() {
    ReadSummary summary = createAlmostSameReadSummary();
    Map<Genotype, Double> llMap = 
        dbn.getReadSummaryLogLikelihood(summary);
    Map<TrioIndividual, ReadSummary> summaryMap = 
        createMapReadSummary(summary, summary, summary);
    InferenceResult result = dbn.performInference(summaryMap);
    
    double answer = llMap.get(AA) * 3 + dbn.getLogLikelihoodFromCPT(CHILD, AA,AA,AA) 
        + dbn.getLogLikelihoodFromCPT(DAD, AA) + dbn.getLogLikelihoodFromCPT(MOM, AA);
    
    assertEquals("ll(AA,AA,AA)", answer, result.maxLikelihood, EPS_SMALL);
  }
}