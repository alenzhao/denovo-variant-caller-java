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

import static com.google.cloud.genomics.denovo.DenovoUtil.Genotypes.AA;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotypes.AC;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotypes.AG;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotypes.AT;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotypes.CC;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotypes.CG;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotypes.CT;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotypes.GG;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotypes.TG;
import static com.google.cloud.genomics.denovo.DenovoUtil.Genotypes.TT;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import com.google.cloud.genomics.denovo.DenovoUtil.Genotypes;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;

import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

/**
 * Tests DenovoBayesNet
 */
public class DenovoBayesNetTest {

  /**
   * Test method for
   * {@link com.google.cloud.genomics.denovo.DenovoBayesNet#DenovoBayesNet(double, double)}.
   */

  private DenovoBayesNet dbn;
  private Map<List<Genotypes>, Double> conditionalProbabilityTable;
  private double EPS = 1e-12;

  @Before
  public void setUp() {
    dbn = new DenovoBayesNet(1e-2, 1e-8);

    conditionalProbabilityTable = new HashMap<>();
    int numGenotypes = Genotypes.values().length;
    for (Genotypes genotype : Genotypes.values()) {
      conditionalProbabilityTable.put(Collections.singletonList(genotype),
          Double.valueOf(1.0 / numGenotypes));
    }

    // makes sure conditionalProbabilityTable is set up properly
    AllTests.assertSumsToOne(conditionalProbabilityTable.values(), EPS);
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
  public void testAddNodeNodeOfTrioIndividualGenotypes() {
    Node<TrioIndividual, Genotypes> dadNode =
        new Node<>(TrioIndividual.DAD, null, conditionalProbabilityTable);

    Node<TrioIndividual, Genotypes> momNode =
        new Node<>(TrioIndividual.MOM, null, conditionalProbabilityTable);

    Node<TrioIndividual, Genotypes> childNode = new Node<>(TrioIndividual.CHILD,
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
    Map<List<Genotypes>, Double> cpt = dbn.createConditionalProbabilityTable(individual);
    int numGenotypes = Genotypes.values().length;

    // check keys
    assertEquals(
        new HashSet<>(
            Arrays.asList(
                Collections.singletonList(Genotypes.AA),
                Collections.singletonList(Genotypes.AC),
                Collections.singletonList(Genotypes.AT),
                Collections.singletonList(Genotypes.AG),
                Collections.singletonList(Genotypes.CC),
                Collections.singletonList(Genotypes.CT),
                Collections.singletonList(Genotypes.CG),
                Collections.singletonList(Genotypes.TT),
                Collections.singletonList(Genotypes.TG),
                Collections.singletonList(Genotypes.GG))),
        cpt.keySet());

    // check values
    double probValue = 1.0 / numGenotypes;
    for (Double value : cpt.values()) {
      assertEquals(probValue, value, EPS);
    }

    // check total probability Values
    AllTests.assertSumsToOne(cpt.values(), EPS);
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
    Map<List<Genotypes>, Double> cpt = dbn.createConditionalProbabilityTable(TrioIndividual.CHILD);

    // check some key values
    assertEquals("TT|AA,AA", 1e-8, cpt.get(Arrays.asList(AA, AA, TT)), 1e-12);
    assertEquals("AT|AA,AA", 1e-8, cpt.get(Arrays.asList(AA, AA, AT)), 1e-12);
    assertEquals("TT|AA,AC", 1e-8, cpt.get(Arrays.asList(AA, AC, TT)), 1e-12);
    assertEquals("CC|AA,AC", 1e-8, cpt.get(Arrays.asList(AA, AC, CC)), 1e-12);
    assertEquals("AT|AA,AC", 1e-8, cpt.get(Arrays.asList(AA, AC, AT)), 1e-12);
    assertEquals("AA|AA,AC", 0.5, cpt.get(Arrays.asList(AA, AC, AA)), 1e-7);
    assertEquals("AC|AA,AC", 0.5, cpt.get(Arrays.asList(AA, AC, AC)), 1e-7);
    assertEquals("GG|AC,TG", 1e-8, cpt.get(Arrays.asList(AC, TG, GG)), 1e-12);
    assertEquals("AA|AC,TG", 1e-8, cpt.get(Arrays.asList(AC, TG, AA)), 1e-12);
    assertEquals("AC|AC,TG", 1e-8, cpt.get(Arrays.asList(AC, TG, AC)), 1e-12);
    assertEquals("TG|AC,TG", 1e-8, cpt.get(Arrays.asList(AC, TG, TG)), 1e-12);
    assertEquals("AT|AC,TG", 0.25, cpt.get(Arrays.asList(AC, TG, AT)), 1e-7);
    assertEquals("CT|AC,TG", 0.25, cpt.get(Arrays.asList(AC, TG, CT)), 1e-7);
    assertEquals("AG|AC,TG", 0.25, cpt.get(Arrays.asList(AC, TG, AG)), 1e-7);
    assertEquals("CG|AC,TG", 0.25, cpt.get(Arrays.asList(AC, TG, CG)), 1e-7);
  }

  @Test
  @Ignore("Known Failure")
  public void testChildCreateConditionalProbabilityTableValuesHarderCases() {
    Map<List<Genotypes>, Double> cpt = dbn.createConditionalProbabilityTable(TrioIndividual.CHILD);

    assertEquals("AA|AC,AC", 0.25, cpt.get(Arrays.asList(AC, AC, AA)), 1e-7);
    assertEquals("AC|AC,AC", 0.5, cpt.get(Arrays.asList(AC, AC, AC)), 1e-7);
    assertEquals("CC|AC,AC", 0.25, cpt.get(Arrays.asList(AC, AC, CC)), 1e-7);

    DenovoBayesNet.printConditionalProbabilityTable(System.out, cpt);
  }

  @Test
  public void testChildCreateConditionalProbabilityTableTotalProbability() {
    Map<List<Genotypes>, Double> cpt = dbn.createConditionalProbabilityTable(TrioIndividual.CHILD);

    // Sanity check - probabilities should add up to 1.0 (almost)
    for (Genotypes genoTypeDad : Genotypes.values()) {
      for (Genotypes genoTypeMom : Genotypes.values()) {
        double totProb = 0.0;
        for (Genotypes genoTypeChild : Genotypes.values()) {
          totProb += cpt.get(Arrays.asList(genoTypeDad, genoTypeMom, genoTypeChild));
        }
        assertEquals(1.0, totProb, EPS);
      }
    }
  }
}
