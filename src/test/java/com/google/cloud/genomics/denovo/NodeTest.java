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

import static org.junit.Assert.assertEquals;

import com.google.cloud.genomics.denovo.DenovoUtil.Genotypes;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Tests Node class
 */
public class NodeTest {

  private Map<List<Genotypes>, Double> conditionalProbabilityTable;
  private double EPS = 1e-12;

  /**
   * Test method for {@link com.google.cloud.genomics.denovo.Node#Node(java.lang.Object,
   * java.util.List, java.util.Map)}.
   */
  @Before
  public void setUp() {
    conditionalProbabilityTable = new HashMap<>();
    int numGenotypes = Genotypes.values().length;
    for (Genotypes genotype : Genotypes.values()) {
      conditionalProbabilityTable.put(Collections.singletonList(genotype),
          Double.valueOf(1.0 / numGenotypes));
    }

    // makes sure conditionalProbabilityTable is set up properly
    AllTests.assertSumsToOne(conditionalProbabilityTable.values(), EPS);
  }

  @After
  public void tearDown() {
    conditionalProbabilityTable = null;
  }

  @Test
  public void testSingleNode() {

    Node<TrioIndividual, Genotypes> dadNode =
        new Node<>(TrioIndividual.DAD, null, conditionalProbabilityTable);
    assertEquals(TrioIndividual.DAD, dadNode.getId());
    assertEquals(null, dadNode.getParents());
    assertEquals(conditionalProbabilityTable, dadNode.getConditionalProbabilityTable());
  }

  @Test
  public void testLinkedNode() {

    Node<TrioIndividual, Genotypes> dadNode =
        new Node<>(TrioIndividual.DAD, null, conditionalProbabilityTable);

    Node<TrioIndividual, Genotypes> childNode = new Node<>(TrioIndividual.CHILD,
        Collections.singletonList(dadNode), conditionalProbabilityTable);

    assertEquals(null, dadNode.getParents());
    assertEquals(1, childNode.getParents().size());
    assertEquals(dadNode, childNode.getParents().get(0));
  }
}
