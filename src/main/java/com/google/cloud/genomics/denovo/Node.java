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

import java.util.List;
import java.util.Map;

/*
 * Individual Node in the Bayes Net
 */
public class Node<T, V> {
  private T id;
  private List<Node<T, V>> parents;
  private Map<List<V>, Double> conditionalProbabilityTable;

  public Node(T individual, List<Node<T, V>> parents, Map<List<V>, Double> map) {
    this.setId(individual);
    this.setParents(parents);
    this.setConditionalProbabilityTable(map);
  }

  /**
   * @return the id
   */
  public T getId() {
    return id;
  }

  /**
   * @param id the id to set
   */
  public void setId(T id) {
    this.id = id;
  }

  /**
   * @return the parents
   */
  public List<Node<T, V>> getParents() {
    return parents;
  }

  /**
   * @param parents the parents to set
   */
  public void setParents(List<Node<T, V>> parents) {
    this.parents = parents;
  }

  /**
   * @return the conditionalProbabilityTable
   */
  public Map<List<V>, Double> getConditionalProbabilityTable() {
    return conditionalProbabilityTable;
  }

  /**
   * @param conditionalProbabilityTable the conditionalProbabilityTable to set
   */
  public void setConditionalProbabilityTable(Map<List<V>, Double> conditionalProbabilityTable) {
    this.conditionalProbabilityTable = conditionalProbabilityTable;
  }
}
