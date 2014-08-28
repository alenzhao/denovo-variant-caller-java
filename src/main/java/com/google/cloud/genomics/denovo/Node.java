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
class Node<T, V> {
  private final T id;
  private final List<Node<T, V>> parents;
  private final Map<List<V>, Double> conditionalProbabilityTable;

  public Node(T individual, List<Node<T, V>> parents, Map<List<V>, Double> cpt) {
    this.id = individual;
    this.parents = parents;
    this.conditionalProbabilityTable = cpt;
  }

  /**
   * @return the id
   */
  T getId() {
    return id;
  }

  /**
   * @return the parents
   */
  List<Node<T, V>> getParents() {
    return parents;
  }

  /**
   * @return the conditionalProbabilityTable
   */
  Map<List<V>, Double> getConditionalProbabilityTable() {
    return conditionalProbabilityTable;
  }
}
