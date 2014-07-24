package com.google.cloud.genomics.denovo;

import java.util.List;
import java.util.Map;

/*
 * Individual Node in the Bayes Net
 */
public class Node<T,V> {
  public T id;
  public List<Node<T,V>> parents;
  public Map<List<V>, Double> conditionalProbabilityTable;

  public Node(T individual, List<Node<T,V>> parents, Map<List<V>, Double> map) {
    this.id = individual;
    this.parents = parents;
    this.conditionalProbabilityTable = map;
  }
}
