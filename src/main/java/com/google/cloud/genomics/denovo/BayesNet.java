package com.google.cloud.genomics.denovo;

import java.util.Map;

/*
 * Bayes net data structure
 */
public abstract class BayesNet<T, V> {
  public Map<T, Node<T, V>> nodeMap;

  public abstract void addNode(Node<T, V> node);
}
