package com.google.cloud.genomics.denovo;

import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;

import java.util.HashMap;
import java.util.Map;

/*
 * Bayes net data structure
 */
public class BayesNet<T,V> {
  public Map<T, Node<T,V>> nodeMap;

  public BayesNet() {
    nodeMap = new HashMap<>();
  }

  public void addNode(Node<T,V> node) {
    nodeMap.put(node.id, node);
  }
}
