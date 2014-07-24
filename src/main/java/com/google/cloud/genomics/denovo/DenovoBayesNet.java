package com.google.cloud.genomics.denovo;

import com.google.cloud.genomics.denovo.DenovoUtil.Genotypes;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;

import java.util.HashMap;


/*
 * DenovoBayesNet implements Generic abstract BayesNet class
 */
public class DenovoBayesNet extends BayesNet<TrioIndividual, Genotypes> {


  public DenovoBayesNet() {
    nodeMap = new HashMap<>();
  }

  @Override
  public void addNode(Node<TrioIndividual, Genotypes> node) {
    nodeMap.put(node.id, node);
  }



}
