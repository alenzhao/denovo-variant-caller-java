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

import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.CHILD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.DAD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.MOM;

import com.google.api.services.genomics.model.Call;
import com.google.api.services.genomics.model.Variant;
import com.google.cloud.genomics.denovo.DenovoUtil.Allele;
import com.google.cloud.genomics.denovo.DenovoUtil.Genotype;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;
import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.base.Optional;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.Range;

import org.javatuples.Pair;

import java.util.Arrays;
import java.util.Deque;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;

/**
 * A Map consisting of queues which hold variant calls from a trio
 */
public class VariantsBuffer {

  private Map<TrioIndividual, Deque<Pair<Variant,Call>>> bufferMap = new HashMap<>();

  public VariantsBuffer() {
    for (TrioIndividual person : TrioIndividual.values()) {
      bufferMap.put(person, new LinkedList<Pair<Variant,Call>>());
    }
  }

  public void push(TrioIndividual person, Pair<Variant,Call> pair) {
    getQueue(person).addLast(pair);
  }

  public Pair<Variant,Call> pop(TrioIndividual person) {
    if (isEmpty(person)) {
      throw new IllegalStateException("Trying to pop from empty queue");
    }
    return getQueue(person).removeFirst();
  }

  /*
   * Checks if first variant in CHILD buffer can be processed
   */
  public boolean canProcess() {
    return (!isEmpty(CHILD) && getEndPosition(MOM) >= getEndPosition(CHILD)
        && getEndPosition(DAD) >= getEndPosition(CHILD)
        && getStartPosition(MOM) <= getStartPosition(CHILD)
        && getStartPosition(DAD) <= getStartPosition(CHILD));
  }

  /*
   * Returns 0 if the buffer is empty for that person otherwise coord position
   */
  public Long getStartPosition(TrioIndividual person) {
    return isEmpty(person) ? 0 : getQueue(person).getFirst().getValue0().getPosition();
  }

  /*
   * Returns 0 if the buffer is empty for that person otherwise coord position
   */
  public Long getEndPosition(TrioIndividual person) {
    return isEmpty(person) ? 0 : getQueue(person).getLast().getValue0().getEnd();
  }

  @Override
  public int hashCode() {
    return bufferMap.hashCode();
  }

  @Override
  public String toString() {
    final Function<Pair<Variant,Call>, String> getStartAndEnd = new Function<Pair<Variant,Call>, String>() {
      @Override
      public String apply(Pair<Variant,Call> pair) {
        // TODO(smoitra): Auto-generated method stub
        return pair.getValue0().getPosition().toString() + "-" 
            + pair.getValue0().getEnd().toString();
      }
    };

    return Joiner.on(", ").join(FluentIterable.from(Arrays.asList(TrioIndividual.values()))
        .transform(new Function<TrioIndividual, String>() {
          @Override
          public String apply(TrioIndividual person) {
            return person.toString() + ":[" + Joiner.on(",").join(
                FluentIterable.from(getQueue(person)).transform(getStartAndEnd)) + "]";
          }
        }));
  }

  public Map<TrioIndividual, Deque<Pair<Variant,Call>>> getBufferMap() {
    return bufferMap;
  }

  public Deque<Pair<Variant,Call>> getQueue(TrioIndividual person) {
    return bufferMap.get(person);
  }

  public boolean isEmpty(TrioIndividual person) {
    return getQueue(person).isEmpty();
  }

  /**
   * Evicts parent variants that are no longer needed
   */
  private void evictParents() {
    if (isEmpty(CHILD)) {
      return;
    }

    for (TrioIndividual parent : Arrays.asList(MOM, DAD)) {
      while (!isEmpty(parent) && getFirst(parent).getValue0().getEnd() < getStartPosition(CHILD)) {
        pop(parent);
      }
    }
  }

  public Pair<Variant,Call> getFirst(TrioIndividual person) {
    return getQueue(person).getFirst();
  }

  /*
   * Check if a call passes filters for queue and then add it
   */
  public boolean checkAndAdd(TrioIndividual person, Pair<Variant, Call> pair) {
    if (person == CHILD && !isSnp(pair)) {
      return false;
    }

    if (callContainsDot(pair.getValue1())) {
      return false;
    }
    
    if (!passesFilter(pair.getValue1())) {
      return false;
    }
    push(person, pair);
    return true;
  }

  /*
   * Returns the next available PositionCall
   */
  public PositionCall retrieveNextCall() {
    evictParents();

    Pair<Variant, Call> childSNP = getNextSNP(CHILD);
    Long snpPosition = childSNP.getValue0().getPosition();
    String referenceBase = childSNP.getValue0().getReferenceBases();
    Map<TrioIndividual, Genotype> genotypeMap = new HashMap<>();

    for (TrioIndividual person : TrioIndividual.values()) {

      if (person == CHILD) {
        genotypeMap.put(CHILD, getGenotypeFromSNP(childSNP));
        continue;
      }
      Pair<Variant, Call> parentVariant = Optional.of(getMatchingPair(person, snpPosition)).get();
      genotypeMap.put(person, isSnp(parentVariant) 
          ? getGenotypeFromSNP(parentVariant) 
          : Genotype.valueOf(referenceBase + referenceBase));
    }
    return new PositionCall(snpPosition, genotypeMap);
  }

  /**
   * Return the genotype corresponding to the SNP
   */
  private Genotype getGenotypeFromSNP(Pair<Variant, Call> pair) {
    Variant variant = pair.getValue0();
    Call call = pair.getValue1();
    Allele[] allelePair = new Allele[2];

    for (int idx = 0; idx < 2; idx++) {
      int gtidx = call.getGenotype().get(idx);
      allelePair[idx] = gtidx > 0 ? Allele.valueOf(variant.getAlternateBases().get(gtidx - 1))
          : Allele.valueOf(variant.getReferenceBases());
    }
    return Genotype.valueOfPairAlleles(allelePair[0], allelePair[1]);
  }

  /*
   * Check if a call has PASS filter annotation
   */
  private boolean passesFilter(Call call) {
    return call.getInfo().containsKey("FILTER") && call.getInfo().get("FILTER").equals("PASS");
  }

  private boolean callContainsDot(Call call) {
    // Check if call is unknown
    return call.getGenotype().contains(-1);
  }
  
  private boolean isSnp(Pair<Variant, Call> pair) {
    Variant variant = pair.getValue0();
    Call call = pair.getValue1();
    if (variant.getEnd() != variant.getPosition() + Long.valueOf(1L)) {
      return false;
    }

    // check for Indels
    if (variant.getReferenceBases().length() != 1) {
      return false;
    }

    for (int gtIdx : call.getGenotype()) {
      if (gtIdx > 0 && variant.getAlternateBases().get(gtIdx - 1).length() != 1) {
        return false;
      }
    }
    return true;
  }

  /**
   * Get Matching variant from the queue
   */
  private Pair<Variant, Call> getMatchingPair(TrioIndividual person, Long snpPosition) {
    for (Pair<Variant, Call> pair : getQueue(person)) {
      Variant variant = pair.getValue0();

      if (!isSnp(pair)) {
        throw new IllegalStateException("Expected SNP ; got : " + pair);
      }

      if (Range.closedOpen(variant.getPosition(), variant.getEnd()).contains(snpPosition)) {
        return pair;
      }
    }
    return null;
  }

  /**
   * Get next available SNP
   */
  private Pair<Variant, Call> getNextSNP(TrioIndividual person) {
    Pair<Variant, Call> pair = getFirst(person);

    if (!isSnp(pair)) { 
      throw new IllegalStateException("Expected SNP : got " + pair);
    }
    return pair;
  }

  /*
   * A data class for storing position
   */
  public static class PositionCall {
    public final Long position;
    public final Map<TrioIndividual, Genotype> genotypeMap;

    public PositionCall(Long snpPosition, Map<TrioIndividual, Genotype> map) {
      this.position = snpPosition;
      this.genotypeMap = map;
    }

    @Override
    public String toString() {
      return "[" + position.toString() + "," + genotypeMap.toString() + "]";
    }
  }
}
