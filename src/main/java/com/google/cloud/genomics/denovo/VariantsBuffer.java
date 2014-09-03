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

import static com.google.cloud.genomics.denovo.DenovoUtil.TrioMember.CHILD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioMember.DAD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioMember.MOM;

import com.google.api.services.genomics.model.Call;
import com.google.api.services.genomics.model.Variant;
import com.google.cloud.genomics.denovo.DenovoUtil.Allele;
import com.google.cloud.genomics.denovo.DenovoUtil.Genotype;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioMember;
import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.base.Optional;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.Range;

import org.javatuples.Pair;

import java.util.Arrays;
import java.util.Deque;
import java.util.LinkedList;
import java.util.Map;
import java.util.TreeMap;

/** Provide Buffering for fetching variants
 */
class VariantsBuffer {
  private Map<TrioMember, Deque<Pair<Variant,Call>>> bufferMap = new TreeMap<>();
  private Map<TrioMember, Long> mostRecentStartPosition = new TreeMap<>();

  VariantsBuffer() {
    for (TrioMember person : TrioMember.values()) {
      bufferMap.put(person, new LinkedList<Pair<Variant,Call>>());
      mostRecentStartPosition.put(person, 0L);
    }
  }

  /** Enqueue element
   * @param person trio member
   * @param pair pair of variant and call
   */
  void push(TrioMember person, Pair<Variant,Call> pair) {
    getQueue(person).addLast(pair);
    mostRecentStartPosition.put(person, pair.getValue0().getPosition());
  }

  /** Pop first element from queue
   * @param person
   * @return first element in queue
   */
  Pair<Variant,Call> pop(TrioMember person) {
    if (isEmpty(person)) {
      throw new IllegalStateException("Trying to pop from empty queue");
    }
    return getQueue(person).removeFirst();
  }

  /** Checks if first variant in CHILD buffer can be processed
   */
  boolean canProcess() {
    return (!isEmpty(CHILD) 
        && getMostRecentStartPosition(MOM) >= getStartPosition(CHILD)
        && getMostRecentStartPosition(DAD) >= getStartPosition(CHILD));
  }

  /**
   * @param person trio member
   * @return the start position of the most recently retrived variant
   */
  private Long getMostRecentStartPosition(TrioMember person) {
    return mostRecentStartPosition.get(person);
  }

  /** start position of first elememt in queue
   * @param person
   * @return 0 if the buffer is empty for that person otherwise coord position
   */
  Long getStartPosition(TrioMember person) {
    return isEmpty(person) ? 0 : getQueue(person).getFirst().getValue0().getPosition();
  }

  /** end position of last element in queue
   * @param person
   * @return 0 if the buffer is empty for that person otherwise coord position
   */
  Long getEndPosition(TrioMember person) {
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
        return pair.getValue0().getPosition().toString() + "-" 
            + pair.getValue0().getEnd().toString();
      }
    };

    return Joiner.on(", ").join(FluentIterable.from(Arrays.asList(TrioMember.values()))
        .transform(new Function<TrioMember, String>() {
          @Override
          public String apply(TrioMember person) {
            return person.toString() + ":[" + Joiner.on(",").join(
                FluentIterable.from(getQueue(person)).transform(getStartAndEnd)) + "]";
          }
        }));
  }

  Map<TrioMember, Deque<Pair<Variant,Call>>> getBufferMap() {
    return bufferMap;
  }

  Deque<Pair<Variant,Call>> getQueue(TrioMember person) {
    return bufferMap.get(person);
  }

  /**
   * @param person trio member
   * @return are all varaints exhausted
   */
  boolean isEmpty(TrioMember person) {
    return getQueue(person).isEmpty();
  }

  /** Evicts parent variants that are no longer needed
   */
  private void evictParents() {
    if (isEmpty(CHILD)) {
      return;
    }

    for (TrioMember parent : Arrays.asList(MOM, DAD)) {
      while (!isEmpty(parent) && getFirst(parent).getValue0().getEnd() < getStartPosition(CHILD)) {
        pop(parent);
      }
    }
  }

  /**
   * @param person
   * @return first element in queue without removing
   */
  Pair<Variant,Call> getFirst(TrioMember person) {
    return getQueue(person).getFirst();
  }

  /** Check if a call passes filters for queue and then add it
   * @param person
   * @param pair
   * @return success
   */
  boolean checkAndAdd(TrioMember person, Pair<Variant, Call> pair) {
    
    Call variant = pair.getValue1();
    if (person == CHILD && !isSnp(pair) 
        ||callContainsDot(variant) 
        || !callIsBiAllelic(variant) 
        || !passesFilter(variant) 
        || isInsertion(pair) 
        || isDeletion(pair)) {
      return false;
    }
    push(person, pair);
    return true;
  }

  /** Does the call contain two alleles
   */
  private boolean callIsBiAllelic(Call call) {
    return call.getGenotype().size() == 2;
  }

  /**
   * Returns the next available PositionCall
   */
  PositionCall retrieveNextCall() {
    evictParents();

    Pair<Variant, Call> childSNP = getNextSNP(CHILD);
    Long snpPosition = childSNP.getValue0().getPosition();
    String referenceBase = childSNP.getValue0().getReferenceBases();
    Map<TrioMember, Genotype> genotypeMap = new TreeMap<>();

    for (TrioMember person : TrioMember.values()) {

      if (person == CHILD) {
        genotypeMap.put(CHILD, getGenotypeFromSNP(childSNP));
        continue;
      }
      
      Optional<Pair<Variant, Call>> parentVariant = 
          Optional.fromNullable(getMatchingPair(person, snpPosition));
      if (!parentVariant.isPresent()) {
        return null;
      }
      
      genotypeMap.put(person, isSnp(parentVariant.get()) 
          ? getGenotypeFromSNP(parentVariant.get()) 
          : Genotype.valueOf(referenceBase + referenceBase));
    }
    return new PositionCall(snpPosition, genotypeMap);
  }

  /** Return the genotype corresponding to the SNP
   * @param pair of variant,call corresponding to snp
   * @return genotype corresponding to the SNP
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

  /** Check if a call has PASS filter annotation
   * @param call
   * @return success
   */
  private boolean passesFilter(Call call) {
    return call.getInfo().containsKey("FILTER") && 
        call.getInfo().get("FILTER").size() == 1 &&
        call.getInfo().get("FILTER").get(0).equals("PASS");
  }

  /** Call contains '.' character 
   * @param call
   * @return succcess
   */
  private boolean callContainsDot(Call call) {
    // Check if call is unknown
    return call.getGenotype().contains(-1);
  }
  
  /** Is variant, call pair a snp
   * @param pair
   * @return success
   */
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
  
  /** Is variant, call a deletion
   * @param pair
   * @return match
   */
  private boolean isDeletion(Pair<Variant, Call> pair) {
    Variant variant = pair.getValue0();
    return variant.getReferenceBases().length() != 1 ;
  }

  /** is variant , call an insertion
   * @param pair
   * @return match
   */
  private boolean isInsertion(Pair<Variant, Call> pair) {
    Variant variant = pair.getValue0();
    Call call = pair.getValue1();
    for (int gtIdx : call.getGenotype()) {
      if (gtIdx > 0 && variant.getAlternateBases().get(gtIdx - 1).length() != 1) {
        return true;
      }
    }
    return false;
  }


  /** Get Matching variant from the queue
   * @param person trio member
   * @param snpPosition
   * @return Get matching variant call pairs at matching position
   */
  private Pair<Variant, Call> getMatchingPair(TrioMember person, Long snpPosition) {
    for (Pair<Variant, Call> pair : getQueue(person)) {
      Variant variant = pair.getValue0();

      if (Range.closedOpen(variant.getPosition(), variant.getEnd()).contains(snpPosition)) {
        return pair;
      }
    }
    return null;
  }

  /** Get next available SNP
   * @param person trio member
   * @return next pair of variant, call corresponding to snp for person
   */
  private Pair<Variant, Call> getNextSNP(TrioMember person) {
    Pair<Variant, Call> pair = getFirst(person);

    if (!isSnp(pair)) { 
      throw new IllegalStateException("Expected SNP : got " + pair);
    }
    return pair;
  }

  /** A data class for storing position
   */
  static class PositionCall {
    private final Long position;
    private final Map<TrioMember, Genotype> genotypeMap;

    PositionCall(Long snpPosition, Map<TrioMember, Genotype> map) {
      this.position = snpPosition;
      this.genotypeMap = map;
    }

    /*
     * Is call denovo
     */
    boolean isDenovo() {
      return DenovoUtil.checkTrioGenoTypeIsDenovo(getGenotypeMap().get(DAD),
            getGenotypeMap().get(MOM),getGenotypeMap().get(CHILD));
    }
    
    @Override
    public String toString() {
      return "[" + getPosition().toString() + "," + getGenotypeMap().toString() + "]";
    }

    /**
     * @return the position
     */
    public Long getPosition() {
      return position;
    }

    /**
     * @return the genotypeMap
     */
    public Map<TrioMember, Genotype> getGenotypeMap() {
      return genotypeMap;
    }
  }
}
