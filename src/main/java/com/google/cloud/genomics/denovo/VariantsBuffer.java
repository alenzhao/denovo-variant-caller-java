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
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;
import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.base.Optional;
import com.google.common.base.Predicate;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

import org.javatuples.Pair;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * A Map consisting of queues which hold variant calls from a trio
 */
public class VariantsBuffer {

  private Map<TrioIndividual, Deque<Variant>> bufferMap = new HashMap<>();
  private final Map<TrioIndividual, String> personToCallsetIdMap;
  
  public VariantsBuffer(ExperimentRunner expRunner) {
    personToCallsetIdMap = expRunner.getPersonToCallsetIdMap();
    for (TrioIndividual person : TrioIndividual.values()) {
      bufferMap.put(person, new LinkedList<Variant>());
    }
  }
  
  public void push(TrioIndividual person, Variant variant) {
    getQueue(person).addLast(variant);
  }

  public Variant pop(TrioIndividual person) {
    if (isEmpty(person)) {
      throw new IllegalStateException("Trying to pop from empty queue");
    }
    return getQueue(person).removeFirst();
  }

  /*
   * Checks if first variant in CHILD buffer can be processed 
   */
  public boolean canProcess() {
    return (!isEmpty(CHILD) 
        && getEndPosition(MOM) >= getEndPosition(CHILD)
        && getEndPosition(DAD) >= getEndPosition(CHILD)
        && getStartPosition(MOM) <= getStartPosition(CHILD)
        && getStartPosition(DAD) <= getStartPosition(CHILD));  
  }
  
  /*
   * Returns 0 if the buffer is empty for that person otherwise coord position
   */
  public Long getStartPosition(TrioIndividual person) {
    return isEmpty(person) ? 0 : getQueue(person).getFirst().getPosition();
  }

  /*
   * Returns 0 if the buffer is empty for that person otherwise coord position
   */
  public Long getEndPosition(TrioIndividual person) {
    return isEmpty(person) ? 0 : getQueue(person).getLast().getEnd();
  }
  
  @Override
  public int hashCode() {
    return bufferMap.hashCode();
  }

  @Override
  public String toString() {
    final Function<Variant, String> getStartAndEnd = new Function<Variant, String>() {
      @Override
      public String apply(Variant input) {
        // TODO(smoitra): Auto-generated method stub
        return input.getPosition().toString() + "-" + input.getEnd().toString();
      }
    };
    
    return Joiner.on(", ").join(FluentIterable
        .from(Arrays.asList(TrioIndividual.values()))
        .transform(new Function<TrioIndividual, String>() {
          @Override
          public String apply(TrioIndividual person) {
            return person.toString() + ":[" + Joiner.on(",").join(
                FluentIterable
                  .from(getQueue(person))
                  .transform(getStartAndEnd)) + "]";
          }
        }));
  }
  
  public Map<TrioIndividual, Deque<Variant>> getBufferMap() {
    return bufferMap;
  }

  public Deque<Variant> getQueue(TrioIndividual person) {
    return bufferMap.get(person);
  }
  
  public boolean isEmpty(TrioIndividual person) {
    return getQueue(person).isEmpty();
  }
  
  /**
   * Evicts parent variants that are no longer needed 
   */
  private void evictParents() {
    if (isEmpty(CHILD)) return;
    
    for (TrioIndividual parent : Arrays.asList(MOM, DAD)) {
      while(!isEmpty(parent) && getEndPosition(parent) <= getStartPosition(CHILD)) {
        pop(parent);
      }
    }
  }

  public Variant getFirst(TrioIndividual person) {
    return getQueue(person).getFirst();
  }

  public Map<TrioIndividual,Iterable<PositionwiseCalls>> retreiveNextCalls() {
    evictParents();
    Variant firstChildVariant = getFirst(CHILD);
    Long startPosition = firstChildVariant.getPosition();
    Long endPosition = firstChildVariant.getEnd();

    Map<TrioIndividual,Iterable<PositionwiseCalls>> callsMap = new HashMap<>();
    for(TrioIndividual person : TrioIndividual.values()) {
      callsMap.put(person, getPositionWiseCallsInRange(startPosition, endPosition, person));
    }
    
    return callsMap;
  }
  
  private Iterable<PositionwiseCalls> getPositionWiseCallsInRange(final Long startPosition,
      final Long endPosition, final TrioIndividual person) {

    Iterable<Variant> variantsInRange = getVariantsInRange(person, startPosition, endPosition);
    Iterable<PositionwiseCalls> positionwiseCalls = Iterables.concat(
        Iterables.transform(variantsInRange, new Function<Variant, Iterable<PositionwiseCalls>>() {
          @Override
          public Iterable<PositionwiseCalls> apply(Variant variant) {
            return convertVariantToPositionwiseCalls(variant, person);
          }
        }));
    return Iterables.filter(positionwiseCalls, new Predicate<PositionwiseCalls>() {
      @Override
      public boolean apply(PositionwiseCalls pcall) {
        return pcall.position >= startPosition && pcall.position <= endPosition;
      }
    });
  }
  
  private Iterable<Variant> getVariantsInRange(TrioIndividual person, final Long startPosition, 
      final Long endPosition) {
    Deque<Variant> queue = getQueue(person);
    return Iterables.filter(queue, new Predicate<Variant>() {
      @Override
      public boolean apply(Variant variant) {
        return !(variant.getEnd() < startPosition || variant.getPosition() > endPosition); 
      }});
  }
  
  private Iterable<PositionwiseCalls> convertVariantToPositionwiseCalls(Variant variant,
      TrioIndividual person) {
    List<PositionwiseCalls> callsList = new LinkedList<>();

    // All calls are conserved
    if (variant.getInfo() != null && variant.getInfo().containsKey("BLOCKAVG_min30p3a")) {
      for (Long pos = variant.getPosition(); pos < variant.getEnd(); pos++) {
        callsList.add(new PositionwiseCalls(pos, Pair.with("0", "0")));
      }
      return callsList;
    }

    // Calls have alternative alleles
    String referenceBases = variant.getReferenceBases();
    List<String> alternateBases = variant.getAlternateBases();
    Call call =
        Optional.of(DenovoUtil.getCallInVariant(variant, person, personToCallsetIdMap)).get();
    List<Integer> genotype = call.getGenotype();

    String[] baseStrings = new String[2];

    // call contains reference sequence or not
    int referenceIndex = genotype.indexOf(Integer.valueOf(0));
    if (referenceIndex != -1) {
      baseStrings[0] = referenceBases;
      if (genotype.get(1 - referenceIndex) == Integer.valueOf(0)) {
        baseStrings[1] = referenceBases;  
      } else {
        baseStrings[1] = alternateBases.get(genotype.get(1 - referenceIndex));
      }
    } else {
      for (int i = 0; i < 2; i++) {
        baseStrings[i] = alternateBases.get(genotype.get(i));
      }
    }

    int idx = 0;
    for (Long pos = variant.getPosition(); pos < variant.getEnd(); pos++) {
      String[] posStrings = new String[2];
      for (int i = 0; i < 2; i++) {
        // Deletion
        if (idx >= baseStrings[i].length()) {
          posStrings[i] = "-";
          continue;
        }
        // Insertion
        if (pos == variant.getEnd() - 1L) {
          posStrings[i] = baseStrings[i].substring(idx);
        }
        // Regular base
        posStrings[i] = String.valueOf(baseStrings[i].charAt(idx));
      }
      callsList.add(new PositionwiseCalls(pos, Pair.with(posStrings[0], posStrings[1])));
      idx++;
    }
    return callsList;
  }

  /*
   * A data class for storing position
   */
  public static class PositionwiseCalls {
    public Long position;
    public Pair<String, String> calls;

    public PositionwiseCalls(Long position, Pair<String, String> calls) {
      this.position = position;
      this.calls = calls;
    }

    @Override
    public String toString() {
      return "[" + position.toString() + "," + calls.toString() + "]";
    }
  }
}