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

import com.google.api.services.genomics.model.Variant;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;
import com.google.common.base.Joiner;
import com.google.common.collect.FluentIterable;

import java.util.Deque;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;

/**
 * A Map consisting of queues which hold variant calls from a trio
 */
public class VariantsBuffer {

  private Map<TrioIndividual, Deque<Variant>> bufferMap = new HashMap<>();
    
  public VariantsBuffer() { 
    for (TrioIndividual person : TrioIndividual.values()) {
      bufferMap.put(person, new LinkedList<Variant>());
    }
  }
  
  public void push(TrioIndividual person, Variant variant) {
    bufferMap.get(person).addLast(variant);
  }

  public void pop(TrioIndividual person) {
    if (bufferMap.get(person).isEmpty()) {
      throw new IllegalStateException("Trying to pop from empty queue");
    }
    bufferMap.get(person).removeFirst();
  }

  /*
   * Returns 0 if the buffer is empty for that person otherwise coord position
   */
  public Long getStartPosition(TrioIndividual person) {
    return bufferMap.get(person).isEmpty() ? 0 : bufferMap.get(person).getFirst().getPosition();
  }

  /*
   * Returns 0 if the buffer is empty for that person otherwise coord position
   */
  public Long getEndPosition(TrioIndividual person) {
    return bufferMap.get(person).isEmpty() ? 0 : bufferMap.get(person).getLast().getEnd();
  }
  
  @Override
  public int hashCode() {
    return bufferMap.hashCode();
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    for (TrioIndividual person : TrioIndividual.values()) {
      sb.append(person);
      sb.append(" : [");
      sb.append(
          Joiner.on(',')
          .join(FluentIterable<String>.from(bufferMap.get(person))));
      sb.append(" : ]\n");
    }
    return sb.toString();
  }
  
}
