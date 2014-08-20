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

import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.DAD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.MOM;

import com.google.api.services.genomics.model.Variant;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;
import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.collect.FluentIterable;

import java.util.Arrays;
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

  public Variant pop(TrioIndividual person) {
    if (bufferMap.get(person).isEmpty()) {
      throw new IllegalStateException("Trying to pop from empty queue");
    }
    return bufferMap.get(person).removeFirst();
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
                  .from(bufferMap.get(person))
                  .transform(getStartAndEnd)) + "]";
          }
        }));
  }
  
  public Map<TrioIndividual, Deque<Variant>> getBufferMap() {
    return bufferMap;
  }
  
  public static void main(String[] args) {
    Variant v = new Variant().setPosition(1L).setEnd(1000L);
    VariantsBuffer vbuf = new VariantsBuffer();
    vbuf.push(DAD, v); vbuf.push(DAD, v); vbuf.push(DAD, v);
    vbuf.push(MOM, v); vbuf.push(MOM, new Variant().setPosition(3L).setEnd(5000L));
    System.out.println(vbuf);
  }
  
}
