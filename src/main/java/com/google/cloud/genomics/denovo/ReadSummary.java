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

import com.google.api.services.genomics.model.Read;
import com.google.cloud.genomics.denovo.DenovoUtil.Allele;

import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * Summary counts of a reads at a particular position
 */
public class ReadSummary {
  private Map<Allele, Integer> count = new TreeMap<>();

  public ReadSummary() {
  }
  
  public ReadSummary(Map<Allele, Integer> count) {
    this.count = count;
  }
  
  public ReadSummary(List<Read> reads, Long candidatePosition) {
    for (Read read : reads) {
      // TODO : Figure out baseAtPos
      String alignedBases = read.getAlignedBases();
      Integer offset = (int) (candidatePosition - read.getPosition());
      String baseAtPos = alignedBases.substring(offset, offset + 1);
      
      if (baseAtPos.equals("-")) {
        continue;
      }
        
      Allele alleleAtPos = Allele.valueOf(baseAtPos);
      count.put(alleleAtPos,
          (count.containsKey(alleleAtPos) ? count.get(alleleAtPos) : 0) + 1);
    }
  }

  @Override
  public String toString() {
    return count.toString();
  }
  
  public Map<Allele, Integer> getCount() {
    return count;
  }

  public ReadSummary setCount(Map<Allele, Integer> count) {
    this.count = count;
    return this;
  }
}
