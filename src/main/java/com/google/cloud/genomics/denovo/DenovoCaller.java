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

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import com.google.api.services.genomics.model.Call;
import com.google.api.services.genomics.model.Variant;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.CHILD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.DAD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.MOM;


public class DenovoCaller {

  private Map<TrioIndividual, Long> lastPosition;
  private Map<TrioIndividual, Call> lastCall;
  private Map<TrioIndividual, String> dictRelationCallsetId;


  public DenovoCaller(Map<TrioIndividual, String> dictRelationCallsetId) {

    this.dictRelationCallsetId = dictRelationCallsetId;

    lastPosition = new HashMap<>();
    lastCall = new HashMap<>();
    // Initialize the data Dicts
    for (TrioIndividual trioType : TrioIndividual.values()) {
      lastPosition.put(trioType, Long.valueOf(0L));
      lastCall.put(trioType, null);
    }
  }

  /*
   * Returns a trio of calls or 'none' if no suitable trio can be found
   */
  private Map<TrioIndividual, Call> getTrioCalls(Long currPosition) {
    Map<TrioIndividual, Call> overlappingCalls = new HashMap<>();
    for (TrioIndividual trioType : TrioIndividual.values()) {

      // Call is not present for this position
      if (lastPosition.get(trioType) < currPosition) {
        return null;
      }

      // Call is not diploid
      Call call = lastCall.get(trioType);
      List<Integer> genotype = DenovoUtil.getGenotype(call);
      if (genotype == null || genotype.size() != 2) {
        return null;
      }
      overlappingCalls.put(trioType, call);

    }
    return overlappingCalls;
  }

  /*
   * Logic for checking Denovo variant calls
   *
   * c1,c2 are child genotypes d1,d2 are dad genotypes m1,m2 are mom genotypes
   *
   * predicate1 = c1 \in {m1,m2} and c2 \in {d1,d2} predicate2 = c2 \in {m1,m2} and c1 \in {d1,d2}
   * predicate = not( predicate1 or predicate2)
   */
  private boolean checkTrioLogic(Map<TrioIndividual, List<Integer>> trioGenoTypes) {
    Iterator<Integer> childIterator = trioGenoTypes.get(CHILD).iterator();
    Integer childAllele1 = childIterator.next();
    Integer childAllele2 = childIterator.next();
    List<Integer> momGenoType = trioGenoTypes.get(MOM);
    List<Integer> dadGenoType = trioGenoTypes.get(DAD);

    boolean predicate1 = momGenoType.contains(childAllele1) & dadGenoType.contains(childAllele2);
    boolean predicate2 = momGenoType.contains(childAllele2) & dadGenoType.contains(childAllele1);
    boolean predicate = !(predicate1 | predicate2);

    return predicate;
  }

  /*
   * Iteration 1 of De novo calling Simply check if a child mutation is different from that of
   * either parents.
   */
  public DenovoResult callDenovoVariantIteration1(Variant variant) {

    // Get all the calls for that variant
    for (Call call : variant.getCalls()) {

      // Reject the call if it doesn't meet quality standards
      if (!call.getInfo().containsKey("GQX") && !call.getInfo().containsKey("QD")
          && !call.getInfo().containsKey("MQ")) {
        continue;
      }

      try {
        if (call.getInfo().containsKey("GQX")
            && Float.parseFloat(call.getInfo().get("GQX").get(0)) < ExperimentRunner.GQX_THRESH) {
          continue;
        }
        if (call.getInfo().containsKey("QD")
            && Float.parseFloat(call.getInfo().get("QD").get(0)) < ExperimentRunner.GQX_THRESH) {
          continue;
        }
        if (call.getInfo().containsKey("MQ")
            && Float.parseFloat(call.getInfo().get("MQ").get(0)) < ExperimentRunner.GQX_THRESH) {
          continue;
        }
      } catch (NumberFormatException e) {
        continue;
      }


      // Update the lastcall and the last position
      for (TrioIndividual trioType : TrioIndividual.values()) {
        if (call.getCallsetId().equals(dictRelationCallsetId.get(trioType))) {
          lastCall.put(trioType, call);
          if (call.getInfo().containsKey("END")) {
            lastPosition.put(trioType, Long.valueOf(call.getInfo().get("END").get(0)));
          } else {
            lastPosition.put(trioType, variant.getPosition());
          }
        }
      }
    }

    // Get the overlapping calls for this variant position
    Map<TrioIndividual, Call> trioCalls = getTrioCalls(variant.getPosition());

    if (trioCalls == null) {
      return null;
    }

    Map<TrioIndividual, List<Integer>> trioGenotypes = new HashMap<>();
    for (TrioIndividual trioType : TrioIndividual.values()) {
      List<Integer> genoTypeList = DenovoUtil.getGenotype(trioCalls.get(trioType));
      trioGenotypes.put(trioType, genoTypeList);
    }

    if (checkTrioLogic(trioGenotypes)) {
      StringBuilder detailsBuilder = new StringBuilder();
      for (TrioIndividual trioType : TrioIndividual.values()) {
        detailsBuilder.append(
            trioType + ":" + trioCalls.get(trioType).getInfo().get("GT").get(0) + ",");
      }
      detailsBuilder.deleteCharAt(detailsBuilder.length() - 1);
      String details = detailsBuilder.toString();

      return new DenovoResult(details);
    }
    return null;
  }
}
