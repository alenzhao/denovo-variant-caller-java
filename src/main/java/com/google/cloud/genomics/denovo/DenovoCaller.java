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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.api.services.genomics.model.Call;
import com.google.api.services.genomics.model.Variant;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual;
import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.base.Optional;
import com.google.common.collect.Iterables;

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
  private Optional<Map<TrioIndividual, Call>> getTrioCalls(Long currPosition) {
    Map<TrioIndividual, Call> overlappingCalls = new HashMap<>();
    for (TrioIndividual trioType : TrioIndividual.values()) {

      // Call is not present for this position
      if (lastPosition.get(trioType) < currPosition) {
        return Optional.absent();
      }

      // Call is not diploid
      Call call = lastCall.get(trioType);
      Optional<List<Integer>> genotypeOption = DenovoUtil.getGenotype(call);
      if (!genotypeOption.isPresent() || genotypeOption.get().size() != 2) {
        return Optional.absent();
      }
      overlappingCalls.put(trioType, call);
    }
    return Optional.of(overlappingCalls);
  }

  /*
   * Logic for checking Denovo variant calls
   *
   * c1,c2 are child genotypes d1,d2 are dad genotypes m1,m2 are mom genotypes
   *
   * predicate1 = c1 \in {m1,m2} and c2 \in {d1,d2} predicate2 = c2 \in {m1,m2} and c1 \in {d1,d2}
   * predicate = not( predicate1 or predicate2)
   */
  private boolean checkTrioLogic(Map<TrioIndividual, DiploidGenotype> trioGenoTypes) {
    DiploidGenotype childGenotype = trioGenoTypes.get(CHILD);
	Integer childAllele1 = childGenotype .getFirstAllele();
    Integer childAllele2 = childGenotype.getSecondAllele();
    
    List<Integer> momGenoType = trioGenoTypes.get(MOM).getAllAlleles();
    List<Integer> dadGenoType = trioGenoTypes.get(DAD).getAllAlleles();

    boolean childAllelesinMomAndDad = momGenoType.contains(childAllele1) && 
        dadGenoType.contains(childAllele2);
    boolean childAllelesinMomAndDadMirrored = momGenoType.contains(childAllele2) && 
        dadGenoType.contains(childAllele1);
    return !(childAllelesinMomAndDad || childAllelesinMomAndDadMirrored);
  }

  /*
   * Iteration 1 of De novo calling Simply check if a child mutation is different from that of
   * either parents.
   */
  public Optional<String> callDenovoFromVarstore(Variant variant) {

    // Get all the calls for that variant
    for (Call call : variant.getCalls()) {
      if (passesAllQualityFilters(call)) {

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
    }

    // Get the overlapping calls for this variant position
    Optional<Map<TrioIndividual, Call>> trioCallsOptional = getTrioCalls(variant.getPosition());
	
    if (!trioCallsOptional.isPresent()) {
      return Optional.absent();
    } 
    
    final Map<TrioIndividual, Call> trioCalls = trioCallsOptional.get();	
    
    Map<TrioIndividual, DiploidGenotype> trioGenotypes = new HashMap<>();
    for (TrioIndividual trioType : TrioIndividual.values()) {
      List<Integer> genoTypeList = DenovoUtil.getGenotype(trioCalls.get(trioType)).get();
      trioGenotypes.put(trioType, new DiploidGenotype(genoTypeList));
    }

    return checkTrioLogic(trioGenotypes) ? Optional.of( Joiner.on(",").join(Iterables.transform(
        Arrays.asList(TrioIndividual.values()),
        new Function<TrioIndividual,String>() {
          @Override
          public String apply(TrioIndividual trioType) {
            return String.format("%s:%s",trioType.name(),
                trioCalls.get(trioType).getInfo().get("GT").get(0));
          }
        }))) :  Optional.<String>absent(); 
  }

	/*
	 * Does the call pass all the quality filters
	 */
	private boolean passesAllQualityFilters(Call call) {
		List<String> qualityKeysPresent = new ArrayList<>();

		for (String qualityKey : Arrays.asList("GQX", "QD", "MQ")) {
			if (call.getInfo().containsKey(qualityKey)) {
				qualityKeysPresent.add(qualityKey);
			}
		}

		// no valid quality keys present
		if (qualityKeysPresent.size() < 1) {
			return false;
		}

		boolean passesFilters = true;
		for (String qualityKey : qualityKeysPresent) {
			if (!passesQualityFilter(call, qualityKey)) {
				passesFilters = false;
				break;
			}
		}
		return passesFilters;
	}

	/*
	 * Checks whether a call passes a particular quality filter
	 */
	private boolean passesQualityFilter(Call call, String qualityKey) {

		Float threshold = ExperimentRunner.qualityThresholdMap.get(qualityKey);

		try {
			Float parseFloat = Float.parseFloat(call.getInfo().get(qualityKey).get(0));
			if (parseFloat < threshold) {
				return false;
			}
		} catch(NumberFormatException e) {
			return false;
		}
		return true;
	}
  
  public static class DiploidGenotype {
	  private List<Integer> genotype;
	  public DiploidGenotype(List<Integer> genotype) {
		  if (genotype.size() != 2) {
			  throw new IllegalStateException("Expected Diploid Genotype ; got"
					  +genotype.toString());
		  }
		  this.genotype = genotype;
	  }
	  /*
	   * Get first or second allele (0/1)
	   */
	  public Integer getFirstAllele() {
		  return genotype.get(0);
	  }
	  public Integer getSecondAllele() {
		  return genotype.get(1);
	  }
	  public List<Integer> getAllAlleles() {
		  return genotype;
	  }
  }
  
}
