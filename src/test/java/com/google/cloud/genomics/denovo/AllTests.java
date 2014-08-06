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

import static org.junit.Assert.assertEquals;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

import java.util.Collection;

@RunWith(Suite.class)
@SuiteClasses({NodeTest.class,
  DenovoBayesNetTest.class,
  DenovoUtilTest.class,
  BayesInferMapTest.class,
  BayesInferBayesTest.class,
  BayesInferLRTTest.class
  })
public class AllTests {

  public static void assertSumsToOne(Collection<Double> collection, double EPS) {
    // makes sure conditionalProbabilityTable is set up properly
    double totProb = 0.0;
    for (Double prob : collection) {
      totProb += prob;
    }
    assertEquals(1.0, totProb, EPS);
  }
}
