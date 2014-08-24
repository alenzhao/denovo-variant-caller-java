/*
Copyright 2014 Google Inc. All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
package com.google.cloud.genomics.denovo;

import static org.junit.Assert.assertEquals;

import com.google.api.services.genomics.Genomics;

import org.junit.After;
import org.junit.Before;
import org.mockito.Mock;
import org.mockito.Mockito;
import org.mockito.MockitoAnnotations;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.Collection;

public abstract class DenovoTest {
  @Mock Genomics genomics;
  @Mock ExperimentRunner expRunner;
  
  @Mock Genomics.Callsets callsets;
  @Mock Genomics.Callsets.Search callsetSearch;

  @Mock Genomics.Reads reads;
  @Mock Genomics.Reads.Search readSearch;

  @Mock Genomics.Readsets readsets;
  @Mock Genomics.Readsets.Search readsetSearch;

  @Mock Genomics.Variants variants;
  @Mock Genomics.Variants.GetSummary variantSummary;
  @Mock Genomics.Variants.Search variantSearch;
  
  ByteArrayOutputStream outContent = new ByteArrayOutputStream();

  @Before
  public void initMocks() {
    MockitoAnnotations.initMocks(this);

    Mockito.when(genomics.callsets()).thenReturn(callsets);
    Mockito.when(genomics.reads()).thenReturn(reads);
    Mockito.when(genomics.readsets()).thenReturn(readsets);
    Mockito.when(genomics.variants()).thenReturn(variants);

    Mockito.when(readsetSearch.setFields(Mockito.anyString())).thenReturn(readsetSearch);
  }

  @Before
  public void setUpStreams() {
    System.setOut(new PrintStream(outContent));
  }

  @After
  public void cleanUpStreams() {
    System.setOut(null);
  }
  
  public static void assertSumsToOne(Collection<Double> collection, double EPS) {
    // makes sure conditionalProbabilityTable is set up properly
    double totProb = 0.0;
    for (Double prob : collection) {
      totProb += prob;
    }
    assertEquals(1.0, totProb, EPS);
  }

}