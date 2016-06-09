/*
 * Copyright 2014 Google Inc. All rights reserved.
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

import static com.google.cloud.genomics.denovo.DenovoUtil.TrioMember.DAD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioMember.MOM;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.api.services.genomics.model.Variant;
import com.google.api.services.genomics.model.VariantCall;
import com.google.cloud.genomics.denovo.DenovoUtil.Genotype;
import com.google.cloud.genomics.denovo.DenovoUtil.TrioMember;

import org.javatuples.Pair;
import org.junit.Before;
import org.junit.Test;

import java.util.HashMap;
import java.util.Map;

/**
 * Test Cases for Variants Buffer
 */
public class VariantsBufferTest extends DenovoTest {

  VariantsBuffer vbuf;
  VariantCall dummyCall;
  Variant variant;

  @Before
  public void setUp() throws Exception {
    vbuf = new VariantsBuffer();
    dummyCall = new VariantCall();
  }

  @Test
  public void testPush() {
    Variant v = new Variant().setStart(1L).setEnd(100001L);
    Pair<Variant, VariantCall> pair = Pair.with(v,dummyCall);

    vbuf.push(DAD, pair);

    assertEquals(1, vbuf.getBufferMap().get(DAD).size());
    assertEquals(pair, vbuf.getBufferMap().get(DAD).getFirst());
  }

  @Test
  public void testPop() {
    Variant v = new Variant().setStart(1L).setEnd(100001L);
    Pair<Variant, VariantCall> pair = Pair.with(v,dummyCall);
    vbuf.push(DAD, pair);
    Pair<Variant,VariantCall> popv = vbuf.pop(DAD);

    assertEquals(0, vbuf.getBufferMap().get(DAD).size());
    assertEquals(pair, popv);
  }

  @Test
  public void testGetStartPosition() {
    vbuf.push(DAD, Pair.with(new Variant().setStart(1L).setEnd(10001L), dummyCall));
    vbuf.push(DAD, Pair.with(new Variant().setStart(10002L).setEnd(10003L), dummyCall));

    assertEquals(Long.valueOf(0L), vbuf.getStartPosition(MOM));
    assertEquals(Long.valueOf(1L), vbuf.getStartPosition(DAD));
  }

  @Test
  public void testGetEndPosition() {
    vbuf.push(DAD, Pair.with(new Variant().setStart(1L).setEnd(10001L), dummyCall));
    vbuf.push(DAD, Pair.with(new Variant().setStart(10002L).setEnd(10003L), dummyCall));

    assertEquals(Long.valueOf(0L), vbuf.getEndPosition(MOM));
    assertEquals(Long.valueOf(10003L), vbuf.getEndPosition(DAD));
  }

  @Test
  public void testToString() {
    Variant v = new Variant().setStart(1L).setEnd(1000L);
    Pair<Variant, VariantCall> pair = Pair.with(v, dummyCall);
    vbuf.push(DAD, pair);
    vbuf.push(DAD, pair);
    vbuf.push(DAD, pair);
    vbuf.push(MOM, pair);
    vbuf.push(MOM, Pair.with(new Variant().setStart(3L).setEnd(5000L), dummyCall));
    assertEquals("CHILD:[], MOM:[1-1000,3-5000], DAD:[1-1000,1-1000,1-1000]", vbuf.toString());
  }

  // Position VariantCall tests
  @Test
  public void testPositionCall_1() {
    Map<TrioMember, Genotype> map = new HashMap<>();
    map.put(TrioMember.DAD, Genotype.AA);
    map.put(TrioMember.MOM, Genotype.AA);
    map.put(TrioMember.CHILD, Genotype.AC);
    VariantsBuffer.PositionCall pc = new VariantsBuffer.PositionCall(10000L, map);
    assertTrue(pc.isDenovo());
    assertEquals(map, pc.getGenotypeMap());
  }
}