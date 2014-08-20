package com.google.cloud.genomics.denovo;

import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.DAD;
import static com.google.cloud.genomics.denovo.DenovoUtil.TrioIndividual.MOM;
import static org.junit.Assert.assertEquals;

import com.google.api.services.genomics.model.Variant;

import org.junit.Before;
import org.junit.Test;

/**
 * Test Cases for Variants Buffer
 */
public class VariantsBufferTest {

  VariantsBuffer vbuf;

  @Before
  public void setUp() throws Exception {
    vbuf = new VariantsBuffer();
  }

  @Test
  public void testPush() {
    Variant v = new Variant().setPosition(1L).setEnd(100001L);
    vbuf.push(DAD, v);

    assertEquals(1, vbuf.getBufferMap().get(DAD).size());
    assertEquals(v, vbuf.getBufferMap().get(DAD).getFirst());
  }

  @Test
  public void testPop() {
    Variant v = new Variant().setPosition(1L).setEnd(100001L);
    vbuf.push(DAD, v);
    Variant popv = vbuf.pop(DAD);

    assertEquals(0, vbuf.getBufferMap().get(DAD).size());
    assertEquals(v, popv);
  }

  @Test
  public void testGetStartPosition() {
    vbuf.push(DAD, new Variant().setPosition(1L).setEnd(10001L));
    vbuf.push(DAD, new Variant().setPosition(10002L).setEnd(10003L));

    assertEquals(Long.valueOf(0L), vbuf.getStartPosition(MOM));
    assertEquals(Long.valueOf(1L), vbuf.getStartPosition(DAD));
  }

  @Test
  public void testGetEndPosition() {
    vbuf.push(DAD, new Variant().setPosition(1L).setEnd(10001L));
    vbuf.push(DAD, new Variant().setPosition(10002L).setEnd(10003L));

    assertEquals(Long.valueOf(0L), vbuf.getEndPosition(MOM));
    assertEquals(Long.valueOf(10003L), vbuf.getEndPosition(DAD));
  }

  @Test
  public void testToString() {
    Variant v = new Variant().setPosition(1L).setEnd(1000L);
    vbuf.push(DAD, v);
    vbuf.push(DAD, v);
    vbuf.push(DAD, v);
    vbuf.push(MOM, v);
    vbuf.push(MOM, new Variant().setPosition(3L).setEnd(5000L));
    assertEquals("DAD:[1-1000,1-1000,1-1000], MOM:[1-1000,3-5000], CHILD:[]", vbuf.toString());
  }
}
