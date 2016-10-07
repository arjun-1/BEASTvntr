package test.beast.evolution.datatype;

import static org.junit.Assert.*;
import org.junit.Test;

import beast.evolution.datatype.FiniteIntegerData;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.Alignment;
import junit.framework.TestCase;;

public class FiniteIntegerDataTest extends TestCase {
  @Test
  public void testFiniteIntegerData() {
    FiniteIntegerData datatype = new FiniteIntegerData();
    datatype.setInputValue("minRepeat", 5);
    datatype.setInputValue("maxRepeat", 17);
    datatype.initAndValidate();

    assertEquals("?", datatype.getCode(-2));
    assertEquals("?", datatype.getCode(-1));
    assertEquals("0", datatype.getCode(0));
    assertEquals("1", datatype.getCode(1));
    assertEquals("10", datatype.getCode(10));
    assertEquals("11", datatype.getCode(11));
    assertEquals("12", datatype.getCode(12));

    assertEquals('?', datatype.getChar(-2));
    assertEquals('?', datatype.getChar(-1));
    assertEquals('0', datatype.getChar(0));
    assertEquals('1', datatype.getChar(1));
    assertEquals('1' , datatype.getChar(1));
    assertEquals('0' + 10, datatype.getChar(10));
    assertEquals('0' + 11, datatype.getChar(11));
    assertEquals('0' + 12, datatype.getChar(12));

    final int[] expectedAll = new int[]{0,1,2,3,4,5,6,7,8,9,10,11,12};
    final int[] expected0 = new int[]{0};
    final int[] expected1 = new int[]{1};
    final int[] expected10 = new int[]{10};
    final int[] expected12 = new int[]{12};

    assertArrayEquals(expectedAll, datatype.getStatesForCode(-2));
    assertArrayEquals(expectedAll, datatype.getStatesForCode(-1));

    assertArrayEquals(expected0, datatype.getStatesForCode(0));
    assertArrayEquals(expected1, datatype.getStatesForCode(1));
    assertArrayEquals(expected10, datatype.getStatesForCode(10));
    assertArrayEquals(expected12, datatype.getStatesForCode(12));

    Sequence sequence = new Sequence();
    sequence.init(13, "myTaxon", "5,6,16,17");
    sequence.setID("mySequence");

    Alignment alignment = new Alignment();
    alignment.sequenceInput.setValue(sequence, alignment);
    alignment.setID("myAlignment");
    alignment.setInputValue("userDataType", datatype);
    alignment.initAndValidate();

    assertEquals(0, alignment.getPattern(0, 0));
    assertEquals(1, alignment.getPattern(0, 1));
    assertEquals(11, alignment.getPattern(0, 2));
    assertEquals(12, alignment.getPattern(0, 3));
  }
}
