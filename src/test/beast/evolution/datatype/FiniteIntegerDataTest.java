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
    final int minRepeat = 5;
    final int maxRepeat = 17;
    final int stateCount = maxRepeat - minRepeat + 1;

    datatype.setInputValue("minRepeat", 5);
    datatype.setInputValue("maxRepeat", 17);
    datatype.initAndValidate();

    assertEquals("?", datatype.getCode(stateCount));
    assertEquals("5", datatype.getCode(0));
    assertEquals("9", datatype.getCode(4));
    assertEquals("10", datatype.getCode(5));
    assertEquals("17", datatype.getCode(12));

    assertEquals('?', datatype.getChar(13));
    assertEquals('5', datatype.getChar(0));
    assertEquals('9', datatype.getChar(4));
    assertEquals('0' + minRepeat + 10, datatype.getChar(10));
    assertEquals('0' + minRepeat + 11, datatype.getChar(11));
    assertEquals('0' + minRepeat + 12, datatype.getChar(12));

    final int[] expectedAll = new int[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    final int[] expected0   = new int[] {0};
    final int[] expected1   = new int[] {1};
    final int[] expected10  = new int[] {10};
    final int[] expected12  = new int[] {12};

    assertArrayEquals(expectedAll, datatype.getStatesForCode(stateCount));

    assertArrayEquals(expected0, datatype.getStatesForCode(0));
    assertArrayEquals(expected1, datatype.getStatesForCode(1));
    assertArrayEquals(expected10, datatype.getStatesForCode(10));
    assertArrayEquals(expected12, datatype.getStatesForCode(12));
  }
}
