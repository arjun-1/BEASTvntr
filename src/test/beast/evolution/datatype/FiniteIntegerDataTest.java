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

		int[] statesTest;
		statesTest = new int[]{0,1,2,3,4,5,6,7,8,9,10,11,12};

		assertArrayEquals(statesTest, datatype.getStatesForCode(-2));
		assertArrayEquals(statesTest, datatype.getStatesForCode(-1));

		statesTest = new int[]{0};
		assertArrayEquals(statesTest, datatype.getStatesForCode(0));

		statesTest = new int[]{1};
		assertArrayEquals(statesTest, datatype.getStatesForCode(1));

		statesTest = new int[]{10};
		assertArrayEquals(statesTest, datatype.getStatesForCode(10));

		statesTest = new int[]{12};
		assertArrayEquals(statesTest, datatype.getStatesForCode(12));

		Sequence sequence = new Sequence();
		sequence.init(13, "myTaxon", "5,6,16,17");
		sequence.setID("mySequence");

		Alignment m_alignment = new Alignment();
		m_alignment.sequenceInput.setValue(sequence, m_alignment);
		m_alignment.setID("myAlignment");
		m_alignment.setInputValue("userDataType", datatype);
		m_alignment.initAndValidate();

		assertEquals(0, m_alignment.getPattern(0, 0));
		assertEquals(1, m_alignment.getPattern(0, 1));
		assertEquals(11, m_alignment.getPattern(0, 2));
		assertEquals(12, m_alignment.getPattern(0, 3));
	}
}
