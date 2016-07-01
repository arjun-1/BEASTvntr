package test.beast.evolution.datatype;

import static org.junit.Assert.*;
import org.junit.Test;

import beast.evolution.datatype.FiniteIntegerData;
import junit.framework.TestCase;;

public class FiniteIntegerDataTest extends TestCase {

	
	@Test
	public void testFiniteIntegerData() {
		//IntegerData datatype = new IntegerData();
		FiniteIntegerData datatype = new FiniteIntegerData();
		datatype.setInputValue("minRepeat", 5);
		datatype.setInputValue("maxRepeat", 17);
		datatype.initAndValidate();

		assertEquals("?", datatype.getCode(-2));
		assertEquals("?", datatype.getCode(-1));
		assertEquals("5", datatype.getCode(0));
		assertEquals("6", datatype.getCode(1));
		assertEquals("15", datatype.getCode(10));

		assertEquals('?', datatype.getChar(-2));
		assertEquals('?', datatype.getChar(-1));
		assertEquals('5', datatype.getChar(0));
		assertEquals('6', datatype.getChar(1));
		assertEquals('0' + 6, datatype.getChar(1));
		assertEquals('0' + 15, datatype.getChar(10));

		int[] statesTest;
		statesTest = new int[]{0,1,2,3,4,5,6,7,8,9,10,11,12};

		assertArrayEquals(statesTest, datatype.getStatesForCode(-2));
		assertArrayEquals(statesTest, datatype.getStatesForCode(-1));

		statesTest = new int[]{0};
		assertArrayEquals(statesTest, datatype.getStatesForCode(5));

		statesTest = new int[]{1};
		assertArrayEquals(statesTest, datatype.getStatesForCode(6));

		statesTest = new int[]{10};
		assertArrayEquals(statesTest, datatype.getStatesForCode(15));

		statesTest = new int[]{12};
		assertArrayEquals(statesTest, datatype.getStatesForCode(17));


	}
}
