/*
* File FiniteIntegerData.java
*
* Copyright (C) 2017 Arjun Dhawan
*
* This file is part of BEASTvntr.
*
* BEASTvntr is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* BEASTvntr is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with BEASTvntr.  If not, see <http://www.gnu.org/licenses/>.
*/

package vntr.evolution.datatype;

import beast.base.core.Description;
import beast.base.evolution.datatype.DataType.Base;
import beast.base.core.Input;

import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
//import java.util.stream.Collectors;
//import java.util.stream.Stream;

@Description("Datatype for finite integer sequences")
public class FiniteIntegerData extends Base {

  final public Input<Integer> maxRepeatInput = new Input<>("maxRepeat", "specifies the highest state");
  final public Input<Integer> minRepeatInput = new Input<>("minRepeat", "specifies the lowest state");

  private Map<String, Integer> mapStringToCode;
  private Map<Integer, String> mapCodeToString;
  private Map<Integer, Character> mapCodeToChar;
  private int [] missing; // missing data state representation

  @Override
  public void initAndValidate() {
    if (maxRepeatInput.get() - minRepeatInput.get() < 0 || minRepeatInput.get() < 0) {
      throw new IllegalArgumentException("Invalid values for maximum repeat: " + maxRepeatInput.get()
              + " or minimum repeat: " + minRepeatInput.get());
    }
    stateCount = maxRepeatInput.get() - minRepeatInput.get() + 1;
    createCodeMapping();
  }

  private void createCodeMapping() {
    // Map the external states (data): {minRepeat, ..., maxRepeat}
    // to the internal states (codes): {0, ..., stateCount - 1}
    mapStringToCode = new HashMap<>();
    mapCodeToString = new HashMap<>();
    mapCodeToChar = new HashMap<>();
    for (int i = 0; i < stateCount; i++) {
      mapStringToCode.put(Integer.toString(i + minRepeatInput.get()), i);
      mapCodeToString.put(i, Integer.toString(i + minRepeatInput.get()));
      mapCodeToChar.put(i, (char)('0' + i + minRepeatInput.get()));
    }
    // Map the ambiguous characters
    mapStringToCode.put(Character.toString(GAP_CHAR), stateCount);
    mapStringToCode.put(Character.toString(MISSING_CHAR), stateCount);
    mapStringToCode.put("-1", stateCount); // backwards compatibility with BEASTvntr 0.1.0
    mapCodeToString.put(stateCount, Character.toString(MISSING_CHAR));
    mapCodeToChar.put(stateCount, MISSING_CHAR);

    // Map the codes to state sets
    mapCodeToStateSet = new int[mapCodeToString.size()][];
    for (int i = 0; i < stateCount; i++) {
      mapCodeToStateSet[i] = new int[] { i };
    }
    missing = new int[stateCount];
    for (int i = 0; i < stateCount; i++) {
    	missing[i] = i;
    }
    // Map the ambiguous character to all states
    int [] stateSetFull = new int[stateCount];
    for (int i = 0; i < stateCount; i++) {
      stateSetFull[i] = i;
    }
    mapCodeToStateSet[stateCount] = stateSetFull;
    
    codeLength = -1;
  }

  /**
   * Maps states from {minRepeat, ..., maxRepeat} -> {0, ..., stateCount - 1}
   */
  @Override
  public List<Integer> string2state(String data) {
    // return Stream.of(data.split(",")).map(this::getState).collect(Collectors.toList());
    List<Integer> states = new ArrayList<>();
    String[] dataSplit = data.split(",");
    for (String str : dataSplit) {
      states.add(getState(str));
    }
    return states;
  }

  /**
   * Maps states from {0, ..., stateCount - 1} -> {minRepeat, ..., maxRepeat}
   */
  @Override
  public String state2string(int[] states) {
    // return Arrays.stream(states).mapToObj(this::getCode).collect(Collectors.joining(","));
    String[] strings = new String[states.length];
    for (int i = 0; i < states.length; i++) {
      strings[i] = getCode(states[i]);
    }
    return String.join(",", strings);
  }

  /**
   * Maps states from {minRepeat, ..., maxRepeat} -> {0, ..., stateCount - 1}
   */
  private int getState(String str) {
    if (!mapStringToCode.containsKey(str)) {
      throw new IllegalArgumentException("Not a legal state: " + str);
    }
    return mapStringToCode.get(str);
  }

  
	@Override
	public List<Integer> stringToEncoding(String data) throws IllegalArgumentException {
		List<Integer> sequence = new ArrayList<>();
		String [] strs = data.split(",");
		for (String str : strs) {
			sequence.add(Integer.parseInt(str.trim()));
		}
		return sequence;
	}
	
  /**
   * This implementation of getCode is defined such that string2state is its inverse,
   * as expected by FilteredAlignment.
   * The 'code' in the name 'getCode' thus does not refer to the 'code' in mapCodeToStateSet or mapCodeToString.
   *
   * Maps states from {0, ..., stateCount - 1} -> {minRepeat, ..., maxRepeat}
   */
  @Override
  public String getCode(int state) {
    if (!mapCodeToString.containsKey(state)) {
      throw new IllegalArgumentException("Not a legal state: " + state);
    }
    return mapCodeToString.get(state);
  }
  
  @Override
  public int[] getStatesForCode(int code) {
		if (code >= 0) {
			return super.getStatesForCode(code);
		}
		return missing; 
  }

  /**
   * Maps states from {0, ..., stateCount - 1} -> char form of {minRepeat, ..., maxRepeat}
   */
  @Override
  public char getChar(int code) {
    if (!mapCodeToChar.containsKey(code)) {
      throw new IllegalArgumentException("Not a legal state: " + code);
    }
    return mapCodeToChar.get(code);
  }

  @Override
  public boolean isAmbiguousState(int state) {
    return (state >= stateCount);
  }

  @Override
  public String getTypeDescription() {
    return "finiteinteger";
  }
  
  // @Override -- it does override method introduced in BEAST v2.6.3
  public boolean hasConstantCodeLength() {
  	return false;
  }
}