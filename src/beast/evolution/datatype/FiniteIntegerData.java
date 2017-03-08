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

package beast.evolution.datatype;

import beast.core.Description;
import beast.evolution.datatype.DataType.Base;
import beast.core.Input;

import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;

@Description("Datatype for finite integer sequences")
public class FiniteIntegerData extends Base {

  final public Input<Integer> maxRepeatInput = new Input<>("maxRepeat", "specifies the highest state");
  final public Input<Integer> minRepeatInput = new Input<>("minRepeat", "specifies the lowest state");

  private Map<String, Integer> mapStringToCode;
  private Map<Integer, String> mapCodeToString;
  private Map<Integer, Character> mapCodeToChar;

  @Override
  public void initAndValidate() {
    if (maxRepeatInput.get() - minRepeatInput.get() >= 0
            && minRepeatInput.get() >= 0) {
      stateCount = maxRepeatInput.get() - minRepeatInput.get() + 1;
    } else {
      throw new IllegalArgumentException("Invalid values for maximum repeat: " + maxRepeatInput.get() + " or minimum repeat: " + minRepeatInput.get());
    }
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
    // Map the ambiguous character to all states
    int [] stateSetFull = new int[stateCount];
    for (int i = 0; i < stateCount; i++) {
      stateSetFull[i] = i;
    }
    mapCodeToStateSet[stateCount] = stateSetFull;
  }

  @Override
  public List<Integer> string2state(String data) {
    List<Integer> sequence;
    sequence = new ArrayList<>();
    String[] strs = data.split(",");
    for (String str : strs) {
      if (mapStringToCode.containsKey(str)) {
          sequence.add(mapStringToCode.get(str));
      } else {
        throw new IllegalArgumentException("Not a legal state: " + str);
      }
    }
    return sequence;
  }

  @Override
  public String state2string(int[] states) {
    StringBuilder sb = new StringBuilder();
    for (int i = 0; i < states.length - 1; i++) {
      if (mapCodeToString.containsKey(states[i])) {
        sb.append(mapCodeToString.get(states[i]));
        sb.append(",");
      } else {
        throw new IllegalArgumentException("Not a legal state: " + states[i]);
      }
    }
    if (mapCodeToString.containsKey(states[states.length - 1])) {
      sb.append(mapCodeToString.get(states[states.length - 1]));
    } else {
      throw new IllegalArgumentException("Not a legal state: " + states[states.length - 1]);
    }
    return sb.toString();
  }

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

  /**
   * This implementation of getCode does not refer to the 'code' in mapCodeToStateSet,
   * but is defined as the inverse of string2state for a single state. This is expected by FilteredAlignment.
   */
  @Override
  public String getCode(int state) {
    return isAmbiguousState(state) ? String.valueOf(getChar(state)) : Integer.toString(getChar(state) - '0');
  }

  @Override
  public String getTypeDescription() {
    return "finiteinteger";
  }
}