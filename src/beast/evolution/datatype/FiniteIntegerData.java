package beast.evolution.datatype;

import beast.core.Description;
import beast.evolution.datatype.DataType.Base;
import beast.core.Input;
import java.util.List;

import java.util.ArrayList;
import java.util.Collections;


// (Almost) copied from StandardModel.java
@Description("Datatype for finite integer sequences")
public class FiniteIntegerData extends Base {

  final public Input<Integer> maxRepeatInput = new Input<>("maxRepeat", "specifies the highest state");
  final public Input<Integer> minRepeatInput = new Input<>("minRepeat", "specifies the lowest state");

  private ArrayList<String> codeMapping;

  @Override
  public void initAndValidate() {
    if (maxRepeatInput.get() - minRepeatInput.get() >= 0 &&
      minRepeatInput.get() >= 0) {
      stateCount = maxRepeatInput.get() - minRepeatInput.get() + 1;
    } else {
      throw new IllegalArgumentException("Bad values for maximum repeat: " + maxRepeatInput.get() + ", minimum repeat: " +minRepeatInput.get());
    }

    mapCodeToStateSet = null;
    codeLength = -1;
    codeMap = null;
    createCodeMapping();
  }

  private void createCodeMapping() {
    codeMapping = new ArrayList<>();
    for (int i=0; i<stateCount; i++) {
      codeMapping.add(Integer.toString(i));
    }
    codeMapping.add(Character.toString(GAP_CHAR));
    codeMapping.add(Character.toString(MISSING_CHAR));

    mapCodeToStateSet = new int[codeMapping.size()][];
    for (int i = 0; i < codeMapping.size() - 2; i++) {
      int [] stateSet = new int[1];
      stateSet[0] = Integer.parseInt(codeMapping.get(i));
      mapCodeToStateSet[i] = stateSet;
    }
    
    // TODO: is this the correct way to deal with stateCount == -1?
    int n = stateCount >= 0 ? stateCount : 10;
    int [] stateSet = new int[n];
    for (int i = 0; i < n; i++) {
      stateSet[i] = i;
    }
    // GAP_CHAR
    mapCodeToStateSet[mapCodeToStateSet.length - 2] = stateSet;
    // MISSING_CHAR
    mapCodeToStateSet[mapCodeToStateSet.length - 1] = stateSet;
  }
  
  @Override
  public int[] getStatesForCode(int state) {
    if (state >= 0) {
      return mapCodeToStateSet[state];
    } else {
      return mapCodeToStateSet[mapCodeToStateSet.length - 1];
    }
  }

  @Override
  public String getTypeDescription() {
    return "finiteinteger";
  }

  @Override
  public char getChar(int state) {
    if (state < 0) {
      return '?';
    }
    return (char)('0'+state);
  }

  @Override
  public String getCode(int state) {
    if (state < 0) {
      return "?";
    }
    return codeMapping.get(state);
  }

  @Override
  public List<Integer> string2state(String data) {
    List<Integer> sequence;
    sequence = new ArrayList<>();
    // remove spaces
    data = data.replaceAll("\\s", "");
    data = data.toUpperCase();

    // assume it is a comma separated string of integers
    String[] strs = data.split(",");
    for (String str : strs) {
      try {                
        if (Integer.parseInt(str)  >= 0 && Integer.parseInt(str)  - minRepeatInput.get() < 0) {
          throw new IllegalArgumentException("Encountered repeat out of model bounds: " + Integer.parseInt(str) + " < " +minRepeatInput.get()); 
        } else if (Integer.parseInt(str)  - maxRepeatInput.get() > 0) {
          throw new IllegalArgumentException("Encountered repeat out of model bounds: " + Integer.parseInt(str) + " > " + maxRepeatInput.get());  
        } else {
          sequence.add(Integer.parseInt(str) - minRepeatInput.get());
        }
      } catch (NumberFormatException e) {
        sequence.add(-1);
      }
    }
    return sequence;
  } // string2state

}