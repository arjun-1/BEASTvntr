package beast.evolution.datatype;

import beast.core.Description;
import beast.evolution.datatype.DataType.Base;
import beast.core.Input;
import java.util.List;

import java.util.ArrayList;
import java.util.Collections;



@Description("Datatype for finite integer sequences")
public class FiniteIntegerData extends Base {

    final public Input<Integer> max_iInput = new Input<>("max_i", "specifies the highest state");
    final public Input<Integer> min_iInput = new Input<>("min_i", "specifies the lowest state");

    private ArrayList<String> codeMapping;

    @Override
    public void initAndValidate() {
        if (max_iInput.get() != null && max_iInput.get() != 0 &&
            min_iInput.get() != null) {
            stateCount = max_iInput.get() - min_iInput.get() + 1;
        } else {
            stateCount = -1;
        }

        mapCodeToStateSet = null;
        codeLength = -1;
        codeMap = null;
        createCodeMapping();
    }

    private void createCodeMapping() {
        codeMapping = new ArrayList<>();
        for (int i=0; i<stateCount; i++) {
            codeMapping.add(Integer.toString(i + min_iInput.get()));
        }
        codeMapping.add(Character.toString(GAP_CHAR));
        codeMapping.add(Character.toString(MISSING_CHAR));

        mapCodeToStateSet = new int[codeMapping.size()][];
        for (int i = 0; i < codeMapping.size() - 2; i++) {
            int [] stateSet = new int[codeMapping.get(i).length()];
            for (int k = 0; k < stateSet.length; k++) {
                stateSet[k] = (codeMapping.get(i).charAt(k) - '0');
            }
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
    public int[] getStatesForCode(int code) {
        if (code >= 0) {
            return mapCodeToStateSet[code - min_iInput.get()];
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
        return (char)('0' + state + min_iInput.get());
    }
    
    @Override
    public String getCode(int state) {
        return codeMapping.get(state);
    }

}