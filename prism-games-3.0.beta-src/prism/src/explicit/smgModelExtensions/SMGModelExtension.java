package explicit.smgModelExtensions;

import explicit.Distribution;
import explicit.STPG;
import explicit.STPGExplicit;
import parser.State;

import java.util.BitSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

public abstract class SMGModelExtension {

    List<State> statesList;
    STPGExplicit stpg;
    BitSet remain;
    boolean useExtension;

    /**
     * Constructor for SMGModelExtensions
     * @param stpg STPG to extend
     * @param remain
     * @param statesList Is necessary, because stateList from STPG is protected in the package
     * @param use Should this extension be used
     */
    public SMGModelExtension(STPGExplicit stpg, BitSet remain, List<State> statesList, boolean use) {
        this.statesList = statesList;
        this.stpg = stpg;
        this.remain = remain;
    }

    protected Distribution getDistribution(STPG stpg, int state, int action) {
        Distribution d = new Distribution();
        for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, action); it.hasNext(); ) {
            Map.Entry<Integer, Double> tr = it.next();
            d.add(tr.getKey(), tr.getValue());
        }
        return d;
    }

    public boolean shouldUseExtension() {
        return useExtension;
    }

    protected void addToStateList(STPGExplicit stpg, int statesFrom, int statesTo) {
        for (int i = statesFrom; i<=statesTo; i++) {
            State s = new State(1);
            s.setValue(0, (Integer) i);
            this.statesList.add(s);
        }
    }

    /**
     *
     * @param oldInitialState the initial state before this extension
     * @return
     */
    public abstract ModelExtensionResult extendSMG(int oldInitialState);
}
