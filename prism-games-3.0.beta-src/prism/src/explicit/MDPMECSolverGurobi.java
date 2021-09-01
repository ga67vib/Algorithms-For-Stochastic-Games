package explicit;

import gurobi.*;
import prism.Pair;
import prism.PrismException;

import java.util.*;

public class MDPMECSolverGurobi {
    protected List<BitSet> mecs;
    protected MDP mdp;
    BitSet known;
    double[] knownValues;
    boolean min;
    private GRBModel lp;
    private GRBVar[] stateVars;
    private HashMap<Integer, Integer> map;

    public MDPMECSolverGurobi(GRBModel lp, GRBVar[] stateVars, List<BitSet> mecs, MDP mdp, boolean min, BitSet known, double[] knownValues) {
        this.mecs=mecs;
        this.mdp=mdp;
        this.known=known;
        this.knownValues=knownValues;
        this.min = min;
        this.lp=lp;
        this.stateVars=stateVars;
        this.map = null;
    }

    public void solve() throws PrismException, GRBException{
        for (BitSet mec : mecs) {
            createStrategyTable(mec);
        }
    }

    protected void createStrategyTable(BitSet mec) throws PrismException, GRBException {
        if (this.min) {
            setMECtoZero(mec);
        }
        else {
            BitSet exits = getMECExits(mec);
            setMECtoBestExit(mec, exits);
        }
    }

    protected void setMECtoZero(BitSet mec) throws GRBException {
        int numECStates=mec.cardinality();
        GRBLinExpr[] lhs = new GRBLinExpr[numECStates];
        double[] rhs = new double[numECStates];
        char[] senses = new char[numECStates];
        Arrays.fill(senses, GRB.EQUAL);
        int state=mec.nextSetBit(0);
        for (int i=0; i<numECStates; i++) {
            lhs[i] = new GRBLinExpr();
            lhs[i].addTerm(1.0, stateVars[state]);
            state=mec.nextSetBit(state+1);
        }
        lp.addConstrs(lhs, senses, rhs, null);
    }

    protected void setMECtoBestExit(BitSet mec, BitSet exits) throws GRBException {
        ArrayList<Pair<Integer, Integer>> exitStateActions = getExitStateActionPairs(mec, exits);
        int sizeExitActions=exitStateActions.size();
        if (sizeExitActions==0) {
            setMECtoZero(mec);
            return;
        }
        GRBVar[] extraVars = createGRBVars(sizeExitActions);

        GRBLinExpr[] lhs = new GRBLinExpr[sizeExitActions];
        double[] rhs = new double[sizeExitActions];
        char[] senses = new char[sizeExitActions];
        Arrays.fill(senses, GRB.EQUAL);

        for (int i=0; i<sizeExitActions; i++) {
            Pair<Integer, Integer> pair = exitStateActions.get(i);
            lhs[i]=new GRBLinExpr();
            lhs[i].addTerm(1.0, extraVars[i]);

            // if we have P(a,b,c) but c leads into EC and a, b are exits then the best exit is
            // then exit probability is a/(a+b) and b/(a+b)
            double divisor=0;
            boolean allExit=true;
            for (Iterator<Map.Entry<Integer, Double>> it = mdp.getTransitionsIterator(pair.getKey(), pair.getValue()); it.hasNext(); ) {
                Map.Entry<Integer, Double> tr = it.next();
                if (exits.get(tr.getKey())) {
                    divisor+=tr.getValue();
                }
                else {
                    allExit=false;
                }
            }
            if (allExit) divisor=1.0; //just to increase precision

            for (Iterator<Map.Entry<Integer, Double>> it = mdp.getTransitionsIterator(pair.getKey(), pair.getValue()); it.hasNext(); ) {
                Map.Entry<Integer, Double> tr = it.next();
                if (exits.get(tr.getKey())) { //Add every leaving Probability
                    if (known.get(tr.getKey())) {
                        lhs[i].addConstant(-1.0*tr.getValue()/divisor * knownValues[tr.getKey()]);
                    }
                    else {
                        lhs[i].addTerm(-1.0*tr.getValue()/divisor, stateVars[tr.getKey()]);
                    }
                }
            }

        }

        lp.addConstrs(lhs, senses, rhs, null);

        GRBVar[] vars=new GRBVar[mec.cardinality()];
        int i=0;
        for (int state = mec.nextSetBit(0); state>=0; state=mec.nextSetBit(state+1)) {
            if (state == Integer.MAX_VALUE) {
                break;
            }
            //lp.addGenConstrMax(stateVars[state], extraVars, 0.0, ""); //Creates Max constraint
            vars[i]=stateVars[state];
            i++;
        }
        createMinMaxConstraint(vars, extraVars, true);
    }

    private void createMinMaxConstraint(GRBVar[] baseVars, GRBVar[] vars, boolean max) throws GRBException{
        GRBVar[] sVars = createGRBVars(vars.length);
        double[] lb=new double[sVars.length];
        double[] ub=new double[sVars.length];
        char[] type=new char[sVars.length];
        Arrays.fill(type, GRB.BINARY);
        Arrays.fill(ub, 1.0);
        GRBVar[] zVars = lp.addVars(lb, ub,null, type, null);

        GRBLinExpr[] lhs=new GRBLinExpr[vars.length];
        double[] rhs=new double[vars.length];
        char[] senses=new char[vars.length];
        Arrays.fill(senses, GRB.EQUAL);

        for (int k=0; k<baseVars.length; k++) {
            for (int i = 0; i < vars.length; i++) {
                lhs[i] = new GRBLinExpr();
                lhs[i].addTerm(1.0, baseVars[k]);
                lhs[i].addTerm(-1.0, vars[i]);
                if (max) {
                    lhs[i].addTerm(-1.0, sVars[i]);
                } else {
                    lhs[i].addTerm(1.0, sVars[i]);
                }
            }
            lp.addConstrs(lhs, senses, rhs, null);
        }

        //only one z is set
        GRBLinExpr oneZ = new GRBLinExpr();
        for (int i=0; i<vars.length; i++) {
            oneZ.addTerm(1.0, zVars[i]);
        }
        lp.addConstr(oneZ, GRB.EQUAL, 1.0, null);
        //s+z<=1 only if z is set to a var where s is 0 this is correct for all
        //s+z-1<=0
        //Is z if set really 1? Check in Documentation
        lhs=new GRBLinExpr[vars.length];
        rhs=new double[vars.length];
        senses=new char[vars.length];
        Arrays.fill(senses, GRB.LESS_EQUAL);
        for (int i=0; i<vars.length; i++) {
            lhs[i]=new GRBLinExpr();
            lhs[i].addTerm(1.0, sVars[i]);
            lhs[i].addTerm(1.0, zVars[i]);
            lhs[i].addConstant(-1.0);
        }
        lp.addConstrs(lhs, senses, rhs, null);
    }

    protected ArrayList<Pair<Integer, Integer>> getExitStateActionPairs(BitSet mec, BitSet exits) {
        ArrayList<Pair<Integer, Integer>> res=new ArrayList<>();
        for (int state = mec.nextSetBit(0); state>=0; state=mec.nextSetBit(state+1)) {
            if (state == Integer.MAX_VALUE) {
                break;
            }
            int numChoices = mdp.getNumChoices(state);
            for (int choice=0; choice<numChoices; choice++) {
                for (Iterator<Map.Entry<Integer, Double>> it = mdp.getTransitionsIterator(state, choice); it.hasNext(); ) {
                    Map.Entry<Integer, Double> tr = it.next();
                    if (exits.get(tr.getKey())) {
                        res.add(new Pair<>(state, choice));
                        break;
                    }
                }
            }

        }

        return res;
    }

    private GRBVar[] createGRBVars(int size) throws  GRBException {
        double[] ub = new double[size];
        double[] lb = new double[size];
        Arrays.fill(ub, 1.0);

        return lp.addVars(lb, ub, null, null, null);
    }

    /**
     * Gets all the Exits the MEC might have
     */
    protected BitSet getMECExits(BitSet mec) {
        BitSet exits = new BitSet();
        for (int state = mec.nextSetBit(0); state>=0; state=mec.nextSetBit(state+1)) {
            if (state==Integer.MAX_VALUE) {
                break;
            }

            int choiceNum = mdp.getNumChoices(state);
            for (int choice = 0; choice < choiceNum; choice++) {
                for (Iterator<Map.Entry<Integer, Double>> it = mdp.getTransitionsIterator(state, choice); it.hasNext(); ) {
                    Map.Entry<Integer, Double> tr = it.next();
                    if (!mec.get(tr.getKey()) && !exits.get(tr.getKey())) {
                        exits.set(tr.getKey());
                    }
                }
            }
        }
        return exits;
    }
}
