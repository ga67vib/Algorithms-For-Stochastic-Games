package explicit;

//import com.sun.xml.internal.bind.v2.runtime.BinderImpl;
import gurobi.*;
import prism.Pair;
import prism.PrismException;
import Jama.*;

import java.lang.reflect.Array;
import java.util.*;

public class SMGMECSolverGRB extends SMGMECSolver {

    private GRBModel qp;
    private GRBVar[] stateVars;
    private HashMap<Integer, Integer> map;

    public SMGMECSolverGRB(GRBModel qp, GRBVar[] stateVars, List<BitSet> mecs, STPG stpg, BitSet known, double[] knownValues) {
        super(mecs, stpg, known, knownValues);
        this.qp=qp;
        this.stateVars=stateVars;
        this.map = null;
    }

    @Override
    public void solve() throws PrismException, GRBException{
        for (BitSet mec : mecs) {
            createStrategyTable(mec);
        }
    }

    public void solve(BitSet mec, HashMap<Integer, Integer> map) throws PrismException, GRBException {
        this.map = map;
        createStrategyTable(mec);
    }

    @Override
    protected void handleMixedMEC(BitSet mec, BitSet exits) throws GRBException, PrismException{

        //pairs of (state, action)
        //Now the second Integer can be substituted by only booleans because every state can only
        //have 2 Actions. -> There is a nicer solution which represents the fixed strategy by one
        //binary number
        //Changed to 4 ArrayLists instead of 2 of Pairs for "indexOf"-Function in solveMarkovChain
        ArrayList<Integer> maxStatesNames = new ArrayList<>();
        ArrayList<Integer> maxStatesActions = new ArrayList<>();
        ArrayList<Integer> minStatesNames = new ArrayList<>();
        ArrayList<Integer> minStatesActions = new ArrayList<>();
        int maxActions = 1;
        int minActions = 1;

        for (int state = mec.nextSetBit(0); state>=0; state=mec.nextSetBit(state+1)) {
            if (state==Integer.MAX_VALUE) {
                break;
            }

            if (stpg.getNumChoices(state) > 1) { //States that don't have choices are taking the same action always anyways
                if (stpg.getPlayer(state) == 1) {
                    maxStatesNames.add(state);
                    maxStatesActions.add(0);
                    maxActions = maxActions * stpg.getNumChoices(state);
                } else {
                    minStatesNames.add(state);
                    minStatesActions.add(0);
                    minActions = minActions * stpg.getNumChoices(state);
                }
            }
            if (maxActions<0) {
                reportTooManyActions(maxStatesNames, true);
            }
            else if (minActions<0) {
                reportTooManyActions(maxStatesNames, false);
            }
        }

        int numMECStates = maxStatesNames.size() + minStatesNames.size();
        //Gurobi and LEqS-Variables
        GRBVar[][] minConstrVars=new GRBVar[maxActions][numMECStates]; //For each maxAction we have a set of minConstrVariables -> In the whole table there
        GRBVar[][] equalVars=new GRBVar[minActions][numMECStates]; // is one maxAction fixed with all equal variables: minActions x numMECStates


        //Shouldn't be String a Set of Linear expressions?
        for (int x = 0; x < maxActions; x++) {
            minConstrVars[x]=createGRBVars(numMECStates);

            for (int i=0; i<minActions; i++) {
                equalVars[i] = createGRBVars(numMECStates);
            }


            for (int y = 0; y < minActions; y++) {
                Matrix exitProbs = solveMarkovChain(maxStatesNames, maxStatesActions, minStatesNames, minStatesActions, exits);
                getEqualConstraints(exitProbs, maxStatesNames.size()+minStatesNames.size(), exits, equalVars[y]);
                incrementAction(minStatesActions, minStatesNames);
            }
            getMinConstraints(minConstrVars[x], equalVars);
            incrementAction(maxStatesActions, maxStatesNames);
        }
        getMaxConstraints(minConstrVars, maxStatesNames, minStatesNames);


    }

    private void getEqualConstraints(Matrix probs, int numECStates, BitSet exits, GRBVar[] equalVars) throws PrismException, GRBException {
        GRBLinExpr[] lhs = new GRBLinExpr[numECStates];
        double[] rhs = new double[numECStates];
        char[] senses = new char[numECStates];
        Arrays.fill(senses, GRB.EQUAL); //minimize over all extraVariables

        for (int i = 0; i < numECStates; i++) {

            lhs[i] = new GRBLinExpr();
            lhs[i].addTerm(1.0, equalVars[i]);

            int exit=exits.nextSetBit(0);
            for (int j = 0; j < probs.getColumnDimension(); j++) {
                double prob = probs.get(i, j);
                if (prob != 0.0) {
                    if (known.get(exit)) {
                        lhs[i].addConstant(-1.0 * prob * knownValues[exit]);
                    }
                    else {
                        lhs[i].addTerm(-1.0 * prob, getStateVar(exit));//stateVars[exit]);//ToDo: Change here accordingly
                    }
                }
                exit=exits.nextSetBit(exit+1);
            }

        }
        qp.addConstrs(lhs, senses, rhs, null);
    }

    /**
     * Takes the Constraints from Probability-Matrix and transforms them into the Min-Constraints
     */
    private void getMinConstraints(GRBVar[] minConstrVars, GRBVar[][] equalVars)
            throws PrismException, GRBException {

        GRBVar[][] equalVarsTranspose = transpose(equalVars);

        for (int i=0; i<minConstrVars.length; i++) {
            //qp.addGenConstrMin(minConstrVars[i], equalVarsTranspose[i], 1.0, "");
            createMinMaxConstraint(minConstrVars[i], equalVarsTranspose[i], false);
        }
    }

    private void getMaxConstraints (GRBVar[][] minConstrVars, ArrayList<Integer> maxStatesNames, ArrayList<Integer> minStatesNames)
            throws PrismException, GRBException {
        GRBVar[][] minConstrVarsTranspose = transpose(minConstrVars);

        int numECStates=maxStatesNames.size()+minStatesNames.size();

        ArrayList<Integer> states;
        int listIndex;

        for (int run=0; run<2; run++) {
            listIndex=0;
            if (run==0) {
                states=maxStatesNames;
            }
            else {
                states=minStatesNames;
            }

            //This loop runs first time vom 0 -> maxStates.size -1 , second time maxStates.size -> numECStates-1
            for (int i=maxStatesNames.size()*run; i<maxStatesNames.size()+minStatesNames.size()*run; i++) {
                //qp.addGenConstrMax(stateVars[states.get(listIndex)], minConstrVarsTranspose[i], 0.0, "");

                createMinMaxConstraint(getStateVar(states.get(listIndex)), minConstrVarsTranspose[i], true);
                listIndex++;
            }
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
        qp.addConstrs(lhs, senses, rhs, null);
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
            for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(pair.getKey(), pair.getValue()); it.hasNext(); ) {
                Map.Entry<Integer, Double> tr = it.next();
                if (exits.get(tr.getKey())) {
                    divisor+=tr.getValue();
                }
                else {
                    allExit=false;
                }
            }
            if (allExit) divisor=1.0; //just to increase precision

            for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(pair.getKey(), pair.getValue()); it.hasNext(); ) {
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

        qp.addConstrs(lhs, senses, rhs, null);

        GRBVar[] vars=new GRBVar[mec.cardinality()];
        int i=0;
        for (int state = mec.nextSetBit(0); state>=0; state=mec.nextSetBit(state+1)) {
            if (state == Integer.MAX_VALUE) {
                break;
            }
            //qp.addGenConstrMax(stateVars[state], extraVars, 0.0, ""); //Creates Max constraint
            vars[i]=stateVars[state];
            i++;
        }
        createMinMaxConstraint(vars, extraVars, true);
    }

    private GRBVar[] createGRBVars(int size) throws  GRBException {
        double[] ub = new double[size];
        double[] lb = new double[size];
        Arrays.fill(ub, 1.0);

        return qp.addVars(lb, ub, null, null, null);
    }

    private GRBVar[][] transpose (GRBVar[][] mat) {
        GRBVar[][] transpose=new GRBVar[mat[0].length][mat.length];
        for (int i=0; i<mat.length; i++) {
            for (int j=0; j<mat[0].length; j++) {
                transpose[j][i]=mat[i][j];
            }
        }
        return transpose;
    }

    private GRBVar getStateVar(int i) {
        if (map == null) {
            return stateVars[i];
        }
        else {
            return stateVars[map.get(i)];
        }
    }

    private void createMinMaxConstraint(GRBVar[] baseVars, GRBVar[] vars, boolean max) throws GRBException{
      /*
       GRBVar[] sVars = createGRBVars(vars.length);
       double[] z=new double[sVars.length];
       Arrays.fill(z, 1.0);
       for (int i=0; i<z.length; i++) {
           z[i]=((double)i)*1.0;
       }

       GRBLinExpr[] lhs=new GRBLinExpr[vars.length];
       double[] rhs=new double[vars.length];
       char[] senses=new char[vars.length];
       Arrays.fill(senses, GRB.EQUAL);

       for (int i=0; i<vars.length; i++) {
           lhs[i]=new GRBLinExpr();
           lhs[i].addTerm(1.0, var);
           lhs[i].addTerm(-1.0, vars[i]);
           if (max) {
               lhs[i].addTerm(-1.0, sVars[i]);
           }
           else {
               lhs[i].addTerm(1.0, sVars[i]);
           }
       }
       qp.addConstrs(lhs, senses, rhs, null);
       qp.addSOS(sVars, z, GRB.SOS_TYPE1);
       */
        GRBVar[] sVars = createGRBVars(vars.length);
        double[] lb=new double[sVars.length];
        double[] ub=new double[sVars.length];
        char[] type=new char[sVars.length];
        Arrays.fill(type, GRB.BINARY);
        Arrays.fill(ub, 1.0);
        GRBVar[] zVars = qp.addVars(lb, ub,null, type, null);

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
            qp.addConstrs(lhs, senses, rhs, null);
        }

        //only one z is set
        GRBLinExpr oneZ = new GRBLinExpr();
        for (int i=0; i<vars.length; i++) {
            oneZ.addTerm(1.0, zVars[i]);
        }
        qp.addConstr(oneZ, GRB.EQUAL, 1.0, null);
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
        qp.addConstrs(lhs, senses, rhs, null);
    }

    private void createMinMaxConstraint(GRBVar var, GRBVar[] vars, boolean max) throws GRBException{
      /*
       GRBVar[] sVars = createGRBVars(vars.length);
       double[] z=new double[sVars.length];
       Arrays.fill(z, 1.0);
       for (int i=0; i<z.length; i++) {
           z[i]=((double)i)*1.0;
       }

       GRBLinExpr[] lhs=new GRBLinExpr[vars.length];
       double[] rhs=new double[vars.length];
       char[] senses=new char[vars.length];
       Arrays.fill(senses, GRB.EQUAL);

       for (int i=0; i<vars.length; i++) {
           lhs[i]=new GRBLinExpr();
           lhs[i].addTerm(1.0, var);
           lhs[i].addTerm(-1.0, vars[i]);
           if (max) {
               lhs[i].addTerm(-1.0, sVars[i]);
           }
           else {
               lhs[i].addTerm(1.0, sVars[i]);
           }
       }
       qp.addConstrs(lhs, senses, rhs, null);
       qp.addSOS(sVars, z, GRB.SOS_TYPE1);
       */
        GRBVar[] sVars = createGRBVars(vars.length);
        double[] lb=new double[sVars.length];
        double[] ub=new double[sVars.length];
        char[] type=new char[sVars.length];
        Arrays.fill(type, GRB.BINARY);
        Arrays.fill(ub, 1.0);
        GRBVar[] zVars = qp.addVars(lb, ub,null, type, null);

        GRBLinExpr[] lhs=new GRBLinExpr[vars.length];
        double[] rhs=new double[vars.length];
        char[] senses=new char[vars.length];
        Arrays.fill(senses, GRB.EQUAL);

        for (int i = 0; i < vars.length; i++) {
            lhs[i] = new GRBLinExpr();
            lhs[i].addTerm(1.0, var);
            lhs[i].addTerm(-1.0, vars[i]);
            if (max) {
                lhs[i].addTerm(-1.0, sVars[i]);
            } else {
                lhs[i].addTerm(1.0, sVars[i]);
            }
        }
        qp.addConstrs(lhs, senses, rhs, null);

        //only one z is set
        GRBLinExpr oneZ = new GRBLinExpr();
        for (int i=0; i<vars.length; i++) {
            oneZ.addTerm(1.0, zVars[i]);
        }
        qp.addConstr(oneZ, GRB.EQUAL, 1.0, null);
        //s+z<=1 only if z is set to a var where s is 0 this is correct for all
        //s+z-1<=0
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
        qp.addConstrs(lhs, senses, rhs, null);
    }


}
