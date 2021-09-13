package explicit;

import Jama.Matrix;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBVar;
import prism.Pair;
import prism.PrismException;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;

public class SMGMECSolverAMPL extends SMGMECSolver {

    BufferedWriter writer;
    long tmpVars;
    String[] stateNames;

    final String BEST_EXIT_PREFIX = "tmpBestExit";
    final String EQUAL_PREFIX = "ec_EQ";
    final String MIN_PREFIX = "ec_MIN";
    final String MAX_PREFIX = "ec_MAX";

    private HashMap<Integer, Integer> map;

    public SMGMECSolverAMPL(List<BitSet> mecs, STPG stpg, BitSet known, double[] knownValues, BufferedWriter writer, String[] stateNames) {
        super(mecs, stpg, known, knownValues);
        this.writer = writer;
        this.tmpVars = 0;
        this.stateNames = stateNames;
        this.map=null;
    }

    //Would be cool if I could remove GRBException here
    @Override
    public void solve() throws PrismException, GRBException {
        for (BitSet mec : mecs) {
            createStrategyTable(mec);
        }
    }

    public void solve(BitSet mec, HashMap<Integer, Integer> map) throws PrismException, GRBException {
        this.map = map;
        createStrategyTable(mec);
    }

    @Override
    protected void setMECtoZero(BitSet mec) {
        String stateName;
        String constraint;

        for (int state = mec.nextSetBit(0); state >= 0; state = mec.nextSetBit(state + 1)) {
            stateName = stateNames[getStateIndex(state)];

            String constraintMin = "subject to const_ZeroMEC_"+stateName+"_Leq: "+stateName+" <= 0;\n";
            String constraintMax = "subject to const_ZeroMEC_"+stateName+"_Geq: "+stateName+" >= 0;\n";

            //constraint = "redeclare param ";
            //constraint += stateName + " = 0;\n";

            append(constraintMin);
            append(constraintMax);
        }

    }

    @Override
    protected void setMECtoBestExit(BitSet mec, BitSet exits) {
        ArrayList<Pair<Integer, Integer>> exitStateActions = getExitStateActionPairs(mec, exits);
        ArrayList<String> tmpVarNames = new ArrayList<>();
        int sizeExitActions = exitStateActions.size();
        if (sizeExitActions == 0) {
            setMECtoZero(mec);
            return;
        }

        String stateName;
        String constraint;

        String tmpVar;
        String rhs;

        //Create tmpVars over that we maximize later
        //Each tmpVar represents one pair V(s,a)
        for (int i = 0; i < sizeExitActions; i++) {
            tmpVar = BEST_EXIT_PREFIX + tmpVars;
            tmpVarNames.add(tmpVar);
            tmpVars++;

            constraint = "param " + tmpVar + " =";

            Pair<Integer, Integer> pair = exitStateActions.get(i);
            // if we have P(a,b,c) but c leads into EC and a, b are exits then the best exit is
            // then exit probability is a/(a+b) and b/(a+b)
            double divisor = 0;
            boolean allExit = true;
            for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(pair.getKey(), pair.getValue()); it.hasNext(); ) {
                Map.Entry<Integer, Double> tr = it.next();
                if (exits.get(tr.getKey())) {
                    divisor += tr.getValue();
                } else {
                    allExit = false;
                }
            }
            if (allExit) divisor = 1.0; //just to increase precision

            rhs = "";
            for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(pair.getKey(), pair.getValue()); it.hasNext(); ) {
                Map.Entry<Integer, Double> tr = it.next();
                if (exits.get(tr.getKey())) { //Add every leaving Probability
                    rhs += " " + tr.getValue() / divisor + " * " + getValueIfKnown(tr.getKey());
                    rhs += it.hasNext() ? " + " : ";\n";
                }
            }
            constraint = createParameter(tmpVar, rhs);
            append(constraint);
        }
        int i = 0;
        for (int state = mec.nextSetBit(0); state >= 0; state = mec.nextSetBit(state + 1)) {
            if (state == Integer.MAX_VALUE) {
                break;
            }
            stateName = stateNames[getStateIndex(state)];
            String constraintMin = "subject to const_"+"BestExit_"+stateName+"_Min: "+stateName+" >=";
            String constraintMax = "subject to const_"+"BestExit_"+stateName+"_Max: "+stateName+" <=";
            //constraint = "redeclare param "+stateName+" =";

            String[] stringArray = new String[0];
            stringArray = tmpVarNames.toArray(stringArray);
            rhs = constructMinMax(stringArray, true);
            constraintMin += " " + rhs + ";\n";
            constraintMax += " " + rhs + ";\n";

            append(constraintMin);
            append(constraintMax);

            i++;
        }
    }

    protected void handleMixedMECTMP(BitSet mec, BitSet exits) {
        ArrayList<Integer> maxStatesNames = new ArrayList<>();
        ArrayList<Integer> maxStatesActions = new ArrayList<>();
        ArrayList<Integer> minStatesNames = new ArrayList<>();
        ArrayList<Integer> minStatesActions = new ArrayList<>();
        int maxActions = 1;
        int minActions = 1;

        for (int state = mec.nextSetBit(0); state >= 0; state = mec.nextSetBit(state + 1)) {
            if (state == Integer.MAX_VALUE) {
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
        }

        int numMECStates = maxStatesNames.size() + minStatesNames.size();

        String[][] minConstrVars = new String[maxActions][numMECStates]; //For each maxAction we have a set of minConstrVariables -> In the whole table there
        String[][] equalVars = new String[minActions][numMECStates]; // is one maxAction fixed with all equal variables: minActions x numMECStates
        initializeVarNames(minConstrVars, MIN_PREFIX);


        for (int x = 0; x < maxActions; x++) {
            initializeVarNames(equalVars, EQUAL_PREFIX);
            //minConstrVars[x]=createGRBVars(numMECStates);
            for (int y = 0; y < minActions; y++) {
                Matrix exitProbs = solveMarkovChain(maxStatesNames, maxStatesActions, minStatesNames, minStatesActions, exits);
                getEqualConstraints(exitProbs, maxStatesNames.size() + minStatesNames.size(), exits, equalVars[y]);
                incrementAction(minStatesActions, minStatesNames);
            }
            getMinConstraints(minConstrVars[x], equalVars);
            incrementAction(maxStatesActions, maxStatesNames);
        }
        getMaxConstraints(minConstrVars, maxStatesNames, minStatesNames);

    }

    @Override
    protected void handleMixedMEC(BitSet mec, BitSet exits) throws PrismException{
        ArrayList<Integer> maxStatesNames = new ArrayList<>();
        ArrayList<Integer> maxStatesActions = new ArrayList<>();
        ArrayList<Integer> minStatesNames = new ArrayList<>();
        ArrayList<Integer> minStatesActions = new ArrayList<>();
        int maxActions = 1;
        int minActions = 1;

        for (int state = mec.nextSetBit(0); state >= 0; state = mec.nextSetBit(state + 1)) {
            if (state == Integer.MAX_VALUE) {
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

        String[][] minConstrVars = new String[maxActions][numMECStates]; //For each maxAction we have a set of minConstrVariables -> In the whole table there
        //initializeVarNames(minConstrVars, MIN_PREFIX);

        String[][][] equalVars = new String[maxActions][minActions][numMECStates];


        for (int x = 0; x < maxActions; x++) {
            //initializeVarNames(equalVars[x], EQUAL_PREFIX);
            //minConstrVars[x]=createGRBVars(numMECStates);
            for (int y = 0; y < minActions; y++) {
                Matrix exitProbs = solveMarkovChain(maxStatesNames, maxStatesActions, minStatesNames, minStatesActions, exits);
                getEqualConstraintsRhs(exitProbs, maxStatesNames.size() + minStatesNames.size(), exits, equalVars[x][y]);
                incrementAction(minStatesActions, minStatesNames);
            }
            getMinConstraintsRhs(minConstrVars[x], equalVars[x]);
            incrementAction(maxStatesActions, maxStatesNames);
        }
        getMaxConstraints(minConstrVars, maxStatesNames, minStatesNames);

    }

    private void getEqualConstraintsRhs(Matrix probs, int numECStates, BitSet exits, String[] eqPlaceholder) {
        String rhs;
        Boolean reachesExit = false;

        for (int i = 0; i < numECStates; i++) {
            rhs = "";
            int exit = exits.nextSetBit(0);
            for (int j = 0; j < probs.getColumnDimension(); j++) {
                double prob = probs.get(i, j);
                if (prob != 0.0) {
                    reachesExit = true;
                    rhs += " " + prob + "*" + getValueIfKnown(exit);
                } else {
                    rhs += " 0";
                }

                if (j < probs.getColumnDimension() - 1) {
                    rhs += " +";
                }
                exit = exits.nextSetBit(exit + 1);
            }
            if (!reachesExit) {
                rhs = " 0";
            }
            eqPlaceholder[i] = rhs;
        }
    }

    private void getEqualConstraints(Matrix probs, int numECStates, BitSet exits, String[] varNames) {
        String constraint;
        String rhs;
        Boolean reachesExit = false;

        for (int i = 0; i < numECStates; i++) {
            constraint = "var " + varNames[i] + ";\n";
            append(constraint);

            for (int run = 0; run < 2; run++) {
                int exit = exits.nextSetBit(0);
                constraint = "subject to " + "EC_CONST_" + varNames[i] + "_" + run + ": " + varNames[i];
                if (run == 0) {
                    constraint += " >= ";
                } else {
                    constraint += " <= ";
                }
                rhs = "";
                for (int j = 0; j < probs.getColumnDimension(); j++) {
                    double prob = probs.get(i, j);
                    if (prob != 0.0) {
                        reachesExit = true;
                        rhs += " " + prob + "*" + getValueIfKnown(exit);
                    } else {
                        rhs += " 0";
                    }

                    if (j < probs.getColumnDimension() - 1) {
                        rhs += " +";
                    }
                    exit = exits.nextSetBit(exit + 1);
                }
                if (!reachesExit) {
                    rhs = " 0";
                }
                constraint += rhs + ";\n";
                append(constraint);

            }
        }
    }

    private void getMinConstraints(String[] minConstrVars, String[][] equalVars) {
        String[][] equalVarsTranspose = transpose(equalVars);

        String constraint;
        String rhs;

        for (int i = 0; i < minConstrVars.length; i++) {
            constraint = "param " + minConstrVars[i] + " =";
            rhs = constructMinMax(equalVarsTranspose[i], false);
            constraint += " " + rhs + ";\n";
            constraint = createParameter(minConstrVars[i], rhs);
            append(constraint);
        }
    }

    private void getMinConstraintsRhs(String[] minConstrVars, String[][] equalVars) {
        String[][] equalVarsTranspose = transpose(equalVars);
        String rhs;

        for (int i = 0; i < minConstrVars.length; i++) {
            rhs = constructMinMax(equalVarsTranspose[i], false);
            minConstrVars[i] = rhs;
        }
    }

    private void getMaxConstraints(String[][] minConstrVars, ArrayList<Integer> maxStatesNames, ArrayList<Integer> minStatesNames) {
        String[][] minConstrVarsTranspose = transpose(minConstrVars);

        int numECStates = maxStatesNames.size() + minStatesNames.size();

        ArrayList<Integer> states;
        int listIndex;

        String constraint;
        String rhs;
        String stateName;

        for (int run = 0; run < 2; run++) {
            listIndex = 0;
            if (run == 0) {
                states = maxStatesNames;
            } else {
                states = minStatesNames;
            }

            //This loop runs first time from 0 -> maxStates.size -1 , second time maxStates.size -> numECStates-1
            for (int i = maxStatesNames.size() * run; i < maxStatesNames.size() + minStatesNames.size() * run; i++) {
                stateName = getValueIfKnown(states.get(listIndex));
                rhs = constructMinMax(minConstrVarsTranspose[i], true);
                for (int run_constraint = 0; run_constraint<2; run_constraint++) {
                    constraint = "subject to "+MAX_PREFIX+"_"+stateName+"_"+run_constraint+": ";
                    constraint += stateName + (run_constraint==0 ? " >=" : " <=");
                    constraint += " " + rhs +";\n";
                    append(constraint);
                }
                //stateName = getValueIfKnown(states.get(listIndex));
                //constraint = "redeclare param ";
                //constraint += stateName + " =";

                //constraint += " " + rhs + ";\n";
                //append(constraint);
                listIndex++;
            }
        }
    }

    private String createParameter(String paramName, String rhs) {
        String paramDefinition = "";

        paramDefinition = "var "+paramName+" >= 0, <= 1;\n";
        paramDefinition += "subject to const_"+paramName+"_LEQ: "+ paramName+" <= " + rhs+"";
        paramDefinition += "subject to const_"+paramName+"_GEQ: "+ paramName+" >= " + rhs+"";

        return paramDefinition;
    }

    private void initializeVarNames(String[][] vars, String prefix) {
        for (int i = 0; i < vars.length; i++) {
            for (int j = 0; j < vars[i].length; j++) {
                vars[i][j] = prefix + tmpVars;
                tmpVars++;
            }
        }
    }

    private String[][] transpose(String[][] mat) {
        String[][] transpose = new String[mat[0].length][mat.length];
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat[0].length; j++) {
                transpose[j][i] = mat[i][j];
            }
        }
        return transpose;
    }

    /**
     * THIS DOES ONLY PROVIDE max(a,b,c,...) / min(a,b,c,...)
     *
     * @param inVars
     * @param max
     * @return
     */
    private String constructMinMax(String[] inVars, boolean max) {
        String rhs = "";
        String constraint = max ? "max" : "min";
        constraint += "(";
        for (int inVar = 0; inVar < inVars.length; inVar++) {
            rhs += " " + inVars[inVar];
            if (inVar < inVars.length - 1) {
                rhs += ",";
            }
        }
        constraint += rhs + ")";
        return constraint;
    }


    private void append(String s) {
        try {
            writer.append(s);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private int getStateIndex(int index) {
        return (map != null) ? map.get(index) : index;
    }

    private String getValueIfKnown(int index) {
        return known.get(index) ?
                ("" +knownValues[index]) :
                (map != null ?
                        stateNames[map.get(index)] :
                        stateNames[index]);
    }
}
