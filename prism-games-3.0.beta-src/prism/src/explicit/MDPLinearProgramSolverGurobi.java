package explicit;

import gurobi.*;
import prism.PrismException;
import prism.PrismSettings;

import java.util.*;

public class MDPLinearProgramSolverGurobi {
    MDPModelChecker modelChecker;
    MDP mdp;
    BitSet no;
    BitSet yes;
    boolean min;
    double init[];
    BitSet known;
    protected PrismSettings settings;
    protected int verbosity;
    BitSet remain;
    BitSet target;

    boolean mecDetected;
    //TOPOLOGICAL
    HashMap<Integer, Integer> mdpIndexToSCCIndex;

    public MDPLinearProgramSolverGurobi(MDPModelChecker modelChecker, MDP mdp, BitSet no, BitSet yes, boolean min, double[] init, BitSet known, PrismSettings settings) {
        this.modelChecker = modelChecker;
        this.mdp = mdp;
        this.no = no;
        this.yes = yes;
        this.min = min;
        this.init = init;
        this.known = known;
        this.settings = settings;
        this.verbosity = 2; // ToDo: Should get fixed one day to better verbosity

        this.mecDetected = false;
        mdpIndexToSCCIndex = null;
    }
    public ModelCheckerResult solve() throws PrismException{
        ModelCheckerResult res = new ModelCheckerResult();
        try {
            res = this.computeReachProbsLinearProgrammingGurobi(mdp, no, yes, this.min, null);
        }
        catch (GRBException e) {
            e.printStackTrace();
        }
        return res;
    }

    protected List<BitSet> getMECs(MDP mdp, BitSet no, BitSet yes)
            throws PrismException
    {

        //compute MECs one time, use the decomposition in every iteration; SECs still have to be recomputed
        explicit.ECComputerDefault ec = (ECComputerDefault) ECComputer.createECComputer(this.modelChecker, mdp);
        BitSet all = new BitSet();
        all.set(0,mdp.getNumStates());
        all.xor(yes);//ignore states that are already finished
        all.xor(no);
        ec.computeMECStates(all);
        List<BitSet> mecs = ec.getMECStates();

        if (!mecs.isEmpty()) {
            mecDetected = true;
        }
        return mecs;
    }

    protected int getIndex(int i) {
        if (mdpIndexToSCCIndex == null) {
            return i;
        }
        else {
            return mdpIndexToSCCIndex.get(i);
        }
    }

    /**
     * Compute the value of a state that builds an SCC alone. We expect every action transition from the state leading
     * to an already computed SCC!
     * If for whatever reason we have here a 1-state-MEC that wasn't found by yes and no,
     * it will also be set to 0 (which is correct)
     * @param mdp
     * @param state
     * @param values
     */
    protected void computeSingletonSCC(MDP mdp, int state, double[] values, boolean minMDP) {
        int numChoices = mdp.getNumChoices(state);

        double max=0.0; //can compute on the fly
        double min=1.0;
        double actualValue;
        for (int choice = 0; choice<numChoices; choice++) {
            actualValue = 0.0;
            for (Iterator<Map.Entry<Integer, Double>> it = mdp.getTransitionsIterator(state, choice); it.hasNext(); ) {
                Map.Entry<Integer, Double> tr = it.next();
                actualValue += tr.getValue() * values[tr.getKey()];
            }
            max = Math.max(actualValue, max);
            min = Math.min(actualValue, min);
        }

        if (!minMDP) {
            values[state] = max;
        }
        else {
            values[state] = min;
        }
    }

    public ModelCheckerResult computeReachProbsLinearProgrammingGurobi (MDP mdp, BitSet no, BitSet yes, boolean min, int strat[]) throws GRBException, PrismException {
        int numStates = mdp.getNumStates();
        GRBEnv grbEnv = new GRBEnv();
        GRBModel linearProgram = new GRBModel(grbEnv);

        List<BitSet> mecs = getMECs(mdp, no, yes);

        //=================CREATE VARIABLES===============================
        double lb[] = new double[numStates]; // also initializes array to 0
        double ub[] = new double[numStates];

        // Initialize all upper bounds to 1
        Arrays.fill(ub, 1.0);

        // Create state variables, type = null implies all variables are of continuous type
        GRBVar[] stateVars = linearProgram.addVars(lb, ub, null, null, null);

        //======================CREATE CONSTRAINTS=======================
        int numChoices = mdp.getNumChoices();
        GRBLinExpr[] lhs = new GRBLinExpr[numChoices];
        double[] rhs = new double[numChoices];
        char[] senses = new char[numChoices];

        int k = 0; //constraint counter
        for (int state = 0; state < numStates; state++) {
            //Skip all MEC-states because their constraints are built differently
            boolean skip = false;
            for (BitSet mec : mecs) {
                if (mec.get(state)) {
                    if (mdp.getNumChoices() > 1) { //We still need to implement equality constraints for the states without choice
                        skip = true;
                    }
                    break;
                }
            }
            if (skip) continue;

            int numStateChoices = mdp.getNumChoices(state);
            for (int choice = 0; choice < numStateChoices; choice++) {
                // Add v(i)
                lhs[k] = new GRBLinExpr();
                lhs[k].addTerm(1.0, stateVars[state]);

                //putting the correct sense. ToDo HAS TO BE ADJUSTED TO min1, min2
                if (numStateChoices == 1) {
                    senses[k] = GRB.EQUAL;
                } //supposing 1 -> maxPlayer 2-> minplayer
                else if (!min) {
                    senses[k] = GRB.GREATER_EQUAL;
                } else if (min) {
                    senses[k] = GRB.LESS_EQUAL;
                }

                //If this is a yes-state, just model it as a state with a self-loop
                if (yes.get(state)) {
                    lhs[k].addConstant(-1.0);
                    k++;
                    continue;
                }

                for (Iterator<Map.Entry<Integer, Double>> it = mdp.getTransitionsIterator(state, choice); it.hasNext(); ) {
                    Map.Entry<Integer, Double> tr = it.next();

                    if (yes.get(tr.getKey())) { //The transition leads to a yes-sink
                        lhs[k].addConstant(-1.0 * tr.getValue());
                    } else if (!no.get(tr.getKey())) { //neither yes/no -> normal vertex
                        lhs[k].addTerm(-1.0 * tr.getValue(), stateVars[tr.getKey()]);
                    }

                }
                k++;
            }
        }

        //All Constraints that are reserved for states in MECs are innecessary
        //We cut them out since there we will set up special constraints for them
        GRBLinExpr[] lhsWithoutMEC;
        double[] rhsWithoutMEC;
        char[] sensesWithoutMEC;

        if (!mecs.isEmpty()) {
            int size = 0;
            for (GRBLinExpr lh : lhs) {
                if (lh != null) {
                    size++;
                }
            }
            lhsWithoutMEC = new GRBLinExpr[size];
            rhsWithoutMEC = new double[size];
            sensesWithoutMEC = new char[size];
            int j = 0;
            for (int i = 0; i < lhs.length; i++) {
                if (lhs[i] != null) {
                    lhsWithoutMEC[j] = lhs[i];
                    rhsWithoutMEC[j] = rhs[i];
                    sensesWithoutMEC[j] = senses[i];
                    j++;
                }
            }
        } else {
            lhsWithoutMEC = lhs;
            rhsWithoutMEC = rhs;
            sensesWithoutMEC = senses;
        }

        // Add non-MEC-constaints to LP
        linearProgram.addConstrs(lhsWithoutMEC, sensesWithoutMEC, rhsWithoutMEC, null);

        //MEC-Solving - also adds the constraints
        MDPMECSolverGurobi mecSolver = new MDPMECSolverGurobi(linearProgram, stateVars, mecs, mdp, min, new BitSet(), null);
        mecSolver.solve(); //we have to treat everything

        //-----------------Objective Function---------------------
        GRBLinExpr objExpr = new GRBLinExpr();
        for (int state = 0; state < numStates; state++) {
            objExpr.addTerm(1.0, stateVars[state]);
        }

        linearProgram.setObjective(objExpr, min ? GRB.MAXIMIZE : GRB.MINIMIZE);

        linearProgram.set(GRB.DoubleParam.FeasibilityTol, 1e-9);
        linearProgram.set(GRB.DoubleParam.IntFeasTol, 1e-9);
        linearProgram.set(GRB.DoubleParam.OptimalityTol, 1e-9);

        linearProgram.set(GRB.DoubleParam.BarQCPConvTol, 0);
        linearProgram.set(GRB.DoubleParam.BarConvTol, 0);
        linearProgram.set(GRB.DoubleParam.MarkowitzTol, 1e-4);
        linearProgram.set(GRB.DoubleParam.MIPGap, 0);
        linearProgram.set(GRB.DoubleParam.MIPGapAbs, 0);
        linearProgram.set(GRB.DoubleParam.PSDTol, 0);

        linearProgram.optimize();

        int optimstatus = linearProgram.get(GRB.IntAttr.Status);

        if (optimstatus == GRB.Status.INFEASIBLE) {
            System.out.println("Model is infeasible");

            // Compute and write out IIS
            linearProgram.computeIIS();
            linearProgram.write("model.ilp");
        }

        // Process result and put it in res
        ModelCheckerResult res = new ModelCheckerResult();
        res.soln = new double[numStates];
        for (int i = 0; i < numStates; i++) {
            res.soln[i] = stateVars[i].get(GRB.DoubleAttr.X);
        }
        return res;
    }
}
