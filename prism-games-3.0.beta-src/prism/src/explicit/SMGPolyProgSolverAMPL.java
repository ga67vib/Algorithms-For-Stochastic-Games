package explicit;

import com.ampl.*;
import common.IntSet;
import gurobi.GRBException;
import prism.PrismException;
import prism.PrismSettings;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class SMGPolyProgSolverAMPL extends SMGPolyProgSolver{

    //AMPL-RELATED STUFF:
    final static String STATE_PREFIX = "s";
    final static String AMPL_PATH = "";
    final static String AMPL_PATH_FROM_SCRIPT= "../prism-games-3.0.beta-src/prism";
    final static String MODEL_NAME = "model.mod";
    final static String MODEL_DIR = "";
    final static String MODEL_PATH = MODEL_DIR+MODEL_NAME;
    String solver;

    Environment env;
    AMPL ampl;

    public SMGPolyProgSolverAMPL(STPGModelChecker modelChecker, STPG stpg, BitSet no, BitSet yes, boolean min1,
                                   boolean min2, double[] init, BitSet known, PrismSettings settings, int verbosity) {
        super(modelChecker, stpg, no, yes, min1, min2, init, known, settings, verbosity);
        try {
            env = new Environment(AMPL_PATH);
            ampl = new AMPL(env);
        }
        catch (RuntimeException e) {
            System.out.println("USING FOLLOWING PATH TO AMPL EXECUTABLE: "+AMPL_PATH_FROM_SCRIPT);
            env = new Environment(AMPL_PATH_FROM_SCRIPT);
            ampl = new AMPL(env);
        }
    }

    public ModelCheckerResult solve(int variant) throws GRBException, PrismException {
        ModelCheckerResult res=null;
        String method = "";

        buildTimeStart = System.nanoTime();
        highestPrecision = true;
        mecDetected = true;

        solver = "minos";
        switch (variant) {
            case 10:
                method="AMPL minos with warm start"; //ToDo
                res = modelChecker.computeReachProbsValIter(stpg, no, yes, min1, min2, init, known, false);
                res = computeReachProbsPolyProgAMPL(stpg,no,yes,min1,min2,init,known,false,res);
                break;
            case 11:
                method="AMPL minos with topological sort";
                res = computeReachProbsPolyProgAMPL(stpg,no,yes,min1,min2,init,known,true,null);
                break;
            default:
                method="AMPL minos";
                res = computeReachProbsPolyProgAMPL(stpg,no,yes,min1,min2,init,known,false,null);
                break;
        }

        printTimes();
        return res;
    }

    private ModelCheckerResult solveTopological(STPG stpg, BitSet no, BitSet yes, boolean min1, boolean min2, double init[], List<BitSet> mecs, BitSet known, ModelCheckerResult preRes)
            throws  PrismException{
        long tranformTime = System.nanoTime() - buildTimeStart;
        long buildAndSolveTime = System.nanoTime();

        BitSet unknown;
        SCCInfo sccs;
        unknown = new BitSet();
        unknown.set(0, stpg.getNumStates()); // what does this do?
        unknown.andNot(yes);
        unknown.andNot(no);

        BitSet statesInSCC;
        known = new BitSet(); //I suppose I know nothing
        double[] values = new double[stpg.getNumStates()];
        initializeYesValues(values, yes);
        known.or(yes);
        known.or(no);

        System.out.println("Getting topologically ordered SCCs...");
        sccs = SCCComputer.computeTopologicalOrdering(modelChecker, stpg, true, unknown::get);

        ModelCheckerResult res = new ModelCheckerResult();
        res.soln = new double[stpg.getNumStates()];
        ModelCheckerResult tmpRes;
        for (int scc=0; scc<sccs.getNumSCCs(); scc++) {
            if (sccs.isSingletonSCC(scc)) {
                int state = sccs.getStatesForSCC(scc).iterator().nextInt();
                computeSingletonSCC(stpg, state, values);
                known.set(state);
            }
            else {
                IntSet sccIntset = sccs.getStatesForSCC(scc);
                statesInSCC = new BitSet();
                for (int state : sccIntset) {
                    statesInSCC.set(state);
                }
                stpgIndexToSCCIndex = new HashMap<>();
                int i=0;
                for (int state=statesInSCC.nextSetBit(0); state>=0; state=statesInSCC.nextSetBit(state+1)) {
                    stpgIndexToSCCIndex.put(state, i);
                    i++;
                }
                tmpRes = computeReachProbsPolyProgAMPL(stpg, no, yes, min1, min2, init, mecs, known, values, statesInSCC, preRes);
                i=0;
                for (int state = statesInSCC.nextSetBit(0); state>=0; state=statesInSCC.nextSetBit(state+1)) {
                    values[state]=tmpRes.soln[i];
                    i++;
                }
                known.or(statesInSCC);
            }
        }

        buildAndSolveTime = System.nanoTime() - buildAndSolveTime;
        res.soln = values;
        res.timeTaken = buildAndSolveTime + tranformTime;
        return res;
    }

    private void constructOptimizationProblem(STPG stpg, BitSet no, BitSet yes, boolean min1, boolean min2, double init[], List<BitSet> mecs, BitSet known, double[] knownValues, BitSet compute, ModelCheckerResult preRes, ArrayList<String> variables) throws PrismException {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(MODEL_PATH, false));
            boolean topological;
            if (compute==null) {
                topological = false;
                compute = new BitSet();
                compute.set(0,stpg.getNumStates());
            }
            else {
                topological = compute.cardinality() != stpg.getNumStates();
            }
            if (known==null) {
                known = new BitSet();
            }
            if (mecs==null) {
                mecs = new ArrayList<>();
            }



            int stateNum = compute.cardinality();
            String[] stateVars = new String[stateNum];
            String stateName;

            String s = "";

            boolean inMec;
            //Create Variables
            s="";
            writer.append("#variable declarations\n");
            int i=0;
            for (int state=compute.nextSetBit(0); state>=0; state = compute.nextSetBit(state+1)) {
                stateName = STATE_PREFIX+state;
                stateVars[i]=(stateName);
                if (yes.get(state)) {
                    s= "param "+stateName+" =1; #Yes";
                }
                else if (no.get(state)) {
                    s= "param "+stateName+" =0; #No";
                }
                //else if (stpg.getNumChoices(i)==1) {
                //    s= "param "+stateName;
                //}
                else {
                    inMec=false;
                    for (BitSet mec : mecs) {
                        if (mec.get(state)) {
                            inMec = true;
                            //s= "param "+stateName+"; #MECState";
                            s= "var "+stateName+" >=0, <=1; #MECState";
                            break;
                        }
                    }
                    if (!inMec) {
                        variables.add(stateName);
                        s = "var "+stateName+" >=0, <=1;";
                    }

                }

                //else {
                //    variables.add(stateName);
                //    s = "var "+stateName+" >=0, <=1;";
                //}
                i++;
                writer.append(s+"\n");
            }

            writer.append("\n");

            //Constraints for non-MEC, non-trivial states
            writer.append("\n#Constraints for non-MEC, non-trivial states\n");
            int choiceNum;
            String sign;
            String term;
            String rhs;
            String lhs;
            boolean repeat=false;
            for (int state=compute.nextSetBit(0); state>=0; state = compute.nextSetBit(state+1)) {
                //skip all the states in MECs because the constraints are built differently
                //Note that for precomputed games yes/no-states don't count as ECs


                boolean skip=false;
                for (BitSet bs : mecs) {
                    if (bs.get(state)) {
                        if (stpg.getNumChoices(state)>1) { //We still need to implement equality constraints for the states without choice
                            skip = true;
                        }
                        break;
                    }
                    if (skip) break;
                }
                if (skip) continue;



                stateName = stateVars[getIndex(state)];

                choiceNum = stpg.getNumChoices(state);

                //We may skip state in theory since it's bounded already by bounds
                if (no.get(state) || yes.get(state)) {
                    continue;
                }
                if (stpg.getNumChoices(state) == 1) {
                    repeat = true;
                }
                for (int choice=0; choice<choiceNum; choice++) {
                    lhs = stateName;
                    rhs = "";

                    if (choiceNum==1) {
                        sign = "=";
                        if (repeat) {
                            sign = ">=";
                        }
                        else {
                            sign = "<=";
                        }
                    }
                    else {
                        sign = (stpg.getPlayer(state) == 1) ? ">=" : "<=";
                    }

                    for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext(); ) {
                        Map.Entry<Integer, Double> tr = it.next();
                        if (yes.get(tr.getKey())) { //if yes sink -1
                            rhs += " "+tr.getValue();
                            if (it.hasNext()) {
                                rhs += " +";
                            }
                        }
                        else if (!no.get(tr.getKey())) { //neither yes/no -> normal vertex
                            if (tr.getValue()==1.0) {
                                rhs += " " + getValueIfKnown(tr.getKey(), known, stateVars, knownValues);
                                if (it.hasNext()) { //This is there in case of very small numbers
                                    rhs += " +";
                                }
                            }
                            else {
                                rhs += " "+tr.getValue()+"*"+getValueIfKnown(tr.getKey(), known, stateVars, knownValues);
                                if (it.hasNext()) {
                                    rhs += " +";
                                }
                            }
                        }
                        else {
                            rhs += " 0";
                            if (it.hasNext()) {
                                rhs += " +";
                            }
                        }
                    }
                    if (choiceNum == 1 && repeat) {
                        choice--;
                        repeat=false;
                    }
                    term = "subject to const_"+stateName+"_a"+(choice+1)+": "+lhs + " " + sign + rhs + ";";
                    writer.append(term+"\n");
                }
            }

            //MEC-Solving
            writer.append("\n#MEC-Constraints\n");
            SMGMECSolverAMPL mecSolver= new SMGMECSolverAMPL(mecs, stpg, known, knownValues, writer, stateVars);
            try {
                if (!topological) {
                    mecSolver.solve(); //we have to treat everything
                }
                else { //We only want to solve mecs that are relevant in this iteration
                    for (BitSet mec : mecs) {
                        //only solve mecs of which at least one states appears in the compute set
                        if (compute.get(mec.nextSetBit(0))) {
                            mecSolver.solve(mec, stpgIndexToSCCIndex);
                        }

                    }
                }
            }
            catch (GRBException e) {
                System.out.println("GUROBI-EXCEPTION IN AMPL ARISED. SHOULDN'T HAPPEN");
                e.printStackTrace();
            }


            //Objective:
            writer.append("\n#Objective function\n");
            writer.append("minimize value: ");
            s="";

            boolean impar;
            for (int state=compute.nextSetBit(0); state>=0; state = compute.nextSetBit(state+1)) {
                stateName = getValueIfKnown(state, known, stateVars, knownValues);
                choiceNum = stpg.getNumChoices(state);

                if (!variables.contains(stateName)) {
                    continue; //in this case it's just a parameter
                }
                //Nothing to optimize for these states
                if (choiceNum == 1) {
                    continue;
                }
                impar = (choiceNum%2 == 1);

                for (int choice=0; choice<choiceNum; choice++) {
                    term = "("+stateName+" - (";

                    for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext(); ) {
                        Map.Entry<Integer, Double> tr = it.next();
                        if (yes.get(tr.getKey())) { //if yes sink -1
                            term += " "+tr.getValue();
                        }
                        else if (no.get(tr.getKey())) {
                            term += " 0"; //I add this in case the state goes only to zero-sinks
                        }
                        else if (!no.get(tr.getKey())) { //neither yes/no -> normal vertex
                            if (tr.getValue() == 1.0) {
                                term += " " + getValueIfKnown(tr.getKey(), known, stateVars, knownValues);
                            } else {
                                term += " " + tr.getValue() + "*" + getValueIfKnown(tr.getKey(), known, stateVars, knownValues);
                            }
                        }

                        if (it.hasNext()) {
                            term += " + ";
                        }

                    }
                    if (impar) {
                        choice--;
                        impar=false;
                    }
                    term += "))";
                    if (choice<choiceNum-1) {
                        term += " * ";
                    }
                    else {
                        term += " + \n";
                    }
                    writer.append(term);

                }
            }

            s +="0;\n"; //not a very elegant way, but I don't want to change the "writer.append appearance"
            writer.append(s);

            System.out.println("finished constructing model");

            writer.close();
        }
        catch (java.io.IOException e) {
            System.out.println("---------------Writing failed!!!-----------------");
            e.printStackTrace();
        }

    }

    private ModelCheckerResult computeReachProbsPolyProgAMPL(STPG stpg, BitSet no, BitSet yes, boolean min1, boolean min2, double init[], BitSet known, boolean topological, ModelCheckerResult preRes)
            throws  PrismException{
        List<BitSet> mecs=getMECs((MDPSimple) stpg, no, yes);
        if (topological) {
            return solveTopological(stpg, no, yes, min1, min2, init, mecs, known, preRes);
        }
        else {
            BitSet compute = new BitSet();
            compute.set(0, stpg.getNumStates());
            return computeReachProbsPolyProgAMPL(stpg, no, yes, min1, min2, init, mecs, known, null, compute, preRes);
        }
    }

    private ModelCheckerResult computeReachProbsPolyProgAMPL(STPG stpg, BitSet no, BitSet yes, boolean min1, boolean min2, double init[], List<BitSet> mecs, BitSet known, double[] knownValues, BitSet compute, ModelCheckerResult preRes)
            throws  PrismException{


        ArrayList<String> variables = new ArrayList<>();
        constructOptimizationProblem(stpg,no,yes,min1,min2,init,mecs,known,knownValues,compute,preRes,variables);

        try {
            // Interpret the two files
            ampl.reset();
            ampl.read(MODEL_PATH);

            buildTimeEnd = System.nanoTime();
            solveTimeStart = System.nanoTime();

            // Solve
            settingsAMPL(ampl, solver);
            System.out.println("Starting AMPL with "+solver+" ...");
            ampl.solve();
            solveTimeEnd = System.nanoTime();

            ModelCheckerResult res = new ModelCheckerResult();

            Objective objective = ampl.getObjective("value");
            try {
                System.out.format("Objective value is: %f%n", objective.value());
            }
            catch (Exception e) {
                System.out.println("Objective value could not be printed out, probably the whole SG/SCC is just made out of MECs...");
            }

            int numStates = compute.cardinality();
            res.soln = new double[numStates];
            String stateName;
            int i=0;
            for (int state=compute.nextSetBit(0); state>=0; state = compute.nextSetBit(state+1)) {
                stateName = "s"+state;
                if (!variables.contains(stateName)) {
                    Parameter y = ampl.getParameter(stateName);
                    try {
                        res.soln[i] = (Double) ampl.getParameter(stateName).get();
                    }
                    catch (NullPointerException e) {
                        res.soln[i] = ampl.getVariable(stateName).value();
                    }
                }

                else {
                    res.soln[i] = ampl.getVariable(stateName).value();
                }
                i++;
            }
            if (verbosity>=2) {
                System.out.println("SOLUTION FOR STATE 0: "+res.soln[0]);
            }
            res.timeTaken = (solveTimeEnd - solveTimeStart) + (buildTimeEnd - buildTimeStart);
            return res;
        }


        catch (IOException e) {
            System.out.println("IO-Exception in AMPL-solving");
            e.printStackTrace();
        }
        return null;
    }

    private void settingsAMPL(AMPL ampl, String solver) {
        ampl.setOption("solver", solver);
        //ampl.setIntOption("times", 1);
        //ampl.setIntOption("gentimes", 1);
        //ampl.setIntOption("show_stats", 1);

        switch (solver) {
            case "gurobi":
                ampl.eval("options gurobi_options $gurobi_options 'nonconvex 2';");
                ampl.eval("options gurobi_optiosn $gurobi_options 'outlev 1';");
                ampl.eval("options gurobi_optiosn $gurobi_options 'threads 1';");
                if (highestPrecision) {
                    ampl.eval("options gurobi_options $gurobi_options 'feastol 1e-9';");
                    ampl.eval("options gurobi_options $gurobi_options 'intfeastol 1e-9';");
                    ampl.eval("options gurobi_options $gurobi_options 'opttol 1e-9';");

                    ampl.eval("options gurobi_options $gurobi_options 'barconvtol 0';");
                    ampl.eval("options gurobi_options $gurobi_options 'mipgap 0';");
                    ampl.eval("options gurobi_options $gurobi_options 'mipgapabs 0';");
                    ampl.eval("options gurobi_options $gurobi_options 'psdtol 0';");
                    ampl.eval("options gurobi_options $gurobi_options 'mipfocus 3';");
                }
                break;
            case "minos":
                ampl.eval("options minos_options " +
                        "'feasibility_tolerance=1e-6 " +
                        "optimality_tolerance=1e-6 " +
                        "superbasic_limit=5000 " +
                        //"solution=yes" +
                        "';");
                break;
            case "cplex":
            default:
        }
    }

    private String getValueIfKnown(int index, BitSet known, String[] stateVars, double[] knownValues) {
        return known.get(index) ?
                ("" +knownValues[index]) :
                (stateVars[
                        stpgIndexToSCCIndex != null ?
                        stpgIndexToSCCIndex.get(index) :
                        index
                        ]
                );
    }
}
