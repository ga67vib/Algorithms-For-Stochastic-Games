package explicit;

import common.IntSet;
import gurobi.*;
/** =============CPLEX-RELATED=============
import ilog.concert.IloException;
import ilog.concert.IloLPMatrix;
import ilog.cplex.IloCplex;
**/
import prism.ModelChecker;
import prism.PrismException;
import prism.PrismSettings;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.util.*;

public class SMGPolyProgSolverGurobi extends SMGPolyProgSolver {

    public SMGPolyProgSolverGurobi(STPGModelChecker modelChecker, STPG stpg, BitSet no, BitSet yes, boolean min1,
                       boolean min2, double[] init, BitSet known, PrismSettings settings, int verbosity) {
        super(modelChecker, stpg, no, yes, min1, min2, init, known, settings, verbosity);
    }

    public ModelCheckerResult solve(int variant) throws GRBException, PrismException {
        ModelCheckerResult res=null;
        String method = "";

        buildTimeStart = System.nanoTime();
        highestPrecision = true;
        mecDetected = true;
        switch (variant) {
            case 1:
                method="Improved quadratic programming with value iteration for start vector";
                res = modelChecker.computeReachProbsValIter(stpg, no, yes, min1, min2, init, known, ProbModelChecker.SolnMethod.VALUE_ITERATION);
                res = computeReachProbsQPGurobiCondonLoose(stpg, no, yes, min1, min2, init, known, res, true, false);
                break;
            case 2:
                if (verbosity>=1) {
                    System.out.println("WARNING: For this option we expect that the model that is to be solved " +
                            "is already provided in the prism folder and is called 'model.lp'");
                }
                method="CPLEX with model already provided by Gurobi";
                /**===============CPLEX-RELATED==========
                 * res = solveQPCplex();
                 */
                res = throwCplexError();
                break;
            case 3:
                method="CPLEX with model created during run by Gurobi";
                res = computeReachProbsQPGurobiCondonLoose(stpg, no, yes, min1, min2, init, known, null, false, false);
                break;
            case 4:
                method="Condon's original quadratic program";
                res = computeReachProbsQPGurobiCondonStrict(stpg, no, yes, min1, min2, init, known, true);
                break;
            case 5:
                method="Improved quadratic program, but MECs are solved with Shapleys approach with m 1/2-probability states";
                highestPrecision = true;
                res = computeReachProbsQPGurobiCondonLooseEps(stpg, no, yes, min1, min2, init, known, true, true);
                break;
            case 6:
                method="Improved quadratic program, but MECs are solved with Shapleys approach";
                highestPrecision = true;
                res = computeReachProbsQPGurobiCondonLooseEps(stpg, no, yes, min1, min2, init, known, false, true);
                break;
            case 7:
                method="Topological QP Gurobi";
                res = computeReachProbsQPGurobiCondonLoose(stpg, no, yes, min1, min2, init, known, null, true, true);
                break;
            case 8:
                method="Topological QP Gurobi with warm start";
                res = modelChecker.computeReachProbsValIter(stpg, no, yes, min1, min2, init, known, ProbModelChecker.SolnMethod.VALUE_ITERATION);
                res = computeReachProbsQPGurobiCondonLoose(stpg, no, yes, min1, min2, init, known, res, true, true);
                break;
            default:
                method="Improved quadratic programming";
                res = computeReachProbsQPGurobiCondonLoose(stpg, no, yes, min1, min2, init, known, null, true, false);
        }

        printTimes();
        return res;
    }


    private void gurobiModelParameters(GRBModel model) throws GRBException {
        //MUST BE CORRECTED WITH min1, min2
        model.set(GRB.DoubleParam.NodefileStart, 0.5);
        model.set(GRB.IntParam.Method, 1);
        model.set(GRB.IntParam.NonConvex, 2);
        model.set(GRB.IntParam.Threads, 1);
        model.set(GRB.StringParam.LogFile, "");

        if (highestPrecision && mecDetected) {
            //Precision HIGH
            model.set(GRB.DoubleParam.FeasibilityTol, 1e-9);
            model.set(GRB.DoubleParam.IntFeasTol, 1e-9);
            model.set(GRB.DoubleParam.OptimalityTol, 1e-9);

            model.set(GRB.DoubleParam.BarQCPConvTol, 0);
            model.set(GRB.DoubleParam.BarConvTol, 0);
            model.set(GRB.DoubleParam.MarkowitzTol, 1e-4);
            model.set(GRB.DoubleParam.MIPGap, 0);
            model.set(GRB.DoubleParam.MIPGapAbs, 0);
            model.set(GRB.DoubleParam.PSDTol, 0);
        }

        if (verbosity < 1) {
            model.set(GRB.IntParam.LogToConsole, 0);
        }
    }

    private ModelCheckerResult throwCplexError() {
        System.out.println("WARNING: WE DISABLED CPLEX FOR NOW!! YOU WILL HAVE TO AQUIRE CPLEX " +
                "AND A LICENSE YOURSELF. FOLLOW THE GUIDELINES PROVIDED AT:\n" +
                "https://github.com/ga67vib/Algorithms-For-Stochastic-Games/blob/master/README.md#cplex");
        return null;
    }

    /**===================CPLEX-RELATED========================
    private void setCplexParams(IloCplex cplex) throws IloException {
        cplex.setParam(IloCplex.Param.OptimalityTarget, IloCplex.OptimalityTarget.FirstOrder);
        cplex.setParam(IloCplex.Param.MIP.Strategy.File, 3);
        cplex.setParam(IloCplex.Param.WorkMem, 40000);
        cplex.setParam(IloCplex.Param.TimeLimit,60*60*2);
        cplex.setParam(IloCplex.Param.Threads, 1);

        if (highestPrecision && mecDetected) {
            cplex.setParam(IloCplex.Param.Barrier.ConvergeTol, 1e-12);
            cplex.setParam(IloCplex.Param.Barrier.QCPConvergeTol, 1e-12);
            cplex.setParam(IloCplex.Param.MIP.Tolerances.AbsMIPGap, 0);
            cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 0);
            cplex.setParam(IloCplex.Param.MIP.Tolerances.Integrality, 0);
            cplex.setParam(IloCplex.Param.MIP.Tolerances.AbsMIPGap, 0);

            cplex.setParam(IloCplex.Param.Simplex.Tolerances.Markowitz, 0.01);
            cplex.setParam(IloCplex.Param.Simplex.Tolerances.Feasibility, 1e-9);
            cplex.setParam(IloCplex.Param.Benders.Tolerances.feasibilitycut, 1e-9);

        }
    }

    private ModelCheckerResult solveQPCplex() {
        ModelCheckerResult res = new ModelCheckerResult();
        try (IloCplex cplex = new IloCplex()) {
            cplex.importModel("model.lp");
            setCplexParams(cplex);

            try {
                FileOutputStream out = new FileOutputStream("cplexOutput");
                if (verbosity == 0 && logging) {
                    cplex.setOut(out);
                }
                else if (verbosity == 0 && !logging) {
                    cplex.setOut(null);
                }
            }
            catch (FileNotFoundException e) {
                log("Couldn't write to file");
            }

            buildTimeEnd = System.nanoTime();
            readWriteTimeEnd = System.nanoTime();
            solveTimeStart = System.nanoTime();
            if ( cplex.solve() ) {
                solveTimeEnd = System.nanoTime();
                IloLPMatrix lp = (IloLPMatrix)cplex.LPMatrixIterator().next();
                res.soln = cplex.getValues(lp);

                if (verbosity >= 1) {
                    System.out.println("CPLEX-OPTIMUM: "+res.soln[0]);
                    System.out.println("Maximum bound violation = " +
                            cplex.getQuality(IloCplex.QualityType.MaxPrimalInfeas));

                    System.out.println("Solution status = " + cplex.getStatus());
                    System.out.println("Solution value  = " + cplex.getObjValue());
                }
                //log("\nCPLEX-OPTIMUM: "+res.soln[0]);

                cplex.exportModel("qpex1.lp");
            }
            else {
                solveTimeEnd = System.nanoTime();
                log("Couldn't solve CPLEX");
            }
        }
        catch (IloException e) {
            log("Concert exception '" + e + "' caught");
            System.err.println("Concert exception '" + e + "' caught");
            System.exit(-1);
        }
        res.timeTaken = (solveTimeEnd - solveTimeStart) + (buildTimeEnd - buildTimeStart);
        return res;
    }
     **/

    private ModelCheckerResult computeReachProbsQPGurobiCondonStrict(STPG stpg, BitSet no, BitSet yes, boolean min1, boolean min2, double init[], BitSet known, boolean solveGRB)
            throws PrismException, GRBException {

        List<BitSet> mecs=getMECs((MDPSimple) stpg, no, yes);
        STPG condonSTPG = CondonNFTransformation.normalize(stpg, no, yes, mecs);

        computeYesAndNo(condonSTPG);

        return computeReachProbs(condonSTPG, no, yes, min1, min2, init,null, null, null, null, null, solveGRB);
    }

    private ModelCheckerResult computeReachProbsQPGurobiCondonLoose(STPG stpg, BitSet no, BitSet yes, boolean min1, boolean min2, double init[], BitSet known, ModelCheckerResult preRes, boolean solveGRB, boolean topological)
            throws PrismException, GRBException {

        STPG condonSTPG = CondonNFTransformation.normalizeLoose(stpg, no, yes);
        computeYesAndNo(condonSTPG);
        List<BitSet> mecs = getMECs((MDPSimple) stpg,no,yes);

        if (topological){
            return solveQPTopological(stpg,no,yes,min1,min2,init,mecs,known,preRes,solveGRB);
        }
        else {
            return computeReachProbs(condonSTPG, no, yes, min1, min2, init,mecs, null, null, null, preRes, solveGRB);
        }
    }

    private ModelCheckerResult computeReachProbsQPGurobiCondonLooseEps(STPG stpg, BitSet no, BitSet yes, boolean min1, boolean min2, double init[], BitSet known, boolean mStates, boolean solveGRB)
            throws PrismException, GRBException {
        List<BitSet> mecs=getMECs((MDPSimple) stpg, no, yes);
        STPGExplicit condonSTPG = (STPGExplicit) CondonNFTransformation.normalizeLooseAndStopping(stpg, no, yes,mecs, mStates);
        computeYesAndNo(condonSTPG);

        return computeReachProbs(condonSTPG, no, yes, min1, min2, init, null,null, null, null, null, solveGRB);
    }

    private ModelCheckerResult solveQPTopological(STPG stpg, BitSet no, BitSet yes, boolean min1, boolean min2, double init[], List<BitSet> mecs, BitSet known, ModelCheckerResult preRes, boolean solveGRB)
            throws PrismException, GRBException{

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
        ModelCheckerResult tmpRes = new ModelCheckerResult();
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
                tmpRes = computeReachProbs(stpg, no, yes, min1, min2, init, mecs, known, values, statesInSCC, preRes, solveGRB);
                i=0;
                for (int state = statesInSCC.nextSetBit(0); state>=0; state=statesInSCC.nextSetBit(state+1)) {
                    values[state]=tmpRes.soln[i];
                    i++;
                }
                known.or(statesInSCC);

                //We can pass on the same STPG, but just care for the interesting set
                //=> For not-topological just pass on everything
            }
        }

        buildAndSolveTime = System.nanoTime() - buildAndSolveTime;
        res.soln = values;
        res.timeTaken = buildAndSolveTime + tranformTime;
        return res;
    }

    /**
     * Sum over i out of Vmax and Vmin: (v(i)-v(k))(v(i)-v(j))
     * Expressed as (v(i)v(i)-v(i)v(j)-v(i)v(k)+v(j)v(k))
     * However, if count == 2 (action leads to AVG-node), you must use 1/2 * (v(j1)+v(j2))
     * so if j is avg node it could be (v(i)v(i)-1/2*v(i)v(j1)-1/2*v(i)v(j2)-v(i)v(k)+1/2*v(j1)v(k)+1/2*v(j2)v(k))
     * @param stpg
     * @param stateVars
     * @return
     */
    private GRBQuadExpr setObjectiveGurobi(STPG stpg, GRBVar[] stateVars, BitSet known, double[] knownValues, BitSet compute) {
        //Set the objective
        GRBQuadExpr objExpr = new GRBQuadExpr();
        int numStates=stpg.getNumStates();
        boolean tr1Known;
        boolean tr2Known;

        for (int state = compute.nextSetBit(0); state>=0; state=compute.nextSetBit(state+1)) {
            int choices = stpg.getNumChoices(state);

            if (choices == 2) {
                //v(i)v(i)
                objExpr.addTerm(1.0, stateVars[getIndex(state)], stateVars[getIndex(state)]);

                //Ugly. Would be prettier if Gorubi would allow to build QuadExpr through LinExpr*LinExpr
                for (Iterator<Map.Entry<Integer, Double>> it1 = stpg.getTransitionsIterator(state, 0); it1.hasNext(); ) {
                    Map.Entry<Integer, Double> tr1 = it1.next();
                    tr1Known = known.get(tr1.getKey());
                    //-1/2*v(i)v(j1)-1/2*v(i)v(j2) - 1/2 just as possible example. Handled by tr.getValue
                    if (tr1Known) {
                        objExpr.addTerm(-1.0 * tr1.getValue() * knownValues[tr1.getKey()], stateVars[getIndex(state)]);
                    }
                    else {
                        objExpr.addTerm(-1.0 * tr1.getValue(), stateVars[getIndex(state)], stateVars[getIndex(tr1.getKey())]);
                    }

                    for (Iterator<Map.Entry<Integer, Double>> it2 = stpg.getTransitionsIterator(state, 1); it2.hasNext(); ) {
                        Map.Entry<Integer, Double> tr2 = it2.next();
                        //+1/2*v(j1)v(k)+1/2*v(j2)v(k)
                        tr2Known = known.get(tr2.getKey());
                        if (!tr1Known && !tr2Known) {
                            objExpr.addTerm(1.0 * tr1.getValue() * tr2.getValue(), stateVars[getIndex(tr1.getKey())], stateVars[getIndex(tr2.getKey())]);
                        }
                        else if (tr1Known && !tr2Known) {
                            objExpr.addTerm(1.0 * tr1.getValue() * tr2.getValue() * knownValues[tr1.getKey()], stateVars[getIndex(tr2.getKey())]);
                        }
                        else if (!tr1Known && tr2Known) {
                            objExpr.addTerm(1.0 * tr1.getValue() * tr2.getValue() * knownValues[tr2.getKey()], stateVars[getIndex(tr1.getKey())]);
                        }
                        else {
                            objExpr.addConstant(1.0 * tr1.getValue() * tr2.getValue() * knownValues[tr1.getKey()] * knownValues[tr2.getKey()]);
                        }
                    }
                }
                for (Iterator<Map.Entry<Integer, Double>> it2 = stpg.getTransitionsIterator(state, 1); it2.hasNext(); ) {
                    Map.Entry<Integer, Double> tr2 = it2.next();
                    //-v(i)v(k)
                    tr2Known = known.get(tr2.getKey());
                    if (tr2Known) {
                        objExpr.addTerm(-1.0 * tr2.getValue() * knownValues[tr2.getKey()], stateVars[getIndex(state)]);
                    }
                    else {
                        objExpr.addTerm(-1.0 * tr2.getValue(), stateVars[getIndex(state)], stateVars[getIndex(tr2.getKey())]);
                    }
                }
            }
        }
        return objExpr;
    }

    /**
     * This function computes the reachability probabilities of given SG stpg
     * NOTE: Can only compute reachability objectives
     * Other probabilities like safety games may lead to problems
     *       Like: prob0 and prob1 give strange results
     * @param stpg
     * @param no
     * @param yes
     * @param min1
     * @param min2
     * @param init
     * @param mecs contains all the MECs
     * ==TOPOLOGICAL PARAMS==
     * @param known
     * @param knownValues
     * @param compute
     * ==WARM START PARAMS==
     * @param preRes
     * ==SOLUTION PARAMS==
     * @param solveGRB
     * @return
     * @throws PrismException
     * @throws GRBException
     */
    protected ModelCheckerResult computeReachProbs(STPG stpg, BitSet no, BitSet yes, boolean min1, boolean min2, double init[], List<BitSet> mecs,
                                                     BitSet known, double[] knownValues, BitSet compute, ModelCheckerResult preRes, boolean solveGRB)
            throws PrismException, GRBException{
        boolean topological=true;
        if (compute==null) {
            topological = false;
            compute = new BitSet();
            compute.set(0,stpg.getNumStates());
        }
        if (known==null) {
            known = new BitSet();
        }
        if (mecs==null) {
            mecs = new ArrayList<>();
        }

        int numStates = compute.cardinality();

        GRBEnv grbEnv = new GRBEnv();
        GRBModel quadraticProgram = new GRBModel(grbEnv);


        //=================CREATE VARIABLES===============================
        double lb[] = new double[numStates]; // also initializes array to 0
        double ub[] = new double[numStates];

        // Initialize all upper bounds to 1
        Arrays.fill(ub, 1.0);

        // Create state variables, type = null implies all variables are of continuous type
        GRBVar[] stateVars = quadraticProgram.addVars(lb, ub, null, null, null);
        //==================ADDING WARM STARTS FOR GUROBI==========================
        if (preRes != null) {
            double[] preResArr = preRes.soln;
            if (preResArr.length > stateVars.length && !topological) {
                throw new PrismException("More states in initial solution than in states");
            }
            else {
                if (preResArr.length < stateVars.length) {
                    //System.out.println("Initial solution has less states than QP model. Additional states won't have initial value");
                }
                int i=0;
                for (int state = compute.nextSetBit(0); state>=0; state=compute.nextSetBit(state+1)) {
                    if (state < preRes.soln.length) {
                        stateVars[i].set(GRB.DoubleAttr.Start, preResArr[state]);
                        i++;
                    }
                    else {
                        break; //most likely 2Act was applied and states is bigger now. Since 2Act appends after all existing states
                    }
                }
            }
        }

        //======================CREATE CONSTRAINTS=======================
        int numChoices = 0;
        if (topological) {
            for (int state=compute.nextSetBit(0); state>=0; state = compute.nextSetBit(state+1)) {
                numChoices+= stpg.getNumChoices(state);
            }
        }
        else {
            numChoices = stpg.getNumChoices();
        }
        GRBLinExpr[] lhs = new GRBLinExpr[numChoices];
        double[] rhs = new double[numChoices];
        char[] senses = new char[numChoices];

        //Initialize Constraint Matrix. state i, choice j, row k
        int k = 0;
        for (int i=compute.nextSetBit(0); i>=0; i = compute.nextSetBit(i+1)) {
            //skip all the states in MECs because the constraints are built differently
            boolean skip=false;
            for (BitSet bs : mecs) {
                if (bs.get(i)) {
                    if (stpg.getNumChoices(i)>1) { //We still need to implement equality constraints for the states without choice
                        skip = true;
                    }
                    break;
                }
            }
            if (skip) continue;

            int choices = stpg.getNumChoices(i);
            for (int j = 0; j < choices; j++) {
                // Add v(i)
                lhs[k] = new GRBLinExpr();
                lhs[k].addTerm(1.0, stateVars[getIndex(i)]);

                //putting the correct sense. ToDo HAS TO BE ADJUSTED TO min1, min2
                if (choices == 1) {
                    senses[k] = GRB.EQUAL;
                } //supposing 1 -> maxPlayer 2-> minplayer
                else if (stpg.getPlayer(i) == 1) {
                    senses[k] = GRB.GREATER_EQUAL;
                } else if (stpg.getPlayer(i) == 2) {
                    senses[k] = GRB.LESS_EQUAL;
                }

                //If this is a yes-state, just model it as a state with a self-loop
                if (yes.get(i)) {
                    lhs[k].addConstant(-1.0);
                    k++;
                    continue;
                }

                for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(i, j); it.hasNext(); ) {
                    Map.Entry<Integer, Double> tr = it.next();

                    //We aleardy know the value for this transition
                    if (known.get(tr.getKey())) {
                        lhs[k].addConstant(-1.0 * tr.getValue() * knownValues[tr.getKey()]);
                    }
                    else if (yes.get(tr.getKey())) { //The transition leads to a yes-sink
                        lhs[k].addConstant(-1.0 * tr.getValue());
                    } else if (!no.get(tr.getKey())) { //neither yes/no -> normal vertex
                        lhs[k].addTerm(-1.0 * tr.getValue(), stateVars[getIndex(tr.getKey())]);
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
            int size=0;
            for (GRBLinExpr lh : lhs) {
                if (lh != null) {
                    size++;
                }
            }
            lhsWithoutMEC = new GRBLinExpr[size];
            rhsWithoutMEC = new double[size];
            sensesWithoutMEC = new char[size];
            int j=0;
            for (int i=0; i<lhs.length; i++) {
                if (lhs[i]!=null) {
                    lhsWithoutMEC[j]=lhs[i];
                    rhsWithoutMEC[j]=rhs[i];
                    sensesWithoutMEC[j]=senses[i];
                    j++;
                }
            }
        }
        else {
            lhsWithoutMEC = lhs;
            rhsWithoutMEC = rhs;
            sensesWithoutMEC = senses;
        }

        // Add non-MEC-constaints to QP
        quadraticProgram.addConstrs(lhsWithoutMEC, sensesWithoutMEC, rhsWithoutMEC, null);

        //MEC-Solving - also adds the constraints
        SMGMECSolverGRB mecSolver= new SMGMECSolverGRB(quadraticProgram, stateVars, mecs, stpg, known, knownValues);
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

        //Set the objective
        GRBQuadExpr objExpr = setObjectiveGurobi(stpg, stateVars, known, knownValues, compute);

        quadraticProgram.setObjective(objExpr, GRB.MINIMIZE);
        gurobiModelParameters(quadraticProgram);

        buildTimeEnd = System.nanoTime();

        quadraticProgram.write("model.lp");

        if (!solveGRB) {
            readWriteTimeStart = System.nanoTime();
            /** ======================CPLEX-RELATED=================
            return solveQPCplex();
            **/
            return throwCplexError();
        }

        solveTimeStart = System.nanoTime();

        // Do the actual solving
        quadraticProgram.optimize();

        solveTimeEnd = System.nanoTime();

        int optimstatus = quadraticProgram.get(GRB.IntAttr.Status);

        if (optimstatus == GRB.Status.INFEASIBLE) {
            System.out.println("Model is infeasible");

            // Compute and write out IIS
            quadraticProgram.computeIIS();
            quadraticProgram.write("model.ilp");
        }

        // Process result and put it in res
        ModelCheckerResult res = new ModelCheckerResult();
        res.soln = new double[numStates];
        for (int i = 0; i < numStates; i++) {
            res.soln[i] = stateVars[i].get(GRB.DoubleAttr.X);
        }
        if (verbosity>=2) {
            System.out.println("SOLUTION FOR STATE 0: "+res.soln[0]);
        }
        res.timeTaken = (solveTimeEnd - solveTimeStart) + (buildTimeEnd - buildTimeStart);
        return res;
    }
}
