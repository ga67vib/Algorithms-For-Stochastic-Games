package explicit;

import common.IntSet;
import gurobi.GRBException;
import prism.PrismException;
import prism.PrismSettings;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.*;

/**
 * This class exists to solve two-player stoachstic reachability games with polynomial programming
 */
public class SMGPolyProgSolver {
    STPGModelChecker modelChecker;
    STPG stpg;
    BitSet no;
    BitSet yes;
    boolean min1;
    boolean min2;
    double init[];
    BitSet known;
    protected PrismSettings settings;
    protected int verbosity;
    BitSet remain;
    BitSet target;
    boolean logging;

    boolean highestPrecision;
    boolean mecDetected;

    long buildTimeStart = -1;
    long buildTimeEnd = -1;

    long solveTimeStart = -1;
    long solveTimeEnd = -1;

    long readWriteTimeStart = -1;
    long readWriteTimeEnd = -1;

    //TOPOLOGICAL
    HashMap<Integer, Integer> stpgIndexToSCCIndex;

    //Which approach solves which variant?
    protected final static int[] amplSolves={9,10,11};
    protected final static int[] gurobiSolves={0,1,2,3,4,5,6,7,8};

    public SMGPolyProgSolver(STPGModelChecker modelChecker, STPG stpg, BitSet no, BitSet yes, boolean min1,
                       boolean min2, double[] init, BitSet known, PrismSettings settings, int verbosity) {
        this.modelChecker = modelChecker;
        this.stpg = stpg;
        this.no = no;
        this.yes = yes;
        this.min1 = min1;
        this.min2 = min2;
        this.init = init;
        this.known = known;
        this.settings = settings;
        this.verbosity = 2; // ToDo: Should get fixed one day to better verbosity
        this.logging = false; //ToDo: Should be passed as a Parameter

        this.highestPrecision = false;
        this.mecDetected = false;
        stpgIndexToSCCIndex = null;
    }

    public ModelCheckerResult solve(int variant) throws GRBException, PrismException {
        ModelCheckerResult res=null;

        if (contains(amplSolves, variant)) {
            SMGPolyProgSolverAMPL solver = new SMGPolyProgSolverAMPL(modelChecker, stpg, no, yes, min1, min2, init, known, settings, verbosity);
            return solver.solve(variant);
        }
        else {
            SMGPolyProgSolverGurobi solver = new SMGPolyProgSolverGurobi(modelChecker, stpg, no, yes, min1, min2, init, known, settings, verbosity);
            return solver.solve(variant);
        }
    }

    private boolean contains(int[] arr, int i) {
        for (int j=0; j<arr.length; j++) {
            if (arr[j]==i) return true;
        }
        return false;
    }

    protected void printTimes() {
        String errormsg = "There occured an error measuring the time";
        String advicemsg = "Check whether both the start and end markers of time measurement were reached";
        long diff;
        double secs;
        if (verbosity >= 1) {
            if (buildTimeEnd>0 && buildTimeStart>0) {
                diff=buildTimeEnd - buildTimeStart;
                secs= diff/(1e9);
                //m= secs/60;
                System.out.println("Time needed to build:");
                System.out.println("ns:\t\t"+diff);
                System.out.println("secs:\t"+secs);
            }
            else {
                System.out.println(errormsg +" needed to build the model");
                System.out.println(advicemsg);
            }
            if (readWriteTimeStart>0 && readWriteTimeEnd>0) {
                diff=readWriteTimeEnd - readWriteTimeStart;
                secs= diff/(1e9);
                //m= secs/60;
                System.out.println("Time needed to read and write model:");
                System.out.println("ns:\t\t"+diff);
                System.out.println("secs:\t"+secs);
            }
            else {
                //For now I don't want to print anything out for this case
            }
            if (solveTimeEnd>0 && solveTimeStart>0) {
                diff=solveTimeEnd - solveTimeStart;
                secs= diff/(1e9);
                //m= secs/60;
                System.out.println("Time needed to solve:");
                System.out.println("ns:\t\t"+diff);
                System.out.println("secs:\t"+secs);
            }
            else {
                System.out.println(errormsg +" needed to solve the model");
                System.out.println(advicemsg);
            }
        }
    }

    protected List<BitSet> getMECs(MDPSimple stpg, BitSet no, BitSet yes)
            throws PrismException
    {

        //compute MECs one time, use the decomposition in every iteration; SECs still have to be recomputed
        explicit.ECComputerDefault ec = (ECComputerDefault) ECComputer.createECComputer(this.modelChecker, stpg);
        BitSet all = new BitSet();
        all.set(0,stpg.getNumStates());
        all.xor(yes);//ignore states that are already finished
        all.xor(no);
        ec.computeMECStates(all);
        List<BitSet> mecs = ec.getMECStates();

        if (!mecs.isEmpty()) {
            mecDetected = true;
        }
        return mecs;
    }

    protected String STPGAnalysis(STPG stpg) {
        String str="\n";
        str+="Number of States:" + stpg.getNumStates();
        str+="\nNumber of Choices: " + stpg.getNumChoices();
        str+="\nNumber of Transition: " + stpg.getNumTransitions();
        return str;
    }

    protected void logBeforeCondon(STPG stpg, BitSet no, BitSet yes) throws PrismException {
        String str="\n\nAnalysis of STPG before CondonNormalForm";
        str+=STPGAnalysis(stpg);
        List<BitSet> mecs=getMECs((MDPSimple) stpg, yes, no);
        str+="\nAmount of MECs: "+mecs.size();
        int cardMax=0;
        int cardMin=0;
        if (mecs.size()>0) {
            cardMin=Integer.MAX_VALUE;
            for (int i=0; i<mecs.size(); i++) {
                if (mecs.get(i).cardinality()>cardMax) cardMax=mecs.get(i).cardinality();
                if (mecs.get(i).cardinality()<cardMin) cardMin=mecs.get(i).cardinality();
            }
            str+=("\nBiggest MEC has size: "+cardMax);
            str+=("\nSmallest MEC has size: "+cardMin);
        }
        int maxChoices=0;
        int maxTransitions=0;
        int choiceNum;
        double minTrans=1.0;
        for (int state=0; state<stpg.getNumStates(); state++) {
            choiceNum=stpg.getNumChoices(state);
            if (choiceNum>maxChoices) maxChoices=choiceNum;
            for (int choice=0; choice<choiceNum; choice++) {
                if (stpg.getNumTransitions(state,choice)>maxTransitions) maxTransitions=stpg.getNumTransitions(state,choice);
                for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext(); ) {
                    Map.Entry<Integer, Double> tr = it.next();
                    if (minTrans > tr.getValue()) {
                        minTrans = tr.getValue();
                    }
                }
            }
        }
        str+="\nMaximum amount of choices: " + maxChoices;
        str+="\nMaximum amount of transitions: " + maxTransitions;
        str+="\nSmallest Transition: " + minTrans;
        log(str);
    }

    protected void logAfterCondon(STPG stpg, BitSet no, BitSet yes) throws PrismException{
        String str="\nAnalysis of STPG after CondonNormalForm";
        str+=STPGAnalysis(stpg);
        List<BitSet> mecs=getMECs((MDPSimple) stpg, yes, no);
        str+="\nAmount of MECs: "+mecs.size();
        int cardMax=0;
        int cardMin=0;
        if (mecs.size()>0) {
            cardMin=Integer.MAX_VALUE;
            for (int i=0; i<mecs.size(); i++) {
                if (mecs.get(i).cardinality()>cardMax) cardMax=mecs.get(i).cardinality();
                if (mecs.get(i).cardinality()<cardMin) cardMin=mecs.get(i).cardinality();
            }
            str+=("\nBiggest MEC has size: "+cardMax);
            str+=("\nSmallest MEC has size: "+cardMin);
        }
        log(str);
    }

    protected void log(String str) {
        try {
            if (logging) {
                BufferedWriter writer = new BufferedWriter(new FileWriter("log", true));
                writer.append(str);

                writer.close();
            }
        }
        catch (java.io.IOException e) {
            System.out.println("---------------Writing failed!!!-----------------");
        }
    }

    protected void computeYesAndNo(STPG stpg) {
        // Are we generating an optimal adversary?
        boolean genAdv = modelChecker.exportAdv || modelChecker.generateStrategy;
        if (modelChecker.precomp && modelChecker.prob0) {
            no = modelChecker.prob0(stpg, null, yes, min1, min2); //AS I DON'T UPDATE REMAIN THIS SHOULD BE SET TO NULL
        }
        if (modelChecker.precomp && modelChecker.prob1 && !genAdv) {
            yes = modelChecker.prob1(stpg, null, yes, min1, min2); //SAME HERE
        }
    }

    protected void initializeYesValues(double[] values, BitSet yes) {
        for (int state = yes.nextSetBit(0); state>=0; state=yes.nextSetBit(state+1)) {
            values[state]=1.0;
        }
        //No just stays actually 0 and is not necessary
    }

    protected int getIndex(int i) {
        if (stpgIndexToSCCIndex == null) {
            return i;
        }
        else {
            return stpgIndexToSCCIndex.get(i);
        }
    }



    /**
     * Compute the value of a state that builds an SCC alone. We expect every action transition from the state leading
     * to an already computed SCC!
     * If for whatever reason we have here a 1-state-MEC that wasn't found by yes and no,
     * it will also be set to 0 (which is correct)
     * @param stpg
     * @param state
     * @param values
     */
    protected void computeSingletonSCC(STPG stpg, int state, double[] values) {
        int numChoices = stpg.getNumChoices(state);

        double max=0.0; //can compute on the fly
        double min=1.0;
        double actualValue;
        for (int choice = 0; choice<numChoices; choice++) {
            actualValue = 0.0;
            for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext(); ) {
                Map.Entry<Integer, Double> tr = it.next();
                actualValue += tr.getValue() * values[tr.getKey()];
            }
            max = Math.max(actualValue, max);
            min = Math.min(actualValue, min);
        }

        if (stpg.getPlayer(state) == 1) {
            values[state] = max;
        }
        else {
            values[state] = min;
        }
    }


}
