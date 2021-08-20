package explicit;

import prism.PrismException;

import java.util.*;

public class STPGModelAnalyser {
    STPGModelChecker modelChecker;

    public STPGModelAnalyser(STPGModelChecker modelChecker) {
        this.modelChecker = modelChecker;
    }

    protected ModelCheckerResult analyse_model(STPG stpg, BitSet no, BitSet yes, boolean min1, boolean min2, double init[], BitSet known)
            throws PrismException {

        System.out.println("You called the model analysis method which doesn't solve the problem, but tells you more about it!");

        // Store num states
        int n = stpg.getNumStates();

        // Determine set of states actually need to compute values for
        BitSet unknown = new BitSet();
        unknown.set(0, n);
        unknown.andNot(yes);
        unknown.andNot(no);

        System.out.println("NumStates, NumTargets, NumSinks, NumUnknown");

        log("Number of States: " + stpg.getNumStates());
        log("Number of Choices: " + stpg.getNumChoices());
        log("Number of Transition: " + stpg.getNumTransitions());
        log("Number of Targets (States with trivial value 1): " + yes.cardinality());
        log("Number of Sinks (States with trivial value 0): " + no.cardinality());
        log("Number of Unknown States: "+unknown.cardinality());

        int maxChoices=0;
        int maxTransitions=0;
        int choiceNum;
        double minTrans=1.0;
        int numChoicesWithBranches = 0;
        for (int state=0; state<stpg.getNumStates(); state++) {
            choiceNum=stpg.getNumChoices(state);
            if (choiceNum>maxChoices) maxChoices=choiceNum;
            for (int choice=0; choice<choiceNum; choice++) {
                if (stpg.getNumTransitions(state, choice) > 1) numChoicesWithBranches++;
                if (stpg.getNumTransitions(state,choice)>maxTransitions) maxTransitions=stpg.getNumTransitions(state,choice);
                for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext(); ) {
                    Map.Entry<Integer, Double> tr = it.next();
                    if (minTrans > tr.getValue()) {
                        minTrans = tr.getValue();
                    }
                }
            }
        }

        log("Number of maximal choices per state: " + maxChoices);
        log("Number of maximal transitions per choice: " + maxTransitions);
        log("Smallest transition probability: " + minTrans);
        log("Number of Choices with probability: " + numChoicesWithBranches);

        List<BitSet> mecs = null;
        ECComputerDefault ec =null;
        System.out.println("Getting MECs...");
        ec = (ECComputerDefault) ECComputer.createECComputer(this.modelChecker, stpg);
        //need a copy of unknown, since EC computation empties the set as a side effect
        BitSet unknownForEC = new BitSet();
        unknownForEC.or(unknown);
        ec.computeMECStates(unknownForEC);
        mecs = ec.getMECStates();
        log("Number of MECs: " + mecs.size());

        Collections.sort(mecs, (Comparator.comparingInt(BitSet::cardinality)));
        long maximalCardinalityMEC=0;
        long minimalCardinalityMEC=Long.MAX_VALUE;
        double avgMecsize = 0;

        double mecSizeMedian = 0;
        int medianIndex = mecs.size()/2;
        if (mecs.size() % 2 == 1) {
            mecSizeMedian = mecs.get(medianIndex).cardinality();
        }
        else {
            mecSizeMedian = (mecs.get(medianIndex).cardinality() + mecs.get(medianIndex+1).cardinality())/2.0;
        }

        if (mecs.size()>0) {
            for (int i=0; i<mecs.size(); i++) {
                avgMecsize += mecs.get(i).cardinality();
                if (mecs.get(i).cardinality()>maximalCardinalityMEC) maximalCardinalityMEC=mecs.get(i).cardinality();
                if (mecs.get(i).cardinality()<minimalCardinalityMEC) minimalCardinalityMEC=mecs.get(i).cardinality();
            }
            avgMecsize/=mecs.size();
            log("Biggest MEC has size: " + maximalCardinalityMEC);
            log("Smallest MEC has size: " + minimalCardinalityMEC);
            log("MEC size on average is: " + avgMecsize);
            log("MEC size median is: " + mecSizeMedian);
        }
        else {
            //Currently only there to make writing and reading csv easier
            log("Biggest MEC has size: " + "");
            log("Smallest MEC has size: " + "");
        }
        System.out.println(("I could also tell you the size of each MEC or more about it. Controlled MEC? SEC?"));

        SCCInfo sccs = null;
        System.out.println("Getting topologically ordered SCCs...");
        sccs = SCCComputer.computeTopologicalOrdering(this.modelChecker, stpg, true, unknown::get);
        System.out.println("Number of SCCs: " + sccs.getNumSCCs());

        System.out.println("I also want to tell you about the min/avg/max length of a path from init to a target.");
        System.out.println("Also about the number of probabilistic loops? For each SCC? The occurring probabilities?");

        // Return result that makes clear that it is not valid
        ModelCheckerResult res = new ModelCheckerResult();
        res.soln = new double[1];
        res.soln[0] = -1;
        res.numIters = -1;
        res.timeTaken = -1;
        return res;

        //TODO:
        /**
         * Remember the following ideas from discussion with Maxi, Muqsit and Tobi on 22.07.21:
         * - Top down: Analyse method, select method based on features
         * - Bottom up: Prepend some subgame/component (exhibiting certain features), check what kind of impact it has
         *
         * 	Afterwards, maybe even prove that a certain component/feature always has a certain impact.
         *
         * 	Regarding OVI: It is bad according to iters-measure, but this might be due to:
         * 		- Iters containing both progress and verification phase
         * 		- after failing a few times (once), the required precision to start verification phase is much smaller (half), which might explain why we need to work for so long. Maybe sth other than half is better? Investigate this hyperparameter. Maybe even make hyperparameter depend on features?
         *
         * 	One more idea: Detect "hard regions" for VI (lots of loops in a few states); solve hard regions by QP/SI, i.e. "inline" them, and then continue VI (or recursively look for more hard regions and inline them as well)
         *
         */
    }
    /**
     * Simple test program.
     */

    protected void log(String s) {
        System.out.println(s);
    }
}
