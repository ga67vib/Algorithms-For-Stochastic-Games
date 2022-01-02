package explicit;

import common.IntSet;
import prism.Pair;
import prism.PrismException;

import java.util.*;

public class STPGModelAnalyser {
    STPGModelChecker modelChecker;

    public STPGModelAnalyser(STPGModelChecker modelChecker) {
        this.modelChecker = modelChecker;
    }

    protected ModelCheckerResult analyse_model(STPG stpg, BitSet no, BitSet yes, boolean min1, boolean min2, double init[],
                                               BitSet known, BitSet target, long timerProb0, long timerProb1)
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
        float avgTransitions=0;
        int choiceNum;
        double minTrans=1.0;
        int numChoicesWithBranches = 0;
        int maximizerStates = 0;
        int minimizerStates = 0;
        int backwardsTransitions = 0;

        int[] incomingTransitions = new int[stpg.getNumStates()];

        for (int state=0; state<stpg.getNumStates(); state++) {
            choiceNum=stpg.getNumChoices(state);
            if (stpg.getPlayer(state) == 1) {
                maximizerStates++;
            }
            else {
                minimizerStates++;
            }

            if (choiceNum>maxChoices) maxChoices=choiceNum;
            for (int choice=0; choice<choiceNum; choice++) {
                avgTransitions+=stpg.getNumTransitions(state, choice);
                if (stpg.getNumTransitions(state, choice) > 1) numChoicesWithBranches++;
                if (stpg.getNumTransitions(state,choice)>maxTransitions) maxTransitions=stpg.getNumTransitions(state,choice);
                for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext(); ) {
                    Map.Entry<Integer, Double> tr = it.next();
                    if (tr.getKey() < state) {
                        backwardsTransitions++;
                    }

                    if (minTrans > tr.getValue()) {
                        minTrans = tr.getValue();
                    }
                    incomingTransitions[tr.getKey()]++;
                }
            }
        }

        log("Prob0 Time in s: "+(timerProb0/1000.0));
        log("Prob1 Time in s: "+(timerProb1/1000.0));

        log("Number of maximal choices per state: " + maxChoices);
        log("Number of maximal transitions per choice: " + maxTransitions);
        log("Average Number of Choices per state: " + (float)stpg.getNumChoices() / (float)stpg.getNumStates());
        log("Average Number of Transitions per state: " + avgTransitions/stpg.getNumStates());
        log("Average Number of Transitions per choice: " + avgTransitions/stpg.getNumChoices());
        log("Smallest transition probability: " + minTrans);
        log("Number of Choices with probability: " + numChoicesWithBranches);
        log("Number of Maximizer States: " + maximizerStates);
        log("Number of Minimizer States: " + minimizerStates);
        log("Percentage of Backwards Transitions: "+ (float) backwardsTransitions / (float) stpg.getNumTransitions());

        float incomingTransitionAvg = 0;
        for (int i = 0; i < incomingTransitions.length; i++) {
            incomingTransitionAvg+=incomingTransitions[i];
        }
        incomingTransitionAvg/=stpg.getNumStates();
        log("Average of Incoming Transitions: " + incomingTransitionAvg);



        try {
        List<BitSet> mecs = null;
        ECComputerDefault ec =null;
        System.out.println("Getting MECs...");
        long timeBeforeMECCompute = System.currentTimeMillis();
        ec = (ECComputerDefault) ECComputer.createECComputer(this.modelChecker, stpg);
        //need a copy of unknown, since EC computation empties the set as a side effect
        BitSet unknownForEC = new BitSet();
        unknownForEC.or(unknown);
        ec.computeMECStates(unknownForEC);
        mecs = ec.getMECStates();
        long timeAfterMECCompute = System.currentTimeMillis();
        log("Number of MECs: " + mecs.size());
        log("Time to compute MECs (s): "+((timeAfterMECCompute-timeBeforeMECCompute)/1000.0));

        Collections.sort(mecs, (Comparator.comparingInt(BitSet::cardinality)));
        long maximalCardinalityMEC=0;
        long minimalCardinalityMEC=Long.MAX_VALUE;
        double avgMecsize = 0;

        double mecSizeMedian = 0;
        int medianIndex = mecs.size()/2;

        if (mecs.size()>0) {
            for (int i=0; i<mecs.size(); i++) {
                avgMecsize += mecs.get(i).cardinality();
                if (mecs.get(i).cardinality()>maximalCardinalityMEC) maximalCardinalityMEC=mecs.get(i).cardinality();
                if (mecs.get(i).cardinality()<minimalCardinalityMEC) minimalCardinalityMEC=mecs.get(i).cardinality();
            }
            avgMecsize/=mecs.size();
            if (mecs.size() % 2 == 1) {
                mecSizeMedian = mecs.get(medianIndex).cardinality();
            }
            else {
                mecSizeMedian = (mecs.get(medianIndex-1).cardinality() + mecs.get(medianIndex).cardinality())/2.0;
            }
            log("Biggest MEC has size: " + maximalCardinalityMEC);
            log("Smallest MEC has size: " + minimalCardinalityMEC);
            log("MEC size on average is: " + avgMecsize);
            log("MEC size median is: " + mecSizeMedian);
        }
        else {
            //Currently only there to make writing and reading csv easier
            log("Biggest MEC has size: " + "");
            log("Smallest MEC has size: " + "");
            log("MEC size on average is: " + "");
            log("MEC size median is: " + "");
        }
        System.out.println(("I could also tell you the size of each MEC or more about it. Controlled MEC? SEC?"));
        }
        catch (Exception e) {
            System.out.println("MEC computation failed");
            log("Number of MECs: " + "?");
            log("Biggest MEC has size: " + "?");
            log("Smallest MEC has size: " + "?");
            log("MEC size on average is: " + "?");
            log("MEC size median is: " + "?");
        }

        SCCInfo sccs = null;
        System.out.println("Getting topologically ordered SCCs...");
        long timeBeforeSCCComp = System.currentTimeMillis();
        sccs = SCCComputer.computeTopologicalOrdering(this.modelChecker, stpg, true, null);
        long timeAfterSCCComp = System.currentTimeMillis();
        System.out.println("Number of SCCs: " + sccs.getNumSCCs());
        System.out.println("SCC computation took (s): " + ((timeAfterSCCComp-timeBeforeSCCComp)/1000.0));
        int numSCCs = sccs.getNumSCCs();
        int numNonSingleton = numSCCs;
        long maximalCardinalitySCC=0;
        long minimalCardinalitySCC=Long.MAX_VALUE;
        long minimalCardinalitySCCNonSingleton=Long.MAX_VALUE;
        double avgSCCsizeNonSingleton = 0.0;
        double avgSCCsize = 0.0;
        if (numSCCs > 0) {
            for (int i = 0; i < sccs.getNumSCCs(); i++) {
                IntSet scc = sccs.getStatesForSCC(i);
                long cardinality = scc.cardinality();
                if (cardinality <= 1) numNonSingleton--;
                else {
                    minimalCardinalitySCCNonSingleton = Math.min(cardinality, minimalCardinalitySCCNonSingleton);
                    avgSCCsizeNonSingleton += cardinality;
                }
                maximalCardinalitySCC = Math.max(cardinality, maximalCardinalitySCC);
                minimalCardinalitySCC = Math.min(cardinality, minimalCardinalitySCC);
                avgSCCsize += cardinality;
            }
            avgSCCsize /= sccs.getNumSCCs();
        }
        else {
            maximalCardinalitySCC = 0;
            minimalCardinalitySCC = 0;
            avgSCCsize = 0;
        }

        if(numNonSingleton > 0) {
            avgSCCsizeNonSingleton /= numNonSingleton;
        }

        log("Biggest SCC has size: "+maximalCardinalitySCC);
        log("Smallest SCC has size: "+minimalCardinalitySCC);
        log("Average SCC has size: "+avgSCCsize);

        log("Number of non-Singleton SCCs: "+(numNonSingleton));
        log("Smallest non-Singleton SCC has size: "+(numNonSingleton > 0 ? minimalCardinalitySCCNonSingleton : 0));
        log("Average non-Singleton SCC has size: "+(avgSCCsizeNonSingleton));

        int maximumSCCDepth = getLongestSCCChain(stpg, sccs);
        log("Longest Chain of SCC has length: "+maximumSCCDepth);

        cycleFreeAnalysis(stpg, yes, no, target);

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

    protected void log(String s) {
        System.out.println(s);
    }

    protected void cycleFreeAnalysis(STPG stpg, BitSet yes, BitSet no, BitSet targets) {
        System.out.println("I also want to tell you about the min/avg/max length of a path from init to a target.");
        int initalState = stpg.getFirstInitialState();
        ArrayList<Integer> sortedTargetCyclefreeDistances = new ArrayList<>();
        long nearestTarget = Integer.MAX_VALUE;
        long furthestTarget = 0;
        double avgCyclefreeDistance = 0.0;
        double medianCyclefreeDistance = 0.0;

        for (int target = targets.nextSetBit(0); target >= 0; target = targets.nextSetBit(target + 1)) {
            BitSet onlyOneTarget = new BitSet();
            onlyOneTarget.set(target);
            int[] attractorDistances = computeAttractor(stpg, onlyOneTarget);
            nearestTarget = Math.min(nearestTarget, attractorDistances[initalState]);

            //Prevent unreachable targets
            if (attractorDistances[initalState] <= stpg.getNumStates()) {
                furthestTarget = Math.max(furthestTarget, attractorDistances[initalState]);
                avgCyclefreeDistance += attractorDistances[initalState];
                sortedTargetCyclefreeDistances.add(attractorDistances[initalState]);
            }
        }

        // Get Median
        int medianCyclefreeDistanceIndex = sortedTargetCyclefreeDistances.size()/2;
        Collections.sort(sortedTargetCyclefreeDistances);
        if (sortedTargetCyclefreeDistances.size() % 2 == 1) {
            medianCyclefreeDistance = sortedTargetCyclefreeDistances.get(medianCyclefreeDistanceIndex);
        }
        else {
            try {
                medianCyclefreeDistance = (sortedTargetCyclefreeDistances.get(medianCyclefreeDistanceIndex) + sortedTargetCyclefreeDistances.get(medianCyclefreeDistanceIndex + 1)) / 2.0;
            }
            catch (Exception e) {
                medianCyclefreeDistance = 0;
            }
        }

        log("Nearest Target from any Initial State: "+nearestTarget);
        log("Furthest Target from any Initial State: "+furthestTarget);
        log("Target-distance Average: "+avgCyclefreeDistance/sortedTargetCyclefreeDistances.size());
        log("Target-distance Median: "+medianCyclefreeDistance);
    }

    /**
     * Computes Attractor as described in appendix A of paper https://arxiv.org/abs/1804.04901
     * Addtionally, may only compute attractor of a subset of states, given that they either lead only into relevantstates, targets or into already computed states.
     * If you want to compute the whole attractor, simply leave @relevantStates, @alreadyComputedStates and @alreadyComputedAttractorDistances empty
     * @param stpg
     * @param yes
     * @return the smallest distance to any target in the attractor. Index i is distance of state i to any target. If state i does not reach any target, then the value will be Integer.MAX_VALUE - 10
     */
    public static int[] computeAttractor(STPG stpg, BitSet yes) {
        int[] attractorDistances = new int[stpg.getNumStates()];
        Arrays.fill(attractorDistances, Integer.MAX_VALUE - 10);

        BitSet lastIteration = new BitSet();
        BitSet currentIteration = new BitSet();

        // Use already obtained results
        for (int target = yes.nextSetBit(0); target >= 0; target = yes.nextSetBit(target + 1)) {
            attractorDistances[target] = 0;
        }

        // Start with targets
        currentIteration.or(yes);

        while(!lastIteration.equals(currentIteration)) {
            // Set last to current
            lastIteration.or(currentIteration);
            currentIteration.clear();

            currentIteration.or(yes);

            for (int state = 0; state < stpg.getNumStates(); state++){
                boolean atLeastOneActionLeadingIntoLastIteration = false;
                boolean allActionsLeadingIntoLastIteration = true;
                boolean isMaximizerState = stpg.getPlayer(state) == 1;

                // Attractor sets are monotonic. We can use this as shortcut
                if (lastIteration.get(state)) {
                    currentIteration.set(state);
                    continue;
                }

                int nearestNextState = state;
                for (int choice = 0; choice < stpg.getNumChoices(state); choice++) {
                    // Can this action reach the last iteration?
                    boolean actionLeadsIntoLastIteration = false;
                    for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext(); ) {
                        Map.Entry<Integer, Double> tr = it.next();
                        if (!(tr.getValue() > 0)) continue;
                        int nextState = tr.getKey();

                        if (lastIteration.get(nextState)) {
                            actionLeadsIntoLastIteration = true;
                            if (attractorDistances[nearestNextState] == Integer.MAX_VALUE - 10 || attractorDistances[nextState] + 1 > attractorDistances[nearestNextState]) {
                                nearestNextState = nextState;
                            }
                            break;
                        }
                    }

                    if (actionLeadsIntoLastIteration) {
                        atLeastOneActionLeadingIntoLastIteration = true;
                    } else {
                        allActionsLeadingIntoLastIteration = false;
                    }
                }
                // Should state get added to iteration?
                if (isMaximizerState && atLeastOneActionLeadingIntoLastIteration) {
                    currentIteration.set(state);
                    if (attractorDistances[nearestNextState] + 1 < attractorDistances[state])
                        attractorDistances[state] = attractorDistances[nearestNextState] + 1;
                }
                else if(!isMaximizerState && allActionsLeadingIntoLastIteration) {
                    currentIteration.set(state);
                    if (attractorDistances[nearestNextState] + 1 < attractorDistances[state])
                        attractorDistances[state] = attractorDistances[nearestNextState] + 1;
                }
            }
        }

        return attractorDistances;
    }

    /**
     * Do a BFS over each SCC from SCC of initialState and see what the maximum from there is
     * @param stpg
     * @param sccInfo
     * @return
     */
    public int getLongestSCCChain(STPG stpg, SCCInfo sccInfo) {
        int maximumSCCChainLength = 1;

        PriorityQueue<int[]> priorityQueue = new PriorityQueue<>(Comparator.comparingInt(a -> a[1]));
        HashMap<Integer, Integer> sccIndexToMaximumLength = new HashMap<>();
        int currentSCCIndex = sccInfo.getSCCIndex(stpg.getFirstInitialState() != -1 ? stpg.getFirstInitialState() : 0);

        sccIndexToMaximumLength.put(currentSCCIndex, 1);

        IntSet currentSCC;

        LinkedList<Integer> sccIndexQueue = new LinkedList<>();
        sccIndexQueue.add(currentSCCIndex);
        while (!sccIndexQueue.isEmpty()) {
            currentSCCIndex = sccIndexQueue.poll();
            currentSCC = sccInfo.getStatesForSCC(currentSCCIndex);
            int chainLengthForCurrentSCC = sccIndexToMaximumLength.get(currentSCCIndex);

            // Go over every state of the SCC and see which other SCC it can reach
            for (Integer state : currentSCC) {
                for (int choice = 0; choice < stpg.getNumChoices(state); choice++) {
                    for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext(); ) {
                        Map.Entry<Integer, Double> tr = it.next();
                        int sccIndexOfNextState = sccInfo.getSCCIndex(tr.getKey());

                        // If we discovered a new SCC, add it to the Queue and set it's length to + 1 of where we are now
                        if (!sccIndexToMaximumLength.containsKey(sccIndexOfNextState)) {
                            sccIndexToMaximumLength.put(sccIndexOfNextState, chainLengthForCurrentSCC + 1);
                            maximumSCCChainLength = Math.max(maximumSCCChainLength, chainLengthForCurrentSCC + 1);
                            sccIndexQueue.add(sccIndexOfNextState);
                        }
                    }
                }
            }
        }

        return maximumSCCChainLength;
    }

    /**
     * Since Edge lengths are only 1, we can use BFS instead of Dijkstra
     * Note that this may be very incorrect since here we a assume that a minimizer would take an action leading straight to a target
     * @param stpg
     * @return
     */
    protected int[] getDistanceFromInitialState(STPG stpg, BitSet yes, BitSet no) {
        int[] distances = new int[stpg.getNumStates()];
        int[] predecessor = new int[stpg.getNumStates()];

        BitSet done = new BitSet();
        // Queue of arrays with two entries [state, distance]
        LinkedList<int[]> queue = new LinkedList<>();
        for (Integer initialState : stpg.getInitialStates()) {
            queue.addLast(new int[]{initialState, 0});
            predecessor[initialState] = -1;
        }

        while(!queue.isEmpty() && done.cardinality() < stpg.getNumStates()) {
            int[] currentState = queue.poll();

            // If visited, no need to visit again
            if (done.get(currentState[0])) {
                continue;
            }
            distances[currentState[0]] = currentState[1];
            done.set(currentState[0]);

            // No need to go on from here
            if (yes.get(currentState[0]) || no.get(currentState[0])) {
                continue;
            }

            for (int choice = 0; choice < stpg.getNumChoices(currentState[0]); choice++) {
                for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(currentState[0], choice); it.hasNext(); ) {
                    Map.Entry<Integer, Double> tr = it.next();
                    //Since edges have only length 1, the done question
                    if (tr.getValue() > 0 && !done.get(tr.getKey())) queue.addLast(new int[]{tr.getKey(), currentState[1] + 1});
                }
            }
        }


        return distances;
    }
}
