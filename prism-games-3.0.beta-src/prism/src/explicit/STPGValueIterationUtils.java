package explicit;

import prism.PrismComponent;
import prism.PrismException;
import prism.PrismUtils;

import java.util.*;

public class STPGValueIterationUtils {

    /**
     * Computes a strategy from the currently given bounds on the values of the states.
     * @param stpg
     * @param yes
     * @param lowerBounds
     * @param upperBounds
     * @param relevantStates Constrains the analysis only to states set in the bitSet. If null, all states will be analysed
     * @param alreadyComputedValues If you want to use topological features: Which states have already been computed. If not, leave empty
     * @param alreadyComputedAttractorDistances If you want to use topological features: What are the states distances to any target in attractor. If not, leave empty
     * @return array [sigma, tau, distancesToAnyTargetInAttractor], where sigma and tau are of int[] and contain store at index i the strategy for state i. The distancesToAnyTargetInAttractor are useful for topological computations
     */
    public static int[][] computeStrategyFromBounds(STPG stpg, BitSet yes, double[] lowerBounds, double[] upperBounds, BitSet relevantStates, BitSet alreadyComputedValues, int[] alreadyComputedAttractorDistances) {
        int[] sigma = new int[stpg.getNumStates()];
        int[] tau = new int[stpg.getNumStates()];

        if (relevantStates == null) {
            relevantStates = new BitSet();
            relevantStates.set(0, stpg.getNumStates()-1);
        }
        ArrayList<Integer>[] allowedMaximizerActions = new ArrayList[stpg.getNumStates()];

        // Find allowd actions for maximizer and strategy for minimizer
        for (int state = relevantStates.nextSetBit(0); state >= 0; state = relevantStates.nextSetBit(state+1)) {
            boolean isMaximizerState = stpg.getPlayer(state) == 1;

            double[] referenceBound;

            // Play conservative
            if (isMaximizerState) {
                referenceBound = lowerBounds;
                allowedMaximizerActions[state] = new ArrayList<>();
            }
            else {
                referenceBound = upperBounds;
            }

            double bestChoiceValueMaximizer = 0;
            double bestChoiceValueMinimizer = 1.0;
            for (int choice = 0; choice < stpg.getNumChoices(state); choice++) {
                double choiceValue = 0;

                for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext(); ) {
                    Map.Entry<Integer, Double> tr = it.next();
                    choiceValue += tr.getValue() * referenceBound[tr.getKey()];
                }

                //Minimizer may take any optimal action
                if (!isMaximizerState && bestChoiceValueMinimizer > choiceValue) {
                    bestChoiceValueMinimizer = choiceValue;
                    tau[state] = choice;
                }

                //If this action is better than the last ones, start considerations from here
                else if (isMaximizerState && bestChoiceValueMaximizer < choiceValue) {
                    bestChoiceValueMaximizer = choiceValue;
                    allowedMaximizerActions[state].clear();
                    allowedMaximizerActions[state].add(choice);
                }
                //If it is equally good as your current actions, add it to your considerations
                else if (isMaximizerState && bestChoiceValueMaximizer == choiceValue) {
                    allowedMaximizerActions[state].add(choice);
                }
            }
        }

        int[] attractorDistances = computeAttractor(stpg, yes,allowedMaximizerActions, relevantStates, alreadyComputedValues, alreadyComputedAttractorDistances);

        // Take any action that truly reduces distance towards the targets
        for (int state = relevantStates.nextSetBit(0); state >= 0; state = relevantStates.nextSetBit(state+1)) {
            // Only Maximizer interesting
            if (stpg.getPlayer(state) != 1) continue;
            for (int choice : allowedMaximizerActions[state]) {
                boolean useThisAction = false;
                for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext(); ) {
                    Map.Entry<Integer, Double> tr = it.next();
                    if (attractorDistances[tr.getKey()] < attractorDistances[state]) {
                        useThisAction = true;
                        break;
                    }
                }

                if (useThisAction) {
                    sigma[state] = choice;
                    break;
                }
            }
        }

        return new int[][]{sigma, tau, attractorDistances};
    }

    /**
     * Computes Attractor as described in appendix A of paper https://arxiv.org/abs/1804.04901
     * Addtionally, may only compute attractor of a subset of states, given that they either lead only into relevantstates, targets or into already computed states.
     * If you want to compute the whole attractor, simply leave @relevantStates, @alreadyComputedStates and @alreadyComputedAttractorDistances empty
     * @param stpg
     * @param yes
     * @param optimalMaximizerActions
     * @param relevantStates
     * @param alreadyComputedStates
     * @param alreadyComputedAttractorDistances
     * @return the smallest distance to any target in the attractor. Index i is distance of state i to any target. If state i does not reach any target, then the value will be Integer.MAX_VALUE - 10
     */
    public static int[] computeAttractor(STPG stpg, BitSet yes, ArrayList<Integer>[] optimalMaximizerActions, BitSet relevantStates, BitSet alreadyComputedStates, int[] alreadyComputedAttractorDistances) {
        int[] attractorDistances = new int[stpg.getNumStates()];
        Arrays.fill(attractorDistances, Integer.MAX_VALUE - 10);

        if (relevantStates == null) {
            relevantStates = new BitSet();
            relevantStates.set(0, stpg.getNumStates()-1);
        }

        BitSet lastIteration = new BitSet();
        BitSet currentIteration = new BitSet();

        if (alreadyComputedStates == null) {
            alreadyComputedStates = new BitSet();
        }

        // Use already obtained results
        if (alreadyComputedAttractorDistances != null) {
            attractorDistances = alreadyComputedAttractorDistances;
        }
        else {
            for (int target = yes.nextSetBit(0); target >= 0; target = yes.nextSetBit(target + 1)) {
                attractorDistances[target] = 0;
            }
        }


        // Start with targets
        currentIteration.or(alreadyComputedStates);
        currentIteration.or(yes);

        while(!lastIteration.equals(currentIteration)) {
            // Set last to current
            lastIteration.or(currentIteration);
            currentIteration.clear();

            currentIteration.or(alreadyComputedStates);
            currentIteration.or(yes);

            for (int state = relevantStates.nextSetBit(0); state >= 0; state = relevantStates.nextSetBit(state+1)){
                boolean atLeastOneActionLeadingIntoLastIteration = false;
                boolean allActionsLeadingIntoLastIteration = true;
                boolean isMaximizerState = stpg.getPlayer(state) == 1;

                int nearestNextState = state;
                for (int choice = 0; choice < stpg.getNumChoices(state); choice++) {
                    if (isMaximizerState && !optimalMaximizerActions[state].contains(choice)) continue;

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

    public static double[] getValueFromStrategies(PrismComponent prismComponent, STPG stpg, int[] sigma, int[] tau, BitSet yes, BitSet no, BitSet relevantStates, BitSet alreadyComputedStates, double[] fixedValues, double precision) throws PrismException {
        /* This creates an whole MDP, but I suppose it's not necessary since we can create Local MDPs
        MDP maxFixedMDP = createMDPFromSTPG(stpg, sigma, true, yes, no, relevantStates, alreadyComputedStates, fixedValues);
        BitSet maxFixedNoBitSet = extendNoSet(no, maxFixedMDP, yes, relevantStates, alreadyComputedStates);
        MDP minFixedMDP = createMDPFromSTPG(stpg, tau, false, yes, no, relevantStates, alreadyComputedStates, fixedValues);
        BitSet minFixedNoBitSet = extendNoSet(no, minFixedMDP, yes, relevantStates, alreadyComputedStates);
         */

        if (fixedValues == null) {
            fixedValues = new double[stpg.getNumStates()];
        }
        LocalMDPResult maxFixedLocalMDP = createLocalMDPFromSTPG(stpg, sigma, true, yes, no, relevantStates, alreadyComputedStates, fixedValues);
        LocalMDPResult minFixedLocalMDP = createLocalMDPFromSTPG(stpg, tau, false, yes, no, relevantStates, alreadyComputedStates, fixedValues);



        MDPModelChecker mdpModelChecker = new MDPModelChecker(prismComponent);
        mdpModelChecker.maxIters = 10000;


        ModelCheckerResult maxFixedResult = mdpModelChecker.computeReachProbsLinearProgrammingGurobi(maxFixedLocalMDP.mdp, maxFixedLocalMDP.sinks, maxFixedLocalMDP.targets, true, null);
        ModelCheckerResult minFixedResult = mdpModelChecker.computeReachProbsLinearProgrammingGurobi(minFixedLocalMDP.mdp, minFixedLocalMDP.sinks, minFixedLocalMDP.targets, false, null);

        //ModelCheckerResult maxFixedResult = mdpModelChecker.computeReachProbsPolIter(maxFixedLocalMDP.mdp, maxFixedLocalMDP.sinks, maxFixedLocalMDP.targets, true, null);
        //ModelCheckerResult minFixedResult = mdpModelChecker.computeReachProbsPolIter(minFixedLocalMDP.mdp, minFixedLocalMDP.sinks, minFixedLocalMDP.targets, false, null);


        for (int state = relevantStates.nextSetBit(0); state >= 0; state = relevantStates.nextSetBit(state + 1)) {
            double maxFixedStateValue = maxFixedResult.soln[maxFixedLocalMDP.stpgStatesToMdpStates.get(state)];
            double minFixedStateValue = minFixedResult.soln[minFixedLocalMDP.stpgStatesToMdpStates.get(state)];
            if (!PrismUtils.doublesAreClose(maxFixedStateValue, minFixedStateValue, precision, true)) {
                throw new PrismException("If Minimizer has free choice "+state+" gets value "+maxFixedStateValue+
                        " but if Maximizer has free choice states gets value "+minFixedStateValue+
                        " ==> Topological value iteration failed");
            }
            fixedValues[state] = maxFixedStateValue;
        }


        return fixedValues;
    }

    /**
     * Creates a MDP that contains only the relevant states, a target and a sink.
     * MUST FULFILL FOLLOWING REQUIREMENT:
     * If a state s in relevantStates has a transition to another state s' then s' must be either in yes, no, relevantStates or alreadyComputedStates
     *
     * @param stpg
     * @param strategy
     * @param maxIsFixed
     * @param yes
     * @param no
     * @param relevantStates
     * @param alreadyComputedStates
     * @param fixedValues
     * @return
     * @throws PrismException
     */
    private static LocalMDPResult createLocalMDPFromSTPG(STPG stpg, int[] strategy, boolean maxIsFixed, BitSet yes, BitSet no, BitSet relevantStates, BitSet alreadyComputedStates, double[] fixedValues) throws PrismException{
        MDPSimple mdp = new MDPSimple();
        HashMap<Integer, Integer> stpgToMdp = new HashMap<>();
        HashMap<Integer, Integer> mdpToStpg = new HashMap<>();

        if (relevantStates == null) {
            relevantStates = new BitSet();
            relevantStates.set(0, stpg.getNumStates()-1);
        }

        for (int state = relevantStates.nextSetBit(0); state >= 0; state = relevantStates.nextSetBit(state + 1)) {
            int mdpState = mdp.addState();
            stpgToMdp.put(state, mdpState);
            mdpToStpg.put(mdpState, state);
        }

        int sink = mdp.addState();
        int target = mdp.addState();

        //Set self-loops for sink and target
        Distribution d = new Distribution();
        d.add(sink, 1.0);
        mdp.addChoice(sink, d);

        d = new Distribution();
        d.add(target, 1.0);
        mdp.addChoice(target, d);

        //Copy actions and fix those that have to be fixed
        for (int state = relevantStates.nextSetBit(0); state >= 0; state = relevantStates.nextSetBit(state + 1)) {
            boolean isMaximizerState = stpg.getPlayer(state) == 1;
            for (int choice = 0; choice < stpg.getNumChoices(state); choice++) {
                if ((isMaximizerState && maxIsFixed) || (!isMaximizerState && !maxIsFixed)) {
                    if (strategy[state] != choice) continue;
                }
                Distribution distribution = new Distribution();
                for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext(); ) {
                    Map.Entry<Integer, Double> tr = it.next();

                    int nextState = tr.getKey();
                    //Redirect action to new target
                    if (yes.get(nextState)) {
                        distribution.add(target, tr.getValue());
                    }
                    //Redirect action to new sink
                    else if (no.get(nextState)) {
                        distribution.add(sink, tr.getValue());
                    }
                    //Stay inside SCC
                    else if (relevantStates.get(nextState)) {
                        distribution.add(stpgToMdp.get(nextState), tr.getValue());
                    }
                    //Moving to a already fixed state -> Should already have a fixed value and thus can be replaced
                    else if (alreadyComputedStates.get(nextState)) {
                        distribution.add(target, fixedValues[nextState] * tr.getValue());
                    }
                    //
                    else {
                        throw new PrismException("State " + state + " has in action " + choice + " a transition to state " + nextState + ", which is not already computed. This method is only suited for topological cases, in which this case should never arrise.");
                    }
                }
                mdp.addChoice(stpgToMdp.get(state), distribution);
            }
        }

        // Take any Initial State, since it does not matter
        mdp.addInitialState(0);

        LocalMDPResult mdpResult = new LocalMDPResult();
        mdpResult.mdp = mdp;
        mdpResult.stpgStatesToMdpStates = stpgToMdp;
        mdpResult.targets = new BitSet();
        mdpResult.targets.set(target);
        mdpResult.sinks = new BitSet();
        mdpResult.sinks.set(sink);
        return mdpResult;
    }

    /**
     * Transforms whole STPG to MDP. States that are not relevant will be converted into sinks
     * States that are already computed will get one actions that yields exactly that value.
     * @param stpg
     * @param strategy
     * @param maxIsFixed
     * @param yes
     * @param no
     * @param relevantStates
     * @param alreadyComputedStates
     * @param fixedValues
     * @return
     */
    private static MDP createMDPFromSTPG(STPG stpg, int[] strategy, boolean maxIsFixed, BitSet yes, BitSet no, BitSet relevantStates, BitSet alreadyComputedStates, double[] fixedValues) {
        MDPSimple mdp = new MDPSimple();
        mdp.addStates(stpg.getNumStates());

        if (relevantStates == null) {
            relevantStates = new BitSet();
            relevantStates.set(0, stpg.getNumStates()-1);
        }

        int sink = mdp.addState();
        int target = yes.nextSetBit(0);

        //Copy actions and fix those that have to be fixed
        for (int state = relevantStates.nextSetBit(0); state >= 0; state = relevantStates.nextSetBit(state + 1)) {
            boolean isMaximizerState = stpg.getPlayer(state) == 1;
            for (int choice = 0; choice < stpg.getNumChoices(state); choice++) {
                if ((isMaximizerState && maxIsFixed) || (!isMaximizerState && !maxIsFixed)) {
                    if (strategy[state] != choice) continue;
                }
                Distribution d = getDistribution(stpg, state, choice);
                mdp.addChoice(state, d);
            }
        }

        //Set values for all states that have already computed values
        for (int state = alreadyComputedStates.nextSetBit(0); state >= 0; state = alreadyComputedStates.nextSetBit(state + 1)) {
            Distribution d = new Distribution();
            d.add(target, fixedValues[state]);
            if (fixedValues[state] < 1.0) d.add(sink, 1.0-fixedValues[state]);
            mdp.addChoice(state, d);
        }

        //Every other state should just have a self-loop since we don't care
        if (relevantStates.cardinality() != mdp.getNumStates()) {
            for (int state = 0; state <= mdp.getNumStates(); state++) {
                if (relevantStates.get(state) || alreadyComputedStates.get(state)) continue;

                Distribution d = new Distribution();
                d.add(state, 1.0);
                mdp.addChoice(state, d);
            }
        }

        for (Integer initialState : stpg.getInitialStates()) {
            mdp.addInitialState(initialState);
        }

        return mdp;
    }

    private static BitSet extendNoSet(BitSet no, MDP mdp, BitSet yes, BitSet relevantStates, BitSet alreadyComputedStates) {
        BitSet extendedBitSet = new BitSet();
        extendedBitSet.set(0, mdp.getNumStates());
        extendedBitSet.andNot(yes);
        extendedBitSet.andNot(relevantStates);
        extendedBitSet.andNot(alreadyComputedStates);
        return extendedBitSet;
    }

    private static Distribution getDistribution(STPG stpg, int state, int action) {
        Distribution d = new Distribution();
        for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, action); it.hasNext(); ) {
            Map.Entry<Integer, Double> tr = it.next();
            d.add(tr.getKey(), tr.getValue());
        }
        return d;
    }

    /**
     * Since Edge lengths are only 1, we can use BFS instead of Dijkstra
     * @param stpg
     * @return
     */
    public static int[] getDistanceToTargets(STPG stpg, BitSet yes, BitSet no) {
        int[] distances = new int[stpg.getNumStates()];
        int[] predecessor = new int[stpg.getNumStates()];

        BitSet done = new BitSet();
        // Queue of arrays with two entries [state, distance]
        LinkedList<int[]> queue = new LinkedList<>();
        for (Integer initialState : stpg.getInitialStates()) {
            queue.addLast(new int[]{initialState, 0});
            predecessor[initialState] = -1;
        }

        while(!queue.isEmpty() || done.cardinality() < stpg.getNumStates()) {
            int[] currentState = queue.poll();

            // If visited, no need to visit again
            if (done.get(currentState[0])) {
                continue;
            }
            distances[currentState[0]] = currentState[1];
            done.set(currentState[0]);

            // No need to go on from here
            if (yes.get(currentState[0])) {
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

    private static class LocalMDPResult {
        MDP mdp;
        BitSet targets;
        BitSet sinks;
        HashMap<Integer, Integer> stpgStatesToMdpStates;
    }
}
