package explicit;

import prism.PrismComponent;
import prism.PrismException;
import prism.PrismUtils;

import java.util.*;

public class STPGValueIterationUtils {

    public static long dtmcBuildingTime = 0;
    public static long dtmcSolvingTime = 0;
    public static long inverseCalcTime = 0;
    public static long getAllowedActionsTime = 0;
    public static long verificationTime = 0;
    public static long assignmentTime = 0;

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
    public static int[][] computeStrategyFromBounds(STPG stpg, BitSet yes, double[] lowerBounds, double[] upperBounds, double precision, BitSet relevantStates, BitSet alreadyComputedValues, int[] alreadyComputedAttractorDistances) {
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
                //If it is equally good as your current actions, add it to your considerations
                else if (isMaximizerState && PrismUtils.doublesAreClose(bestChoiceValueMaximizer, choiceValue, precision, true)) {
                    allowedMaximizerActions[state].add(choice);
                }
                //If this action is better than the last ones, start considerations from here
                else if (isMaximizerState && bestChoiceValueMaximizer < choiceValue) {
                    bestChoiceValueMaximizer = choiceValue;
                    allowedMaximizerActions[state].clear();
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

    public static double[] getValueFromStrategiesWithoutAttractorWithLP(PrismComponent prismComponent, STPG stpg, BitSet yes, BitSet no, BitSet relevantStates, BitSet alreadyComputedStates, double[] fixedValues, double precision, double[] lowerBounds, double[] upperBounds) throws PrismException {
        // Compute all Action Possibilities
        if (relevantStates == null) {
            relevantStates = new BitSet();
            relevantStates.set(0, stpg.getNumStates()-1);
        }

        long t1;
        long t2;

        t1 = System.currentTimeMillis();
        ArrayList<Integer>[] allowedMaximizerActions = new ArrayList[stpg.getNumStates()];
        ArrayList<Integer>[] allowedMinimizerActions = new ArrayList[stpg.getNumStates()];
        //int[] allowedMinimizerActions = new int[stpg.getNumStates()];

        // Important hack for VI and OVI
        if (upperBounds == null) {
            upperBounds = lowerBounds;
        }
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
                allowedMinimizerActions[state] = new ArrayList<>();
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

                    allowedMinimizerActions[state].clear();
                    allowedMinimizerActions[state].add(choice);
                }
                //If it is equally good as your current actions, add it to your considerations
                else if (isMaximizerState && PrismUtils.doublesAreClose(bestChoiceValueMaximizer, choiceValue, precision, true)) {
                    allowedMaximizerActions[state].add(choice);
                }
                //If this action is better than the last ones, start considerations from here
                else if (isMaximizerState && bestChoiceValueMaximizer < choiceValue) {
                    bestChoiceValueMaximizer = choiceValue;
                    allowedMaximizerActions[state].clear();
                    allowedMaximizerActions[state].add(choice);
                }
            }
        }
        t2 = System.currentTimeMillis();
        getAllowedActionsTime += (t2-t1);

        if (fixedValues == null) {
            fixedValues = new double[stpg.getNumStates()];
        }
        t1 = System.currentTimeMillis();
        LocalMDPResult maxFixedLocalMDP = createLocalMDPFromSTPG(stpg, allowedMaximizerActions, true, yes, no, relevantStates, alreadyComputedStates, fixedValues, lowerBounds);
        LocalMDPResult minFixedLocalMDP = createLocalMDPFromSTPG(stpg, allowedMinimizerActions, false, yes, no, relevantStates, alreadyComputedStates, fixedValues, upperBounds);
        t2 = System.currentTimeMillis();
        dtmcBuildingTime += (t2-t1);

        MDPModelChecker mdpModelChecker = new MDPModelChecker(prismComponent);

        ModelCheckerResult maxFixedResult = mdpModelChecker.computeReachProbsLinearProgrammingGurobi(maxFixedLocalMDP.mdp, maxFixedLocalMDP.sinks, maxFixedLocalMDP.targets, true, null, maxFixedLocalMDP.initVector);
        ModelCheckerResult minFixedResult = mdpModelChecker.computeReachProbsLinearProgrammingGurobi(minFixedLocalMDP.mdp, minFixedLocalMDP.sinks, minFixedLocalMDP.targets, false, null, minFixedLocalMDP.initVector);

        for (int state = relevantStates.nextSetBit(0); state >= 0; state = relevantStates.nextSetBit(state+1)) {
            double valInMaxFixed = maxFixedResult.soln[maxFixedLocalMDP.stpgStatesToMdpStates.get(state)];
            double valInMinFixed = minFixedResult.soln[minFixedLocalMDP.stpgStatesToMdpStates.get(state)];
            if (!PrismUtils.doublesAreClose(valInMaxFixed, valInMinFixed, precision, true)) {
                throw new PrismException("State " + state + " would get in the MaximizerMDP value " + valInMinFixed + " but " + valInMaxFixed + " in the MinimizerMDP." +
                        "Thus, topological failed. Bound of the state in the STPG: [" + lowerBounds[state] + ", " + upperBounds[state] + "]");
            }
            fixedValues[state] = valInMaxFixed;
        }

        return fixedValues;
    }


    public static double[] getValueFromStrategiesWithoutAttractor(PrismComponent prismComponent, STPG stpg, BitSet yes, BitSet no, BitSet relevantStates, BitSet alreadyComputedStates, double[] fixedValues, double precision, double[] lowerBounds, double[] upperBounds) throws PrismException {
        // Compute all Action Possibilities
        if (relevantStates == null) {
            relevantStates = new BitSet();
            relevantStates.set(0, stpg.getNumStates()-1);
        }

        long t1;
        long t2;

        t1 = System.currentTimeMillis();
        ArrayList<Integer>[] allowedMaximizerActions = new ArrayList[stpg.getNumStates()];
        int[] allowedMinimizerActions = new int[stpg.getNumStates()];

        // Important hack for VI and OVI
        if (upperBounds == null) {
            upperBounds = lowerBounds;
        }
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

                    allowedMinimizerActions[state] = (choice);
                }
                //If it is equally good as your current actions, add it to your considerations
                else if (isMaximizerState && PrismUtils.doublesAreClose(bestChoiceValueMaximizer, choiceValue, precision, true)) {
                    allowedMaximizerActions[state].add(choice);
                }
                //If this action is better than the last ones, start considerations from here
                else if (isMaximizerState && bestChoiceValueMaximizer < choiceValue) {
                    bestChoiceValueMaximizer = choiceValue;
                    allowedMaximizerActions[state].clear();
                    allowedMaximizerActions[state].add(choice);
                }
            }
        }
        t2 = System.currentTimeMillis();
        getAllowedActionsTime += (t2-t1);

        if (fixedValues == null) {
            fixedValues = new double[stpg.getNumStates()];
        }
        t1 = System.currentTimeMillis();
        LocalDTMCTransformation suggestedLocalDTMC = createLocalDTMCFromStpg(
                stpg,
                allowedMaximizerActions, allowedMinimizerActions,
                yes, no,
                relevantStates,
                alreadyComputedStates, fixedValues,
                lowerBounds, upperBounds);
        t2 = System.currentTimeMillis();
        dtmcBuildingTime += (t2-t1);

        // Solve the created MarkovChain
        t1 = System.currentTimeMillis();
        DTMCNonIterativeSolutionMethods dtmcNonIterativeSolutionMethods = new DTMCNonIterativeSolutionMethods();
        ModelCheckerResult dtmcSolution = dtmcNonIterativeSolutionMethods.solveMarkovChain(suggestedLocalDTMC.dtmc, suggestedLocalDTMC.targets, suggestedLocalDTMC.upperboundsInDTMC, precision);
        t2 = System.currentTimeMillis();
        inverseCalcTime = DTMCNonIterativeSolutionMethods.inverseCalcTime;
        dtmcSolvingTime += (t2-t1);

        // We now have supposedly correct values for every state in relevantStates - Try to confirm them
        t1 = System.currentTimeMillis();
        if (!isDTMCResultConsistentWithSTPG(
                stpg,
                yes, no,
                relevantStates,
                alreadyComputedStates, fixedValues,
                suggestedLocalDTMC.stpgStatesToDtmcStates, suggestedLocalDTMC.dtmcStatesToStpgStates, dtmcSolution.soln,
                allowedMaximizerActions, allowedMinimizerActions,
                precision //The precision must be very close to the suggested Value
            )
        ) {
            throw new PrismException("Not consistent!!");
        }
        t2 = System.currentTimeMillis();
        verificationTime += (t2-t1);

        t1 = System.currentTimeMillis();
        // Set Values to the ones computed in the DTMC, since it's the values that should correspond to the optimal strategy
        for (int state = relevantStates.nextSetBit(0); state >= 0; state = relevantStates.nextSetBit(state+1)) {
            fixedValues[state] = dtmcSolution.soln[suggestedLocalDTMC.stpgStatesToDtmcStates.get(state)];
        }
        t2 = System.currentTimeMillis();
        assignmentTime += (t2-t1);

        return fixedValues;
    }

    private static boolean isDTMCResultConsistentWithSTPG(STPG stpg,
                                                          BitSet yes, BitSet no,
                                                          BitSet relevantStates,
                                                          BitSet alreadyComputedStates, double[] fixedValues,
                                                          HashMap<Integer, Integer> stpgToDtmc, HashMap<Integer, Integer> dtmcToStpg, double[] suggestedSGValue,
                                                          ArrayList<Integer>[] suggestedMaximizerActions, int[] suggestedMinimizerActions,
                                                          double precision)
    throws PrismException {

        for (int state = relevantStates.nextSetBit(0); state >= 0; state = relevantStates.nextSetBit(state+1)) {
            boolean isMaximizerState = stpg.getPlayer(state) == 1;
            // Start at worst-possible outcome as value
            double bestChoiceValue = 0;
            if (!isMaximizerState) bestChoiceValue = 1.0f;

            int bestChoice = -1;

            ArrayList<Double> choiceValues = new ArrayList<>();

            for (int choice = 0; choice < stpg.getNumChoices(state); choice++) {
                double choiceValue = 0;

                for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext(); ) {
                    Map.Entry<Integer, Double> tr = it.next();
                    int nextState = tr.getKey();
                    if (yes.get(nextState)) {
                        choiceValue+=tr.getValue();
                    }
                    else if (no.get(nextState)) {
                        // Nothing - value is 0 anyways
                    }
                    else if (relevantStates.get(nextState)) {
                        // Add value that state gets in DTMC - this should be the true value
                        choiceValue+= tr.getValue() * suggestedSGValue[stpgToDtmc.get(nextState)];
                    }
                    // Not part of the SCC -> must be already fixed
                    else if (alreadyComputedStates.get(nextState)){
                        choiceValue += tr.getValue() * fixedValues[nextState];
                    }
                    else {
                        throw new PrismException("State "+state+" has choice "+choice+ " that leads to state "+nextState+" which is neither computed already nor part of " +
                                "the current SCC. This method requires topological computation, and this this case should never happen.");
                    }
                }

                if (
                        (isMaximizerState && bestChoiceValue < choiceValue) ||
                        (!isMaximizerState && bestChoiceValue > choiceValue) ||
                                bestChoice == -1
                ) {
                    bestChoiceValue = choiceValue;
                    bestChoice = choice;
                }
                choiceValues.add(choiceValue);
            }

            ArrayList<Integer> argMaxChoices = new ArrayList<>();
            for (int choice = 0; choice < stpg.getNumChoices(state); choice++) {
                if (choice == bestChoice || PrismUtils.doublesAreClose(choiceValues.get(choice), bestChoiceValue, precision, true)) {
                    argMaxChoices.add(choice);
                }
            }

            if (isMaximizerState) {
                for (int suggestedAction : suggestedMaximizerActions[state]) {
                    if (!argMaxChoices.contains(suggestedAction)) {
                        throw new PrismException("State "+state+" was suspected to get value "+bestChoiceValue+" by playing "+
                                (isMaximizerState ? "an action from set "+suggestedMaximizerActions[state].toString() : "action "+suggestedMinimizerActions[state])+
                                " but "+suggestedAction+" would only yield value "+choiceValues.get(suggestedAction)+", which is not optimal." +
                                " Thus, the TOP heuristic failed. Try another solution method..");
                    }
                }
            }
        }
        return true;
    }

    /*
    Uses half-baked STPG -> MDP -> DTMC approach
    public static double[] getValueFromStrategiesWithoutAttractor(PrismComponent prismComponent, STPG stpg, BitSet yes, BitSet no, BitSet relevantStates, BitSet alreadyComputedStates, double[] fixedValues, double precision, double[] lowerBounds, double[] upperBounds) throws PrismException {
        // Compute all Action Possibilities
        if (relevantStates == null) {
            relevantStates = new BitSet();
            relevantStates.set(0, stpg.getNumStates()-1);
        }
        ArrayList<Integer>[] allowedMaximizerActions = new ArrayList[stpg.getNumStates()];
        int[] allowedMinimizerActions = new int[stpg.getNumStates()];

        // Find allowd actions for maximizer and strategy for minimizer
        for (int state = relevantStates.nextSetBit(0); state >= 0; state = relevantStates.nextSetBit(state+1)) {
            if (state == 55) {
                System.out.println("State 55 here!");
            }
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

                    allowedMinimizerActions[state] = (choice);
                }
                //If it is equally good as your current actions, add it to your considerations
                else if (isMaximizerState && PrismUtils.doublesAreClose(bestChoiceValueMaximizer, choiceValue, precision, true)) {
                    allowedMaximizerActions[state].add(choice);
                }
                //If this action is better than the last ones, start considerations from here
                else if (isMaximizerState && bestChoiceValueMaximizer < choiceValue) {
                    bestChoiceValueMaximizer = choiceValue;
                    allowedMaximizerActions[state].clear();
                    allowedMaximizerActions[state].add(choice);
                }
            }
        }

        if (fixedValues == null) {
            fixedValues = new double[stpg.getNumStates()];
        }
        LocalMDPResult maxFixedLocalMDP = createLocalMDPFromSTPG(stpg, allowedMaximizerActions, true, yes, no, relevantStates, alreadyComputedStates, fixedValues, lowerBounds);
        LocalMDPResult minFixedLocalMDP = createLocalMDPFromSTPG(stpg, allowedMinimizerActions, false, yes, no, relevantStates, alreadyComputedStates, fixedValues, upperBounds);

        MDPModelChecker mdpModelChecker = new MDPModelChecker(prismComponent);
        mdpModelChecker.termCritParam = precision;
        mdpModelChecker.maxIters = 10000;

        boolean useLP = false;

        ModelCheckerResult maxFixedResult;
        ModelCheckerResult minFixedResult;

        if (useLP) {
            maxFixedResult = mdpModelChecker.computeReachProbsLinearProgrammingGurobi(maxFixedLocalMDP.mdp, maxFixedLocalMDP.sinks, maxFixedLocalMDP.targets, true, null, maxFixedLocalMDP.initVector);
            minFixedResult = mdpModelChecker.computeReachProbsLinearProgrammingGurobi(minFixedLocalMDP.mdp, minFixedLocalMDP.sinks, minFixedLocalMDP.targets, false, null, minFixedLocalMDP.initVector);

        }
        else {
            int[] sigmaMDPRecommendation = new int[minFixedLocalMDP.mdp.getNumStates()];
            int[] tauMDPRecommendation = new int[maxFixedLocalMDP.mdp.getNumStates()];

            for (int state = relevantStates.nextSetBit(0); state >= 0; state = relevantStates.nextSetBit(state + 1)) {
                // Max player
                if (stpg.getPlayer(state) == 1) {
                    // Always assign first allowed action since by copying and merging actions there may be less now, but the one at index 0
                    // will also have the same index in the reduced MDP
                    sigmaMDPRecommendation[minFixedLocalMDP.stpgStatesToMdpStates.get(state)] = allowedMaximizerActions[state].get(0);
                }
                else {
                    tauMDPRecommendation[maxFixedLocalMDP.stpgStatesToMdpStates.get(state)] = allowedMinimizerActions[state];
                }
            }

            maxFixedResult = mdpModelChecker.computeReachProbsPolIter(maxFixedLocalMDP.mdp, maxFixedLocalMDP.sinks, maxFixedLocalMDP.targets, true, tauMDPRecommendation, false, maxFixedLocalMDP.initVector, minFixedLocalMDP.initVector);
            minFixedResult = mdpModelChecker.computeReachProbsPolIter(minFixedLocalMDP.mdp, minFixedLocalMDP.sinks, minFixedLocalMDP.targets, false, sigmaMDPRecommendation, false, minFixedLocalMDP.initVector, minFixedLocalMDP.initVector);
        }

        for (int state = relevantStates.nextSetBit(0); state >= 0; state = relevantStates.nextSetBit(state + 1)) {
            double maxFixedStateValue = maxFixedResult.soln[maxFixedLocalMDP.stpgStatesToMdpStates.get(state)];
            double minFixedStateValue = minFixedResult.soln[minFixedLocalMDP.stpgStatesToMdpStates.get(state)];
            if (!PrismUtils.doublesAreClose(maxFixedStateValue, minFixedStateValue, precision, true)) {
                throw new PrismException("If Minimizer has free choice then state "+state+" gets value "+maxFixedStateValue+
                        " but if Maximizer has free choice states gets value "+minFixedStateValue+
                        " ==> Topological value iteration failed");
            }
            fixedValues[state] = maxFixedStateValue;
        }


        return fixedValues;
    }
     */

    public static double[] getValueFromStrategies(PrismComponent prismComponent, STPG stpg, int[] sigma, int[] tau, BitSet yes, BitSet no, BitSet relevantStates, BitSet alreadyComputedStates, double[] fixedValues, double precision, double[] lowerbound, double[] upperBound) throws PrismException {
        /* This creates an whole MDP, but I suppose it's not necessary since we can create Local MDPs
        MDP maxFixedMDP = createMDPFromSTPG(stpg, sigma, true, yes, no, relevantStates, alreadyComputedStates, fixedValues);
        BitSet maxFixedNoBitSet = extendNoSet(no, maxFixedMDP, yes, relevantStates, alreadyComputedStates);
        MDP minFixedMDP = createMDPFromSTPG(stpg, tau, false, yes, no, relevantStates, alreadyComputedStates, fixedValues);
        BitSet minFixedNoBitSet = extendNoSet(no, minFixedMDP, yes, relevantStates, alreadyComputedStates);
         */

        if (fixedValues == null) {
            fixedValues = new double[stpg.getNumStates()];
        }
        LocalMDPResult maxFixedLocalMDP = createLocalMDPFromSTPG(stpg, sigma, true, yes, no, relevantStates, alreadyComputedStates, fixedValues, lowerbound);
        LocalMDPResult minFixedLocalMDP = createLocalMDPFromSTPG(stpg, tau, false, yes, no, relevantStates, alreadyComputedStates, fixedValues, upperBound);


        MDPModelChecker mdpModelChecker = new MDPModelChecker(prismComponent);
        mdpModelChecker.maxIters = 10000;

        boolean useLP = false;

        ModelCheckerResult maxFixedResult;
        ModelCheckerResult minFixedResult;

        if (useLP) {
            maxFixedResult = mdpModelChecker.computeReachProbsLinearProgrammingGurobi(maxFixedLocalMDP.mdp, maxFixedLocalMDP.sinks, maxFixedLocalMDP.targets, true, null, maxFixedLocalMDP.initVector);
            minFixedResult = mdpModelChecker.computeReachProbsLinearProgrammingGurobi(minFixedLocalMDP.mdp, minFixedLocalMDP.sinks, minFixedLocalMDP.targets, false, null, minFixedLocalMDP.initVector);

        }
        else {
            int[] sigmaMDPRecommendation = new int[minFixedLocalMDP.mdp.getNumStates()];
            int[] tauMDPRecommendation = new int[maxFixedLocalMDP.mdp.getNumStates()];

            for (int state = relevantStates.nextSetBit(0); state >= 0; state = relevantStates.nextSetBit(state + 1)) {
                // Max player
                if (stpg.getPlayer(state) == 1) {
                    sigmaMDPRecommendation[minFixedLocalMDP.stpgStatesToMdpStates.get(state)] = sigma[state];
                }
                else {
                    tauMDPRecommendation[maxFixedLocalMDP.stpgStatesToMdpStates.get(state)] = tau[state];
                }
            }

            maxFixedResult = mdpModelChecker.computeReachProbsPolIter(maxFixedLocalMDP.mdp, maxFixedLocalMDP.sinks, maxFixedLocalMDP.targets, true, tauMDPRecommendation, false, maxFixedLocalMDP.initVector, minFixedLocalMDP.initVector);
            minFixedResult = mdpModelChecker.computeReachProbsPolIter(minFixedLocalMDP.mdp, minFixedLocalMDP.sinks, minFixedLocalMDP.targets, false, sigmaMDPRecommendation, false, minFixedLocalMDP.initVector, minFixedLocalMDP.initVector);
        }

        for (int state = relevantStates.nextSetBit(0); state >= 0; state = relevantStates.nextSetBit(state + 1)) {
            double maxFixedStateValue = maxFixedResult.soln[maxFixedLocalMDP.stpgStatesToMdpStates.get(state)];
            double minFixedStateValue = minFixedResult.soln[minFixedLocalMDP.stpgStatesToMdpStates.get(state)];
            if (!PrismUtils.doublesAreClose(maxFixedStateValue, minFixedStateValue, precision, true)) {
                throw new PrismException("If Minimizer has free choice then state "+state+" gets value "+maxFixedStateValue+
                        " but if Maximizer has free choice states gets value "+minFixedStateValue+
                        " ==> Topological value iteration failed");
            }
            fixedValues[state] = maxFixedStateValue;
        }


        return fixedValues;
    }

    private static LocalDTMCTransformation createLocalDTMCFromStpg(STPG stpg, ArrayList<Integer>[] suggestedMaximizerActions, int[] suggestedMinimizerAction,
                                                BitSet yes, BitSet no,
                                                BitSet relevantStates,
                                                BitSet alreadyComputedStates, double[] fixedValues,
                                                double[] lowerBounds, double[] upperBounds) throws PrismException{
        DTMCSimple dtmc = new DTMCSimple();
        HashMap<Integer, Integer> stpgToDtmc = new HashMap<>();
        HashMap<Integer, Integer> dtmcToStpg = new HashMap<>();

        if (relevantStates == null) {
            relevantStates = new BitSet();
            relevantStates.set(0, stpg.getNumStates()-1);
        }

        for (int state = relevantStates.nextSetBit(0); state >= 0; state = relevantStates.nextSetBit(state + 1)) {
            int mdpState = dtmc.addState();
            stpgToDtmc.put(state, mdpState);
            dtmcToStpg.put(mdpState, state);
        }

        int sink = dtmc.addState();
        int target = dtmc.addState();

        double[] upperboundsInDTMC = new double[dtmc.getNumStates()];
        upperboundsInDTMC[sink] = 0;
        upperboundsInDTMC[target] = 1;

        //Set self-loops for sink and target
        dtmc.addToProbability(sink, sink, 1.0);
        dtmc.addToProbability(target, target, 1.0);

        //Copy the suggested actions (merge all suggested actions of maximizer into one
        for (int state = relevantStates.nextSetBit(0); state >= 0; state = relevantStates.nextSetBit(state + 1)) {
            boolean isMaximizerState = stpg.getPlayer(state) == 1;
            double probabilityScaling = 1.0;

            upperboundsInDTMC[stpgToDtmc.get(state)] = upperBounds[state];

            // One Distribution for all actions that will be contracted into one
            for (int choice = 0; choice < stpg.getNumChoices(state); choice++) {
                Distribution distribution;
                // Discard all actions that we are not interested in
                if ((isMaximizerState && !suggestedMaximizerActions[state].contains(choice)) || (!isMaximizerState && suggestedMinimizerAction[state] != choice)) {
                    continue;
                }

                if (isMaximizerState) {
                    // Distribute uniformly over all suggested Actions
                    probabilityScaling = 1.0 / suggestedMaximizerActions[state].size();
                }

                for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext(); ) {
                    Map.Entry<Integer, Double> tr = it.next();

                    int nextState = tr.getKey();
                    //Redirect action to new target
                    if (yes.get(nextState)) {
                        dtmc.addToProbability(stpgToDtmc.get(state), target, tr.getValue() * probabilityScaling);
                    }
                    //Redirect action to new sink
                    else if (no.get(nextState)) {
                        dtmc.addToProbability(stpgToDtmc.get(state), sink, tr.getValue() * probabilityScaling);
                    }
                    //Stay inside SCC
                    else if (relevantStates.get(nextState)) {
                        dtmc.addToProbability(stpgToDtmc.get(state), stpgToDtmc.get(nextState), tr.getValue() * probabilityScaling);
                    }
                    //Moving to a already fixed state -> Should already have a fixed value and thus can be replaced
                    else if (alreadyComputedStates.get(nextState)) {
                        double targetReachProb = fixedValues[nextState] * tr.getValue();
                        dtmc.addToProbability(stpgToDtmc.get(state), target, targetReachProb * probabilityScaling);
                        if (targetReachProb < 1) {
                            dtmc.addToProbability(stpgToDtmc.get(state), sink, (1.0f - fixedValues[nextState]) * tr.getValue() * probabilityScaling);
                        }
                    } else {
                        throw new PrismException("State " + state + " has in action " + choice + " a transition to state " + nextState + ", which is not already computed nor in the current SCC. This method is only suited for topological cases, in which this case should never arrise.");
                    }
                }
            }
        }

        boolean sanityCheck = true;
        if (sanityCheck) {
            for (int state = 0; state < dtmc.getNumStates(); state++) {
                double probability = 0;
                for (Iterator<Map.Entry<Integer, Double>> it = dtmc.getTransitionsIterator(state); it.hasNext(); ) {
                    Map.Entry<Integer, Double> tr = it.next();
                    probability += tr.getValue();
                }
                if (probability > 1.01 || probability < 0.99) {
                    throw new PrismException("For State " + dtmcToStpg.get(state) + " probability of all actions is " + probability + " and does not add up to 1 :(");
                }
            }
        }

        LocalDTMCTransformation transformationResult = new LocalDTMCTransformation();
        transformationResult.upperboundsInDTMC = upperboundsInDTMC;
        transformationResult.dtmc = dtmc;
        transformationResult.dtmcStatesToStpgStates = dtmcToStpg;
        transformationResult.stpgStatesToDtmcStates = stpgToDtmc;
        transformationResult.targets = new BitSet();
        transformationResult.targets.set(target);
        transformationResult.sinks = new BitSet();
        transformationResult.sinks.set(sink);
        return transformationResult;
    }

    /**
     * Creates a MDP that contains only the relevant states, a target and a sink.
     * MUST FULFILL FOLLOWING REQUIREMENT:
     * If a state s in relevantStates has a transition to another state s' then s' must be either in yes, no, relevantStates or alreadyComputedStates
     *
     * @param stpg
     * @param maxIsFixed
     * @param yes
     * @param no
     * @param relevantStates
     * @param alreadyComputedStates
     * @param fixedValues
     * @return
     * @throws PrismException
     */
    private static LocalMDPResult createLocalMDPFromSTPG(STPG stpg, ArrayList<Integer>[] allowedActions, boolean maxIsFixed, BitSet yes, BitSet no, BitSet relevantStates, BitSet alreadyComputedStates, double[] fixedValues, double[] bound) throws PrismException{
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

        double[] initVector = new double[mdp.getNumStates()];
        initVector[sink] = 0;
        initVector[target] = 1;


        //Set self-loops for sink and target
        Distribution d = new Distribution();
        d.add(sink, 1.0);
        mdp.addChoice(sink, d);

        d = new Distribution();
        d.add(target, 1.0);
        mdp.addChoice(target, d);

        //Copy actions and fix those that have to be fixed
        //Copy the suggested actions (merge all suggested actions of maximizer into one
        for (int state = relevantStates.nextSetBit(0); state >= 0; state = relevantStates.nextSetBit(state + 1)) {
            boolean isMaximizerState = stpg.getPlayer(state) == 1;
            double probabilityScaling = 1.0;

            if ((isMaximizerState && maxIsFixed) || (!isMaximizerState && !maxIsFixed)) {
                Distribution interActionDistribution = new Distribution();
                probabilityScaling = 1.0/((double) allowedActions[state].size());
                if (state == 6501) {
                    System.out.println("6501!! Probscaling: "+probabilityScaling);

                }
                double probability = 0.0;
                for (int choice : allowedActions[state]) {
                    if (state == 6501) {
                        System.out.println("Action: "+choice);

                    }
                    for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext(); ) {
                        Map.Entry<Integer, Double> tr = it.next();
                        int nextState = tr.getKey();
                        //Redirect action to new target
                        if (yes.get(nextState)) {
                            if (state == 6501) {
                                System.out.println("TARGET: "+tr.getValue()+", * prob -> "+(tr.getValue()*probabilityScaling));
                            }
                            probability+=tr.getValue() * probabilityScaling;
                            interActionDistribution.add(target, tr.getValue() * probabilityScaling);
                        }
                        //Redirect action to new sink
                        else if (no.get(nextState)) {
                            probability+=tr.getValue() * probabilityScaling;
                            interActionDistribution.add(sink, tr.getValue() * probabilityScaling);
                            if (state == 6501) {
                                System.out.println("SINK: "+tr.getValue()+", * prob -> "+(tr.getValue()*probabilityScaling));
                            }
                        }
                        //Stay inside SCC
                        else if (relevantStates.get(nextState)) {
                            probability+=tr.getValue() * probabilityScaling;
                            interActionDistribution.add(stpgToMdp.get(nextState), tr.getValue() * probabilityScaling);
                            if (state == 6501) {
                                System.out.println("IN SCC: "+tr.getValue()+", * prob -> "+(tr.getValue()*probabilityScaling));
                            }
                        }
                        //Moving to a already fixed state -> Should already have a fixed value and thus can be replaced
                        else if (alreadyComputedStates.get(nextState)) {
                            double targetReachProb = fixedValues[nextState] * tr.getValue();
                            interActionDistribution.add(target, targetReachProb * probabilityScaling);
                            if (targetReachProb < 1) {
                                interActionDistribution.add(sink, (1.0f - fixedValues[nextState]) * tr.getValue() * probabilityScaling);
                            }
                            probability+= targetReachProb * probabilityScaling;
                            probability+= (1.0f - fixedValues[nextState]) * tr.getValue() * probabilityScaling;
                            if (state == 6501) {
                                System.out.println("OUT OF SCC: Probability: "+tr.getValue()+", Value of Outside state: "+(fixedValues[nextState]));
                                System.out.println("ROUTED FROM 6501 TO TARGET: "+targetReachProb * probabilityScaling+", TO SINK: "+((1.0f - fixedValues[nextState]) * tr.getValue() * probabilityScaling));
                            }
                        } else {
                            throw new PrismException("State " + state + " has in action " + choice + " a transition to state " + nextState + ", which is not already computed nor in the current SCC. This method is only suited for topological cases, in which this case should never arrise.");
                        }
                    }
                }
                if (probability > 1.01 || probability < 0.99) {
                    throw new PrismException("For state "+state+" and and allowed actions "+allowedActions[state]+" the summed probability is "+probability);
                }
                mdp.addChoice(stpgToMdp.get(state), interActionDistribution);
            }
            else {
                for (int choice = 0; choice < stpg.getNumChoices(state); choice++) {
                    double probability = 0.0;
                    Distribution distribution = new Distribution();
                    for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext(); ) {
                        Map.Entry<Integer, Double> tr = it.next();
                        int nextState = tr.getKey();
                        //Redirect action to new target
                        if (yes.get(nextState)) {
                            probability+=tr.getValue();
                            distribution.add(target, tr.getValue());
                        }
                        //Redirect action to new sink
                        else if (no.get(nextState)) {
                            distribution.add(sink, tr.getValue());
                            probability+=tr.getValue();
                        }
                        //Stay inside SCC
                        else if (relevantStates.get(nextState)) {
                            distribution.add(stpgToMdp.get(nextState), tr.getValue());
                            probability+=tr.getValue();
                        }
                        //Moving to a already fixed state -> Should already have a fixed value and thus can be replaced
                        else if (alreadyComputedStates.get(nextState)) {
                            double targetReachProb = fixedValues[nextState] * tr.getValue();
                            distribution.add(target, targetReachProb);
                            if (targetReachProb < 1) {
                                distribution.add(sink, (1.0f - fixedValues[nextState]) * tr.getValue());
                            }
                            probability+=fixedValues[nextState] * tr.getValue();
                            probability+=(1.0f - fixedValues[nextState]) * tr.getValue();

                        } else {
                            throw new PrismException("State " + state + " has in action " + choice + " a transition to state " + nextState + ", which is not already computed nor in the current SCC. This method is only suited for topological cases, in which this case should never arrise.");
                        }
                    }
                    if (probability > 1.01 || probability < 0.99) {
                        throw new PrismException("For state "+state+" and and allowed actions "+allowedActions[state]+" the summed probability is "+probability);
                    }
                    mdp.addChoice(stpgToMdp.get(state), distribution);
                }
            }
        }

        // Take any Initial State, since it does not matter
        mdp.addInitialState(0);

        LocalMDPResult mdpResult = new LocalMDPResult();
        mdpResult.mdp = mdp;
        mdpResult.stpgStatesToMdpStates = stpgToMdp;
        mdpResult.mdpStatesToStpgStates = mdpToStpg;
        mdpResult.targets = new BitSet();
        mdpResult.targets.set(target);
        mdpResult.sinks = new BitSet();
        mdpResult.sinks.set(sink);
        mdpResult.initVector = initVector;
        return mdpResult;
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
    private static LocalMDPResult createLocalMDPFromSTPG(STPG stpg, int[] strategy, boolean maxIsFixed, BitSet yes, BitSet no, BitSet relevantStates, BitSet alreadyComputedStates, double[] fixedValues, double[] bound) throws PrismException{
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

        double[] initVector = new double[mdp.getNumStates()];
        initVector[sink] = 0;
        initVector[target] = 1;


        //Set self-loops for sink and target
        Distribution d = new Distribution();
        d.add(sink, 1.0);
        mdp.addChoice(sink, d);

        d = new Distribution();
        d.add(target, 1.0);
        mdp.addChoice(target, d);

        //Copy actions and fix those that have to be fixed
        for (int state = relevantStates.nextSetBit(0); state >= 0; state = relevantStates.nextSetBit(state + 1)) {
            initVector[stpgToMdp.get(state)] = bound[state];
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
                        if (fixedValues[nextState] * tr.getValue() < 1) {
                            distribution.add(sink, fixedValues[nextState] * tr.getValue());
                        }
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
        mdpResult.mdpStatesToStpgStates = mdpToStpg;
        mdpResult.targets = new BitSet();
        mdpResult.targets.set(target);
        mdpResult.sinks = new BitSet();
        mdpResult.sinks.set(sink);
        mdpResult.initVector = initVector;
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

    /**
     * Implements Widest Path Deflating as presented in Widest Paths and Global Propagation in Bounded Value Iteration for Stochastic Games
     * (https://link.springer.com/chapter/10.1007%2F978-3-030-53291-8_19)
     * Put this method BEFORE swapping lowerBoundsNew and lowerBounds and BEFORE swapping upperBounds and upperBoundsNew
     * 
     * NOTE: CURRENTLY DOES NOT (!!!) SUPPORT GAUSS-SEIDEL OPTIONS OR TOPOLOGICAL!
     * @param stpg
     * @param iters
     * @param numberOfIterationsBeforeApplication : How many iterations to wait before applying deflating (In paper by default 5)
     * @param yes
     * @param lowerBounds
     * @param lowerBoundsNew
     * @param upperBoundsNew
     * @return {lowerBoundsNew, lowerBounds, upperBoundsNew}
     */
    public static double[][] widestPathDeflating(STPGExplicit stpg, int iters, int numberOfIterationsBeforeApplication, BitSet yes, double[] lowerBounds, double[] lowerBoundsNew, double[] upperBoundsNew) {
        int n = stpg.getNumStates();
        // New Algo from WP
        if (iters % numberOfIterationsBeforeApplication == 0) {

            // create new graph
            LinkedList<Pair<Double, Integer>> G[] = new LinkedList[n];
            int visit[] = new int[n];
            Heap heap = new Heap(n);

            for (int s = 0; s < n; s++) {
                G[s] = new LinkedList<>();
                heap.pointer[s] = -1;
            }

            for (int s = 0; s < n; s++) {
                // Maximizer
                if (stpg.getPlayer(s) == 1) {
                    for (int a = 0; a < stpg.getNumChoices(s); a++) {
                        Distribution distr = ((SMG) stpg).trans.get(s).get(a);
                        double up_val = 0.0;
                        for (Map.Entry<Integer, Double> e : distr) {
                            int t = e.getKey();
                            double prob = e.getValue();
                            up_val += prob * upperBoundsNew[t];
                        }
                        for (Map.Entry<Integer, Double> e : distr) {
                            int t = e.getKey();
                            G[t].add(new Pair(up_val, s));
                        }
                    }
                }
                // Minimizer
                else {
                    for (int a = 0; a < stpg.getNumChoices(s); a++) {
                        Distribution distr = ((SMG) stpg).trans.get(s).get(a);
                        double up_val = 0.0, low_val = 0.0;
                        for (Map.Entry<Integer, Double> e : distr) {
                            int t = e.getKey();
                            double prob = e.getValue();
                            up_val += prob * upperBoundsNew[t];
                            low_val += prob * lowerBounds[t];
                        }
                        if (low_val <= lowerBoundsNew[s]) {
                            for (Map.Entry<Integer, Double> e : distr) {
                                int t = e.getKey();
                                G[t].add(new Pair(up_val, s));
                            }
                        }
                    }
                }
            }

            // mainLog.println("finish creating graph in itr "+iters);

            for (int s = yes.nextSetBit(0); s >= 0; s = yes.nextSetBit(s + 1))
                heap.append(1.0, s);

            while (heap.heap_size > 0) {
                Pair<Double, Integer> top = heap.pop();
                double v = top.x;
                int t = top.y;
                visit[t] = 1;
                upperBoundsNew[t] = v;

                for (Pair<Double, Integer> p : G[t]) {
                    double w = p.x;
                    int s = p.y;
                    if (visit[s] == 0) {
                        heap.update(min(v, w, upperBoundsNew[s]), s);
                    }
                }
            }

            for (int s = 0; s < n; s++) {
                if (visit[s] == 0)
                    upperBoundsNew[s] = 0;
            }
        }

        return new double[][]{lowerBounds, lowerBoundsNew, upperBoundsNew};
    }
    private static class LocalMDPResult {
        MDP mdp;
        BitSet targets;
        BitSet sinks;
        HashMap<Integer, Integer> stpgStatesToMdpStates;
        HashMap<Integer, Integer> mdpStatesToStpgStates;
        double[] initVector;
    }

    private static class LocalDTMCTransformation {
        DTMC dtmc;
        BitSet targets;
        BitSet sinks;
        HashMap<Integer, Integer> stpgStatesToDtmcStates;
        HashMap<Integer, Integer> dtmcStatesToStpgStates;
        double[] upperboundsInDTMC;
    }

    // Widest Path Stuff

    public static class Pair<X, Y> {
        public X x;
        public Y y;
        public Pair(X x, Y y) {
            this.x = x;
            this.y = y;
        }
    }

    public static class Heap {
        public Pair<Double, Integer>[] heap;
        public int[] pointer;
        public int heap_size;

        public Heap(int n) {
            heap = new Pair[n];
            pointer = new int[n];
            heap_size = 0;
        }

        public void append(double p, int s) {
            heap[heap_size] = new Pair(p,s);
            pointer[s] = heap_size;
            heap_size++;
        }

        public Pair<Double, Integer> pop() {
            Pair<Double, Integer> top = heap[0];
            heap[0] = heap[heap_size-1];
            heap_size--;
            down(0);
            return top;
        }

        public void update(double p, int s) {
            if(pointer[s] == -1) {
                append(p,s);
                up(heap_size-1);
            }
            else if(heap[pointer[s]].x < p){
                heap[pointer[s]].x = p;
                up(pointer[s]);
            }
        }

        private void swap(int i,int j) {
            pointer[heap[i].y] = j;
            pointer[heap[j].y] = i;

            Pair tmp;
            tmp = heap[i];
            heap[i] = heap[j];
            heap[j] = tmp;
        }

        private void up(int i) {
            while(i > 0) {
                int p = (i-1)/2;
                double v_i = heap[i].x;
                double v_p = heap[p].x;
                if(v_i > v_p) {
                    swap(i,p);
                    i = p;
                }
                else
                    break;
            }
        }

        private void down(int i) {
            while(true) {
                int left = 2*i+1;
                int right = 2*i+2;

                if(left >= heap_size)
                    break;

                if(right >= heap_size) {
                    double v_left = heap[left].x;
                    double v_i = heap[i].x;
                    if(v_i < v_left) {
                        swap(i,left);
                        i = left;
                    }
                    else
                        break;
                }
                else {
                    double v_left = heap[left].x;
                    double v_right = heap[right].x;
                    double v_i = heap[i].x;
                    if (v_left >= v_right && v_i < v_left) {
                        swap(i, left);
                        i = left;
                    }
                    else if (v_right > v_left && v_i < v_right) {
                        swap(i, right);
                        i = right;
                    }
                    else
                        break;
                }
            }
        }
    }

    private static double min(double a,double b,double c) {
        if(a <= b && a <= c) return a;
        if(b <= a && b <= c) return b;
        return c;
    }
}
