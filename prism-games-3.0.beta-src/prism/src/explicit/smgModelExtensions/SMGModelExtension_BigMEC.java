package explicit.smgModelExtensions;

import explicit.Distribution;
import explicit.STPGExplicit;
import parser.State;

import java.util.BitSet;
import java.util.List;

public class SMGModelExtension_BigMEC extends SMGModelExtension {
    int numChainedMECs;
    int chainLengthInMEC;
    double sinkProbabilityPerMEC;
    int sinkState;
    double probabilityInMEC;

    /**
     *
     * @param stpg STPG to extend
     * @param remain
     * @param statesList Is necessary, because stateList from STPG is protected in the package
     * @param useThisExtension Should this extension be used
     * @param numChainedMECs number of BigMEC-components chained one after the other
     * @param chainLengthInMEC length of chain per BigMEC
     * @param sinkProbabilityPerMEC probability leading from the final action of a chain in a BigMEC to the sink
     * @param sinkState sink state
     * @param probabilityInMEC probability of an action inside the chain to jump to the parallel state of the same chain
     */
    public SMGModelExtension_BigMEC(STPGExplicit stpg, BitSet remain, List<State> statesList, boolean useThisExtension, int numChainedMECs, int chainLengthInMEC,
                                    double sinkProbabilityPerMEC, int sinkState, double probabilityInMEC) {
        super(stpg, remain, statesList, useThisExtension);
        this.numChainedMECs = numChainedMECs;
        this.chainLengthInMEC = chainLengthInMEC;
        this.sinkProbabilityPerMEC = sinkProbabilityPerMEC;
        this.sinkState = sinkState;
        this.probabilityInMEC = probabilityInMEC;

    }

    @Override
    public ModelExtensionResult extendSMG(int oldInitialState) {
        if (numChainedMECs <= 0) {
            throw new IllegalArgumentException("[Model Extension]: The number of added MECs per initial state must be at least 1");
        }
        else if (chainLengthInMEC <= 0) {
            throw new IllegalArgumentException("[Model Extension]: The length of added MEC-Chains per initial state must be at least 1");
        }
        else if (sinkProbabilityPerMEC < 0 || sinkProbabilityPerMEC > 1) {
            throw new IllegalArgumentException("[Model Extension]: Probability of reaching a sink was set to "+sinkProbabilityPerMEC+" which is not in [0,1]");
        }
        else if (probabilityInMEC < 0 || probabilityInMEC > 1) {
            throw new IllegalArgumentException("[Model Extension]: Probability to branch in chain was set to "+probabilityInMEC+" which is not in [0,1]");
        }

        int numChains = 2;
        int numNewStatesPerMEC = 2 * numChainedMECs+1;
        int numNewStates = numChainedMECs*numNewStatesPerMEC;
        int addedSTPGInitialState = stpg.getNumStates();
        int currentInitialState = addedSTPGInitialState;
        int nextInitialState;

        int currentState;

        stpg.addStates(numNewStates);
        if (remain != null) remain.set(stpg.getNumStates()-numNewStates, stpg.getNumStates());
        addToStateList(stpg, stpg.getNumStates()-numNewStates, stpg.getNumStates());

        for (int chainedMEC = 0; chainedMEC < numChainedMECs; chainedMEC++) {
            // Make sure to know which will be the next initialState
            nextInitialState = currentInitialState + numNewStatesPerMEC;

            // Index out of bounds -> reached the end of prepended MEC. Next state the MEC leads into must be old initial state
            if (nextInitialState == stpg.getNumStates()) {
                nextInitialState = oldInitialState;
            }

            stpg.setPlayer(currentInitialState, 2); //Min-Player gets to decide the beginning of MEC

            // Add transitions from chains initial state to both chains
            Distribution choiceDistribution = new Distribution();
            choiceDistribution.add(currentInitialState+1, 1.0);
            stpg.addChoice(currentInitialState, choiceDistribution);

            choiceDistribution = new Distribution();
            choiceDistribution.add(currentInitialState+2, 1.0);
            stpg.addChoice(currentInitialState, choiceDistribution);

            // Loop for the chain-complexes
            for (int chainMember = 0; chainMember < chainLengthInMEC; chainMember++) {
                // Loop for upper and lower chain
                for (int chain = 1; chain < numChains + 1; chain++) {
                    currentState = currentInitialState + chainMember * numChains + chain;


                    stpg.setPlayer(currentState, 1);
                    int previousState = currentState - numChains;
                    int nextState = currentState + numChains;

                    if (previousState <= currentInitialState) {
                        previousState = currentInitialState;
                    }
                    if (nextState >= currentInitialState + numNewStatesPerMEC) {
                        nextState = nextInitialState;
                    }

                    // add Choice to go back
                    choiceDistribution = new Distribution();
                    choiceDistribution.add(previousState, 1);
                    stpg.addChoice(currentState, choiceDistribution);

                    // add Choice to go to next
                    if (nextState == nextInitialState) {
                        // exit choices of the MEC
                        double probabilityNotToSink = 1.0 - sinkProbabilityPerMEC;
                        double probabilityToReachNextMEC = probabilityNotToSink * 0.5;
                        double probabilityToGoBackToCurrentInitialState = probabilityNotToSink * 0.5;
                        if (chain == 1) {
                            probabilityToReachNextMEC = probabilityNotToSink * 0.4;
                            probabilityToGoBackToCurrentInitialState = probabilityNotToSink * 0.6;
                        }
                        choiceDistribution = new Distribution();
                        if (sinkProbabilityPerMEC > 0) choiceDistribution.add(sinkState, sinkProbabilityPerMEC);
                        if (sinkProbabilityPerMEC < 1) {
                            choiceDistribution.add(nextInitialState, probabilityToReachNextMEC);
                            choiceDistribution.add(currentInitialState, probabilityToGoBackToCurrentInitialState);
                        }
                    }
                    else {
                        // If not exit, one can move forward or jump to other chain
                        int parallelState = (chain == 1) ? (currentState + 1) : (currentState - 1);
                        choiceDistribution = new Distribution();
                        choiceDistribution.add(nextState, 1.0 - probabilityInMEC);
                        if (probabilityInMEC > 0) choiceDistribution.add(parallelState, probabilityInMEC);
                    }
                    stpg.addChoice(currentState, choiceDistribution);
                }
            }
            currentInitialState = nextInitialState;
        }

        ModelExtensionResult extensionResult = new ModelExtensionResult(addedSTPGInitialState, stpg);
        return extensionResult;
    }
}
