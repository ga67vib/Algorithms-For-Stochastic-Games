package explicit.smgModelExtensions;

import explicit.Distribution;
import explicit.STPGExplicit;
import parser.State;

import java.util.BitSet;
import java.util.List;
import java.util.Random;

public class SMGModelExtension_ActionTree extends SMGModelExtension{
    int oldInitialState;
    int numComponents;
    int numStates;
    int treeBranchingFactor;
    double probabilityToReachComponentsInitialState;
    double probabilityToReachSink;
    int sinkState;

    /**
     *
     * @param stpg STPG to extend
     * @param remain
     * @param statesList Is necessary, because stateList from STPG is protected in the package
     * @param useThisExtension Should this extension be used
     * @param numComponents how many components should there be?
     * @param numStates how many states should the tree have
     * @param treeBranchingFactor how many actions does each node of the tree have
     * @param probabilityToReachComponentsInitialState what is the probability to reach the initial state while using an action?
     * @param probabilityToReachSink what is the probability to reach the sink while using an action. The sum of probabilityToReachSink and probabilityToReachComponentsInitialState must be in [0, 1]
     * @param sinkState sink
     */
    public SMGModelExtension_ActionTree(STPGExplicit stpg, BitSet remain, List<State> statesList, boolean useThisExtension, int numComponents, int numStates, int treeBranchingFactor,
                                        double probabilityToReachComponentsInitialState, double probabilityToReachSink, int sinkState) {
        super(stpg, remain, statesList, useThisExtension);
        this.sinkState = sinkState;
        this.numComponents = numComponents;
        this.numStates = numStates;
        this.treeBranchingFactor = treeBranchingFactor;
        this.probabilityToReachComponentsInitialState = probabilityToReachComponentsInitialState;
        this.probabilityToReachSink = probabilityToReachSink;

    }

    @Override
    public ModelExtensionResult extendSMG(int oldInitialState) {
        if (numComponents <= 0) {
            throw new IllegalArgumentException("[Model Extension]: The number of added Components per initial state must be at least 1");
        }
        else if (numStates <= 0) {
            throw new IllegalArgumentException("[Model Extension]: There must be at least one state");
        }
        else if (treeBranchingFactor <= 0) {
            throw new IllegalArgumentException("[Model Extension]: The treeBranchingFactor must be at least 1");
        }
        else if (probabilityToReachSink < 0 || probabilityToReachSink > 1) {
            throw new IllegalArgumentException("[Model Extension]: Probability of reaching a sink was set to "+probabilityToReachSink+" which is not in [0,1]");
        }
        else if (probabilityToReachComponentsInitialState < 0 || probabilityToReachComponentsInitialState > 1) {
            throw new IllegalArgumentException("[Model Extension]: Probability to branch in tree was set to "+probabilityToReachComponentsInitialState+" which is not in [0,1]");
        }
        else if (probabilityToReachSink + probabilityToReachComponentsInitialState > 1) {
            throw new IllegalArgumentException("[Model Extension]: The sum of the probabilities of reaching a sink and reaching the initial states is "+(probabilityToReachSink+probabilityToReachComponentsInitialState)+" which is not in [0,1]");
        }

        Random random = new Random();
        random.setSeed(100); //Set deterministic Seed because we do want deterministic models but don't care about who the states belong

        long numNewStatesPerComponentLong = numStates;

        if (numNewStatesPerComponentLong * numComponents > Integer.MAX_VALUE || numNewStatesPerComponentLong * numComponents < 0) {
            throw new IllegalArgumentException("[Model Extension]: The number of new states is too big. You exceed Integer.MAX_VALUE. Try reducing tree depth or the number of copmonents");
        }

        int numNewStatesPerComponent = (int) numNewStatesPerComponentLong;
        int numNewStates = numNewStatesPerComponent*numComponents;


        int addedSTPGInitialState = stpg.getNumStates();
        int currentInitialState = addedSTPGInitialState;
        int nextInitialState;

        int currentState;

        stpg.addStates(numNewStates);
        if (remain != null) remain.set(stpg.getNumStates()-numNewStates, stpg.getNumStates());
        addToStateList(stpg, stpg.getNumStates()-numNewStates, stpg.getNumStates());

        for (int component = 0; component < numComponents; component++) {
            // Make sure to know which will be the next initialState
            nextInitialState = currentInitialState + numNewStatesPerComponent;

            // Index out of bounds -> reached the end of prepended MEC. Next state the MEC leads into must be old initial state
            if (nextInitialState == stpg.getNumStates()) {
                nextInitialState = oldInitialState;
            }

            for (currentState = currentInitialState; currentState < currentInitialState + numNewStatesPerComponent; currentState++) {
                stpg.setPlayer(currentState, random.nextInt(2)+1); //Set player to either 1 or 2
                double actionCounter = 0;
                int childrenOffset = (currentState - currentInitialState) * treeBranchingFactor + currentInitialState; //Indexing in Trees
                for (int nextState = childrenOffset + 1; nextState <= childrenOffset + treeBranchingFactor; nextState++) {
                    int trueNextState = nextState;
                    if (nextState >= currentInitialState + numNewStatesPerComponent) {
                        trueNextState = nextInitialState;
                    }
                    Distribution transitionDistribution = new Distribution();

                    if (probabilityToReachSink > 0) {
                        transitionDistribution.add(sinkState, probabilityToReachSink);
                    }
                    if (probabilityToReachSink < 1) {
                        double probToShare = 1.0 - probabilityToReachSink - probabilityToReachComponentsInitialState;
                        double probToSelf = probToShare - (1.0 * actionCounter / (this.treeBranchingFactor+1));
                        double probToNext = probToShare - probToSelf;
                        if (probabilityToReachComponentsInitialState > 0) transitionDistribution.add(currentInitialState, probabilityToReachComponentsInitialState);
                        transitionDistribution.add(trueNextState, (probToNext));
                        transitionDistribution.add(currentState, probToSelf);
                    }
                    stpg.addChoice(currentState, transitionDistribution);
                    actionCounter++;
                }
            }
            currentInitialState = nextInitialState;
        }

        ModelExtensionResult extensionResult = new ModelExtensionResult(addedSTPGInitialState, stpg);
        return extensionResult;
    }
}
