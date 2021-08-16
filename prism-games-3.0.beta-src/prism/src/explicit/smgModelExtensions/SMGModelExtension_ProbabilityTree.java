package explicit.smgModelExtensions;

import explicit.Distribution;
import explicit.SMGModelExtender;
import explicit.STPGExplicit;
import parser.State;

import java.util.BitSet;
import java.util.List;
import java.util.Random;

public class SMGModelExtension_ProbabilityTree extends SMGModelExtension{
    int oldInitialState;
    int numComponents;
    int treeDepth;
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
     * @param treeDepth how many levels should the tree have (root level is 0)
     * @param treeBranchingFactor how many actions does each node of the tree have
     * @param probabilityToReachComponentsInitialState what is the probability to reach the initial state while using an action?
     * @param probabilityToReachSink what is the probability to reach the sink while using an action. The sum of probabilityToReachSink and probabilityToReachComponentsInitialState must be in [0, 1]
     * @param sinkState sink
     */
    public SMGModelExtension_ProbabilityTree(STPGExplicit stpg, BitSet remain, List<State> statesList, boolean useThisExtension, int numComponents, int treeDepth, int treeBranchingFactor,
                                             double probabilityToReachComponentsInitialState, double probabilityToReachSink, int sinkState) {
        super(stpg, remain, statesList, useThisExtension);
        this.sinkState = sinkState;
        this.numComponents = numComponents;
        this.treeDepth = treeDepth;
        this.treeBranchingFactor = treeBranchingFactor;
        this.probabilityToReachComponentsInitialState = probabilityToReachComponentsInitialState;
        this.probabilityToReachSink = probabilityToReachSink;

    }

    @Override
    public ModelExtensionResult extendSMG(int oldInitialState) {
        if (numComponents <= 0) {
            throw new IllegalArgumentException("[Model Extension]: The number of added Components per initial state must be at least 1");
        }
        else if (treeDepth < 0) {
            throw new IllegalArgumentException("[Model Extension]: The depth of added ProbabilityTrees per initial state must be at least 0");
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

        int numNewStatesPerComponent = 1;
        for (int depth = 1; depth <= treeDepth; depth++) {
            numNewStatesPerComponent += Math.pow(treeBranchingFactor, depth);
        }
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
                int childrenOffset = (currentState - currentInitialState) * treeBranchingFactor + currentInitialState; //Indexing in Trees
                for (int nextState = childrenOffset + 1; nextState <= childrenOffset + treeBranchingFactor; nextState++) {
                    if (nextState >= currentInitialState + numNewStatesPerComponent) {
                        nextState = nextInitialState;
                    }
                    Distribution transitionDistribution = new Distribution();

                    if (probabilityToReachSink > 0) {
                        transitionDistribution.add(sinkState, probabilityToReachSink);
                    }
                    if (probabilityToReachSink < 1) {
                        transitionDistribution.add(currentInitialState, probabilityToReachComponentsInitialState);
                        transitionDistribution.add(nextState, (1-probabilityToReachComponentsInitialState-probabilityToReachSink));
                    }
                    stpg.addChoice(currentState, transitionDistribution);

                    //Without this statement the for-loop may never end
                    if (nextState == nextInitialState) {
                        break;
                    }
                }
            }
            currentInitialState = nextInitialState;
        }

        ModelExtensionResult extensionResult = new ModelExtensionResult(addedSTPGInitialState, stpg);
        return extensionResult;
    }
}
