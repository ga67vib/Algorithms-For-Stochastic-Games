package explicit;

import prism.Pair;
import prism.PrismException;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.*;

public class CondonNFTransformation {

    /**
     *
     * @param stpg
     * @param yes
     * @param no
     * @return
     * @throws PrismException
     */
    public static STPG normalize(STPG stpg, BitSet no, BitSet yes, List<BitSet> mecs) throws PrismException{
        STPG newSTPG=create2SuccForm(stpg, no, yes);
        newSTPG=createHalfProbsForm(newSTPG);
        double n = newSTPG.getNumStates();
        int m = 2 * ((int)n);
        //m = (int)(Math.ceil(Math.log(1.0/(1- Math.pow(1.0-(1.0/Math.pow(4.0,n)),(1.0/n))))) + 1.0);
        if (!mecs.isEmpty()) {
            newSTPG = toStopping(newSTPG, m, no, yes);
        }
        addUselessAction(((STPGExplicit)newSTPG), no, yes);
        return newSTPG;
    }

    public static STPG normalizeLoose(STPG stpg, BitSet no, BitSet yes) throws PrismException{
        return create2SuccForm(stpg, no, yes);
    }

    public static STPG normalizeLooseAndStopping(STPG stpg, BitSet no, BitSet yes, List<BitSet> mecs, boolean mStates) throws PrismException{
        STPG newSTPG = create2SuccForm(stpg, no, yes);
        double p = 0.5;
        newSTPG = toMECStopping(newSTPG, mecs, p, no, mStates);
        return newSTPG;
    }

    /**
     * Checks in a more loose form:
     * Are there 1 or 2 actions for each node?
     *
     * @param stpg
     * @param yes
     * @param no
     * @return
     * @throws PrismException
     */
    public static boolean checkCondonNormalformLoose(STPG stpg, BitSet no, BitSet yes) throws PrismException {
        int numStates = stpg.getNumStates();
        for (int i = 0; i < numStates; i++) {
            int choices = stpg.getNumChoices(i);

            if (choices != 2 && choices != 1) {
                throw new PrismException("State " + i + " has neither one nor two choices");
            }
        }
        return true;
    }

    /**
     * Checks, whether given stpg is in Condon Normalform. The method is checking for:
     * <p>
     * Are there 1 or 2 actions for each node?
     * Every state with 1 Action must be in yes/no (this is more general than real CondonNF because there every state with 1 Action must be a sink)
     * There can only be one yes-sink and one no-sink, which have to be the n-th and the n-1th nodes of the graph
     * Each action may have 1 or 2 transitions
     *
     * @param stpg
     * @param yes
     * @param no
     * @return
     * @throws PrismException
     */
    public static boolean checkCondonNormalformExplicit(STPG stpg, BitSet no, BitSet yes) throws PrismException {
        int numStates = stpg.getNumStates();
        for (int i = 0; i < numStates; i++) {
            int choices = stpg.getNumChoices(i);


            if (choices == 1) { //must be sink state
                //Vertex n must be 1-sink, Vertex n-1 must be 0-sink
                if (!(yes.get(i) && i == numStates - 1) && !(no.get(i) && i == numStates - 2)) {
                    throw new PrismException("State " + i + " has only one action which doesn't lead into itself!");
                }
            } else if (choices == 2) {
                for (int j = 0; j < choices; j++) {

                    int count = 0; //checking how many transitions an action has
                    for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(i, j); it.hasNext(); it.next()) {
                        count++;
                    }

                    if (count == 1) {
                        //should be nothing to do
                    } else if (count == 2) {
                        for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(i, j); it.hasNext(); ) {
                            Map.Entry<Integer, Double> tr = it.next();

                            if (tr.getValue() != 0.5) {
                                throw new PrismException("State " + i + " with choice " + j + " has two transitions but at least one with probability not equal to 0.5");
                            }

                        }
                    } else {
                        throw new PrismException("State " + i + " with choice " + j + " has neither one nor two transitions!");
                    }
                }
            } else {
                throw new PrismException("State " + i + " has neither one nor two choices");
            }
        }
        return true;
    }

    /**
     * Blows up the graph by replacing all nodes with more than 2 actions by binary-node-trees
     * Doesn't delete any nodes for now.
     * <p>
     * !!IMPORTANT!! FUNCTION ADDS THE TREE-NODES INDEX-WISE AFTER THE OLD GRAPH
     * THEREFOR IT THE NEW STPG DEPENDS ON THE yes AND no BITSET!
     * IF NODES ARE DELETED THERE MAY BE PROBLEMS WITH THIS BECAUSE THERE ARE NO NEW BITSETS COMPUTED!!
     * <p>
     * ALSO THIS METHOD EXPECTS THE yes-BITSET TO BE NON-EMPTY (Though this should be always given)
     *
     * @param stpg
     * @param no
     * @param yes
     * @return
     * @throws PrismException
     */
    private static STPG create2SuccForm(STPG stpg, BitSet no, BitSet yes) throws PrismException {
        //Don't delete vertices in this Function because the function depends on having the same
        //amout of vertices as the stpg from the parameter
        STPGExplicit newSTPG = new STPGExplicit();
        int stateNum = stpg.getNumStates();
        int choiceNum;
        Distribution d;

        int maxActions=0;
        for (int i=0; i<stpg.getNumStates(); i++) {
            if (stpg.getNumChoices(i) > maxActions) maxActions = stpg.getNumChoices(i);
        }



        for (int state = 0; state < stateNum; state++) {
            newSTPG.addState(stpg.getPlayer(state));
        }
        for (Integer initstate : stpg.getInitialStates()) {
            newSTPG.initialStates.add(initstate);
        }

        for (int state = 0; state < stateNum; state++) {
            choiceNum = stpg.getNumChoices(state);

            if (choiceNum == 0) {
                throw new PrismException("State " + state + " has no Transition");
            }
            else if (choiceNum == 1) {
                newSTPG.addChoice(state, getDistribution(stpg, state, 0));

            } else if (choiceNum == 2) {
                newSTPG.addChoice(state, getDistribution(stpg, state, 0));
                newSTPG.addChoice(state, getDistribution(stpg, state, 1));
            }

            //transition has >2 states
            else if (choiceNum > 2) {

                //Iterator<Entry<Integer, Double>> it = stpg.getTransitionsIterator(0, 0);
                //Entry<Integer, Double> tr = it.next();

                int transState;
                LinkedList<Pair<Integer, Integer>> queue = new LinkedList<>(); //(state, action) Tuples
                for (int action = 0; action < choiceNum; action++) {
                    queue.add(new Pair<>(state, action));
                }

                Pair<Integer, Integer> pair;

                while (queue.size() > 1) {
                    if (queue.size() == 2) { //So the root of the tree is still the actual vertex from the original graph
                        transState = state;
                    } else {
                        transState = newSTPG.addState(stpg.getPlayer(state));
                    }

                    for (int index = 0; index < 2; index++) {
                        pair = queue.poll();

                        //If we talk about our root, we need to assign his action to a tree-node
                        if (pair.getKey() == state) {
                            d = getDistribution(stpg, pair.getKey(), pair.getValue());
                        }
                        //Else we are talking about about assigning a tree-node to another tree-node(or to the root)
                        else {
                            d = new Distribution();
                            d.add(pair.getKey(), 1.0);
                        }
                        //New node should point to the nodes in pair / should take over old action
                        newSTPG.addChoice(transState, d);
                    }
                    queue.add(new Pair<>(transState, 0)); //action doesn't matter here

                }
            }
        }
        return newSTPG;
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
     * @param stpg
     * @return
     */
    private static STPG splitProbabilities(STPG stpg) {
        STPGExplicit newSTPG = new STPGExplicit();
        int stateNum = stpg.getNumStates();
        int choiceNum;
        int transNum;
        Distribution d;

        for (int state = 0; state < stateNum; state++) {
            newSTPG.addState(stpg.getPlayer(state));
        }
        for (Integer initstate : stpg.getInitialStates()) {
            newSTPG.initialStates.add(initstate);
        }

        for (int state = 0; state < stateNum; state++) {
            choiceNum = stpg.getNumChoices(state);
            for (int choice = 0; choice < choiceNum; choice++) {
                transNum=stpg.getNumTransitions(state, choice);
                d=getDistribution(stpg, state, choice);
                if (transNum<=2) { //If there are no changes to be made, just add the choice to the state
                    d=getDistribution(stpg, state, choice);
                    newSTPG.addChoice(state, d);
                }
                else {
                    {
                        int tmpState;
                        double[] probabilities=new double[transNum];
                        int[] destinations=new int[transNum];
                        int i=0;
                        for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext();) {
                            Map.Entry<Integer, Double> tr = it.next();
                            destinations[i]=tr.getKey();
                            probabilities[i]=tr.getValue();
                            i++;
                        }
                        LinkedList<ArrayList<Integer>> vertexHoldsTransIndices = new LinkedList<>(); //which transitions does this vertex lead to?
                        LinkedList<Pair<Integer, Integer>> queue = new LinkedList<>(); //(state, vertexHoldsTransIndex)


                        for (int trans = 0; trans < transNum; trans++) {
                            ArrayList<Integer> arr=new ArrayList<>();
                            arr.add(trans);
                            vertexHoldsTransIndices.add(arr);
                            queue.add(new Pair<Integer, Integer>(state, trans));
                        }

                        Pair<Integer, Integer> pair;
                        Pair<Integer, Integer> pair1;
                        Pair<Integer, Integer> pair2;
                        int tmpStateIndex;
                        ArrayList<Integer> tmpStateHoldsTransIndices;
                        int destination1;
                        int destination2;

                        while (queue.size() > 1) {
                            if (queue.size() == 2) { //So the root of the tree is still the actual vertex from the original graph
                                tmpState = state;
                            } else {
                                tmpState = newSTPG.addState(stpg.getPlayer(state));
                            }

                            tmpStateHoldsTransIndices=new ArrayList<>();
                            vertexHoldsTransIndices.add(tmpStateHoldsTransIndices);
                            queue.add(new Pair<>(tmpState, vertexHoldsTransIndices.size()-1));

                            d = new Distribution();
                            pair1 = queue.poll();
                            pair2 = queue.poll();
                            if (pair1.getKey() == state) {
                                //this only works because the index is also the index of destinations for the initial state (see trans above)
                                destination1=destinations[pair1.getValue()];
                            }
                            else {
                                destination1=pair1.getKey();
                            }

                            if (pair2.getKey() == state) {
                                destination2=destinations[pair2.getValue()];
                            }
                            else {
                                destination2=pair2.getKey();
                            }

                            double d1Sum=0;
                            double d2Sum=0;
                            double dSum;
                            for (Integer index : vertexHoldsTransIndices.get(pair1.getValue())) {
                                tmpStateHoldsTransIndices.add(index);
                                d1Sum+=probabilities[index];
                            }
                            for (Integer index : vertexHoldsTransIndices.get(pair2.getValue())) {
                                tmpStateHoldsTransIndices.add(index);
                                d2Sum+=probabilities[index];
                            }

                            dSum=d1Sum+d2Sum;
                            d.add(destination1, d1Sum/dSum);
                            d.add(destination2, d2Sum/dSum);

                            newSTPG.addChoice(tmpState, d);
                        }
                    }
                }
            }
        }


        return newSTPG;
    }

    /**
     * Precision: 4 Digits of the probability get considered
     * @param stpg
     * @return
     */
    private static STPG createHalfProbsForm(STPG stpg) {

        STPGExplicit newSTPG = new STPGExplicit();
        int stateNum = stpg.getNumStates();
        int choiceNum;
        int transNum;
        Distribution d;

        for (int state=0; state<stateNum; state++) {
            newSTPG.addState(stpg.getPlayer(state));
        }
        for (Integer initstate : stpg.getInitialStates()) {
            newSTPG.initialStates.add(initstate);
        }
        for (int state=0; state<stateNum; state++) {
            choiceNum = stpg.getNumChoices(state);
            for (int choice = 0; choice < choiceNum; choice++) {
                transNum = stpg.getNumTransitions(state, choice);
                if (transNum == 1) {
                    for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext(); ) {
                        Map.Entry<Integer, Double> tr = it.next();
                        d = new Distribution();
                        d.add(tr.getKey(), tr.getValue());
                        newSTPG.addChoice(state, d);
                    }
                }
                else {
                    if (transNum == 2) {
                        double prob1;
                        double prob2;
                        int dest1;
                        int dest2;
                        Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice);
                        Map.Entry<Integer, Double> tr = it.next();
                        dest1 = tr.getKey();
                        prob1 = tr.getValue();
                        tr = it.next();
                        dest2 = tr.getKey();
                        prob2 = tr.getValue();
                        if (prob1==0.5) {
                            d=new Distribution();
                            d.add(dest1, prob1);
                            d.add(dest2, prob2);
                            newSTPG.addChoice(state, d);
                            continue;
                        }
                    }
                    ArrayList<Double> probs = new ArrayList<>();
                    ArrayList<Integer> dest = new ArrayList<>();

                    for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext();) {
                        Map.Entry<Integer, Double> tr = it.next();
                        dest.add(tr.getKey());
                        probs.add(tr.getValue());
                    }

                    int tmpState;
                    LinkedList<Integer> queueOfStates = new LinkedList<>(); //used to create tree
                        /*
                        Idea: We need to find common denominator q' of all the transitions.
                        To do so we would need to first to represent the floats as fractions
                        This is barely possible with repeating decimals (1/3 = 0.333333...)
                        So we have to make a cut somewhere. Even if q is 1000 (3 digits)
                        This would mean that for EVERY action where we don't have 1/2 1/2-transition
                        we would create 1024 leafs (and 1022 inner nodes). So this is a lot.
                        We would need to find a good way to represent floats as rationals and
                        where the smallest denominator possible is taken
                         */

                    //for now take first 3 digits
                    ArrayList<Integer> bounds = new ArrayList<>();
                    int precision = 1000;
                    for (Double prob : probs) {
                        bounds.add((int) (precision * prob));
                    }
                    int divisor = bounds.get(0);
                    for (int bound : bounds) {
                        divisor = gcd(divisor, bound); //get gcd of all states
                    }

                    int tmpNum;
                    int boundSum=0;
                    for (int i=0; i<bounds.size(); i++) {
                        tmpNum=bounds.get(i)/divisor;
                        boundSum+=tmpNum;
                        bounds.set(i,tmpNum);
                    }

                    int q = 1;
                    while (q< boundSum) {
                        q = q*2;
                    }

                    for (int i=0; i<bounds.size(); i++) {
                        int bound = bounds.get(i);
                        for (int j=0; j<bound; j++) {
                            queueOfStates.add(dest.get(i));
                        }
                    }
                    int backtrans = q - boundSum;
                    for (int i = 0; i < backtrans; i++) {
                        queueOfStates.add(state);
                    }
                    int dest1;
                    int dest2;
                    while (queueOfStates.size() > 1) {
                        if (queueOfStates.size() == 2) {
                            tmpState = state;
                        } else {
                            tmpState = newSTPG.addState();
                        }

                        dest1=queueOfStates.poll();
                        dest2=queueOfStates.poll();
                        d=new Distribution();
                        d.add(dest1, 0.5);
                        d.add(dest2, 0.5);

                        newSTPG.addChoice(tmpState, d);
                        queueOfStates.add(tmpState);
                    }

                }
            }
        }

        return newSTPG;
    }

    /**
     * Adds for EVERY transition the possibility to get to a 0-sink so that this is a 1/2^m - stopping game
     * @param stpg
     * @param m
     * @param no
     * @param yes
     * @return
     * @throws PrismException
     */
    private static STPG toStopping(STPG stpg, int m, BitSet no, BitSet yes) throws PrismException{
        //If there is no 0-sink, one must be created!! Set it then also in no

        STPGExplicit newSTPG=new STPGExplicit();
        int stateNum=stpg.getNumStates();
        int choiceNum;
        int transNum;
        int goal=0; //shouldn't be initialized
        int[] mStates;
        Distribution d;

        for (int state=0; state<stateNum; state++) {
             newSTPG.addState(stpg.getPlayer(state));
        }

        int yesSink = yes.nextSetBit(0);
        int noSink = no.nextSetBit(0);
        if (noSink == -1) {
            noSink = newSTPG.addState(2);
            d = new Distribution();
            d.add(noSink, 1.0);
            newSTPG.addChoice(noSink, d);
            no.set(noSink);
        }

        for (Integer initstate : stpg.getInitialStates()) {
            newSTPG.initialStates.add(initstate);
        }
        for (int state=0; state<stateNum; state++) {
            choiceNum = stpg.getNumChoices(state);
            if (yes.get(state)) {
                d=new Distribution();
                d.add(yesSink, 1.0);
                newSTPG.addChoice(state,d);
                continue;
            }
            else if (no.get(state)) {
                d=new Distribution();
                d.add(noSink, 1.0);
                newSTPG.addChoice(state,d);
                continue;
            }
            for (int choice=0; choice<choiceNum; choice++) {
                addMStatesPerAction(stpg, newSTPG, m, state, choice, noSink, 0.5);
            }
        }
        return newSTPG;
    }


    /**
     * No return necessary
     * @param stpg
     * @param newSTPG
     * @param m
     * @param state
     * @param choice
     * @param noSink
     */
    private static void addMStatesPerAction(STPG stpg, STPGExplicit newSTPG, int m, int state, int choice, int noSink, double probGoIntoChain) {
        int[] mStates;
        int transNum;
        int goal=0;
        Distribution d;



        mStates = new int[m];
        transNum = stpg.getNumTransitions(state, choice);

        for (int i = 0; i<m; i++) {
            mStates[i]=newSTPG.addState(stpg.getPlayer(state)); //it doesn't matter which player, as the vertex has no choice
        }

        if (transNum>1) {
            goal=newSTPG.addState(stpg.getPlayer(state));
            d=getDistribution(stpg, state, choice);
            newSTPG.addChoice(goal, d); //bring in extra vertex to maintain 1/2 probabilities
        }
        else {
            for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext(); ) {
                Map.Entry<Integer, Double> tr = it.next();
                goal = tr.getKey();
            }
        }

        addMChain(newSTPG, mStates, goal, state, noSink, probGoIntoChain);
    }

    /**
     * We can add AFTER the transition probability a state with 1/2 to goal and 1/2 m states
     * @param stpg
     * @param newSTPG
     * @param m
     * @param state
     * @param choice
     * @param noSink
     */
    private static void addMStatesPerTransition(STPG stpg, STPGExplicit newSTPG, int m, int state, int choice, int noSink, double probGoIntoChain) {
        int[] mStates;
        int transNum;
        int goal;



        mStates = new int[m];
        transNum = stpg.getNumTransitions(state, choice);

        Distribution distStateWithToken = new Distribution();
        for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext();) {
            Map.Entry<Integer, Double> tr = it.next();
            goal = tr.getKey();
            for (int i = 0; i<m; i++) {
                mStates[i]=newSTPG.addState(stpg.getPlayer(state)); //it doesn't matter which player, as the vertex has no choice
            };
            if (transNum>1) {
                //Always go to condonState and from there on have m-transitions
                int condonState = newSTPG.addState(stpg.getPlayer(state));
                distStateWithToken.add(condonState, tr.getValue());
                addMChain(newSTPG, mStates, goal, condonState, noSink, probGoIntoChain);

            }
            else {
                addMChain(newSTPG, mStates, goal, state, noSink, probGoIntoChain);
            }
        }
        //if == 1 then handled in addMChain
        if (transNum>1) {
            newSTPG.addChoice(state, distStateWithToken);
        }
    }

    /**
     *
     * @param newSTPG
     * @param mStates
     * @param goal
     * @param stateToChain
     * @param noSink
     * @param probGoingIntoChain This is the probability that you will NOT reach the goal in one step (in Condons approach 1/2)
     */
    private static void addMChain(STPGExplicit newSTPG, int[] mStates, int goal, int stateToChain, int noSink, double probGoingIntoChain) {
        double failure = probGoingIntoChain;
        double success = 1.0-failure;
        Distribution d=new Distribution();
        d.add(mStates[0], failure);
        d.add(goal, success);
        newSTPG.addChoice(stateToChain, d);

        d=new Distribution();
        d.add(noSink, failure);
        d.add(goal, success);
        newSTPG.addChoice(mStates[mStates.length-1], d);

        for (int i=0; i<mStates.length-1; i++) {
            d=new Distribution();
            d.add(mStates[i+1], failure);
            d.add(goal, success);
            newSTPG.addChoice(mStates[i],d);
        }
    }

    private static void addEpsToTransition(STPG stpg, STPGExplicit newSTPG, double eps, int state, int choice, int noSink) {
        double succ;
        Distribution d=new Distribution();
        for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext();) {
            Map.Entry<Integer, Double> tr = it.next();
            succ = tr.getValue() - eps;
            d.add(tr.getKey(), succ);
            d.add(noSink, eps);
        }
        newSTPG.addChoice(state, d);

    }


    private static STPGExplicit toMECStopping(STPG stpg, List<BitSet> mecs, double goIntoChain, BitSet no, boolean mStates) {
        STPGExplicit newSTPG = new STPGExplicit();
        int stateNum = stpg.getNumStates();
        int noSink = no.nextSetBit(0);
        int choiceNum;
        int m;
        boolean inMec;
        boolean allToMec;

        for (int state=0; state<stateNum; state++) {
            newSTPG.addState(stpg.getPlayer(state));
        }
        for (Integer initstate : stpg.getInitialStates()) {
            newSTPG.initialStates.add(initstate);
        }

        for (int state=0; state<stateNum; state++) {
            inMec = false;
            choiceNum = stpg.getNumChoices(state);
            for (BitSet mec : mecs) {
                if (mec.get(state)) {
                    inMec=true;
                    for (int choice=0; choice<choiceNum; choice++) {
                        allToMec=true;
                        for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext(); ) {
                            Map.Entry<Integer, Double> tr = it.next();
                            if (!mec.get(tr.getKey())) {
                                allToMec=false;
                                break; //action not part of MEC
                            }
                        }
                        if (allToMec) { //MEC-Action -> include EPS
                            m = 2 * mec.cardinality() + 1;
                            //m = 2 * stpg.getNumStates() + 1;
                            if (mStates) {
                                addMStatesPerTransition(stpg, newSTPG, m, state, choice, noSink, goIntoChain);
                            }
                            else {
                                addEpsToTransition(stpg, newSTPG, Math.pow(goIntoChain, m), state, choice, noSink);
                            }
                        }
                        else { //Not MEC-Action -> leave it as it is
                            Distribution d = getDistribution(stpg, state, choice);
                            newSTPG.addChoice(state, d);
                        }

                    }
                    break; //a state belongs to atmost one MEC
                }
            }
            if (!inMec) { //just add his action and leave the state be
                for (int choice=0; choice<choiceNum; choice++) {
                    Distribution d = getDistribution(stpg, state, choice);
                    newSTPG.addChoice(state, d);
                }
            }
        }


        return newSTPG;
    }

    private static void addUselessAction(STPGExplicit stpg, BitSet no, BitSet yes) throws PrismException {
        int numState=stpg.getNumStates();
        int numChoices;
        for (int state=0; state<numState; state++) {
            numChoices=stpg.getNumChoices(state);
            if (numChoices==1) {
                int yesSink = yes.nextSetBit(0);
                int noSink = no.nextSetBit(0);
                Distribution distr = new Distribution();
                if (!no.get(state) && !yes.get(state)) {

                    if (stpg.getPlayer(state) == 1) {
                        //If there are no no-sinks the value of this node must be 1
                        //Because the no-set also computes min-player-only-MECs
                        if (noSink == -1) {
                            distr.add(yesSink, 1.0);
                        } else {
                            distr.add(noSink, 1.0);
                        }
                    } else {
                        distr.add(yesSink, 1.0);
                    }
                    stpg.addChoice(state, distr);
                }
            }
        }
    }

    // non-recursive implementation of Eucledean GCD
    private static int gcd(int p, int q) {
        while (q != 0) {
            int temp = q;
            q = p % q;
            p = temp;
        }
        return p;
    }
 }
