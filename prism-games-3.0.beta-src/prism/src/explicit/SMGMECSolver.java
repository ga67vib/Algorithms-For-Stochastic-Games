package explicit;

//import com.sun.xml.internal.bind.v2.runtime.BinderImpl;
import gurobi.*;
import prism.Pair;
import prism.PrismException;
import Jama.*;

import java.lang.reflect.Array;
import java.util.*;

public abstract class SMGMECSolver {

    protected List<BitSet> mecs;
    protected STPG stpg;
    BitSet known;
    double[] knownValues;

    public SMGMECSolver(List<BitSet> mecs, STPG stpg, BitSet known, double[] knownValues) {
        this.mecs=mecs;
        this.stpg=stpg;
        this.known=known;
        this.knownValues=knownValues;
    }

    //
    public abstract void solve () throws PrismException, GRBException;

    protected void createStrategyTable(BitSet mec) throws PrismException, GRBException {
        boolean foundMin = false;
        boolean foundMax = false;
        boolean foundProbs = false; //are there any transitions that have probability less than 1?
        int choiceNum;
        for (int state = mec.nextSetBit(0); state >= 0; state = mec.nextSetBit(state + 1)) {
            if (state == Integer.MAX_VALUE) {
                break;
            }
            if (stpg.getPlayer(state) == 1) {
                foundMax = true;
            } else {
                foundMin = true;
            }

            if (foundMin && foundMax) {
                break;
            }

        }

        BitSet exits = getMECExits(mec);

        if (foundMin && !foundMax) {
            setMECtoZero(mec);
            return;
        } else if (foundMax && !foundMin) {
            setMECtoBestExit(mec, exits);
            return;
        }

        handleMixedMEC(mec, exits);
    }

    protected abstract void handleMixedMEC(BitSet mec, BitSet exits) throws GRBException, PrismException;

    protected abstract void setMECtoZero(BitSet mec) throws GRBException, PrismException;

    protected abstract void setMECtoBestExit(BitSet mec, BitSet exits) throws GRBException, PrismException;

    // print matrix to standard output
    protected String matrixToString(Matrix M) {
        String s="";
        for (int i = 0; i < M.getRowDimension(); i++) {
            for (int j = 0; j < M.getColumnDimension(); j++) {
                String num = Double.toString(M.get(i,j));
                if (num.length()>4) {
                    num = num.substring(0, 4);
                }
                s += num + "\t";
            }
            s+="\n";
        }
        return s;
    }

    /**
     * Gets all the Exits the MEC might have
     */
    protected BitSet getMECExits(BitSet mec) {
        BitSet exits = new BitSet();
        for (int state = mec.nextSetBit(0); state>=0; state=mec.nextSetBit(state+1)) {
            if (state==Integer.MAX_VALUE) {
                break;
            }

            int choiceNum = stpg.getNumChoices(state);
            for (int choice = 0; choice < choiceNum; choice++) {
                for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext(); ) {
                    Map.Entry<Integer, Double> tr = it.next();
                    if (!mec.get(tr.getKey()) && !exits.get(tr.getKey())) {
                        exits.set(tr.getKey());
                    }
                }
            }
        }
        return exits;
    }

    protected ArrayList<Integer> getRemovableRows(Matrix A, Matrix P) {
        int sizeA=A.getColumnDimension();
        int sizeP=P.getColumnDimension();
        int size=sizeA+sizeP;

        ArrayList<Integer> delete=new ArrayList<>();

        double[][] reachabilityMatrix=floydWarshall(A, P);


        boolean reachesExit;
        for (int i=0; i<A.getRowDimension(); i++) {
            reachesExit=false;
            for (int j=sizeA; j<size; j++) {
                if (reachabilityMatrix[i][j]!=Double.POSITIVE_INFINITY) {
                    reachesExit=true;
                    break;
                }
            }
            if (!reachesExit) {
                delete.add(i);
            }
        }

        return delete;
    }

    /**
     * We want to know whether the vertices in A can even reach P. Because if they don't They can be cut out
     * @param A
     * @param P
     */
    protected double[][] floydWarshall(Matrix A, Matrix P)
    {
        int sizeA=A.getColumnDimension();
        int sizeP=P.getColumnDimension();
        int size=sizeA+sizeP;
        double[][] dist = new double[size][size];

        int i, j, k;

        /* Initialize the solution matrix same as input graph matrix.
           Or we can say the initial values of shortest distances
           are based on shortest paths considering no intermediate
           vertex. */
        Matrix tmp;
        int index;

        for (i = 0; i<A.getRowDimension(); i++) {
            for (int offset = 0; offset<2; offset++) {
                if (offset==0) {
                    tmp=A;
                }
                else {
                    tmp=P;
                }
                index=0;
                for (j=sizeA*offset; j<sizeA + sizeP*offset; j++) {
                    double transProb=tmp.get(i,index);
                    if (transProb==0.0) {
                        dist[i][j]= Double.POSITIVE_INFINITY;
                    }
                    else {
                        dist[i][j] = tmp.get(i, index);
                    }
                    index++;
                }
            }
        }
        for (i=A.getRowDimension(); i<size; i++) {
            for (j=0; j<size; j++) {
                if (i==j) {
                    dist[i][i]=1.0; //the last rows only consist of the sinks having self loops
                }
                else {
                    dist[i][j]=Double.POSITIVE_INFINITY;
                }
            }
        }

        /* Add all vertices one by one to the set of intermediate
           vertices.
          ---> Before start of an iteration, we have shortest
               distances between all pairs of vertices such that
               the shortest distances consider only the vertices in
               set {0, 1, 2, .. k-1} as intermediate vertices.
          ----> After the end of an iteration, vertex no. k is added
                to the set of intermediate vertices and the set
                becomes {0, 1, 2, .. k} */
        for (k = 0; k < size; k++)
        {
            // Pick all vertices as source one by one
            for (i = 0; i < size; i++)
            {
                // Pick all vertices as destination for the
                // above picked source
                for (j = 0; j < size; j++)
                {
                    // If vertex k is on the shortest path from
                    // i to j, then update the value of dist[i][j]
                    if (dist[i][k] + dist[k][j] < dist[i][j])
                        dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
        return dist;
    }


    protected ArrayList<Pair<Integer, Integer>> getExitStateActionPairs(BitSet mec, BitSet exits) {
        ArrayList<Pair<Integer, Integer>> res=new ArrayList<>();
        for (int state = mec.nextSetBit(0); state>=0; state=mec.nextSetBit(state+1)) {
            if (state == Integer.MAX_VALUE) {
                break;
            }
            int numChoices = stpg.getNumChoices(state);
            for (int choice=0; choice<numChoices; choice++) {
                for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, choice); it.hasNext(); ) {
                    Map.Entry<Integer, Double> tr = it.next();
                    if (exits.get(tr.getKey())) {
                        res.add(new Pair<>(state, choice));
                        break;
                    }
                }
            }

        }

        return res;
    }

    protected void incrementAction(ArrayList<Integer> actions, ArrayList<Integer> states) {
        int i = 0;
        boolean done = false;
        while (!done && i < actions.size()) {
            int action = actions.get(i);
            if (stpg.getNumChoices(states.get(i)) - 1 == action) {
                actions.set(i, 0);
                i++;
            } else {
                actions.set(i, action + 1);
                done = true;
            }
        }
    }

    protected void reportTooManyActions(ArrayList<Integer> states, boolean max) throws  PrismException {
        String player = max ? "Maximizer" : "Minimizer";
        int choices = 0;
        for (Integer state : states) {
            choices += stpg.getNumChoices(state);
        }
        throw new PrismException(player + " has "+choices+" choices in MEC, which results in 2^("+choices+") possibilities to be tested." +
                "This would be most likely too big to be handled and causes furthermore an Integer Overflow. Try other solution methods.");
    }

    /**
     * Returns the Matrix of leaving probabilities of states, where entry ij is the probability that
     * state i will take exit j.
     * The Matrix is organized in the way that the first rows are taken up by the maxStates
     * (should be in ascending order) and then come the minStates
     */
    protected Matrix solveMarkovChain(ArrayList<Integer> maxStatesNames, ArrayList<Integer> maxStatesActions, ArrayList<Integer> minStatesNames, ArrayList<Integer> minStatesActions, BitSet exits) {

        int stateNum = maxStatesNames.size() + minStatesNames.size();
        double[][] transInEC = new double[stateNum][stateNum];
        double[][] transToExit = new double[stateNum][exits.cardinality()];
        int minOffset = maxStatesNames.size();

        ArrayList<Integer> states;
        ArrayList<Integer> actions;
        int index;
        int arrayIndex;
        for (int runs = 0; runs < 2; runs++) {
            if (runs == 0) {
                states = maxStatesNames;
                actions = maxStatesActions;
            } else {
                states = minStatesNames;
                actions = minStatesActions;
            }

            arrayIndex=0;
            for (index = maxStatesNames.size() * runs; index < maxStatesNames.size() + minStatesNames.size()*runs; index++) {
                int state = states.get(arrayIndex);
                int action = actions.get(arrayIndex);
                //System.out.println("state: "+state+", action: "+action);
                for (Iterator<Map.Entry<Integer, Double>> it = stpg.getTransitionsIterator(state, action); it.hasNext(); ) {
                    Map.Entry<Integer, Double> tr = it.next();
                    int to = tr.getKey();
                    double prob = tr.getValue();
                    int j = maxStatesNames.indexOf(to); //as every state name appears once, this should be no problem
                    if (j != -1) {
                        transInEC[index][j] = prob;
                    } else {
                        j = minStatesNames.indexOf(to);
                        if (j != -1) {
                            transInEC[index][j + minOffset] = prob;
                        } else { //this transition goes out of the MEC (so that it's an exit)
                            int exit=exits.nextSetBit(0);
                            for (j = 0; j < exits.cardinality(); j++) {
                                if (exit == to) {
                                    transToExit[index][j] = prob;
                                    break;
                                }
                                exit=exits.nextSetBit(exit+1);
                            }
                        }
                    }
                }
                arrayIndex++;
            }
        }
        Matrix A = new Matrix(transInEC); //n MecStates, n MecStates -> nxn
        Matrix P = new Matrix(transToExit); //n MecStates, m Exits -> nxm
        boolean unequal=false;
        for (int i=0; i<P.getRowDimension(); i++) {
            for (int j=0; j<P.getColumnDimension(); j++) {
                if (P.get(i,j)!=0.0) {
                    unequal=true;
                    break;
                }
            }
            if (unequal) {
                break;
            }
        }
        if (!unequal) return P; //If P is 0 Matrix then

        //System.out.println("---------------------------------");
        //System.out.println(maxStatesNames + "" + minStatesNames);
        //System.out.println(maxStatesActions + "" + minStatesActions);
        //System.out.println("---------------------------------");

        //System.out.println("Matrix A:\n"+matrixToString(A));
        //System.out.println("Matrix P:\n"+matrixToString(P));


        ArrayList<Integer> delete=getRemovableRows(A, P);
        Matrix Ashort=shortenMatrix(A, delete);
        Matrix M = Matrix.identity(Ashort.getRowDimension(), Ashort.getColumnDimension()).minus(Ashort);
        Matrix MInverse=M.inverse();
        MInverse=extendMatrix(MInverse, delete);


        Matrix res = MInverse.times(P); // nxn * nxm -> nxm

        //System.out.println("Resulting Matrix:\n"+matrixToString(res));

        return res;
    }

    protected boolean matrixEquals(Matrix A, Matrix B) {
        if (A.getRowDimension()!=B.getRowDimension()) return false;
        if (A.getColumnDimension()!=B.getColumnDimension()) return false;
        for (int i=0; i<A.getRowDimension(); i++) {
            for (int j=0; j<B.getColumnDimension(); j++) {
                if (A.get(i,j)!=B.get(i,j)) return false;
            }
        }
        return true;
    }

    protected Matrix shortenMatrix(Matrix m, ArrayList<Integer> rowsToDelete) {
        int size = m.getColumnDimension()-rowsToDelete.size();
        Matrix a=new Matrix(size, size);

        int iNew=0;
        int jNew;

        for (int i=0; i<m.getRowDimension(); i++) {
            jNew=0;
            if (!rowsToDelete.contains(i)) {
                for (int j=0; j<m.getColumnDimension(); j++) {
                    if (!rowsToDelete.contains(j)) {
                        a.set(iNew, jNew, m.get(i,j));

                        jNew++;
                    }
                }
                iNew++;
            }
        }

        //System.out.println("SHORTEND MATRIX U:\n"+matrixToString(a));

        return a;
    }

    protected Matrix extendMatrix(Matrix m, ArrayList<Integer> rowsToInsert) {
        int size = m.getColumnDimension()+rowsToInsert.size();
        Matrix a=new Matrix(size, size);

        int iNew=0;
        int jNew;

        for (int i=0; i<a.getRowDimension(); i++) {
            jNew=0;
            if (!rowsToInsert.contains(i)) {
                for (int j=0; j<a.getColumnDimension(); j++) {
                    if (!rowsToInsert.contains(j)) {
                        a.set(i, j, m.get(iNew,jNew));

                        jNew++;
                    }
                }
                iNew++;
            }
        }

        return a;
    }


}
