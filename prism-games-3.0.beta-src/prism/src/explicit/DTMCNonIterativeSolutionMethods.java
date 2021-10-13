package explicit;

import Jama.Matrix;
import explicit.rewards.MCRewards;
import prism.Pair;
import prism.PrismException;

import java.lang.reflect.Array;
import java.util.*;

public class DTMCNonIterativeSolutionMethods {

    public static long inverseCalcTime = 0;
    /**
     * Returns the Matrix of leaving probabilities of states, where entry ij is the probability that
     * state i will take exit j.
     *
     *
     * This is interpreted in a ModelCheckerResult and the transition probs are returned
     */
    protected ModelCheckerResult solveMarkovChain(DTMC dtmc, BitSet targets, double[] upperbounds, double precision) {

        ModelCheckerResult result = new ModelCheckerResult();
        result.soln = new double[dtmc.getNumStates()];
        result.numIters = 0;

        for (int target = targets.nextSetBit(0); target >= 0; target = targets.nextSetBit(target +1)) {
            result.soln[target] = 1;
        }

        /*
        Create a (n+m) x (n+m) Transition Probability Matrix M, where n+m the number of state in dtmc and m the number of targets.
        Thus, n is the rest of states
         */

        // Nasty preprocessing step. Maybe can be done better
        HashMap<Integer, Integer> normalStateToMatrixColumnIndex = new HashMap<>();
        HashMap<Integer, Integer> targetStateToMatrixColumnIndex = new HashMap<>();

        for (int state = 0; state < dtmc.getNumStates(); state++) {
            if (targets.get(state)) {
                targetStateToMatrixColumnIndex.put(state, targetStateToMatrixColumnIndex.size());
            }
            else {
                normalStateToMatrixColumnIndex.put(state, normalStateToMatrixColumnIndex.size());
            }
        }


        int stateNum = dtmc.getNumStates() - targets.cardinality();
        int n = dtmc.getNumStates() - targets.cardinality();
        int m = targets.cardinality();
        double[][] transInDTMC = new double[n][n];
        double[][] transToTarget = new double[n][m];

        for (int state = 0; state < dtmc.getNumStates(); state++) {
            // We don't have for targets transition probabilities
            if (targets.get(state)) {
                continue;
            }

            int matrixRowIndex = normalStateToMatrixColumnIndex.get(state);

            for (Iterator<Map.Entry<Integer, Double>> it = dtmc.getTransitionsIterator(state); it.hasNext(); ) {
                Map.Entry<Integer, Double> transition = it.next();
                // Belongs to Target-part
                if (targets.get(transition.getKey())) {
                    transToTarget[matrixRowIndex][targetStateToMatrixColumnIndex.get(transition.getKey())] = transition.getValue();
                }
                // Leads back into DTMC
                else {
                    transInDTMC[matrixRowIndex][normalStateToMatrixColumnIndex.get(transition.getKey())] = transition.getValue();
                }
            }
        }

        Matrix A = new Matrix(transInDTMC); //n States, n States -> nxn
        Matrix P = new Matrix(transToTarget); //n States, m Targets -> nxm

        //Check whether P isn't simply the 0-matrix -> every target is unreachable
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
        if (!unequal) {
            // No need to modify solutions since 0 anyways
            return result;
        }

        ArrayList<Integer> delete;
        if (upperbounds != null) {
            delete = getRemovableRowsFast(normalStateToMatrixColumnIndex, upperbounds, precision);
        }
        else {
            delete = getRemovableRows(A, P);
        }
        Matrix Ashort=shortenMatrix(A, delete);

        /*
        ArrayList<Integer> delete1 = getRemovableRowsFast(normalStateToMatrixColumnIndex, upperbounds, precision);
        ArrayList<Integer> delete2 = getRemovableRows(A, P);

        System.out.println("Equal? "+(delete1.equals(delete2)));
        if (!delete1.equals(delete2)) {
            System.out.println("Delete 1: "+Arrays.toString(delete1.toArray()));
            System.out.println("Delete 2: "+Arrays.toString(delete2.toArray()));
            System.out.println("Reach Target: "+getReaching(targets.nextSetBit(0), dtmc));
            delete2 = getRemovableRows(A, P);
            getProbabilities(0, dtmc);
            whereLeads(0, dtmc);
        }
        */

        Matrix M = Matrix.identity(Ashort.getRowDimension(), Ashort.getColumnDimension()).minus(Ashort);
        long t1 = System.currentTimeMillis();
        Matrix MInverse=M.inverse();
        long t2 = System.currentTimeMillis();
        inverseCalcTime += (t2-t1);
        MInverse=extendMatrix(MInverse, delete);


        Matrix res = MInverse.times(P); // nxn * nxm -> nxm

        //System.out.println("Resulting Matrix:\n"+matrixToString(res));

        for (int state = 0; state < dtmc.getNumStates(); state++) {
            if (targets.get(state)) continue;

            int row = normalStateToMatrixColumnIndex.get(state);
            for (int target = targets.nextSetBit(0); target >= 0; target = targets.nextSetBit(target +1)) {
                int column = targetStateToMatrixColumnIndex.get(target);

                result.soln[state] += res.get(row, column);
            }

            // Normalize Value
            result.soln[state] = Math.min(result.soln[state], 1.0);
            result.soln[state] = Math.max(result.soln[state], 0.0);
        }

        return result;
    }

    /**
     * Note that due to Floating Point Arithmetics, this is also not 100% exact.
     * @return
     */
    public ModelCheckerResult computeReachProbsExact(DTMC dtmc, MCRewards mcRewards, BitSet target, BitSet inf, double init[], BitSet known, double roundFrom) throws PrismException {
        return solveMarkovChain(dtmc, target, init, roundFrom);
    }

    protected ArrayList<Integer> getRemovableRowsFast(HashMap<Integer, Integer> normalStateToMatrixColumnIndex, double[] upperBounds, double roundFrom) {
        // Assume that every
        ArrayList<Integer> delete=new ArrayList<>();
        for (int state = 0; state < upperBounds.length; state++) {
            // This number may be set up a little higher due to imprecision
            if (upperBounds[state] <= roundFrom) {
                if (normalStateToMatrixColumnIndex.containsKey(state)) delete.add(normalStateToMatrixColumnIndex.get(state));
            }
        }

        return delete;
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

    public ArrayList<Integer> getPredecessor(int state, DTMC dtmc) {
        ArrayList<Integer> pred = new ArrayList<>();
        for (int s = 0; s < dtmc.getNumStates(); s++) {
            for (Iterator<Map.Entry<Integer, Double>> it = dtmc.getTransitionsIterator(s); it.hasNext(); ) {
                Map.Entry<Integer, Double> t = it.next();
                if (t.getKey() == state && t.getValue() > 0) {
                    pred.add(s);
                }
            }
        }
        return pred;
    }

    public void whereLeads(int state, DTMC dtmc) {
        int currentState = state;
        LinkedList<Integer> nextStates = new LinkedList<>();
        nextStates.add(state);
        BitSet seen = new BitSet();
        seen.set(state);
        while (seen.cardinality() < dtmc.getNumStates() && !nextStates.isEmpty()) {
            currentState = nextStates.poll();
            System.out.println("See state: "+currentState);
            seen.set(currentState);

            for (Iterator<Map.Entry<Integer, Double>> it = dtmc.getTransitionsIterator(currentState); it.hasNext(); ) {
                Map.Entry<Integer, Double> t = it.next();
                if (t.getValue() > 0 && !nextStates.contains(t.getKey()) && !seen.get(t.getKey())) {
                    nextStates.addLast(t.getKey());
                }
            }
        }
    }

    public HashMap<Integer, Double> getProbabilities(int state, DTMC dtmc) {
        HashMap<Integer, Double> probs = new HashMap<Integer, Double>();
        for (Iterator<Map.Entry<Integer, Double>> it = dtmc.getTransitionsIterator(state); it.hasNext(); ) {
            Map.Entry<Integer, Double> t = it.next();
            if (!probs.containsKey(t.getKey())) {
                probs.put(t.getKey(), t.getValue());
            }
            else {
                probs.put(t.getKey(), probs.get(t.getKey()) + t.getValue());
            }
        }
        for (Map.Entry<Integer, Double> pair : probs.entrySet()) {
            System.out.println("Key: "+pair.getKey()+", Value: "+pair.getValue());
        }
        return probs;
    }

    public List<Integer> getReaching(int state, DTMC dtmc) {
        ArrayList<Integer> pred = new ArrayList<>();
        ArrayList<Integer> nextIter = new ArrayList<>();
        pred.add(state);
        while (true) {
            for (int s : pred) {
                ArrayList<Integer> sPreds = getPredecessor(s, dtmc);
                for (int x : sPreds) {
                    if (!nextIter.contains(x)) {
                        nextIter.add(x);
                    }
                }
            }
            if (nextIter.equals(pred)) {
                break;
            }
            else {
                pred.clear();
                pred.addAll(nextIter);
                nextIter.clear();
            }
        }
        Collections.sort(pred);
        return pred;
    }
}
