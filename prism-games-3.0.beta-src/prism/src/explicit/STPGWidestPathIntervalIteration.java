package explicit;

import common.IntSet;
import prism.PrismException;
import prism.PrismFileLog;
import prism.PrismLog;
import prism.PrismUtils;

import java.util.*;

public class STPGWidestPathIntervalIteration {
    STPGModelChecker modelChecker;

    boolean debug = false;

    long timeSpentDeflating = 0;

    public STPGWidestPathIntervalIteration(STPGModelChecker modelChecker) {
        this.modelChecker = modelChecker;
    }

    public ModelCheckerResult computeReachProbsValIter(STPG stpg, BitSet no, BitSet yes, boolean min1, boolean min2, double init[], BitSet known, ProbModelChecker.SolnMethod variant)
            throws PrismException {
        System.out.println("Doing value iteration variant " + variant + " with solnMethodOptions=" + modelChecker.solnMethodOptions + ", modelChecker.maxIters=" + modelChecker.maxIters + ", epsilon=" + modelChecker.termCritParam + " and topological=" + modelChecker.getDoTopologicalValueIteration());

        //solnMethodOptions is used to decide whether we use Gauss Seidel or not (%2==1 implies Gauss-Seidel) and for one optimization of OVI. See the iterateOnSubset-method.
        //The idea of solnMethodOptions is: first bit encodes GS. Second bit encodes OVI. Further bits are still free for further opts.


        ModelCheckerResult res = null;
        BitSet unknown;
        int s, n, iters;
        iters=0;
        // Note: For SVI, lowerBounds are stepBoundReach-values (x) and upperBounds are the stepBoundStay-values (y)
        double lowerBounds[], lowerBounds2[], upperBounds[], upperBounds2[], initVal;
        int adv[] = null;
        boolean genAdv;
        long timer;

        // Are we generating an optimal adversary?
        genAdv = modelChecker.exportAdv || modelChecker.generateStrategy;

        // Start value iteration
        timer = System.currentTimeMillis();
        if (modelChecker.verbosity >= 1)
            System.out.println("Starting value iteration (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")...");

        // Store num states
        n = stpg.getNumStates();

        // Determine set of states actually need to compute values for
        unknown = new BitSet();
        unknown.set(0, n);
        unknown.andNot(yes);
        unknown.andNot(no);
        // Maxi changed this on 22.06.20: known can be any lower bound we know, need not be precise; so we have to compute known again
//		if (known != null)
//			unknown.andNot(known);

        boolean needsUpperBounds = true;

        // Create solution vector(s)
        lowerBounds = new double[n];
        lowerBounds2 = new double[n];
        upperBounds = new double[n];
        upperBounds2 = new double[n];

        // Initialise solution vectors. Use (where available) the following in order of preference:
        // (1) exact answer, if already known; (2) 1.0/0.0 if in yes/no; (3) passed in initial value; (4) initVal
        // where initVal is 0.0 or 1.0, depending on whether we converge from below/above.
        // If needsUpperBounds, initVal is set to 0 for lower and 1 for upper
        // The valIterDir thing is only relevant for unguaranteed VI, which might also come from above
        initVal = 0;
        if (init != null) {
            if (known != null) {
                for (s = 0; s < n; s++)
                    lowerBounds[s] = lowerBounds2[s] = known.get(s) ? init[s] : yes.get(s) ? 1.0 : no.get(s) ? 0.0 : init[s];
            } else {
                for (s = 0; s < n; s++)
                    lowerBounds[s] = lowerBounds2[s] = yes.get(s) ? 1.0 : no.get(s) ? 0.0 : init[s];
            }
        } else {
            for (s = 0; s < n; s++)
                lowerBounds[s] = lowerBounds2[s] = yes.get(s) ? 1.0 : no.get(s) ? 0.0 : initVal;
        }

        for (s = 0; s < n; s++)
            upperBounds[s] = upperBounds2[s] = no.get(s) ? 0.0 : 1.0;


        SCCInfo sccs = null;

        if (modelChecker.getDoTopologicalValueIteration()){
            System.out.println("Getting topologically ordered SCCs...");
            sccs = SCCComputer.computeTopologicalOrdering(modelChecker, stpg, true, unknown::get);
            System.out.println("Done; got " + sccs.getNumSCCs() + " SCCs.");
        }


        // Create/initialise adversary storage
        if (genAdv) {
            adv = new int[n];
            for (s = 0; s < n; s++) {
                adv[s] = -1;
            }

            for (s = 0; s < no.length(); s++) {
                s = no.nextSetBit(s);
                for (int c = 0; c < stpg.getNumChoices(s); c++) {
                    if (stpg.allSuccessorsInSet(s, c, no)) {
                        adv[s] = c;
                        break;
                    }
                }
            }
        }

        // Need initstate to determine whether done in bounded case
        int initialState = stpg.getFirstInitialState();


        // Start iterations
        if (modelChecker.getDoTopologicalValueIteration()) {

            boolean useAttractorToGetDTMC = false;
            boolean useTOP = modelChecker.solnMethodOptions != 2;

            BitSet haveAlreadyFixedValues = new BitSet();
            haveAlreadyFixedValues.or(yes);
            // [Small optimization:] We don't care for sinks for the attractors so we won't set them
            int[] alreadyComputedAttractorDistances = new int[stpg.getNumStates()];
            Arrays.fill(alreadyComputedAttractorDistances, Integer.MAX_VALUE -10);
            for (int target = yes.nextSetBit(0); target >= 0; target = yes.nextSetBit(target + 1)) {
                alreadyComputedAttractorDistances[target] = 0;
            }

            for (int scc = 0; scc < sccs.getNumSCCs(); scc++) {
                if (sccs.isSingletonSCC(scc)) {
                    // get the single state in this SCC and finish it. Trivial.
                    int state = sccs.getStatesForSCC(scc).iterator().nextInt();
                    // finish doing that state in all vectors
                    lowerBounds[state] = stpg.mvMultJacMinMaxSingle(state, lowerBounds, min1, min2);
                    lowerBounds2[state] = stpg.mvMultJacMinMaxSingle(state, lowerBounds2, min1, min2);
                    if (upperBounds != null) {
                        upperBounds[state] = stpg.mvMultJacMinMaxSingle(state, upperBounds, min1, min2);
                        upperBounds2[state] = stpg.mvMultJacMinMaxSingle(state, upperBounds2, min1, min2);
                    }

                    BitSet statesForSCC = new BitSet();
                    statesForSCC.set(state);
                    double[] values;
                    if (useTOP) {

                        if (useAttractorToGetDTMC) {
                            int[][] startegyComputationResult = STPGValueIterationUtils.computeStrategyFromBounds(stpg, yes, lowerBounds, upperBounds, modelChecker.termCritParam, statesForSCC, haveAlreadyFixedValues, alreadyComputedAttractorDistances);
                            int[] sigma = startegyComputationResult[0];
                            int[] tau = startegyComputationResult[1];

                            values = STPGValueIterationUtils.getValueFromStrategies(modelChecker, stpg, sigma, tau, yes, no, statesForSCC, haveAlreadyFixedValues, lowerBounds, modelChecker.termCritParam, lowerBounds, upperBounds);
                            alreadyComputedAttractorDistances = startegyComputationResult[2];
                        } else {
                            values = STPGValueIterationUtils.getValueFromStrategiesWithoutAttractor(modelChecker, stpg, yes, no, statesForSCC, haveAlreadyFixedValues, lowerBounds, modelChecker.termCritParam, lowerBounds, upperBounds);
                        }

                        // Fix value of states
                        lowerBounds[state] = values[state];
                        lowerBounds2[state] = values[state];
                        if (upperBounds != null && upperBounds2 != null) {
                            upperBounds[state] = values[state];
                            upperBounds2[state] = values[state];
                        }
                    }

                    iters++;
                    IterationMethod.intervalIterationCheckForProblems(lowerBounds, upperBounds, IntSet.asIntSet(state).iterator());

                    if (useTOP) haveAlreadyFixedValues.set(state);
                } else {
                    // complex SCC: do VI
                    int itersInSCC = 0;

                    IntSet statesForSCCIntSet = sccs.getStatesForSCC(scc);
                    BitSet statesForSCC = new BitSet(n);
                    for (int state : statesForSCCIntSet){
                        statesForSCC.set(state);
                    }
                    //Solve the SCC until all states are close (initialState argument is -1 to ensure all states are solved, not just initial)
                    double[][] subres = iterateOnSubset((STPGExplicit) stpg, min1, min2, upperBounds, upperBounds2, lowerBounds2, lowerBounds, genAdv, adv, itersInSCC, variant, statesForSCC, statesForSCCIntSet, -1, yes);
                    itersInSCC = (int) subres[2][0];
                    lowerBounds = subres[0];
                    upperBounds = subres[1];

                    double[] values;
                    if (useTOP) {
                        if (useAttractorToGetDTMC) {
                            int[][] startegyComputationResult = STPGValueIterationUtils.computeStrategyFromBounds(stpg, yes, lowerBounds, upperBounds, modelChecker.termCritParam, statesForSCC, haveAlreadyFixedValues, alreadyComputedAttractorDistances);
                            int[] sigma = startegyComputationResult[0];
                            int[] tau = startegyComputationResult[1];

                            values = STPGValueIterationUtils.getValueFromStrategies(modelChecker, stpg, sigma, tau, yes, no, statesForSCC, haveAlreadyFixedValues, lowerBounds, modelChecker.termCritParam, lowerBounds, upperBounds);
                            alreadyComputedAttractorDistances = startegyComputationResult[2];
                        } else {
                            values = STPGValueIterationUtils.getValueFromStrategiesWithoutAttractor(modelChecker, stpg, yes, no, statesForSCC, haveAlreadyFixedValues, lowerBounds, modelChecker.termCritParam, lowerBounds, upperBounds);
                        }
                        // Fix value of states
                        for (int state = statesForSCC.nextSetBit(0); state >= 0; state = statesForSCC.nextSetBit(state + 1)) {
                            lowerBounds[state] = values[state];
                            lowerBounds2[state] = values[state];
                            if (upperBounds != null && upperBounds2 != null) {
                                upperBounds[state] = values[state];
                                upperBounds2[state] = values[state];
                            }
                        }
                    }

                    IterationMethod.intervalIterationCheckForProblems(lowerBounds, upperBounds, statesForSCCIntSet.iterator());
//					System.out.println("Non-trivial SCC done in " + itersInSCC + " many iterations");
                    iters+=itersInSCC;

                    if (useTOP) haveAlreadyFixedValues.or(statesForSCC);
                }
//				System.out.println("SCC done. Precision: [" + lowerBounds[sccs.getStatesForSCC(scc).iterator().nextInt()] + "," + upperBounds[sccs.getStatesForSCC(scc).iterator().nextInt()] + "]. \nIters: " + iters);
            }
        }
        else{
            long t1 = System.currentTimeMillis();
            double[][] subres = iterateOnSubset((STPGExplicit) stpg, min1, min2, upperBounds, upperBounds2, lowerBounds2, lowerBounds, genAdv, adv, iters, variant, unknown, null, initialState, yes);
            long t2 = System.currentTimeMillis();
            log("Time Spent on Iterate on Subset: "+((t2-t1)/1000.0));
            iters = (int) subres[2][0];
            lowerBounds = subres[0];
            upperBounds = subres[1];
        }


        // Finished value iteration
        timer = System.currentTimeMillis() - timer;
        if (modelChecker.verbosity >= 1) {
            System.out.print("Value iteration variant "+ variant + "(" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")");
            System.out.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");
            if (modelChecker.doTopologicalValueIteration) {
                System.out.println("--TOP Stats--");
                System.out.println("Get Suggested Actions: "+STPGValueIterationUtils.getAllowedActionsTime);
                System.out.println("Build DTMC: "+STPGValueIterationUtils.dtmcBuildingTime);
                System.out.println("Solve DTMC: "+STPGValueIterationUtils.dtmcSolvingTime);
                System.out.println("Inverse Calc Time: "+STPGValueIterationUtils.inverseCalcTime);
                System.out.println("Verify Value: "+STPGValueIterationUtils.verificationTime);
                System.out.println("Assign SCC Result : "+STPGValueIterationUtils.assignmentTime);
            }
        }

        // Non-convergence is an error (usually)
//		if (errorOnNonConverge) {
//			String msg = "Iterative method did not converge within " + iters + " iterations.";
//			msg += "\nConsider using a different numerical method or increasing the maximum number of iterations";
//			throw new PrismException(msg);
//		}

        // Store results/strategy
        res = new ModelCheckerResult();
        res.soln = lowerBounds;
        res.numIters = iters;
        res.timeTaken = timer / 1000.0;

        if (modelChecker.generateStrategy) {
            BitSet allStates = new BitSet();
            allStates.set(0, stpg.getNumStates()-1);
            int[][] startegyComputationResult = STPGValueIterationUtils.computeStrategyFromBounds(stpg, yes, lowerBounds, lowerBounds, modelChecker.termCritParam, allStates, null, null);
            int[] sigma = startegyComputationResult[0];
            int[] tau = startegyComputationResult[1];
			/*for (int state = 0; state < stpg.getNumStates(); state++) {
				System.out.println("State "+state+" choice: "+((stpg.getPlayer(state) == 1) ? sigma[state] : tau[state]));
			}*/
        }

        // Print adversary
        if (genAdv) {
            PrismLog out = new PrismFileLog(modelChecker.exportAdvFilename);
            if (modelChecker.exportAdvFilename.lastIndexOf('.') != -1 && modelChecker.exportAdvFilename.substring(modelChecker.exportAdvFilename.lastIndexOf('.') + 1).equals("dot")) {
                stpg.exportToDotFileWithStrat(out, null, adv);
            } else {
                for (s = 0; s < n; s++) {
                    out.println(s + " " + (adv[s] != -1 ? stpg.getAction(s, adv[s]) : "-"));
                }
                out.println();
            }
            out.close();
        }

        return res;
    }

    /**
     * Gets the whole context of VI + 2 things:
     * @param subset Set to iterate on; should be unknown if working on whole thing at once and SCC if doing topological VI
     * @param initialState used for checking done; if we only care about initstate, set to its index. If we need everything to be close (topological VI), then set to -1.
     * @return number of iterations; upper and lower bounds are modified as side effect: Changed by Maxi on 07.07.21, since side effect not working as expected
     * Sometimes it returns the previous iteration, thus resulting in result which is not eps-precise
     * Thus we return an array of double arrays: LowerBounds, UpperBounds and an array containing only iters
     */
    private double[][] iterateOnSubset(STPGExplicit stpg, boolean min1, boolean min2, double[] upperBounds, double[] upperBoundsNew, double[] lowerBoundsNew, double[] lowerBounds,
                                       boolean genAdv, int[] adv, int iters, ProbModelChecker.SolnMethod variant,
                                       BitSet subset, IntSet subsetAsIntSet, int initialState, BitSet yes) throws PrismException{
        iters = 0;
        boolean done = false;
        double tmpsoln[];

        long totalDeflating = 0;
        long t1;
        // Gauss-Seidel
        boolean doGS = ((modelChecker.solnMethodOptions%2)==1);
        if (doGS){System.out.println("Doing Gauss Seidel variant");}

        while (!done) {
            iters++;

            if (doGS) {
                // For Gauss Seidel, let the new objects have the most up-to-date values, then only work on the new ones
                lowerBoundsNew = lowerBounds.clone();
                upperBoundsNew = upperBounds != null ? upperBounds.clone() : null;
            }
            //Bellman Updates

            // Update lower and upper bounds for classical Interval Iteration
            if (!doGS) {
                //Non-Gauss Seidel - Save the result in distinct array
                stpg.mvMultMinMax(lowerBounds, min1, min2, lowerBoundsNew, subset, false, genAdv ? adv : null);
                stpg.mvMultMinMax(upperBounds, min1, min2, upperBoundsNew, subset, false, genAdv ? adv : null);

            } else {
                //Gauss Seidel - Do operations in-place
                stpg.mvMultGSMinMax(lowerBoundsNew, min1, min2, subset, false, modelChecker.termCrit == ProbModelChecker.TermCrit.ABSOLUTE);
                stpg.mvMultGSMinMax(upperBoundsNew, min1, min2, subset, false, modelChecker.termCrit == ProbModelChecker.TermCrit.ABSOLUTE);
            }

            t1 = System.currentTimeMillis();
            double[][] wpDeflatingResult = STPGValueIterationUtils.widestPathDeflating(stpg, iters, 5, yes, lowerBounds, lowerBoundsNew, upperBoundsNew);
            totalDeflating += (System.currentTimeMillis() - t1);
            lowerBounds = wpDeflatingResult[0];
            lowerBoundsNew = wpDeflatingResult[1];
            upperBoundsNew = wpDeflatingResult[2];

            // Swap vectors for next iter
            // Now lowerBounds is the most up-to-date approximation, while the lowerBoundsNew contains the previous iteration
            // If lowerBounds were not changed no need to change anything
            tmpsoln = lowerBounds;
            lowerBounds = lowerBoundsNew;
            lowerBoundsNew = tmpsoln;



            tmpsoln = upperBounds;
            upperBounds = upperBoundsNew;
            upperBoundsNew = tmpsoln;

            /**
             * Check termination
             */
            // Interval Iteration Criterion
            done = (initialState != -1) ?
                    upperBounds[initialState] - lowerBounds[initialState] < this.modelChecker.termCritParam :
                    PrismUtils.doublesAreClose(lowerBounds, upperBounds, subsetAsIntSet.iterator(), modelChecker.termCritParam, modelChecker.termCrit == ProbModelChecker.TermCrit.ABSOLUTE);

            if (done) {
                System.out.println("Upper and Lower bounds are eps-close -> BVI iteration criterion met");
                break; // End the while loop
            }
        }



        if(initialState!=-1)
            System.out.println("Resulting interval: ["+lowerBounds[initialState]+","+((upperBounds!=null) ? upperBounds[initialState] : "none")+"]");

        timeSpentDeflating = totalDeflating;
        System.out.println("Time spent Deflating: "+(((double)timeSpentDeflating) / 1000.0)+"s");
        return new double[][]{lowerBounds,upperBounds,{iters}};
    }

    /**
     * KKKW18 deflate operation.
     */
    private double[][] deflate(STPGExplicit stpg, boolean min1, boolean min2, double[] lowerBounds, double[] upperBounds, BitSet mec, ECComputerDefault ec) throws PrismException {

        //TODO: I might turn on repeated adjustment again
        //TODO: I might optimize for MDP handling again (i.e. not recompute SECs and stuff)

        // Find all SECs in given MEC
        List<BitSet> SECs = ec.getSECs(mec, lowerBounds,min1,min2);
        for (int j = 0; j < SECs.size(); j++) {
            // Get best Maximizer exit from SEC
            BitSet sec = SECs.get(j);
            int maxPlayer = min1 ? (min2 ? -1 : 2) : (min2 ? 1 : 3);
            double bestLeavingUpperBound = getBestLeavingValue(sec, stpg, upperBounds, maxPlayer);
            // And deflate all states in SEC
            for (int s = sec.nextSetBit(0); s >= 0; s = SECs.get(j).nextSetBit(s+1)) {
                double formerValue = upperBounds[s];
                if(formerValue>bestLeavingUpperBound) {//monotonic: only decrease if new value smaller
                    upperBounds[s]=bestLeavingUpperBound;
                }
            }
        }
        return new double[][]{upperBounds,lowerBounds};
    }

    /**
     *
     * @param sec
     * @param stpg
     * @param vector
     * @param maxPlayer -1: Both want to minimize => return 0. 1 or 2: Respective player. 3: Both
     * @return
     */
    private static double getBestLeavingValue(BitSet sec, STPGExplicit stpg, double[] vector, int maxPlayer){
        double bestUpperBoundSoFar = 0;
        //find best outgoing upper bound belonging to maxPlayer
        for (int s = sec.nextSetBit(0); s >= 0; s = sec.nextSetBit(s+1)) {
//			//if we find target in simBCEC, all states in there should have value 1.
//			//should not happen since prob1 takes care of this, but let's be safe
//			if (vector[s]==1) {
//				bestUpperBoundSoFar=1;
//				break;
//			}
            if (stpg.getPlayer(s)==maxPlayer || maxPlayer == 3) {
                //System.out.println("Searching for best leaving value; state belongs to maximizer");
                for (int i = 0; i < stpg.getNumChoices(s); i++) {
                    boolean all = stpg.allSuccessorsInSet(s, i, sec);
                    //System.out.println("Action " + i + " all succ in set? " + all);
                    if (!all) {
                        double upperBound = 0;
                        for (int succ : stpg.getChoice(s,i).keySet()){
                            upperBound += stpg.getChoice(s,i).get(succ) * vector[succ];
                        }
                        if (upperBound>bestUpperBoundSoFar){
                            bestUpperBoundSoFar = upperBound;
                        }
                    }
                }
            }
            //TODO: We might want to deflate minimizer stuff to 0. Should be handled by prob0 though.
        }

        return bestUpperBoundSoFar;
    }

    private void log(String s) {
        if (debug) {
            System.out.println(s);
        }
    }
}
