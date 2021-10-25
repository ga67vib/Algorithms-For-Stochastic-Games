package explicit;

import common.IntSet;
import prism.PrismException;
import prism.PrismFileLog;
import prism.PrismLog;
import prism.PrismUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;

public class STPGOptimisticIntervalIteration {
    STPGModelChecker modelChecker;

    int lowerJumps = 0;
    int upperJumps = 0;
    int penaltyForUnsuccessfulVerfication = 1000;

    boolean debug = true;

    double precision = 1e-12;
    boolean absolute = true;

    public STPGOptimisticIntervalIteration(STPGModelChecker modelChecker) {
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




        // For guaranteed VI, we need MECs and an EC computer to find SECs
        List<BitSet> mecs = null;
        explicit.ECComputerDefault ec =null;

        System.out.println("Getting MECs...");
        //compute MECs one time, use the decomposition in every iteration; SECs still have to be recomputed
        ec = (ECComputerDefault) ECComputer.createECComputer(modelChecker, stpg);
        //need a copy of unknown, since EC computation empties the set as a side effect
        BitSet unknownForEC = new BitSet();
        unknownForEC.or(unknown);
        ec.computeMECStates(unknownForEC);
        mecs = ec.getMECStates();
        System.out.println("Number of MECs: " + mecs.size());

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
                    double[][] subres = iterateOnSubset((STPGExplicit) stpg, min1, min2, upperBounds, upperBounds2, lowerBounds2, lowerBounds, genAdv, adv, itersInSCC, variant, mecs, ec, statesForSCC, statesForSCCIntSet, -1);
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
            double[][] subres = iterateOnSubset((STPGExplicit) stpg, min1, min2, upperBounds, upperBounds2, lowerBounds2, lowerBounds, genAdv, adv, iters, variant, mecs, ec, unknown, null, initialState);
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
                                       boolean genAdv, int[] adv, int iters, ProbModelChecker.SolnMethod variant, List<BitSet> mecs, explicit.ECComputerDefault ec,
                                       BitSet subset, IntSet subsetAsIntSet, int initialState) throws PrismException{
        iters = 0;
        boolean done = false;
        double tmpsoln[];

        // Helper variables needed for OVI
        double epsPrimeLower = this.modelChecker.termCritParam; //precision that OVI gives the normal VI phase. Will decrease over time.
        double epsPrimeUpper = this.modelChecker.termCritParam; //precision that OVI gives the normal VI phase. Will decrease over time.


        boolean OVI_L_in_verification_phase = ((this.modelChecker.solnMethodOptions&2) == 2); // Continue to work on L in verification phase?
        int verifIters = 0; //counts how many iterations we made in verification phases
        boolean verifPhase = false; //indicates whether we are in a verification phase right now
        List<BitSet> OVI_SECs = null;

        // OIVI verification holders

        List<BitSet> OVI_Lower_SECs = null;
        List<BitSet> OVI_Upper_SECs = null;

        boolean lowerVerificationPhase = false;
        boolean upperVerificationPhase = false;
        int counterUpperboundsDidNotImproveMuch = 0;

        int counterUntilAllowLowerVerification = 0;
        int counterUntilAllowUpperVerification = 0;

        boolean lowerVerificationTerminated = false;
        boolean upperVerificationTerminated = false;

        double[] lowerVerificationGuessBounds = new double[stpg.getNumStates()];
        double[] lowerVerificationGuessBoundsAfterIteration = new double[stpg.getNumStates()];

        double[] upperVerificationGuessBounds = new double[stpg.getNumStates()];
        double[] upperVerificationGuessBoundsAfterIteration = new double[stpg.getNumStates()];

        Arrays.fill(lowerVerificationGuessBounds, 1.0);
        Arrays.fill(lowerVerificationGuessBoundsAfterIteration, 1.0);

        for (int state = 0; state < stpg.getNumStates(); state++) {
            if (upperBounds[state] == 0) {
                lowerVerificationGuessBounds[state] = 0;
                lowerVerificationGuessBoundsAfterIteration[state] = 0;
            }
            if (lowerBounds[state] == 1) {
                upperVerificationGuessBounds[state] = 1;
                upperVerificationGuessBounds[state] = 1;
            }
        }

        int lowerVerifIters = 0;
        int upperVerifIters = 0;

        // Gauss-Seidel
        boolean doGS = ((modelChecker.solnMethodOptions%2)==1);
        if (doGS){System.out.println("Doing Gauss Seidel variant");}

        while (!done) {
            iters++;


            // Enable
            counterUntilAllowLowerVerification--;
            counterUntilAllowUpperVerification--;
            if (lowerVerificationTerminated && counterUntilAllowLowerVerification < 0) {
                lowerVerificationTerminated = false;
            }
            if (upperVerificationTerminated && counterUntilAllowUpperVerification < 0) {
                upperVerificationTerminated = false;
            }


            if (doGS) {
                // For Gauss Seidel, let the new objects have the most up-to-date values, then only work on the new ones
                lowerBoundsNew = lowerBounds.clone();
                upperBoundsNew = upperBounds != null ? upperBounds.clone() : null;
            }
            //Bellman Updates

            // Update lower and upper bounds for classical Interval Iteration
            // TODO Sasha: Since verification can create a better bound, we should stop normal iterations of the bounds while in verification
            if (!lowerVerificationPhase) {
                if (!doGS) {
                    //Non-Gauss Seidel - Save the result in distinct array
                    stpg.mvMultMinMax(lowerBounds, min1, min2, lowerBoundsNew, subset, false, genAdv ? adv : null);
                } else {
                    //Gauss Seidel - Do operations in-place
                    stpg.mvMultGSMinMax(lowerBoundsNew, min1, min2, subset, false, modelChecker.termCrit == ProbModelChecker.TermCrit.ABSOLUTE);
                }
            }
            else {
                if (!doGS) {
                    //Non-Gauss Seidel - Save the result in distinct array
                    stpg.mvMultMinMax(lowerVerificationGuessBounds, min1, min2, lowerVerificationGuessBoundsAfterIteration, subset, false, genAdv ? adv : null);
                } else {
                    //Gauss Seidel - Do operations in-place
                    stpg.mvMultGSMinMax(lowerVerificationGuessBoundsAfterIteration, min1, min2, subset, false, modelChecker.termCrit == ProbModelChecker.TermCrit.ABSOLUTE);
                }
            }

            if (!upperVerificationPhase) {
                    if (!doGS) {
                    //Non-Gauss Seidel - Save the result in distinct array
                    stpg.mvMultGSMinMax(upperBoundsNew, min1, min2, subset, false, modelChecker.termCrit == ProbModelChecker.TermCrit.ABSOLUTE);
                } else {
                    //Gauss Seidel - Do operations in-place
                    stpg.mvMultGSMinMax(upperBoundsNew, min1, min2, subset, false, modelChecker.termCrit == ProbModelChecker.TermCrit.ABSOLUTE);
                }
            }
            else {
                if (!doGS) {
                    //Non-Gauss Seidel - Save the result in distinct array
                    stpg.mvMultMinMax(upperVerificationGuessBounds, min1, min2, upperVerificationGuessBoundsAfterIteration, subset, false, genAdv ? adv : null);
                } else {
                    //Gauss Seidel - Do operations in-place
                    stpg.mvMultGSMinMax(upperVerificationGuessBoundsAfterIteration, min1, min2, subset, false, modelChecker.termCrit == ProbModelChecker.TermCrit.ABSOLUTE);
                }
            }

            /**
             * DEFLATING
             * for II as in KKKW18. For SVI a bit special. For OVI, use the fixed SECs we found before.
             */

            // Compute deflated upperbounds as in II (KKKW18)
            for (BitSet mec : mecs) {
                if (subset.intersects(mec)) {
                    upperBoundsNew = deflate(stpg, min1, min2, lowerBoundsNew, upperBoundsNew, mec, ec)[0];
                }
            }

            //For OVI: If in verification phase, deflate using the precomputed set of SECs
            //We *must* deflate every iteration, because otherwise we might have some mixed thing (smaller everywhere but in some best exit), and without deflating we conclude inductive lower bound. With deflating, we realize that the best exit became smaller and it is not inductive lower bound.
            //For now, always deflate
            if (lowerVerificationPhase) {
                for (BitSet sec : OVI_Lower_SECs) {
                    lowerVerificationGuessBoundsAfterIteration = deflate(stpg, min1, min2, lowerBoundsNew, lowerVerificationGuessBoundsAfterIteration, sec, ec)[0];
                }
            }
            if (false && upperVerificationPhase) { // Lower bounds don't need deflating and we induce a lower bound from the upper bound
                for (BitSet sec : OVI_Upper_SECs) {
                    upperVerificationGuessBoundsAfterIteration = deflate(stpg, min1, min2, upperBoundsNew, upperVerificationGuessBoundsAfterIteration, sec, ec)[0];
                }
            }

            // Swap vectors for next iter
            // Now lowerBounds is the most up-to-date approximation, while the lowerBoundsNew contains the previous iteration
            // If lowerBounds were not changed no need to change anything
            if (!lowerVerificationPhase) {
                tmpsoln = lowerBounds;
                lowerBounds = lowerBoundsNew;
                lowerBoundsNew = tmpsoln;
            }

            if (!upperVerificationPhase) {
                tmpsoln = upperBounds;
                upperBounds = upperBoundsNew;
                upperBoundsNew = tmpsoln;
            }

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


            // Optimistic Iteration Criterion

            // ------- VERIFICATION FROM THE LOWER BOUND -------
            if (!lowerVerificationPhase) {
                boolean lowerBoundEpsClose = PrismUtils.doublesAreClose(lowerBounds, lowerBoundsNew, epsPrimeLower,
                        modelChecker.termCrit == ProbModelChecker.TermCrit.ABSOLUTE);

                if (lowerBoundEpsClose) {
                    System.out.println("Starting Lower Verification Phase in Iteration " + iters);

                    // Guess the bounds that are above lowerBounds and hopefully above the solution
                    lowerVerificationPhase = true;
                    lowerVerificationGuessBounds = diffPlus(lowerBounds, lowerVerificationGuessBounds, subset);

                    //Also precompute the SECs according to the current lowerBound, which will then be used for deflating
                    OVI_Lower_SECs = new ArrayList<BitSet>();
                    for (BitSet mec : mecs) {
                        if (subset.intersects(mec)) {
                            List<BitSet> SECs = ec.getSECs(mec, lowerBounds, min1, min2);
                            OVI_Lower_SECs.addAll(SECs);
                        }
                    }
                }
            }
            else if (lowerVerificationPhase) {
                // Try to find whether we
                boolean allUp = true;
                boolean allDown = true;
                boolean atLeastOneImprovement = false;

                for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
                    boolean wentDown = lowerVerificationGuessBoundsAfterIteration[s] < lowerVerificationGuessBounds[s];
                    boolean wentUp = lowerVerificationGuessBoundsAfterIteration[s] > lowerVerificationGuessBounds[s];
                    boolean improvedAtLeastSomewhere = allUp && lowerVerificationGuessBoundsAfterIteration[s] > lowerBounds[s];

                    // We should only consider counterexamples that are not due to some fiddly floating-point arithmetics
                    //boolean farEnoughApart = !PrismUtils.doublesAreClose(lowerVerificationGuessBoundsAfterIteration[s], lowerVerificationGuessBounds[s],
                    //        precision, absolute);
                    //wentUp = wentUp && farEnoughApart; // Only count something as an improvment if it REALLY is
                    //wentDown = wentDown && farEnoughApart;



                    if (wentDown) {
                        if (allUp) {
                            log("[L-Verification] Not all go up. In state "+s+" value before iter is: "+lowerVerificationGuessBounds[s]+", after iteration: "+lowerVerificationGuessBoundsAfterIteration[s]);
                        }
                        allUp = false;
                    }
                    if (wentUp) {
                        if (allDown) {
                            log("[L-Verification] Not all go down. In state "+s+" value before iter is: "+lowerVerificationGuessBounds[s]+", after iteration: "+lowerVerificationGuessBoundsAfterIteration[s]);
                        }
                        allDown = false;
                    }
                    // To be a better L, the new bound needs to be in between current L and U
                    if (improvedAtLeastSomewhere && lowerVerificationGuessBoundsAfterIteration[s] <= upperBounds[s]) {
                        //System.out.println("In State "+s+" improvement from "+upperBounds[s]+" to "+upperVerificationGuessBoundsAfterIteration[s]);
                        atLeastOneImprovement = true;
                    }
                }

                if (allDown) {
                    // bound is inductive counter-bound, everything stayed or went down
                    done = true;
                    System.out.println("The induced bound from L was an upper bound so we are finished. Happened in iteration " + iters + " after " + verifIters + " iterations in the verification phase.");
                    // Set result
                    log("Lower: "+lowerBounds[initialState]);
                    log("Guess: "+lowerVerificationGuessBounds[initialState]);
                    log("G+Ite: "+lowerVerificationGuessBoundsAfterIteration[initialState]);
                    for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
                        upperBounds[s] = (lowerBounds[s] + lowerVerificationGuessBoundsAfterIteration[s])/2;
                    }
                } else if (allUp && atLeastOneImprovement) {
                    // bound is improvment to current bound. Go on with verification, but use this bound
                    lowerBounds = lowerVerificationGuessBoundsAfterIteration;
                    lowerJumps++;
                    System.out.println("[Iteration "+iters+"]: The induced bound from L was a lower bound, so we use it now");
                }
                else {
                    System.out.println("Lower bound induction happened in iteration " + iters + " but was unsuccessful");
                    lowerVerificationPhase = false;
                    lowerVerificationTerminated = true;
                    counterUntilAllowLowerVerification = penaltyForUnsuccessfulVerfication;
                    epsPrimeLower/=2.0;
                }
            }

            // ------- VERIFICATION FROM THE UPPER BOUND -------
            if (!done && !upperVerificationPhase) {
                boolean upperBoundEpsClose = PrismUtils.doublesAreClose(upperBounds, upperBoundsNew, epsPrimeUpper, modelChecker.termCrit == ProbModelChecker.TermCrit.ABSOLUTE);

                if (upperBoundEpsClose) {
                    counterUpperboundsDidNotImproveMuch++;
                }
                else {
                    // Reset counter
                    counterUpperboundsDidNotImproveMuch = 0;
                }

                // Only start upper verification phase if there wasn't any big change to upper bound for n times
                // This ensures that it's not because the better value hasn't reached any state so far
                // furthermore, waiting n iterations too long is probably not too horrible
                if (upperBoundEpsClose && counterUpperboundsDidNotImproveMuch >= stpg.getNumStates()) {
                    System.out.println("Starting Upper Verification Phase in Iteration " + iters);

                    // Guess the bounds that are above upperBounds and hopefully above the solution
                    upperVerificationPhase = true;
                    upperVerificationGuessBounds = diffMinus(upperBounds, upperVerificationGuessBounds, subset);
                    log("Upperbounds: "+Arrays.toString(upperBounds));
                    log("UpperGuess : "+Arrays.toString(upperVerificationGuessBounds));
                }
            }
            else if (!done && upperVerificationPhase) {
                // Try to find whether we
                boolean allUp = true;
                boolean allDown = true;
                boolean atLeastOneImprovement = false;
                boolean betterThanLowerBound = false;

                for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
                    boolean wentDown = upperVerificationGuessBoundsAfterIteration[s] < upperVerificationGuessBounds[s];
                    boolean wentUp = upperVerificationGuessBoundsAfterIteration[s] > upperVerificationGuessBounds[s];
                    boolean improvedAtLeastSomewhere = allDown && upperVerificationGuessBoundsAfterIteration[s] < upperBounds[s];

                    // We should only consider counterexamples that are not due to some fiddly floating-point arithmetics
                    boolean farEnoughApart = !PrismUtils.doublesAreClose(upperVerificationGuessBoundsAfterIteration[s], upperVerificationGuessBounds[s],
                            precision, absolute);
                    wentUp = wentUp && farEnoughApart;
                    wentDown = wentDown && farEnoughApart;

                    if (wentDown) {
                        if (allUp) {
                            log("[U-Verification] Not all go up. In state "+s+" value before iter is: "+upperVerificationGuessBounds[s]+", after iteration: "+upperVerificationGuessBoundsAfterIteration[s]);
                        }
                        allUp = false;
                    }
                    if (wentUp) {
                        if (allDown) {
                            log("[U-Verification] Not all go down. In state "+s+" value before iter is: "+upperVerificationGuessBounds[s]+", after iteration: "+upperVerificationGuessBoundsAfterIteration[s]);
                        }
                        allDown = false;
                    }
                    // To be a better U, the new bound needs to be in between current L and U
                    if (improvedAtLeastSomewhere) {
                        //System.out.println("In State "+s+" improvement from "+upperBounds[s]+" to "+upperVerificationGuessBoundsAfterIteration[s]);
                        atLeastOneImprovement = true;
                        if (upperVerificationGuessBoundsAfterIteration[s] >= lowerBounds[s]) {
                            betterThanLowerBound = true;
                        }
                        else {
                            log("Guessed: "+upperVerificationGuessBoundsAfterIteration[s]+", Lower: "+lowerBounds[s]+", Upper: "+upperBounds[s]);
                        }
                    }
                }

                if (allDown && atLeastOneImprovement && betterThanLowerBound) {
                    // bound is improvment to current bound. Go on with verification, but use this bound
                    upperBounds = upperVerificationGuessBoundsAfterIteration;
                    System.out.println("[Iteration "+iters+"]: The induced bound from U was a upper bound, so we use it now");
                    upperJumps++;
                } else if (allUp) {
                    // bound is inductive counter-bound, everything stayed or went down
                    done = true;
                    System.out.println("The induced bound from U was a lower bound so we are finished. Happened in iteration " + iters + " after " + verifIters + " iterations in the verification phase.");
                    // Set result
                    for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
                        lowerBounds[s] = (upperBounds[s] + upperVerificationGuessBoundsAfterIteration[s])/2;
                    }
                } else {
                    System.out.println("Upper bound induction happened in iteration " + iters + " but was unsuccessful.");
                    log("Was inductive Bound improving everywhere? "+allDown);
                    log("Was inductive Bound lower than current upper bound in at least one state? "+atLeastOneImprovement);
                    log("Was inductive Bound higher than current lower bound? "+betterThanLowerBound);
                    upperVerificationPhase = false;
                    upperVerificationTerminated = true;
                    counterUntilAllowUpperVerification = penaltyForUnsuccessfulVerfication;
                    epsPrimeUpper/=2.0;
                }
            }

            /*
            System.out.println("Iteration: "+iters);
            System.out.println("LOWER: \t\t"+Arrays.toString(lowerBounds));
            System.out.println("LOWER OLD: \t"+Arrays.toString(lowerBoundsNew));

            System.out.println("UPPER: \t\t"+Arrays.toString(upperBounds));
            System.out.println("UPPER_OLD: \t"+Arrays.toString(upperBoundsNew));
            */
        }



        if(initialState!=-1)
            System.out.println("Resulting interval: ["+lowerBounds[initialState]+","+((upperBounds!=null) ? upperBounds[initialState] : "none")+"]");

        System.out.println("Lowerjumps: "+lowerJumps+", upperJumps: "+upperJumps);
        System.out.println("Eps' Lower: "+epsPrimeLower);
        System.out.println("Eps' Upper: "+epsPrimeUpper);

        return new double[][]{lowerBounds,upperBounds,{iters}};
    }

    /**
     * KKKW18 deflate operation.
     */
    private double[][] deflate(STPGExplicit stpg, boolean min1, boolean min2, double[] lowerBounds, double[] upperBounds, BitSet mec, explicit.ECComputerDefault ec) throws PrismException {

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

    /**
     * Method for guessing the upper bound for OVI from L. See definition of diff^+ in original OVI paper or the new one.
     * @param lowerBounds The lower bound which is the basis of the guess
     * @param upperBoundsGuess The upper bounds which is filled on the considered subset and kept unchanged outside the subset
     * @param subset The subset that we are currently working on (S? or the current SCC)
     * @return upperBounds, since I don't want to rely on side effects anymore. Sometimes they don't work.
     * On subset, UpperBounds are now set to sth that is epsilon (termcritparam) greater than vector, relatively or absolutely.
     * Corner case: 0 stays 0, nothing greater than 1
     */
    private double[] diffPlus(double[] lowerBounds, double[] upperBoundsGuess, BitSet subset) {
        upperBoundsGuess = lowerBounds.clone();
        for(int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s+1)){
            if(modelChecker.termCrit == ProbModelChecker.TermCrit.ABSOLUTE){
                upperBoundsGuess[s] = Math.min(1,lowerBounds[s] == 0 ? 0 : lowerBounds[s] + this.modelChecker.termCritParam);
            }
            else{ //if(modelChecker.termCrit == TermCrit.RELATIVE){
                upperBoundsGuess[s] = Math.min(1,lowerBounds[s] * (1+modelChecker.termCritParam));
            }
        }
        return upperBoundsGuess;
    }

    /**
     * Method for guessing the lower bound for OVI from U. See definition of diff^+ in original OVI paper or the new one.
     * @param upperBounds The upper bound which is the basis of the guess
     * @param lowerBoundsGuess The lower bounds which is filled on the considered subset and kept unchanged outside the subset
     * @param subset The subset that we are currently working on (S? or the current SCC)
     * @return lowerBounds, since I don't want to rely on side effects anymore. Sometimes they don't work.
     * On subset, UpperBounds are now set to sth that is epsilon (termcritparam) greater than vector, relatively or absolutely.
     * Corner case: 0 stays 0, nothing greater than 1
     */
    private double[] diffMinus(double[] upperBounds, double[] lowerBoundsGuess, BitSet subset){
        lowerBoundsGuess = upperBounds.clone(); //Necessary so that states out of the subset get correct Values for Actions
        for(int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s+1)){
            if(modelChecker.termCrit == ProbModelChecker.TermCrit.ABSOLUTE){
                lowerBoundsGuess[s] = Math.max(0,upperBounds[s] == 1 ? 1 : upperBounds[s] - this.modelChecker.termCritParam);
            }
            else{ //if(modelChecker.termCrit == TermCrit.RELATIVE){
                lowerBoundsGuess[s] = Math.max(0,upperBounds[s] * (1-modelChecker.termCritParam));
            }
        }
        return lowerBoundsGuess;
    }

    private void log(String s) {
        if (debug) {
            System.out.println(s);
        }
    }
}
