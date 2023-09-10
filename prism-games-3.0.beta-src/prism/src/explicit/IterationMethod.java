//==============================================================================
//	
//	Copyright (c) 2016-
//	Authors:
//	* Joachim Klein <klein@tcs.inf.tu-dresden.de> (TU Dresden)
//	
//------------------------------------------------------------------------------
//	
//	This file is part of PRISM.
//	
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//	
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//	
//==============================================================================


package explicit;

import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.PrimitiveIterator;

import common.IntSet;
import common.PeriodicTimer;
import explicit.rewards.MCRewards;
import explicit.rewards.MDPRewards;
import java.util.Set;
import prism.OptionsIntervalIteration;
import prism.PrismException;
import prism.PrismUtils;

/**
 * Abstract class that encapsulates the functionality for the different iteration methods
 * (e.g., Power, Jacobi, Gauss-Seidel, ...).
 * <p>
 * Provides methods as well to do the actual work in a (topological) value or interval iteration.
 */
public abstract class IterationMethod {

	/**
	 * Interface for an object that provides the basic steps for a value iteration.
	 */
	public interface IterationValIter {
		/** Initialise the value iteration with the given solution vector */
		public void init(double[] soln);
		/** Get the current solution vector */
		public double[] getSolnVector();

		/** Perform one iteration (over the set of states) and return true if convergence has been detected. */
		public boolean iterateAndCheckConvergence(IntSet states) throws PrismException;

		/**
		 * Notify that the given states are done (e.g., because the given SCC is finished
		 * during a topological iteration).
		 * <br>
		 * This allows the two-vector iteration methods to store the current values for these
		 * states into the second vector, so that switching the vector does not return to
		 * previous values.
		 */
		public void doneWith(IntSet states);

		/**
		 * Solve for a given singleton SCC consisting of {@code state} using {@code solver},
		 * store the result in the solution vector(s).
		 */
		public void solveSingletonSCC(int state, SingletonSCCSolver solver);

		/** Return the underlying model */
		public Model getModel();
	}

	/**
	 * Interface for an object that provides the atomic steps for a value iteration
	 * in the context of interval iteration, i.e., for the interval iteration
	 * there are two IterationIntervalIter objects, one from below or from above.
	 */
	public interface IterationIntervalIter {
		/** Initialise the value iteration with the given solution vector */
		public void init(double[] soln);
		/** Get the current solution vector */
		public double[] getSolnVector();

		/** Perform one iteration (over the set of states) */
		public void iterate(IntSet states) throws PrismException;

		/**
		 * Notify that the given states are done (e.g., because the given SCC is finished
		 * during a topological iteration).
		 * <br>
		 * This allows the two-vector iteration methods to store the current values for these
		 * states into the second vector, so that switching the vector does not return to
		 * previous values.
		 */
		public void doneWith(IntSet states);

		/**
		 * Solve for a given singleton SCC consisting of {@code state} using {@code solver},
		 * store the result in the solution vector(s).
		 */
		public void solveSingletonSCC(int s, SingletonSCCSolver solver);

		/** Return the underlying model */
		public Model getModel();
	}

	/** Storage for a single solution vector */
	public class IterationBasic {
		protected final Model model;
		protected double[] soln;

		public IterationBasic(Model model)
		{
			this.model = model;
		}

		public void init(double[] soln)
		{
			this.soln = soln;
		}

		public double[] getSolnVector()
		{
			return soln;
		}

		/* see IterationValIter.solveSingletonSCC() */
		public void solveSingletonSCC(int state, SingletonSCCSolver solver)
		{
			solver.solveFor(state, soln);
		}

		/* see IterationValIter.doneWith() */
		public void doneWith(IntSet states)
		{
			// single vector, nothing to do
		}

		public Model getModel()
		{
			return model;
		}
	}

	/** Abstract base class for an IterationValIter with a single solution vector */
	protected abstract class SingleVectorIterationValIter extends IterationBasic implements IterationValIter
	{
		public SingleVectorIterationValIter(Model model)
		{
			super(model);
		}
	}

	/** Abstract base class for an IterationIntervalIter with a single solution vector */
	protected abstract class SingleVectorIterationIntervalIter extends IterationBasic implements IterationIntervalIter
	{
		public SingleVectorIterationIntervalIter(Model model)
		{
			super(model);
		}
	}

	/**
	 * Functional interface for a post-processing step after an iteration that involves
	 * a pair of solution vectors.
	 * <br>
	 * This method may modify solnNew.
	 *
	 * @param solnOld the previous solution vector
	 * @param solnNew the new solution vector
	 * @param states the set of states that are the focus of the current iteration
	 */
	@FunctionalInterface
	interface IterationPostProcessor {
		void apply(double[] solnOld, double[] solnNew, IntSet states) throws PrismException;
	}

	/**
	 * Abstract base class for an IterationValIter / IterationIntervalIter that
	 * requires two solution vectors.
	 * Optionally, a post processing step is performed after each iteration.
	 */
	protected abstract class TwoVectorIteration extends IterationBasic implements IterationValIter, IterationIntervalIter {
		/** The solution vector that serves as the target vector in the iteration step */
		protected double[] soln2;
		/** Post processing, may be null */
		protected final IterationPostProcessor postProcessor;

		/** Constructor */
		protected TwoVectorIteration(Model model, IterationMethod.IterationPostProcessor postProcessor)
		{
			super(model);
			this.postProcessor = postProcessor;
		}

		@Override
		public void init(double[] soln)
		{
			super.init(soln);

			// create and initialise the second solution vector
			soln2 = new double[soln.length];
			System.arraycopy(soln, 0, soln2, 0, soln.length);
		}

		/** Perform one iteration */
		public abstract void doIterate(IntSet states) throws PrismException;

		@Override
		public void iterate(IntSet states) throws PrismException
		{
			// do the iteration
			doIterate(states);
			// optionally, post processing
			if (postProcessor != null) {
				postProcessor.apply(soln, soln2, states);
			}

			// switch vectors
			double[] tmp = soln;
			soln = soln2;
			soln2 = tmp;
		}

		@Override
		public boolean iterateAndCheckConvergence(IntSet states) throws PrismException
		{
			// do the iteration
			doIterate(states);
			// optionally, post processing
			if (postProcessor != null) {
				postProcessor.apply(soln, soln2, states);
			}
			// check convergence (on the set of states)
			boolean done = PrismUtils.doublesAreClose(soln, soln2, states.iterator(), termCritParam, absolute);

			// switch vectors
			double[] tmp = soln;
			soln = soln2;
			soln2 = tmp;

			return done;
		}

		@Override
		public void doneWith(IntSet states)
		{
			// we copy the values for the given states to the
			// second vector, so that switching between vectors
			// does not change their values
			PrimitiveIterator.OfInt it = states.iterator();
			while (it.hasNext()) {
				int state = it.nextInt();
				soln2[state] = soln[state];
			}
		}

		@Override
		public void solveSingletonSCC(int state, SingletonSCCSolver solver)
		{
			// solve and store result in soln vector
			super.solveSingletonSCC(state, solver);
			// copy result to soln2 vector as well
			soln2[state] = soln[state];
		}

	}

	/**
	 * Functional interface for a method that allows to
	 * determine the value for a singleton SCC in the model,
	 * given that all successor value have already been computed.
	 */
	@FunctionalInterface
	public interface SingletonSCCSolver {
		/**
		 * Compute the value for state {@code state}, under the assumption
		 * that it constitutes a (trivial or non-trivial) singleton SCC
		 * and that all successor values have already been computed in {@code soln}.
		 * Stores the result in {@code soln[state]}.
		 */
		public void solveFor(int state, double[] soln);
	}

	/** Convergence check: absolute or relative? */
	protected final boolean absolute;
	/** Convergence check: epsilon value */
	protected final double termCritParam;

	/**
	 * Constructor.
	 * @param absolute For convergence check, perform absolute comparison?
	 * @param termCritParam For convergence check, the epsilon value to use
	 */
	protected IterationMethod(boolean absolute, double termCritParam)
	{
		this.absolute = absolute;
		this.termCritParam = termCritParam;
	}

	// ------------ Abstract DTMC methods ----------------------------

	/** Obtain an Iteration object using mvMult (matrix-vector multiplication) in a DTMC */
	public abstract IterationValIter forMvMult(DTMC dtmc) throws PrismException;

	/**
	 * Obtain an Iteration object (for interval iteration) using mvMult
	 * (matrix-vector multiplication) in a DTMC.
	 * @param fromBelow for interval iteration from below?
	 * @param enforceMonotonic enforce element-wise monotonicity of the solution vector
	 * @param checkMonotonic check the element-wise monotonicity of the solution vector, throw exception if violated
	 */
	public abstract IterationIntervalIter forMvMultInterval(DTMC dtmc, boolean fromBelow, boolean enforceMonotonicity, boolean checkMonotonicity) throws PrismException;

	/** Obtain an Iteration object using mvMultRew (matrix-vector multiplication with rewards) in a DTMC */
	public abstract IterationValIter forMvMultRew(DTMC dtmc, MCRewards rew) throws PrismException;

	/**
	 * Obtain an Iteration object (for interval iteration) using mvMultRew
	 * (matrix-vector multiplication with rewards) in a DTMC.
	 * @param fromBelow for interval iteration from below?
	 * @param enforceMonotonic enforce element-wise monotonicity of the solution vector
	 * @param checkMonotonic check the element-wise monotonicity of the solution vector, throw exception if violated
	 */
	public abstract IterationIntervalIter forMvMultRewInterval(DTMC dtmc, MCRewards rew, boolean fromBelow, boolean enforceMonotonicity, boolean checkMonotonicity) throws PrismException;

	// ------------ Abstract MDP methods ----------------------------

	/**
	 * Obtain an Iteration object using mvMultMinMax (matrix-vector multiplication, followed by min/max)
	 * in an MDP.
	 * @param mdp the MDP
	 * @param min do min?
	 * @param strat optional, storage for strategy, ignored if null
	 */
	public abstract IterationValIter forMvMultMinMax(MDP mdp, boolean min, int[] strat) throws PrismException;

	/**
	 * Obtain an Iteration object using mvMultMinMax (matrix-vector multiplication, followed by min/max)
	 * in an MDP, for interval iteration.
	 * @param mdp the MDP
	 * @param min do min?
	 * @param strat optional, storage for strategy, ignored if null
	 * @param fromBelow for interval iteration from below?
	 * @param enforceMonotonic enforce element-wise monotonicity of the solution vector
	 * @param checkMonotonic check the element-wise monotonicity of the solution vector, throw exception if violated
	 */
	public abstract IterationIntervalIter forMvMultMinMaxInterval(MDP mdp, boolean min, int[] strat, boolean fromBelow, boolean enforceMonotonicity, boolean checkMonotonicity) throws PrismException;

	/**
	 * Obtain an Iteration object using mvMultRewMinMax (matrix-vector multiplication with rewards, followed by min/max)
	 * in an MDP.
	 * @param mdp the MDP
	 * @param rewards the reward structure
	 * @param min do min?
	 * @param strat optional, storage for strategy, ignored if null
	 */
	public abstract IterationValIter forMvMultRewMinMax(MDP mdp, MDPRewards rewards, boolean min, int[] strat) throws PrismException;

	/**
	 * Obtain an Iteration object using mvMultRewMinMax (matrix-vector multiplication with rewards, followed by min/max)
	 * in an MDP, for interval iteration.
	 * @param mdp the MDP
	 * @param rewards the reward structure
	 * @param min do min?
	 * @param strat optional, storage for strategy, ignored if null
	 * @param fromBelow for interval iteration from below?
	 * @param enforceMonotonic enforce element-wise monotonicity of the solution vector
	 * @param checkMonotonic check the element-wise monotonicity of the solution vector, throw exception if violated
	 */
	public abstract IterationIntervalIter forMvMultRewMinMaxInterval(MDP mdp, MDPRewards rewards, boolean min, int[] strat, boolean fromBelow, boolean enforceMonotonicity, boolean checkMonotonicity) throws PrismException;


	// ------------ Abstract generic methods ----------------------------

	/**
	 * Return a description of this iteration method for display.
	 */
	public abstract String getDescriptionShort();


	// ------------ Value iteration implementations ----------------------------

	/**
	 * Perform the actual work of a value iteration, i.e., iterate until convergence or abort.
	 * @param mc ProbModelChecker (for log and settings)
	 * @param description (for logging)
	 * @param iteration The iteration object
	 * @param unknownStates The set of unknown states, i.e., whose value should be determined
	 * @param startTime The start time (for logging purposes, obtained from a call to System.currentTimeMillis())
	 * @param iterationsExport an ExportIterations object (optional, ignored if null)
	 * @return a ModelChecker result with the solution vector and statistics
	 * @throws PrismException on non-convergence (if mc.errorOnNonConverge is set)
	 */
	public ModelCheckerResult doValueIteration(ProbModelChecker mc, String description, IterationValIter iteration, IntSet unknownStates, long startTime, ExportIterations iterationsExport) throws PrismException
	{
		int iters = 0;
		final int maxIters = mc.maxIters;
		boolean done = false;

		PeriodicTimer updatesTimer = new PeriodicTimer(ProbModelChecker.UPDATE_DELAY);
		updatesTimer.start();

		while (!done && iters < maxIters) {
			iters++;
			// do iteration step
			done = iteration.iterateAndCheckConvergence(unknownStates);

			if (iterationsExport != null)
				iterationsExport.exportVector(iteration.getSolnVector(), 0);

			if (!done && updatesTimer.triggered()) {
				mc.getLog().print("Iteration " + iters + ": ");
				mc.getLog().println(PrismUtils.formatDouble2dp(updatesTimer.elapsedMillisTotal() / 1000.0) + " sec so far");
			}
		}

		// Finished value iteration
		long mvCount = iters * countTransitions(iteration.getModel(), unknownStates);
		long timer = System.currentTimeMillis() - startTime;
//		mc.getLog().print("Value iteration (" + description + ")");
//		mc.getLog().print(" took " + iters + " iterations, ");
//		mc.getLog().print(mvCount + " multiplications");
//		mc.getLog().println(" and " + timer / 1000.0 + " seconds.");

		if (iterationsExport != null)
			iterationsExport.close();

		// Non-convergence is an error (usually)
		if (!done && mc.errorOnNonConverge) {
			String msg = "Iterative method did not converge within " + iters + " iterations.";
			msg += "\nConsider using a different numerical method or increasing the maximum number of iterations";
			throw new PrismException(msg);
		}

		// Return results
		ModelCheckerResult res = new ModelCheckerResult();
		res.soln = iteration.getSolnVector();
		res.numIters = iters;
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Perform the actual work of a topological value iteration, i.e., iterate until convergence or abort.
	 *
	 * @param mc ProbModelChecker (for log and settings)
	 * @param description (for logging)
	 * @param sccs The information about the SCCs and topological order
	 * @param iteration The iteration object
	 * @param singletonSCCSolver The solver for singleton SCCs
	 * @param startTime The start time (for logging purposes, obtained from a call to System.currentTimeMillis())
	 * @param iterationsExport an ExportIterations object (optional, ignored if null)
	 * @return a ModelChecker result with the solution vector and statistics
	 * @throws PrismException on non-convergence (if mc.errorOnNonConverge is set)
	 */
	public ModelCheckerResult doTopologicalValueIteration(ProbModelChecker mc, String description, SCCInfo sccs, IterationMethod.IterationValIter iterator, SingletonSCCSolver singletonSCCSolver, long startTime, ExportIterations iterationsExport) throws PrismException
	{
		// Start iterations
		int iters = 0;
		long mvCount = 0;
		final int maxIters = mc.maxIters;

		int numSCCs = sccs.getNumSCCs();
		int numNonSingletonSCCs = sccs.countNonSingletonSCCs();
		int finishedNonSingletonSCCs = 0;

		PeriodicTimer updatesTimer = new PeriodicTimer(ProbModelChecker.UPDATE_DELAY);
		updatesTimer.start();

		boolean done = true;
		for (int scc = 0; scc < numSCCs; scc++) {
			boolean doneSCC;

			if (sccs.isSingletonSCC(scc)) {
				// get the single state in this SCC
				int state = sccs.getStatesForSCC(scc).iterator().nextInt();
				iterator.solveSingletonSCC(state, singletonSCCSolver);

				// no need to call doneWith(...), as solveSingletonSCC updates
				// both vectors for two-iteration methods

				mvCount += countTransitions(iterator.getModel(), IntSet.asIntSet(state));

				iters++;
				if (iterationsExport != null)
					iterationsExport.exportVector(iterator.getSolnVector(), 0);

				doneSCC = true;
			} else {
				// complex SCC: do VI
				doneSCC = false;
				IntSet statesForSCC = sccs.getStatesForSCC(scc);
				int itersInSCC = 0;
				// abort on convergence or if iterations *in this SCC* are above maxIters
				while (!doneSCC && itersInSCC < maxIters) {
					iters++;
					itersInSCC++;
					// do iteration step
					doneSCC = iterator.iterateAndCheckConvergence(statesForSCC);

					if (iterationsExport != null)
						iterationsExport.exportVector(iterator.getSolnVector(), 0);

					if (!doneSCC && updatesTimer.triggered()) {
						mc.getLog().print("Iteration " + iters + ": ");
						mc.getLog().print("Iteration " + itersInSCC + " in SCC " + (finishedNonSingletonSCCs+1) + " of " + numNonSingletonSCCs);
						mc.getLog().println(", " + PrismUtils.formatDouble2dp(updatesTimer.elapsedMillisTotal() / 1000.0) + " sec so far");
					}

				}

				// notify the iterator that the states are done so that
				// their values can be copied to the second vector in a two-vector
				// iterator
				iterator.doneWith(statesForSCC);

				mvCount += itersInSCC * countTransitions(iterator.getModel(), statesForSCC);
			}

			if (!doneSCC) {
				done = false;
				break;
			}
		}

		// Finished value iteration
		long timer = System.currentTimeMillis() - startTime;
		mc.getLog().print("Value iteration (" + description + ", with " + numNonSingletonSCCs + " non-singleton SCCs)");
		mc.getLog().print(" took " + iters + " iterations, ");
		mc.getLog().print(mvCount + " multiplications");
		mc.getLog().println(" and " + timer / 1000.0 + " seconds.");

		if (iterationsExport != null)
			iterationsExport.close();

		// Non-convergence is an error (usually)
		if (!done && mc.errorOnNonConverge) {
			String msg = "Iterative method did not converge within " + iters + " iterations.";
			msg += "\nConsider using a different numerical method or increasing the maximum number of iterations";
			throw new PrismException(msg);
		}

		// Return results
		ModelCheckerResult res = new ModelCheckerResult();
		res.soln = iterator.getSolnVector();
		res.numIters = iters;
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Perform the actual work of an interval iteration, i.e., iterate until convergence or abort.
	 *
	 * @param mc ProbModelChecker (for log and settings)
	 * @param description Description (for logging)
	 * @param below The iteration object for the iteration from below
	 * @param above The iteration object for the iteration from above
	 * @param unknownStates The set of unknown states, i.e., whose value should be determined
	 * @param startTime The start time (for logging purposes, obtained from a call to System.currentTimeMillis())
	 * @param iterationsExport an ExportIterations object (optional, ignored if null)
	 * @return a ModelChecker result with the solution vector and statistics
	 * @throws PrismException on non-convergence (if mc.errorOnNonConverge is set)
	 */
	public ModelCheckerResult doIntervalIteration(ProbModelChecker mc, String description, IterationIntervalIter below, IterationIntervalIter above, IntSet unknownStates, long timer, ExportIterations iterationsExport) throws PrismException {
		try {
			// Start iterations
			int iters = 0;
			final int maxIters = mc.maxIters;
			boolean done = false;

			PeriodicTimer updatesTimer = new PeriodicTimer(ProbModelChecker.UPDATE_DELAY);
			updatesTimer.start();

			while (!done && iters < maxIters) {
				iters++;
				// Matrix-vector multiply
				below.iterate(unknownStates);
				above.iterate(unknownStates);

				if (iterationsExport != null) {
					iterationsExport.exportVector(below.getSolnVector(), 0);
					iterationsExport.exportVector(above.getSolnVector(), 1);
				}

				intervalIterationCheckForProblems(below.getSolnVector(), above.getSolnVector(), unknownStates.iterator());

				// Check termination
				done = PrismUtils.doublesAreClose(below.getSolnVector(), above.getSolnVector(), termCritParam, absolute);

				if (done) {
					double diff = PrismUtils.measureSupNormInterval(below.getSolnVector(), above.getSolnVector(), absolute);
//					mc.getLog().println("Max " + (!absolute ? "relative ": "") +
//							"diff between upper and lower bound on convergence: " + PrismUtils.formatDouble(diff));
					done = true;
				}

				if (!done && updatesTimer.triggered()) {
					double diff = PrismUtils.measureSupNormInterval(below.getSolnVector(), above.getSolnVector(), absolute);
//					mc.getLog().print("Iteration " + iters + ": ");
//					mc.getLog().print("max " + (absolute ? "" : "relative ") + "diff=" + PrismUtils.formatDouble(diff));
//					mc.getLog().println(", " + PrismUtils.formatDouble2dp(updatesTimer.elapsedMillisTotal() / 1000.0) + " sec so far");
				}
			}

			// Finished value iteration
			long mvCount = 2 * iters * countTransitions(below.getModel(), unknownStates);
			timer = System.currentTimeMillis() - timer;
//			mc.getLog().print("Interval iteration (" + description + ")");
//			mc.getLog().print(" took " + iters + " iterations, ");
//			mc.getLog().print(mvCount + " multiplications");
//			mc.getLog().println(" and " + timer / 1000.0 + " seconds.");

			if (done && OptionsIntervalIteration.from(mc.getSettings()).isSelectMidpointForResult()) {
				PrismUtils.selectMidpoint(below.getSolnVector(), above.getSolnVector());

				if (iterationsExport != null) {
					// export midpoint
					iterationsExport.exportVector(below.getSolnVector(), 0);
					iterationsExport.exportVector(below.getSolnVector(), 1);
				}
			}

			// Non-convergence is an error (usually)
			if (!done && mc.errorOnNonConverge) {
				String msg = "Iterative method (interval iteration) did not converge within " + iters + " iterations.";
				msg += "\nConsider using a different numerical method or increasing the maximum number of iterations";
				throw new PrismException(msg);
			}

			// Return results
			ModelCheckerResult res = new ModelCheckerResult();
			res.soln = below.getSolnVector();
			res.numIters = iters;
			res.timeTaken = timer / 1000.0;
			System.out.println("Number of Iterations:" + iters);
			return res;
		} finally {
			if (iterationsExport != null)
				iterationsExport.close();
		}
	}

  /**
   * Perform the actual work of an sound value iteration, i.e., iterate until convergence or abort.
   *
   * @param mc ProbModelChecker (for log and settings)
   * @param description Description (for logging)
   * @param below The iteration object for the iteration from below
   * @param above The iteration object for the iteration from above
   * @param unknownStates The set of unknown states, i.e., whose value should be determined
   * @param startTime The start time (for logging purposes, obtained from a call to System.currentTimeMillis())
   * @param iterationsExport an ExportIterations object (optional, ignored if null)
   * @return a ModelChecker result with the solution vector and statistics
   * @throws PrismException on non-convergence (if mc.errorOnNonConverge is set)
   */
  public double[][] doSoundValueIteration(MDP mdp, boolean min, double[] stepBoundReach, double[] stepBoundReachNew, double[] stepBoundStay, double[] stepBoundStayNew,
      int iters, List<BitSet> mecs, ECComputerDefault ec, BitSet subset, IntSet subsetAsIntSet, int initialState, BitSet yes) throws  PrismException{
    iters = 0;
    boolean done = false;
    double tmpsoln[];

    // Helper variables needed for SVI
    double decisionValue, decisionValueNew, lowerBound, lowerBoundNew, upperBound, upperBoundNew;
    lowerBound = lowerBoundNew = 0;
    upperBound = upperBoundNew = 1;
    decisionValue = min? 1 : 0;

    while (!done) {
      iters++;
      System.out.println("ITERATION:" + iters);

      /**
       * BELLMAN UPDATES
       * (Special for SVI, others just normal Bellmann.)
       **/
      // For SVI, we need the special find action
      // Recall that lower bounds are stepBoundReach and upperBounds are stepBoundStay. So in fact, upperBounds are not upperBounds but sth completely different. We just use the name, because this code is for four different kinds of VI at once
      boolean update = true;

      for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
        // find if the player is min or max
        int a;
        a = findAction(mdp, s, stepBoundReach, stepBoundStay, min, lowerBound, upperBound);
        if(iters==240 && s == 185)
        {
          System.out.println("here");
        }
        decisionValueNew = computeDecisionValue(mdp, stepBoundReach, stepBoundStay, s, a, min);

        decisionValue = min? Math.min(decisionValue, decisionValueNew) : Math.max(decisionValue, decisionValueNew);
        double reachVal = 0.0;
        double stayVal = 0.0;

        for (int succ : mdp.getChoice(s, a).keySet()) {
          reachVal += mdp.getChoice(s, a).get(succ) * (stepBoundReach[succ]);
          stayVal += mdp.getChoice(s, a).get(succ) * (stepBoundStay[succ]);
        }
        stepBoundReachNew[s] = reachVal;
        stepBoundStayNew[s] = stayVal;

      }
      /** to update global bounds if all state less than 1 */
      boolean allStatesStayValLessThan1 = true;
      for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
        if (stepBoundStayNew[s] == 1) {
          allStatesStayValLessThan1 = false;
          break;
        }
      }

      /** Do smart SVI stuff: Compute upper and lower bound. This is used for checking termination.*/
      // When we terminate, the vectors lowerBounds and upperBounds are updated to contain the smarter values, i.e. best lower and upper bound SVI can give us right now
      if (allStatesStayValLessThan1) {
        //  if (allStatesStayValLessThan1 & update) {
        //double lower_val=Double.POSITIVE_INFINITY; //would be for rewards

        double lowerboundCandidate = 1.0;
        for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
          lowerboundCandidate = Math
                .min(lowerboundCandidate, stepBoundReachNew[s] / (1 - stepBoundStayNew[s]));
        }
        if(min){
          if(decisionValue < lowerboundCandidate) System.out.println("Need of DECISION VALUE for LB in iteration " + iters + ". DecVal: "+decisionValue + ", approx_lower: "+ lowerboundCandidate + ", oldlowerBound: "+ lowerBoundNew);
            lowerboundCandidate = Math.min(decisionValue, lowerboundCandidate);
        }
        lowerBound = Math.max(lowerBoundNew, lowerboundCandidate);
        lowerBoundNew = lowerBound; //remember this for next iteration

        System.out.println("lowerBound: "+ lowerBound);

        //double upper_val=Double.NEGATIVE_INFINITY;
        double upperBoundCandidate = 0.0;
        for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
          upperBoundCandidate = Math.max(upperBoundCandidate, stepBoundReachNew[s] / (1 - stepBoundStayNew[s]));
        }
        if(!min) {
          if (decisionValue > upperBoundCandidate) {
            System.out.println(
                "Need of DECISION VALUE for UBin iteration " + iters + ". DecVal: " + decisionValue
                    + ", approx_upper: " + upperBoundCandidate + ", oldupperBound: " + upperBoundNew);

            double[] underApproximationforVI  = new double[stepBoundReach.length];
            double[] overApproximationforVI  = new double[stepBoundReach.length];

            for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
              underApproximationforVI[s] = stepBoundReach[s] + stepBoundStay[s] * lowerBound;
              overApproximationforVI[s] = stepBoundReach[s] + stepBoundStay[s] * upperBound;
            }
            //TODO: call BVI here on current approximation
          }
          upperBoundCandidate = Math.max(decisionValue, upperBoundCandidate);
        }
        upperBound = Math.min(upperBoundNew, upperBoundCandidate);
        upperBoundNew = upperBound; //remember this for next iteration

        System.out.println("upperBound: "+upperBound);
      }

      // Swap vectors for next iter
      // Now lowerBounds is the most up-to-date approximation, while the lowerBoundsNew contains the previous iteration
      tmpsoln = stepBoundReach;
      stepBoundReach = stepBoundReachNew;
      stepBoundReachNew = tmpsoln;
      //System.out.println("Reach: " + Arrays.toString(stepBoundReach));

      tmpsoln = stepBoundStay;
      stepBoundStay = stepBoundStayNew;
      stepBoundStayNew = tmpsoln;
      //System.out.println("Stay: " + Arrays.toString(stepBoundStay));
      double [] overApproximation = new double[stepBoundReach.length];
      double [] underApproximation = new double[stepBoundReach.length];
      for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
        underApproximation[s] = stepBoundReach[s] + stepBoundStay[s] * lowerBound;
        overApproximation[s] = stepBoundReach[s] + stepBoundStay[s] * upperBound;
      }

      //System.out.println("Under-approximation: " + Arrays.toString(underApproximation));
      //System.out.println("Over-approximation: " + Arrays.toString(overApproximation));


      /**
       * Check termination
       */
      double relevantStayVal = 0;
      if (initialState != -1) {
        relevantStayVal = stepBoundStay[initialState];
      } else {
        for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
          relevantStayVal = Math.max(relevantStayVal, stepBoundStay[s]);
        }
      }
      /**modifed termination criterion for the termination testing**/
      done = (upperBound - lowerBound) * relevantStayVal < this.termCritParam;
      //done = relevantStayVal < 2 * this.termCritParam;
      if (done) {
//			print_smg(stpg,lowerBounds,upperBounds);
        //When we are done, we have to insert the smarter values
        if (initialState != -1) {
          //Only in initial state
          //System.out.println("ReachFin: " + Arrays.toString(stepBoundReach) + "\nStay: " + Arrays.toString(stepBoundStay) + "\nLB: " + lowerBound + "\nUB: " + upperBound);
          double lb_i = stepBoundReach[initialState] + stepBoundStay[initialState] * lowerBound;
          double ub_i = stepBoundReach[initialState] + stepBoundStay[initialState] * upperBound;
        } else {
          //in all states, for topological VI. Need the second thing as temp, since meaning of content switches from reach/stayVal to actual lower/upper bound
          for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
//            stepBoundReachNew[s] = stepBoundReach[s] + stepBoundStay[s] * lowerBound;
//            stepBoundStayNew[s] = stepBoundReach[s] + stepBoundStay[s] * upperBound;
            stepBoundReachNew[s] = stepBoundReach[s];
            stepBoundStayNew[s] = stepBoundStay[s];
          }
          stepBoundReach = stepBoundReachNew.clone();
          stepBoundStay = stepBoundStayNew.clone();
        }
      } else{
        //System.out.println("Reach: " + Arrays.toString(stepBoundReach) + "\nStay: " + Arrays.toString(stepBoundStay) + "\nLB: " + lowerBound + "\nUB: " + upperBound);
        for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
          double lbprint = stepBoundReach[s] + stepBoundStay[s] * lowerBound;
          double ubprint =stepBoundReach[s] + stepBoundStay[s] * upperBound;
          ///System.out.println("bounds for state "+ s + " are: [" + lbprint +", "+ ubprint +"]");
        }
      }
    }
    if(initialState!=-1){
      double lbprint = stepBoundReach[initialState] + stepBoundStay[initialState] * lowerBound;
      double ubprint =stepBoundReach[initialState] + stepBoundStay[initialState] * upperBound;
      System.out.println("Resulting interval: ["+lbprint+","+ ubprint + "]");
    }
    for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
      // storing the approximated value in stepBoundReachNew
      stepBoundReachNew[s] = stepBoundReach[s] + stepBoundStay[s] * (lowerBound+upperBound)/2;
    }

    return new double[][]{stepBoundReachNew, stepBoundReach,stepBoundStay,{iters},{lowerBound},{upperBound}};
  }



  public double[][] doSoundValueIterationUntildVal(MDP mdp, boolean min, double[] stepBoundReach, double[] stepBoundReachNew, double[] stepBoundStay, double[] stepBoundStayNew,
      int iters, List<BitSet> mecs, ECComputerDefault ec, BitSet subset, IntSet subsetAsIntSet, int initialState, BitSet yes) throws  PrismException{
    iters = 0;
    boolean done = false;
    double tmpsoln[];

    // Helper variables needed for SVI
    double decisionValue, decisionValueNew, lowerBound, lowerBoundNew, upperBound, upperBoundNew;
    lowerBound = lowerBoundNew = 0;
    upperBound = upperBoundNew = 1;
    decisionValue = min? 1 : 0;

    while (!done) {
      iters++;
      System.out.println("ITERATION:" + iters);

      /**
       * BELLMAN UPDATES
       * (Special for SVI, others just normal Bellmann.)
       **/
      // For SVI, we need the special find action
      // Recall that lower bounds are stepBoundReach and upperBounds are stepBoundStay. So in fact, upperBounds are not upperBounds but sth completely different. We just use the name, because this code is for four different kinds of VI at once
      boolean update = true;

      for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
        // find if the player is min or max
        int a;
        a = findAction(mdp, s, stepBoundReach, stepBoundStay, min, lowerBound, upperBound);
        if(iters==240 && s == 185)
        {
          System.out.println("here");
        }
        decisionValueNew = computeDecisionValue(mdp, stepBoundReach, stepBoundStay, s, a, min);

        decisionValue = min? Math.min(decisionValue, decisionValueNew) : Math.max(decisionValue, decisionValueNew);
        double reachVal = 0.0;
        double stayVal = 0.0;

        for (int succ : mdp.getChoice(s, a).keySet()) {
          reachVal += mdp.getChoice(s, a).get(succ) * (stepBoundReach[succ]);
          stayVal += mdp.getChoice(s, a).get(succ) * (stepBoundStay[succ]);
        }
        stepBoundReachNew[s] = reachVal;
        stepBoundStayNew[s] = stayVal;

      }
      /** to update global bounds if all state less than 1 */
      boolean allStatesStayValLessThan1 = true;
      for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
        if (stepBoundStayNew[s] == 1) {
          allStatesStayValLessThan1 = false;
          break;
        }
      }

      /** Do smart SVI stuff: Compute upper and lower bound. This is used for checking termination.*/
      // When we terminate, the vectors lowerBounds and upperBounds are updated to contain the smarter values, i.e. best lower and upper bound SVI can give us right now
      if (allStatesStayValLessThan1) {
        //  if (allStatesStayValLessThan1 & update) {
        //double lower_val=Double.POSITIVE_INFINITY; //would be for rewards

        double lowerboundCandidate = 1.0;
        for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
          lowerboundCandidate = Math
              .min(lowerboundCandidate, stepBoundReachNew[s] / (1 - stepBoundStayNew[s]));
        }
        if(min){
          if(decisionValue < lowerboundCandidate) {
            System.out.println(
                "Need of DECISION VALUE for LB in iteration " + iters + ". DecVal: " + decisionValue
                    + ", approx_lower: " + lowerboundCandidate + ", oldlowerBound: "
                    + lowerBoundNew);
            double[] underApproximationforVI  = new double[stepBoundReach.length];
            double[] overApproximationforVI  = new double[stepBoundReach.length];

            for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
              underApproximationforVI[s] = stepBoundReach[s] + stepBoundStay[s] * lowerBound;
              overApproximationforVI[s] = stepBoundReach[s] + stepBoundStay[s] * upperBound;
              if(underApproximationforVI[s] > overApproximationforVI[s]){
                System.out.println("ERRORRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR...");
              }
            }
            return new double[][]{underApproximationforVI, overApproximationforVI};
          }
          lowerboundCandidate = Math.min(decisionValue, lowerboundCandidate);


        }
        lowerBound = Math.max(lowerBoundNew, lowerboundCandidate);
        lowerBoundNew = lowerBound; //remember this for next iteration

        System.out.println("lowerBound: "+ lowerBound);

        //double upper_val=Double.NEGATIVE_INFINITY;
        double upperBoundCandidate = 0.0;
        for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
          upperBoundCandidate = Math.max(upperBoundCandidate, stepBoundReachNew[s] / (1 - stepBoundStayNew[s]));
        }
        if(!min) {
          if (decisionValue > upperBoundCandidate) {
            System.out.println(
                "Need of DECISION VALUE for UBin iteration " + iters + ". DecVal: " + decisionValue
                    + ", approx_upper: " + upperBoundCandidate + ", oldupperBound: " + upperBoundNew);

            double[] underApproximationforVI  = new double[stepBoundReach.length];
            double[] overApproximationforVI  = new double[stepBoundReach.length];

            for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
              underApproximationforVI[s] = stepBoundReach[s] + stepBoundStay[s] * lowerBound;
              overApproximationforVI[s] = stepBoundReach[s] + stepBoundStay[s] * upperBound;
              if(underApproximationforVI[s] > overApproximationforVI[s]){
                System.out.println("ERRORRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR...");
              }
            }
            //TODO: call BVI here on current approximation
            return new double[][]{underApproximationforVI, overApproximationforVI};
          }
          upperBoundCandidate = Math.max(decisionValue, upperBoundCandidate);
        }
        upperBound = Math.min(upperBoundNew, upperBoundCandidate);
        upperBoundNew = upperBound; //remember this for next iteration

        System.out.println("upperBound: "+upperBound);
      }

      // Swap vectors for next iter
      // Now lowerBounds is the most up-to-date approximation, while the lowerBoundsNew contains the previous iteration
      tmpsoln = stepBoundReach;
      stepBoundReach = stepBoundReachNew;
      stepBoundReachNew = tmpsoln;
      //System.out.println("Reach: " + Arrays.toString(stepBoundReach));

      tmpsoln = stepBoundStay;
      stepBoundStay = stepBoundStayNew;
      stepBoundStayNew = tmpsoln;
      //System.out.println("Stay: " + Arrays.toString(stepBoundStay));
      double [] overApproximation = new double[stepBoundReach.length];
      double [] underApproximation = new double[stepBoundReach.length];
      for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
        underApproximation[s] = stepBoundReach[s] + stepBoundStay[s] * lowerBound;
        overApproximation[s] = stepBoundReach[s] + stepBoundStay[s] * upperBound;
      }

      //System.out.println("Under-approximation: " + Arrays.toString(underApproximation));
      //System.out.println("Over-approximation: " + Arrays.toString(overApproximation));


      /**
       * Check termination
       */
      double relevantStayVal = 0;
      if (initialState != -1) {
        relevantStayVal = stepBoundStay[initialState];
      } else {
        for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
          relevantStayVal = Math.max(relevantStayVal, stepBoundStay[s]);
        }
      }
      /**modifed termination criterion for the termination testing**/
      done = (upperBound - lowerBound) * relevantStayVal < this.termCritParam;
      //done = relevantStayVal < 2 * this.termCritParam;
      if (done) {
//			print_smg(stpg,lowerBounds,upperBounds);
        //When we are done, we have to insert the smarter values
        if (initialState != -1) {
          //Only in initial state
          //System.out.println("ReachFin: " + Arrays.toString(stepBoundReach) + "\nStay: " + Arrays.toString(stepBoundStay) + "\nLB: " + lowerBound + "\nUB: " + upperBound);
          double lb_i = stepBoundReach[initialState] + stepBoundStay[initialState] * lowerBound;
          double ub_i = stepBoundReach[initialState] + stepBoundStay[initialState] * upperBound;
        } else {
          //in all states, for topological VI. Need the second thing as temp, since meaning of content switches from reach/stayVal to actual lower/upper bound
          for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
//            stepBoundReachNew[s] = stepBoundReach[s] + stepBoundStay[s] * lowerBound;
//            stepBoundStayNew[s] = stepBoundReach[s] + stepBoundStay[s] * upperBound;
            stepBoundReachNew[s] = stepBoundReach[s];
            stepBoundStayNew[s] = stepBoundStay[s];
          }
          stepBoundReach = stepBoundReachNew.clone();
          stepBoundStay = stepBoundStayNew.clone();
        }
      } else{
        //System.out.println("Reach: " + Arrays.toString(stepBoundReach) + "\nStay: " + Arrays.toString(stepBoundStay) + "\nLB: " + lowerBound + "\nUB: " + upperBound);
        for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
          double lbprint = stepBoundReach[s] + stepBoundStay[s] * lowerBound;
          double ubprint =stepBoundReach[s] + stepBoundStay[s] * upperBound;
          ///System.out.println("bounds for state "+ s + " are: [" + lbprint +", "+ ubprint +"]");
        }
      }
    }
    if(initialState!=-1){
      double lbprint = stepBoundReach[initialState] + stepBoundStay[initialState] * lowerBound;
      double ubprint =stepBoundReach[initialState] + stepBoundStay[initialState] * upperBound;
      System.out.println("Resulting interval: ["+lbprint+","+ ubprint + "]");
    }
    for (int s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
      // storing the approximated value in stepBoundReachNew
      stepBoundReachNew[s] = stepBoundReach[s] + stepBoundStay[s] * (lowerBound+upperBound)/2;
    }

    return new double[][]{stepBoundReachNew, stepBoundReach,stepBoundStay,{iters},{lowerBound},{upperBound}};
  }



  private static int findAction(MDP mdp, int s, double[] stepBoundReach, double[] stepBoundStay, boolean min, double lowerbound, double upperbound){
    // best choice -1 means it does not exit
    int bestChoice=-1;
    int numChoices = mdp.getNumChoices(s);

    // TODO: optimize (return action 0, which is 0) in case of only one action
    if(numChoices==1){
      return 0;
    }
    if (!min) {
      double bestValueSoFar = 0;
      //mainLog.println("Searching for best leaving value; state belongs to maximizer");
      //for each state the choices are from 0 to NumChoices-1
      for (int curr_action = 0; curr_action < numChoices; curr_action++) {
        Distribution d = mdp.getChoice(s,curr_action);
        //Set<Integer> successors = d.keySet();
        double currentValue=0;
        for(int succ : d.keySet()){
          currentValue += d.get(succ) * (stepBoundReach[succ]+stepBoundStay[succ]*upperbound);
        }
        if (currentValue>=bestValueSoFar){
          bestValueSoFar = currentValue;
          bestChoice = curr_action;
        }
      }
    }else{
      double bestValueSoFar = 1; //init in case of minimizer
      //mainLog.println("Searching for best leaving value; state belongs to minimizer");
      for (int curr_action = 0; curr_action < numChoices; curr_action++) {
        Distribution d = mdp.getChoice(s,curr_action);
        double currentValue=0;
        for(int succ : d.keySet()){
          currentValue += d.get(succ) * (stepBoundReach[succ]+stepBoundStay[succ]*lowerbound);
        }
        if (currentValue<=bestValueSoFar){
          bestValueSoFar = currentValue;
          bestChoice = curr_action;
        }
      }
    }
    //System.out.println("best choice for state "+ s + " is: "+ bestChoice);
    return bestChoice;
  }



  private double computeDecisionValue(MDP mdp, double[] stepBoundReach, double[] stepBoundStay, int s, int bestChoice, boolean min){
    double decisionValue = min? 1: 0;
    // double decisionValue = max? Double.NEGATIVE_INFINITY: Double.POSITIVE_INFINITY;
    for (int curr_action = 0; curr_action < mdp.getNumChoices(s); curr_action++) {
      if (curr_action != bestChoice) {
        Distribution d_other = mdp.getChoice(s, curr_action);
        Distribution d_choice = mdp.getChoice(s, bestChoice);
        double y_delta = 0;
        Set<Integer> possible_successors = new HashSet<>(d_choice.keySet());
        possible_successors.addAll(d_other.keySet());
        for(int j : possible_successors){
          y_delta += (d_choice.get(j) - d_other.get(j)) * stepBoundStay[j];
        }
        if (y_delta > 0) {
          double x_delta = 0;

          for(int j : possible_successors){
            x_delta += (d_other.get(j) - d_choice.get(j)) * stepBoundReach[j];
          }
          if(min){
            decisionValue = Math.min(decisionValue, x_delta / y_delta);
          }
          else{
            decisionValue = Math.max(decisionValue, x_delta / y_delta);
          }
        }
      }
    }
    return decisionValue;
  }




  /**
	 * Perform the actual work of a topological interval iteration, i.e., iterate until convergence or abort.
	 *
	 * @param mc ProbModelChecker (for log and settings)
	 * @param description Description (for logging)
	 * @param sccs The information about the SCCs and topological order
	 * @param below The iteration object for the value iteration from below
	 * @param above The iteration object for the value iteration from above
	 * @param singletonSCCSolver The solver for singleton SCCs
	 * @param startTime The start time (for logging purposes, obtained from a call to System.currentTimeMillis())
	 * @param iterationsExport an ExportIterations object (optional, ignored if null)
	 * @return a ModelChecker result with the solution vector and statistics
	 * @throws PrismException on non-convergence (if mc.errorOnNonConverge is set)
	 */
	public ModelCheckerResult doTopologicalIntervalIteration(ProbModelChecker mc, String description, SCCInfo sccs, IterationIntervalIter below, IterationIntervalIter above, SingletonSCCSolver singletonSCCSolver, long timer, ExportIterations iterationsExport) throws PrismException {
		try {
			// Start iterations
			int iters = 0;
			long mvCount = 0;
			final int maxIters = mc.maxIters;

			PeriodicTimer updatesTimer = new PeriodicTimer(ProbModelChecker.UPDATE_DELAY);
			updatesTimer.start();

			int numSCCs = sccs.getNumSCCs();
			int numNonSingletonSCCs = sccs.countNonSingletonSCCs();
			int finishedNonSingletonSCCs = 0;

			boolean done = true;
			for (int scc = 0; scc < numSCCs; scc++) {
				boolean doneSCC;

				if (sccs.isSingletonSCC(scc)) {
					// get the single state in this SCC
					int state = sccs.getStatesForSCC(scc).iterator().nextInt();
					below.solveSingletonSCC(state, singletonSCCSolver);
					above.solveSingletonSCC(state, singletonSCCSolver);

					// no need to call doneWith(...), as solveSingletonSCC updates
					// both vectors for two-iteration methods

					iters++;
					mvCount += 2 * countTransitions(below.getModel(), IntSet.asIntSet(state));

					if (iterationsExport != null) {
						iterationsExport.exportVector(below.getSolnVector(), 0);
						iterationsExport.exportVector(above.getSolnVector(), 1);
					}

					intervalIterationCheckForProblems(below.getSolnVector(), above.getSolnVector(), IntSet.asIntSet(state).iterator());

					doneSCC = true;
				} else {
					// complex SCC: do VI
					doneSCC = false;
					int itersInSCC = 0;

					IntSet statesForSCC = sccs.getStatesForSCC(scc);

					// Adjust upper bound by adding 2*epsilon,
					// adding 1*epsilon would be fine, but we are a bit more conservative.
					// TODO: We also don't really need to do adjustment for bottom SCCs...
					PrimitiveIterator.OfInt it = statesForSCC.iterator();
					final double[] solnAbove = above.getSolnVector();
					final double adjustment = 2*termCritParam;
					while (it.hasNext()) {
						solnAbove[it.nextInt()] += adjustment;
					}

					// abort on convergence or if iterations *in this SCC* are above maxIters
					while (!doneSCC && itersInSCC < maxIters) {
						iters++;
						itersInSCC++;

						// do iteration step
						below.iterate(statesForSCC);
						above.iterate(statesForSCC);

						if (iterationsExport != null) {
							iterationsExport.exportVector(below.getSolnVector(), 0);
							iterationsExport.exportVector(above.getSolnVector(), 1);
						}

						intervalIterationCheckForProblems(below.getSolnVector(), above.getSolnVector(), statesForSCC.iterator());

						// Check termination (inside SCC)
						doneSCC = PrismUtils.doublesAreClose(below.getSolnVector(), above.getSolnVector(), statesForSCC.iterator(), termCritParam, absolute);

						if (!doneSCC && updatesTimer.triggered()) {
							double diff = PrismUtils.measureSupNormInterval(below.getSolnVector(), above.getSolnVector(), absolute, statesForSCC.iterator());
							mc.getLog().print("Iteration " + iters + ": ");
							mc.getLog().print("max " + (absolute ? "" : "relative ") + "diff (for iteration " + itersInSCC + " in current SCC " + (finishedNonSingletonSCCs+1) + " of " + numNonSingletonSCCs + ") = " + PrismUtils.formatDouble(diff));
							mc.getLog().println(", " + PrismUtils.formatDouble2dp(updatesTimer.elapsedMillisTotal() / 1000.0) + " sec so far");
						}
					}

					// notify the iterators that the states are done so that
					// their values can be copied to the second vector in a two-vector
					// iterator
					below.doneWith(statesForSCC);
					above.doneWith(statesForSCC);

					mvCount += 2 * itersInSCC * countTransitions(below.getModel(), statesForSCC);
					finishedNonSingletonSCCs++;
				}

				if (!doneSCC) {
					done = false;
					break;
				}
			}

			if (done) {
				double diff = PrismUtils.measureSupNormInterval(below.getSolnVector(), above.getSolnVector(), absolute);
				mc.getLog().println("Max " + (absolute ? "" : "relative ") +
						"diff between upper and lower bound on convergence: " + PrismUtils.formatDouble(diff));
				done = true;
			}

			// Finished value iteration
			timer = System.currentTimeMillis() - timer;
//			mc.getLog().print("Interval iteration (" + description + ", with " + numNonSingletonSCCs + " non-singleton SCCs)");
//			mc.getLog().print(" took " + iters + " iterations, ");
//			mc.getLog().print(mvCount + " multiplications");
//			mc.getLog().println(" and " + timer / 1000.0 + " seconds.");

			if (done && OptionsIntervalIteration.from(mc.getSettings()).isSelectMidpointForResult()) {
				PrismUtils.selectMidpoint(below.getSolnVector(), above.getSolnVector());

				if (iterationsExport != null) {
					// export midpoint
					iterationsExport.exportVector(below.getSolnVector(), 0);
					iterationsExport.exportVector(below.getSolnVector(), 1);
				}
			}

			if (iterationsExport != null)
				iterationsExport.close();

			// Non-convergence is an error (usually)
			if (!done && mc.errorOnNonConverge) {
				String msg = "Iterative method (interval iteration) did not converge within " + iters + " iterations.";
				msg += "\nConsider using a different numerical method or increasing the maximum number of iterations";
				throw new PrismException(msg);
			}

			// Return results
			ModelCheckerResult res = new ModelCheckerResult();
			res.soln = below.getSolnVector();
			res.numIters = iters;
			res.timeTaken = timer / 1000.0;
			return res;
		} finally {
			if (iterationsExport != null)
				iterationsExport.close();
		}
	}

	/**
	 * Compares the current lower and upper solution vectors in an interval iteration
	 * and throws an exception if lower bound values are larger than upper bound values,
	 * as this indicates problems.
	 * @param lower the current lower iteration solution vector
	 * @param upper the current upper iteration solution vector
	 * @param states iterator over the states in question
	 */
	public static void intervalIterationCheckForProblems(double[] lower, double[] upper, PrimitiveIterator.OfInt states) throws PrismException
	{
		while (states.hasNext()) {
			int s = states.nextInt();
			if (upper != null && lower[s] > upper[s] && !PrismUtils.doublesAreClose(lower[s], upper[s], 1E-12, true)) {
				throw new PrismException("In interval iteration, the lower value (" + lower[s] + ") is larger than the upper value (" + upper[s] + ").\n"
						+ "This indicates either problems with numerical stability (rounding, precision of the floating-point representation) or that the initial bounds (for reward computations) are incorrect");
			}
		}
	}

	/**
	 * Perform a post-processing for a two-vector value iteration.
	 * @param solnOld the previous solution vector
	 * @param solnNew the newly computed solution vector
	 * @param states the relevant set of states
	 * @param fromBelow are we iterating from below?
	 * @param enforceMonotonicity if true, enforces monotonicity
	 * @param checkMonotonicity if true, checks for monotonicity (and throws error when non-monotonic)
	 */
	public static void twoVectorPostProcessing(double[] solnOld, double[] solnNew, IntSet states, boolean fromBelow, boolean enforceMonotonicity, boolean checkMonotonicity) throws PrismException
	{
		// TODO: use IntSet states
		if (enforceMonotonicity)
			if (fromBelow) {
				PrismUtils.ensureMonotonicityFromBelow(solnOld, solnNew);
			} else {
				PrismUtils.ensureMonotonicityFromAbove(solnOld, solnNew);
			}

		if (checkMonotonicity) {
			PrismUtils.checkMonotonicity(solnOld, solnNew, !fromBelow);
		}
	}

	protected long countTransitions(Model model, IntSet unknownStates)
	{
		if (model instanceof DTMC) {
			return ((DTMC)model).getNumTransitions(unknownStates.iterator());
		} else if (model instanceof MDP) {
			return ((MDP)model).getNumTransitions(unknownStates.iterator());
		} else {
			throw new IllegalArgumentException("Can only count transitions for DTMCs and MDPs, not for " + model.getModelType());
		}
	}

}