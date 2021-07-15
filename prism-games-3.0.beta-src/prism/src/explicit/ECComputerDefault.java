//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham/Oxford)
//	* Mateusz Ujma <mateusz.ujma@cs.ox.ac.uk> (University of Oxford)
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

import java.util.*;

import prism.PrismComponent;
import prism.PrismException;

/**
 * Explicit maximal end component computer for a nondeterministic model such as an MDP.
 * Implements the algorithm from p.48 of:
 * Luca de Alfaro. Formal Verification of Probabilistic Systems. Ph.D. thesis, Stanford University (1997)
 */
public class ECComputerDefault extends ECComputer
{
	/** The model to compute (M)ECs for **/
	private NondetModel model;

	/** Computed list of MECs **/
	private List<BitSet> mecs = new ArrayList<BitSet>();

	/**
	 * Build (M)EC computer for a given model.
	 */
	public ECComputerDefault(PrismComponent parent, NondetModel model) throws PrismException
	{
		super(parent);
		this.model = model;
	}

	// Methods for ECComputer interface

	@Override
	public void computeMECStates() throws PrismException
	{
		mecs = findEndComponents(null, null);
	}

	@Override
	public void computeMECStates(BitSet restrict) throws PrismException
	{
		mecs = findEndComponents(restrict, null);
	}

	@Override
	public void computeMECStates(BitSet restrict, BitSet accept) throws PrismException
	{
		mecs = findEndComponents(restrict, accept);
	}

	@Override
	public List<BitSet> getMECStates()
	{
		return mecs;
	}

	// Computation
	
	/**
	 * Find all accepting maximal end components (MECs) in the submodel obtained
	 * by restricting this one to the set of states {@code restrict},
	 * where acceptance is defined as those which intersect with {@code accept}.
	 * If {@code restrict} is null, we look at the whole model, not a submodel.
	 * If {@code accept} is null, the acceptance condition is trivially satisfied.
	 * @param restrict BitSet for the set of states to restrict to
	 * @param accept BitSet for the set of accepting states
	 * @return a list of BitSets representing the MECs
	 */
	private List<BitSet> findEndComponents(BitSet restrict, BitSet accept) throws PrismException
	{
		// If restrict is null, look within set of all reachable states
		if (restrict == null) {
			restrict = new BitSet();
			restrict.set(0, model.getNumStates());
		}
		// Initialise L with set of all states to look in (if non-empty)
		List<BitSet> L = new ArrayList<BitSet>();
		if (restrict.isEmpty())
			return L;
		L.add(restrict);
		// Find MECs
		boolean changed = true;
		while (changed) {
			changed = false;
			BitSet E = L.remove(0);
			SubNondetModel submodel = restrict(model, E);
			List<BitSet> sccs = translateStates(submodel, computeSCCs(submodel));
			L = replaceEWithSCCs(L, E, sccs);
			changed = canLBeChanged(L, E);
		}
		// Filter and return those that contain a state in accept
		if (accept != null) {
			int i = 0;
			while (i < L.size()) {
				if (!L.get(i).intersects(accept)) {
					L.remove(i);
				} else {
					i++;
				}
			}
		}
		return L;
	}

	private Set<BitSet> processedSCCs = new HashSet<BitSet>();

	private boolean canLBeChanged(List<BitSet> L, BitSet E)
	{
		processedSCCs.add(E);
		for (int i = 0; i < L.size(); i++) {
			if (!processedSCCs.contains(L.get(i))) {
				return true;
			}
		}
		return false;
	}

	private List<BitSet> replaceEWithSCCs(List<BitSet> L, BitSet E, List<BitSet> sccs)
	{
		if (sccs.size() > 0) {
			List<BitSet> toAdd = new ArrayList<BitSet>();
			for (int i = 0; i < sccs.size(); i++) {
				if (!L.contains(sccs.get(i))) {
					toAdd.add(sccs.get(i));
				}
			}
			if (toAdd.size() > 0) {
				L.addAll(toAdd);
			}
		}
		return L;
	}

	private SubNondetModel restrict(NondetModel model, BitSet states)
	{
		Map<Integer, BitSet> actions = new HashMap<Integer, BitSet>();
		BitSet initialStates = new BitSet();
		initialStates.set(states.nextSetBit(0));

		boolean changed = true;
		while (changed) {
			changed = false;
			actions.clear();
			for (int i = 0; i < model.getNumStates(); i++) {
				BitSet act = new BitSet();
				if (states.get(i)) {
					for (int j = 0; j < model.getNumChoices(i); j++) {
						if (model.allSuccessorsInSet(i, j, states)) {
							act.set(j);
						}
					}
					if (act.isEmpty()) {
						states.clear(i);
						changed = true;
					}
					actions.put(i, act);
				}
			}
		}

		return new SubNondetModel(model, states, actions, initialStates);
	}

	private List<BitSet> computeSCCs(NondetModel model) throws PrismException
	{
		SCCConsumerStore sccs = new SCCConsumerStore();
		SCCComputer sccc = SCCComputer.createSCCComputer(this, model, sccs);
		sccc.computeSCCs();
		return sccs.getSCCs();
	}

	private List<BitSet> translateStates(SubNondetModel model, List<BitSet> sccs)
	{
		List<BitSet> r = new ArrayList<BitSet>();
		for (int i = 0; i < sccs.size(); i++) {
			BitSet set = sccs.get(i);
			BitSet set2 = new BitSet();
			r.add(set2);
			for (int j = set.nextSetBit(0); j >= 0; j = set.nextSetBit(j + 1)) {
				set2.set(model.translateState(j));

			}
		}
		return r;
	}

	private boolean isMEC(BitSet b)
	{
		if (b.isEmpty())
			return false;

		int state = b.nextSetBit(0);
		while (state != -1) {
			boolean atLeastOneAction = false;
			for (int i = 0; i < model.getNumChoices(state); i++) {
				if (model.allSuccessorsInSet(state, i, b)) {
					atLeastOneAction = true;
				}
			}
			if (!atLeastOneAction) {
				return false;
			}
			state = b.nextSetBit(state + 1);
		}

		return true;
	}

	/**
	 * Described in KKKW18, SEC
	 * Finds simple end components in the MEC given to it, fixing the decisions of the minimizer according to L
	 */
	public List<BitSet> getSECs(BitSet mec, double[] L, boolean min1, boolean min2) throws PrismException{
		//TODO: Investigate why caching gPrime for every mec was slower

		if (! (model instanceof STPGExplicit)){
			throw new PrismException("getSECs only works on STPGs.");
		}

		//Catch a corner case
		if (mec.size()==1){
			List<BitSet> result = new ArrayList<>();
			result.add(mec);
			return result;
		}

		// Construct a new STPG gPrime for the MEC, in which we then compute the SECs; add dummy sink for stuff outside MEC
		// It has numStates = |MEC| and we keep a HashMap to transfer the state numbers
		// We only add optimal choices for the Minimizer

		STPGExplicit gPrime = new STPGExplicit(0);
		HashMap<Integer,Integer> g2prime = new HashMap<>();
		HashMap<Integer,Integer> prime2g = new HashMap<>();

		//Traverse MEC once to set both HashMaps (and initialize states with players in gPrime)
		for (int s = mec.nextSetBit(0); s >= 0; s = mec.nextSetBit(s+1)){
			int prime_s = gPrime.addState(((STPGExplicit) model).getPlayer(s));
			g2prime.put(s,prime_s);
			prime2g.put(prime_s,s);
		}

		//Add sink for transitions leading outside of MEC
		int sink = gPrime.addState(1);
		Distribution sinkLoop = new Distribution();
		sinkLoop.add(sink,1);
		gPrime.addChoice(sink,sinkLoop);

		//Traverse a second time to add choices (only possibe after both maps were initialized)
		for (int sPrime = 0; sPrime < gPrime.numStates; sPrime++){
		    if(sPrime==sink){
		        continue;
            }
		    int s = prime2g.get(sPrime);
			// For each state, set the player and the successors.
			int player = gPrime.getPlayer(sPrime);
			// If minPlayer, have to keep only optimal action. Else (if maxPlayer), keep all action.
			// Also: Only keep actions that reach stuff in the MEC
			boolean minPlayer = (player == 1 && min1) || (player == 2 && min2);
			if (minPlayer){
				//Find action with minimum lower bound L(s,a)
				double bestLsa = 1.0;
				Set<Integer> bestActions = new HashSet<>();
				for (int a = 0; a<model.getNumChoices(s);a++){
					Distribution distr = ((STPGExplicit) model).trans.get(s).get(a);
					double Lsa = 0.0;
					for (Map.Entry<Integer, Double> e : distr) {
						int succ = (Integer) e.getKey();
						double prob = (Double) e.getValue();
						Lsa += prob * L[succ];
					}
					if (Lsa < bestLsa){
						bestLsa = Lsa;
						bestActions.clear();
						bestActions.add(a);
					}
					if (Lsa == bestLsa){
						bestActions.add(a);
					}
				}
				for (Integer a : bestActions) {
					gPrime.addChoice(sPrime, makePrimeDistribution(((STPGExplicit) model).getChoice(s, a), g2prime, sink));
				}
			}
			else{
				for (int a = 0; a<model.getNumChoices(s);a++){
					gPrime.addChoice(sPrime, makePrimeDistribution(((STPGExplicit) model).getChoice(s,a),g2prime,sink));
				}
			}
		}

		// MECs in gPrime are SECs in input G
		try {
			explicit.ECComputerDefault ec = (ECComputerDefault) ECComputer.createECComputer(this, gPrime);
			ec.computeMECStates();
			List<BitSet> secsPrime =  ec.getMECStates();
			List<BitSet> secs = new ArrayList<>();
			// Need to transfer indices back
			for (BitSet secPrime : secsPrime) {
				if (secPrime.nextSetBit(0) == sink){
					continue;
				}
				BitSet sec = new BitSet(0);
				for (int s = secPrime.nextSetBit(0); s>=0; s = secPrime.nextSetBit(s+1)) {
					sec.set(prime2g.get(s));
				}
				secs.add(sec);
			}
			return secs;
		} catch (PrismException e) {
			System.err.println("Exception in SEC computing in ECComputerDefault."); //Bad style, but short on time
			return null;
		}
	}

	private Distribution makePrimeDistribution(Distribution d, HashMap<Integer,Integer> g2prime, int sink){
		Distribution distrNew = new Distribution();
		Iterator<Map.Entry<Integer, Double>> i = d.iterator();
		while (i.hasNext()) {
			Map.Entry<Integer, Double> e = i.next();
			int newKey = g2prime.containsKey(e.getKey()) ? g2prime.get(e.getKey()) : sink;
			distrNew.add(newKey, e.getValue());
		}
		return distrNew;
	}

	/**
	 * Same method as above, but with caching. Was slower in the one time I tried, so I was disappointed and commented it out.
	 * @param mec
	 * @param L
	 * @param min1
	 * @param min2
	 * @return
	 * @throws PrismException
	 */
//	public List<BitSet> getSECs(BitSet mec, double[] L, boolean min1, boolean min2) throws PrismException{
//		if (! (model instanceof STPGExplicit)){
//			throw new PrismException("getSECs only works on STPGs.");
//		}
//
//		// Construct a new STPG gPrime for the MEC, in which we then compute the SECs; add dummy sink for stuff outside MEC
//		// It has numStates = |MEC| and we keep a HashMap to transfer the state numbers
//		// We only add optimal choices for the Minimizer
//		//We cache the gPrime, as the only thing that changes in them is the Minimizer choices, and as they should be small (except if some MEC is very large)
//
//		STPGExplicit gPrime;
//		HashMap<Integer,Integer> g2prime;
//		HashMap<Integer,Integer> prime2g;
//		int sink; // dummy sink state for everything leading outside the MEC
//
//		if(mec2gPrime.containsKey(mec)){
//			gPrime = mec2gPrime.get(mec);
//			g2prime = mec2_g2p.get(mec);
//			prime2g = mec2_p2g.get(mec);
//			sink = gPrime.numStates-1;
//		}
//		else{
//			gPrime = new STPGExplicit(0);
//			g2prime = new HashMap<Integer,Integer>();
//			prime2g = new HashMap<Integer,Integer>();
//
//			//Traverse MEC once to set both HashMaps (and initialize states with players in gPrime)
//			for (int s = mec.nextSetBit(0); s >= 0; s = mec.nextSetBit(s+1)){
//				int prime_s = gPrime.addState(((STPGExplicit) model).getPlayer(s));
//				g2prime.put(s,prime_s);
//				prime2g.put(prime_s,s);
//			}
//
//			//Add sink for transitions leading outside of MEC
//			sink = gPrime.addState(1);
//			Distribution sinkLoop = new Distribution();
//			sinkLoop.add(sink,1);
//			gPrime.addChoice(sink,sinkLoop);
//
//			//Traverse a second time to add maximizer choices (only possibe after both maps were initialized)
//			for (int sPrime = 0; sPrime < gPrime.numStates; sPrime++){
//				if(sPrime==sink){
//					continue;
//				}
//				int s = prime2g.get(sPrime);
//				int player = gPrime.getPlayer(sPrime);
//				// If minPlayer, have to keep only optimal action; we do this outside the if(cached), since it has to be done every time. Else (if maxPlayer), keep all action.
//				boolean minPlayer = (player == 1 && min1) || (player == 2 && min2);
//				if (minPlayer){
//					continue;
//				}
//				else{
//					for (int a = 0; a<model.getNumChoices(sPrime);a++){
//						gPrime.addChoice(sPrime, makePrimeDistribution(((STPGExplicit) model).getChoice(s,a),g2prime,sink));
//					}
//				}
//			}
//			mec2gPrime.put(mec,gPrime);
//			mec2_g2p.put(mec,g2prime);
//			mec2_p2g.put(mec,prime2g);
//		}
//
//
//		//Traverse MEC to add Minimzer choices (clear states beforehand)
//		for (int sPrime = 0; sPrime < gPrime.numStates; sPrime++){
//			if(sPrime==sink){
//				continue;
//			}
//			int s = prime2g.get(sPrime);
//			int player = gPrime.getPlayer(sPrime);
//			// If minPlayer, have to keep only optimal action. Else (if maxPlayer), was initialized before
//			boolean minPlayer = (player == 1 && min1) || (player == 2 && min2);
//			if (minPlayer){
//				gPrime.clearState(sPrime);
//				for (int a = 0; a<model.getNumChoices(s);a++){
//					Distribution distr = ((STPGExplicit) model).trans.get(s).get(a);
//					double Lsa = 0.0;
//					for (Map.Entry<Integer, Double> e : distr) {
//						int succ = (Integer) e.getKey();
//						double prob = (Double) e.getValue();
//						Lsa += prob * L[succ];
//					}
//					if (Lsa == L[s]){
//						gPrime.addChoice(sPrime, makePrimeDistribution(((STPGExplicit) model).getChoice(s,a),g2prime,sink));
//					}
//				}
//			}
//		}
//
//		// MECs in gPrime are SECs in input G
//		try {
//			explicit.ECComputerDefault ec = (ECComputerDefault) ECComputer.createECComputer(this, gPrime);
//			ec.computeMECStates();
//			List<BitSet> secsPrime =  ec.getMECStates();
//			List<BitSet> secs = new ArrayList<>();
//			// Need to transfer indices back
//			for (BitSet secPrime : secsPrime) {
//				if (secPrime.nextSetBit(0) == sink){
//					continue;
//				}
//				BitSet sec = new BitSet(0);
//				for (int s = secPrime.nextSetBit(0); s>=0; s = secPrime.nextSetBit(s+1)) {
//					sec.set(prime2g.get(s));
//				}
//				secs.add(sec);
//			}
//			return secs;
//		} catch (PrismException e) {
//			System.err.println("Exception in SEC computing in ECComputerDefault."); //Bad style, but short on time
//			return null;
//		}
//	}

}
