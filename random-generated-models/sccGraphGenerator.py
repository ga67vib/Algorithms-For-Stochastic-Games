import random
from typing import Any, Dict, List
import graphGenerator
import dataclasses
import graphGeneratorTargetTrans
from graphGenParams import GraphGenerationParameters
from randomStateGetter import Permutation_AllStatesPossible
import tarjansAlgorithm
import numpy as np

class SCCGraphGenerator(graphGenerator.GeneratedGraph):
    """
    Respects parameters:
        probability_to_branch
        probability_for_backwards_actions
        minimum_outgoing_edges
        num_states
        force_unknown
        denominator_range

    Respects roughly:
        probability_for_backwards_action: These probabilities only apply to leaf nodes

    """
    def __init__(self) -> None:
        super().__init__()

    def generateGraph(self, params: GraphGenerationParameters):

        self._initVars(params)

        # ==== PARAMETERS THAT SHOULD BE ADJUSTABLE ====
        min_SCC_size = 9990
        max_SCC_size = max(
            min_SCC_size, 
            min(
                self.params.num_states,
                9990
                )
            )
        
        min_incoming_actions_per_SCC = 200
        max_incoming_actions_per_SCC = 400

        # If we decide that 
        min_actions_per_SCC_connection = 1
        max_actions_per_SCC_connection = 1
        # ============== PARAMETER END ================

        numberOfStatesLeftToGenerate = self.params.num_states
        state_offset = 0

        while(numberOfStatesLeftToGenerate > 0):
            # Generate a random Graph that will be forced to be a SCC
            sccSize = min(random.randint(min_SCC_size, max_SCC_size),numberOfStatesLeftToGenerate)
            sccParameters = dataclasses.replace(self.params, num_states = sccSize, force_unknown = False, verbose = False)
            sccParameters.choice_permutator = Permutation_AllStatesPossible()
            sccParameters.transition_permutator = Permutation_AllStatesPossible()
            sccParameters.denominator_range = self.params.denominator_range
            sccParameters.maximum_incoming_edges = self.params.maximum_incoming_edges
            sccParameters.minimum_outgoing_edges = self.params.minimum_incoming_edges
            sccParameters.maximum_outgoing_edges = self.params.maximum_outgoing_edges
            sccParameters.minimum_outgoing_edges = self.params.minimum_outgoing_edges
            sccParameters.probability_to_branch = self.params.probability_to_branch

            graph = graphGenerator.GeneratedGraph()
            graph.generateGraph(
                sccParameters
            )

            # Make sure that the generated graph is truly a SCC
            self._forceSCC(graph, min_actions_per_SCC_connection, max_actions_per_SCC_connection)
            
            [offset_states_of_player_1, offset_states_of_player_2, offset_actions_map] = self._offset_graph(graph, state_offset)
            
            # Merge the Graph to the rest
            self.states_of_player1 += offset_states_of_player_1
            self.states_of_player2 += offset_states_of_player_2
            self.actions_map = {**self.actions_map , **offset_actions_map}

            offset_states = list(offset_actions_map.keys())
            # Connect previous states to this SCC
            if state_offset > 0:
                for _ in range(random.randint(min_incoming_actions_per_SCC, max_incoming_actions_per_SCC)):
                    state_from_previous_sccs = random.randint(0, state_offset-1)
                    state_from_newly_added_scc = offset_states[random.randint(0, len(offset_states)-1)]+state_offset
                    choice = self._generateChoice(state_from_previous_sccs, state_from_newly_added_scc, 0) # Set branching probability to 0 to ensure the SCC
                    self.actions_map[state_from_previous_sccs].append(choice)

            # Prepare for next loop
            state_offset += graph.params.num_states
            numberOfStatesLeftToGenerate -= graph.params.num_states

        if (params.force_unknown):
            self._reduceTrivialStates()

        self._ensureDeadlockFreedom()

        if (params.force_unknown):
            self._reduceNoStates()

        self._computeMaxActionsPerPlayer()

        # Go to target - Redirect last transition
        for i in range(self.params.num_states-2):
            actions = self.actions_map[i]
            for action in actions:
                if ((self.params.num_states-1) in action.keys()):
                    continue

                transition_goal_key = list(action.keys())[0]

                transition_goals = [self.params.num_states-2, self.params.num_states-1, i, i+1]
                not_included = len(transition_goals)+1
                for s in transition_goals:
                    if s in action.keys():
                        not_included-=1

                factor = not_included / len(transition_goals)+1

                val = action[transition_goal_key]
                if (val == "1"):
                    val = "1/1"
                split = val.split("/")
                split[0] = str(int(split[0]))
                new_denom = str(int(split[1])*5)

                action[transition_goal_key] = split[0]+"/"+new_denom
                self._addIfNotThere(action, self.params.num_states-2, split[0]+"/"+new_denom)
                self._addIfNotThere(action, self.params.num_states-1, split[0]+"/"+new_denom)
                self._addIfNotThere(action, i, split[0]+"/"+new_denom)
                self._addIfNotThere(action, i+1, split[0]+"/"+new_denom)

    def _addIfNotThere(self, dictionary, goal, prob ):
        if goal in dictionary.keys():
            dictionary[goal] = f'({prob} + {dictionary[goal]})'
        else:
            dictionary[goal] = prob

    def _forceSCC(self, graph : graphGenerator.GeneratedGraph, min_actions_per_SCC_connect, max_actions_per_SCC_connect):
        # Get SCCs
        tarjan = tarjansAlgorithm.TarjansAlgorithm()
        sccs = tarjan.tarjans_algorithm(graph.actions_map, graph.params.num_states)
        
        # Connect last SCC to first one
        scc_keys : List[Any] = list(sccs.keys())
        last_scc_key = scc_keys[-1]
        first_scc_key = scc_keys[0]

        for _ in range(random.randint(min_actions_per_SCC_connect, max_actions_per_SCC_connect)):
            state_index_in_scc = random.randint(0, len(sccs[last_scc_key])-1)
            state_of_last_scc = sccs[last_scc_key][state_index_in_scc]

            state_index_in_scc = random.randint(0, len(sccs[first_scc_key])-1)
            state_of_first_scc = sccs[first_scc_key][state_index_in_scc]
            choice = self._generateChoice(state_of_last_scc, state_of_first_scc, graph.params.probability_to_branch, upper_bound_state_index=graph.params.num_states)
            if (choice not in graph.actions_map[state_of_last_scc]):
                graph.actions_map[state_of_last_scc].append(choice)

        previous_scc_key = scc_keys.pop()
        while len(scc_keys):
            # Connect current SCC to last SCC
            current_scc_key = scc_keys.pop()

            for _ in range(random.randint(min_actions_per_SCC_connect, max_actions_per_SCC_connect)):
                state_index_in_scc = random.randint(0, len(sccs[previous_scc_key])-1)
                state_of_previous_scc = sccs[previous_scc_key][state_index_in_scc]

                state_index_in_scc = random.randint(0, len(sccs[current_scc_key])-1)
                state_of_current_scc = sccs[current_scc_key][state_index_in_scc]
                choice = self._generateChoice(state_of_current_scc, state_of_previous_scc, graph.params.probability_to_branch, upper_bound_state_index=graph.params.num_states)
                if (choice not in graph.actions_map[state_of_current_scc]):
                    graph.actions_map[state_of_current_scc].append(choice)

            previous_scc_key = current_scc_key
            
    def _offset_graph(self, graph : graphGenerator.GeneratedGraph, offset : int):
        states_of_player_1_new = (np.array(graph.states_of_player1) + offset).tolist()
        states_of_player_2_new = (np.array(graph.states_of_player2) + offset).tolist()
        actions_map_new : Dict[int, List[Dict[int, str]]] = dict()

        for state in graph.actions_map.keys():
            state_actions = graph.actions_map[state]
            new_state_actions = []

            for action in state_actions:
                new_action : Dict[int, str] = dict()
                for target in action:
                    new_action[target+offset] = action[target]

                new_state_actions.append(new_action)

            actions_map_new[state+offset] = new_state_actions

        return [states_of_player_1_new, states_of_player_2_new, actions_map_new]


    def _generateChoices(self):
        choices_per_state = random.randint(1, self.params.minimum_outgoing_edges)
        for state in range(self.params.num_states):
            leftest_child = state*choices_per_state+1

            for next_state in range(leftest_child, leftest_child+choices_per_state):
                if (next_state >= self.params.num_states):
                    self._possiblyAddReturnAction(state)

                else:
                    distribution = self._generateChoice(state, next_state, self.params.probability_to_branch) 
                    self.actions_map[state].append(distribution)
    
    def _reduceNoStates(self):
        """Add to leave nodes action with small target probability"""
        choices_per_state = self.params.minimum_outgoing_edges
        target = self.params.num_states -1
        sink = target -1

        for state in range(self.params.num_states):
            leftest_child = state*choices_per_state+1
            for next_state in range(leftest_child, leftest_child+choices_per_state):
                if (next_state >= self.params.num_states):
                    distribution = dict()
                    distribution[target] = self._probabilityToFractionString(1, 100)
                    distribution[sink] = self._probabilityToFractionString(99, 100)
                    self.actions_map[state].append(distribution)
                    break

    def _possiblyAddReturnAction(self, state):
        if (self._meetsThreshhold(self.params.probability_for_backwards_action)):
            randomState = random.randint(0, state-1)
            distribution = self._generateChoice(state, randomState, 0.2)
            if (distribution not in self.actions_map[state]):
                self.actions_map[state].append(distribution)


def test():
    params : GraphGenerationParameters = GraphGenerationParameters(
        num_states=15,
        minimum_incoming_edges=1,
        maximum_incoming_edges=3,
        probability_to_branch=1,
        probability_to_be_maximizer_state=0.5,
        probability_for_backwards_action=1,
        minimum_outgoing_edges=1,
        force_unknown=True
    )

    graph = SCCGraphGenerator()
    graph.generateGraph(params)
    print('done')

#test()