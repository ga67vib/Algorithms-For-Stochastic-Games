import random
import graphGenerator
from graphGenParams import GraphGenerationParameters
import dataclasses

class MECGraphGenerator(graphGenerator.GeneratedGraph):
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

        generatedSubgraphs = []

        # ==== PARAMETERS THAT SHOULD BE ADJUSTABLE ====
        min_MEC_size = 3
        max_MEC_size = max(
            min_MEC_size, 
            min(
                self.params.num_states,
                8 
                )
            )
        
        min_incoming_actions_per_MEC = 3
        max_incoming_actions_per_MEC = 5

        # If we decide that 
        max_actions_per_MEC_connection = 3
        min_actions_per_MEC_connection = 1
        # ============== PARAMETER END ================

        numberOfStatesLeftToGenerate = self.params.num_states
        state_offset = 0

        while(numberOfStatesLeftToGenerate > 0):
            # Generate a random Graph that will be forced to be a MEC
            mecSize = random.randint(min_MEC_size, max_MEC_size)
            mecParameters = dataclasses.replace(self.params, num_states = mecSize, force_unknown = False)

            graph = graphGenerator.GeneratedGraph()
            graph.generateGraph(
                mecParameters
            )

            # Make sure that the generated graph is truly a MEC
            




            numberOfStatesLeftToGenerate-= mecParameters.num_states


        
        
        self._generateStates()

        self._generateChoices()

        if (params.force_unknown):
            self._reduceTrivialStates()

        self._ensureDeadlockFreedom()

        if (params.force_unknown):
            self._reduceNoStates()

        self._computeMaxActionsPerPlayer()

    def _forceMEC(self, graph : graphGenerator.GeneratedGraph, graph_params : GraphGenerationParameters):
        # Add an action from target to the initial state
        target = graph_params.num_states
        targetIsMaximizer = target in graph.states_of_player1

        choice = self._generateChoice(target, 0, 0)
        graph.actions_map[target].append(choice)

        # Make sure the rest of the graph is a MEC too

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
            distribution = self._generateChoice(state, randomState, 0)
            if (distribution not in self.actions_map[state]):
                self.actions_map[state].append(distribution)