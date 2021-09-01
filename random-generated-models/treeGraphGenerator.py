import random
import graphGenerator
from graphGenParams import GraphGenerationParameters

class TreeGraphGenerator(graphGenerator.GeneratedGraph):
    """
    Respects parameters:
        probability_to_branch
        probability_for_backwards_actions
        minimum_outgoing_edges
        num_states

    Respects roughly:
        probability_for_backwards_action: These probabilities only apply to leaf nodes

    """
    def __init__(self) -> None:
        super().__init__()

    def generateGraph(self, params: GraphGenerationParameters):

        self._initVars(params)
        
        self._generateStates()

        self._generateChoices()

        self._reduceTrivialStates()

        self._ensureDeadlockFreedom()

        self._computeMaxActionsPerPlayer()

    def _generateChoices(self):
        choices_per_state = self.params.minimum_outgoing_edges
        for state in range(self.params.num_states):
            leftest_child = state*choices_per_state+1

            for next_state in range(leftest_child, leftest_child+choices_per_state):
                if (next_state >= self.params.num_states):
                    self._possiblyAddReturnAction(state)

                else:
                    distribution = self._generateChoice(state, next_state, self.params.probability_to_branch) 
                    self.actions_map[state].append(distribution)

    def _possiblyAddReturnAction(self, state):
        if (self._meetsThreshhold(self.params.probability_for_backwards_action)):
            randomState = random.randint(0, state-1)
            distribution = self._generateChoice(state, randomState, 0)
            if (distribution not in self.actions_map[state]):
                self.actions_map[state].append(distribution)