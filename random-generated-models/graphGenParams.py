from randomStateGetter import *

class GraphGenerationParameters:
    def __init__(self,
    num_states : int,
    minimum_incoming_edges : int = 1, 
    maximum_incoming_edges : int = 10,
    probability_to_branch : float = 0.2, 
    probability_for_backwards_action : float = 0.05,
    probability_to_be_maximizer_state : float = 0.5,
    guaranteed_sinks : int = 1,
    minimum_outgoing_edges : int = 1,
    transition_permutator : Permutation = Permutation(),
    choice_permutator : Permutation = Permutation(),
    ) -> None:
        self.num_states = num_states
        self.minimum_incoming_edges = minimum_incoming_edges 
        self.maximum_incoming_edges = maximum_incoming_edges 
        self.probability_to_branch = probability_to_branch 
        self.probability_for_backwards_action = probability_for_backwards_action
        self.probability_to_be_maximizer_state = probability_to_be_maximizer_state
        self.guaranteed_sinks = guaranteed_sinks
        self.minimum_outgoing_edges = minimum_outgoing_edges
        self.transition_permutator = transition_permutator
        self.choice_permutator = choice_permutator