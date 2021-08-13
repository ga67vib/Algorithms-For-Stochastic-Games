import random
class GeneratedGraph:
    def generateGraph(
        self,
        num_states : int,
        minimum_incoming_edges : int = 1, 
        maximum_incoming_edges : int = 10, 
        probability_to_branch : float = 0.2, 
        probability_for_backwards_action : float = 0.05,
        probability_to_be_maximizer_state : float = 0.5,
        guaranteed_sinks : int = 1
        ):
        """
        Generates a Graph

            Parameters:
                num_states: how many states should there be?

                minimum_incoming_edges: how many actions should at least lead into this state

                maximum_incoming_edges: how many actions should atmost leat into this state

                    NOTE: The actual incoming edges may violate these parameters sometimes. 
                        Every state except for the initial state has always at least one action leading into it

        """

        self.states_of_player1 = []
        self.states_of_player2 = []
        self.actions_map = dict() #Each state has a List of dicts, where each dict is mapping of states to probabilities 

        for state in range(num_states):
            # Assign State to a player
            if self._meetsThreshhold(probability_to_be_maximizer_state):
                self.states_of_player1.append(state)
            else:
                self.states_of_player2.append(state)
            self.actions_map[state] = []

            # Connect backwards to states
            choice_targets_permutation = Permutation()

            # Cannot have more actions leading into a state than overall available states
            num_incoming_actions = random.randint(min(state, minimum_incoming_edges), min(state, maximum_incoming_edges))
            already_included_states = []
            for _ in range(num_incoming_actions):
                outgoing_state = choice_targets_permutation.next(state, exclude = already_included_states)

                #Add choice
                already_included_states.append(outgoing_state)
                distribution = self._generateChoice(outgoing_state, state, probability_to_branch)
                self.actions_map[outgoing_state].append(distribution)

        #Same procedure backwards to add MECs:
        count_down = list(range(num_states-1))
        count_down.reverse()
        for state in count_down:
            if (self._meetsThreshhold(probability_for_backwards_action)):
                # TODO Currently per state atmost one backwards action
                target_state = random.randint(0, state)
                distribution = self._generateChoice(state, target_state, probability_to_branch)
                self.actions_map[state].append(distribution)

        #Every state that has no action should have self-loop
        for state in range(num_states):
            if (len(self.actions_map[state]) == 0):
                distribution = dict()
                distribution[state] = self._probabilityToFractionString(1,1)
                self.actions_map[state].append(distribution)


        #Compute maximal number of actions used
        self.max_player_1_actions = 0
        self.max_player_2_actions = 0
        for state in self.actions_map:
            if (state in self.states_of_player1):
                self.max_player_1_actions = max(len(self.actions_map[state]),self.max_player_1_actions)
            else:
                self.max_player_2_actions = max(len(self.actions_map[state]),self.max_player_2_actions)


    def _generateChoice(self, source_state, target_state, probability_to_branch):
        distribution = dict()
        if (self._meetsThreshhold(probability_to_branch)):
            distribution[target_state] = self._probabilityToFractionString(1,1)
        else:
            distribution = self._getProbabilityDistributionForBranchingAction(target_state)
        return distribution
                

    def _getProbabilityDistributionForBranchingAction(self, target_state):
        """TODO"""
        transition_targets_permutation = Permutation()
        distribution = dict()
        other_state = transition_targets_permutation.next(target_state)
        if (target_state != other_state):
            distribution[target_state] = self._probabilityToFractionString(1,2)
            distribution[other_state] = self._probabilityToFractionString(1,2)
        else:
            distribution[target_state] = self._probabilityToFractionString(1,1)
        return distribution

    def _meetsThreshhold(self, probability):
        return float(random.randint(1, 100))/100.0 <= probability

    def _probabilityToFractionString(self, numerator, denominator):
        if (numerator == denominator):
            return "1"
        else:
            return str(numerator)+"/"+str(denominator)

class Permutation:
    def __init__(self) -> None:
        self.currentPermutation = []

    def _next(self, lowerBound, upperBound):
        """Gets a random state in [0, upperbound-1]"""
        if (upperBound == 0):
            return 0
        if (self.currentPermutation == []):
            self.currentPermutation = list(range(lowerBound,upperBound)) #Only use last 30% of states
            random.shuffle(self.currentPermutation)
        return self.currentPermutation.pop()

    def next(self, upperBound, exclude = []):
        """Gets a random state in [0, upperbound-1] that is not included in exclude-list"""
        outgoing_state = self._next(int(0.8*upperBound), upperBound)

        counter = 0
        while(outgoing_state in exclude):
            # Has to succeed at some point because num_incoming_actions <= state-1
            lowerBound = max(0,int((0.8-(counter//10)*0.1)*upperBound))
            outgoing_state = self._next(lowerBound, upperBound)
            counter+=1
        return outgoing_state