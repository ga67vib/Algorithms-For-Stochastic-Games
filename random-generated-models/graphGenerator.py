import random

from seaborn import distributions
from graphGenParams import GraphGenerationParameters
class GeneratedGraph:
    def generateGraph(
        self,
        params : GraphGenerationParameters
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

        self._initVars(params)

        self._generateStates()

        self._generateChoices()

        if (params.force_unknown):
            self._reduceTrivialStates()

        self._ensureDeadlockFreedom()

        if (params.force_unknown):
            self._reduceNoStates()

        self._computeMaxActionsPerPlayer()

    def _initVars(self, params : GraphGenerationParameters):
        self.states_of_player1 = []
        self.states_of_player2 = []
        self.actions_map = dict() #Each state has a List of dicts, where each dict is mapping of states to probabilities 

        self.params : GraphGenerationParameters = params

    def _computeMaxActionsPerPlayer(self):
        #Compute maximal number of actions used
        self.max_player_1_actions = 0
        self.max_player_2_actions = 0
        for state in self.actions_map:
            if (state in self.states_of_player1):
                self.max_player_1_actions = max(len(self.actions_map[state]),self.max_player_1_actions)
            else:
                self.max_player_2_actions = max(len(self.actions_map[state]),self.max_player_2_actions)

    def _reduceTrivialStates(self):
        """ Add after the old target a new target with 90% chance to go to a sink """
        newSink = self.params.num_states
        newTarget = self.params.num_states+1
        oldTarget = self.params.num_states-1

        self.params.num_states = self.params.num_states + 2
        self.states_of_player1.append(newTarget)
        self.states_of_player2.append(newSink)

        self.actions_map[newSink] = []
        self.actions_map[newTarget] = []

        distribution = dict()
        distribution[newTarget] = self._probabilityToFractionString(9,10)
        distribution[newSink] = self._probabilityToFractionString(1,10)

        self.actions_map[oldTarget].append(distribution)

    def _reduceNoStates(self):
        pass

    def _generateStates(self):
        for state in range(self.params.num_states):
            # Assign State to a player
            if self._meetsThreshhold(self.params.probability_to_be_maximizer_state):
                self.states_of_player1.append(state)
            else:
                self.states_of_player2.append(state)
            self.actions_map[state] = []
    
    def _ensureDeadlockFreedom(self):
        #Every state that has no action should have self-loop
        for state in range(self.params.num_states):
            if (len(self.actions_map[state]) == 0):
                distribution = dict()
                distribution[state] = self._probabilityToFractionString(1,1)
                self.actions_map[state].append(distribution)

    def _generateChoices(self):
        for state in range(self.params.num_states):
            # Connect backwards to states
            choice_targets_permutation = self.params.choice_permutator

            # Cannot have more actions leading into a state than overall available states
            num_incoming_actions = random.randint(min(state, self.params.minimum_incoming_edges), min(state, self.params.maximum_incoming_edges))
            already_included_states = []
            for _ in range(num_incoming_actions):
                outgoing_state = choice_targets_permutation.next(state, exclude = already_included_states)

                #Add choice
                already_included_states.append(outgoing_state)
                #Try to force MECs by setting branching probability low
                distribution = self._generateChoice(outgoing_state, state, 0.1)
                self.actions_map[outgoing_state].append(distribution)

        #Same procedure backwards to add MECs:
        count_down = list(range(self.params.num_states-1))
        count_down.reverse()
        for state in count_down:
            if (self._meetsThreshhold(self.params.probability_for_backwards_action)):
                # TODO Currently per state atmost one backwards action
                target_state = random.randint(0, state)
                distribution = self._generateChoice(state, target_state, self.params.probability_to_branch)
                self.actions_map[state].append(distribution)

    def _generateChoice(self, source_state, target_state, probability_to_branch):
        distribution = dict()
        if (self._meetsThreshhold(probability_to_branch)):
            distribution = self._getProbabilityDistributionForBranchingAction(target_state)
        else:
            distribution[target_state] = self._probabilityToFractionString(1,1)
        return distribution
                

    def _getProbabilityDistributionForBranchingAction(self, target_state):
        counter = 1
        transition_targets_permutation = self.params.transition_permutator
        distribution = dict()
        intDistribution = dict()

        denominator = self.params.denominator_range
        maximal_branch_count = 10

        #Assign a sure branch to the target state
        probability = random.randint(1, denominator)
        intDistribution[target_state] = probability
        totalProbability = probability

        while (counter < maximal_branch_count and self._meetsThreshhold(self.params.probability_to_branch) and totalProbability < denominator):
            other_state = transition_targets_permutation.next(target_state)
            probability = random.randint(1, denominator-totalProbability)
            if not other_state in intDistribution.keys():
                intDistribution[other_state] = probability
            else:
                intDistribution[other_state] += probability

            counter+=1
            totalProbability += probability

        #If probabilities don't add up to one, increase target probability
        if totalProbability != denominator:
            intDistribution[target_state] += denominator-totalProbability

        if (len(intDistribution) == 1):
            distribution[target_state] = self._probabilityToFractionString(1,1) 
        else:
            for state in intDistribution.keys():
                distribution[state] = self._probabilityToFractionString(intDistribution[state], denominator)
        return distribution

    def _meetsThreshhold(self, probability):
        return float(random.randint(1, 100))/100.0 <= probability

    def _probabilityToFractionString(self, numerator, denominator):
        if (numerator == denominator):
            return "1"
        else:
            return str(numerator)+"/"+str(denominator)