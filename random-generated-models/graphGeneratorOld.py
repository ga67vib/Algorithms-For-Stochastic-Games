import random

from seaborn import distributions
from graphGenParams import GraphGenerationParameters
from typing import Dict, List
class GeneratedGraph:
    def generateGraph(
        self,
        params : GraphGenerationParameters
        ):
        """
        Generates a Graph

        Respects parameters:
            num_states
            force_unknown
            denominator_range
            probability_to_branch
            probability_for_backwards_actions
            minimum_outgoing_edges

        Respects roughly:
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
        self.actions_map : Dict[int, List[Dict[int, str]]] = dict() #Each state has a List of dicts, where each dict is mapping of states to probabilities 

        self.params : GraphGenerationParameters = params

        # If several graphs are created with the same permutators, we might run into out of bounds if we don't reset them
        self.params.choice_permutator.reset_permutation()
        self.params.transition_permutator.reset_permutation()

    def _computeMaxActionsPerPlayer(self):
        #Compute maximal number of actions used
        self.max_player_1_actions = 0
        self.max_player_2_actions = 0

        counter = 0

        for state in self.states_of_player1:

            if (counter % 100000 == 0 and self.params.verbose >= 1):
                print(f"Count the maximum number of actions. Processed: {counter} / {self.params.num_states }")
            counter+=1

            self.max_player_1_actions = max(len(self.actions_map[state]), self.max_player_1_actions)

        for state in self.states_of_player2:

            if (counter % 100000 == 0 and self.params.verbose >= 1):
                print(f"Count the maximum number of actions. Processed: {counter} / {self.params.num_states }")
            counter+=1

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

            if (state % 100000 == 0 and self.params.verbose >= 1):
                print(f"Ensure Deadlock Freedom for State {state} / {self.params.num_states }")

            if (len(self.actions_map[state]) == 0):
                distribution = dict()
                distribution[state] = self._probabilityToFractionString(1,1)
                self.actions_map[state].append(distribution)

    def _generateChoices(self):
        for state in range(self.params.num_states):
            if (state % 100000 == 0 and self.params.verbose >= 1):
                print(f"Generate Choices for State {state} / {self.params.num_states }")
            # Random-number generator
            choice_targets_permutation = self.params.choice_permutator

            # Cannot have more actions leading into a state than overall available states
            num_incoming_actions = random.randint(min(state, self.params.minimum_incoming_edges), min(state, self.params.maximum_incoming_edges))
            if num_incoming_actions > 0:
                num_incoming_actions = random.randint(1, num_incoming_actions)
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

    def _generateChoice(self, source_state, target_state, probability_to_branch) -> Dict[int, str]:
        """ 
        Generates a Choice from source_state to target_state with probability_to_branch.
        Returns a mapping State -> Probability to reach as fraction 
        """
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

        # Create a transition to the target state to ensure that this action can reach it
        probability = random.randint(1, denominator)
        intDistribution[target_state] = probability
        totalProbability = probability

        # Decide if you want to have more transitions in this action leading to random states
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

    def __str__(self):
        s = ""
        for state in range(self.params.num_states):
            s+=(f'State {state}:\n')
            index = 0
            for action in self.actions_map[state]:
                s+=(f'\tAction {index}:\n')
                for key in action:
                    s+=(f'\t\t{key} : {action[key]}\n')
                index += 1
        return s