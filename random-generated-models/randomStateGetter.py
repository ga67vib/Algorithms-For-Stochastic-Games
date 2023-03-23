import random

class Permutation:
    """
    Implements a Backwards Permutation
    """

    def __init__(self, fast_transitions = False) -> None:
        self.fast_transitions = fast_transitions
        self.reset_permutation()

    def reset_permutation(self):
        self.currentPermutation = []

    def _next_fast(self, lowerBound, upperBound, exclude = []):
        """Gets a random state in [lowerBound, upperbound-1]"""
        if (upperBound == 0):
            return 0
        if (self.currentPermutation == []):
            self.currentPermutation = list(range(lowerBound,upperBound))
            random.shuffle(self.currentPermutation)
        return self.currentPermutation.pop()

    def _next(self, lowerBound, upperBound, exclude = []):
        """Gets a random state in [lowerBound, upperbound-1]"""
        if (upperBound == 0):
            return 0

        pick_from = list(set(range(lowerBound, upperBound)) - set(exclude))
        picked_index = random.randint(0, len(pick_from)-1)
        return pick_from[picked_index]

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

class Permutation_AllStatesPossible(Permutation):
    """
    Every state is eligible as next state
    """
    def __init__(self, fast_transitions = False) -> None:
        super().__init__(fast_transitions=fast_transitions)

    def next(self, upperBound, exclude = []):
        """Gets a random state in [0, upperbound-1] that is not included in exclude-list"""
        if (self.fast_transitions):
            outgoing_state = self._next_fast(0, upperBound, exclude=exclude)
        else:    
            outgoing_state = self._next(0, upperBound, exclude=exclude)
        return outgoing_state