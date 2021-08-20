import random

class Permutation:
    """
    Implements a Backwards Permutation
    """
    def __init__(self) -> None:
        self.currentPermutation = []

    def _next(self, lowerBound, upperBound):
        """Gets a random state in [lowerBound, upperbound-1]"""
        if (upperBound == 0):
            return 0
        if (self.currentPermutation == []):
            self.currentPermutation = list(range(lowerBound,upperBound))
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

class Permutation_AllStatesPossible(Permutation):
    """
    Every state is eligible as next state
    """
    def __init__(self) -> None:
        super().__init__()

    def next(self, upperBound, exclude = []):
        """Gets a random state in [0, upperbound-1] that is not included in exclude-list"""
        outgoing_state = self._next(0, upperBound)
        lowerBound = 0

        counter = 0
        while(outgoing_state in exclude):
            # Has to succeed at some point because num_incoming_actions <= state-1
            outgoing_state = self._next(lowerBound, upperBound)
            counter+=1
        return outgoing_state