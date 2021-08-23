import os

models_per_conf = 50
types = ["tree"]
state_sizes = [100, 1000, 10000]
num_min_actions = [1, 5, 20]
branchingProbs = [0.8]#[0.1, 0.8]
backwardsProbs = [0.8]
smallestProbs = [0.01]



for type in types:
    for branchingProb in branchingProbs:
        for backwardsProb in backwardsProbs:
            for num_min_action in num_min_actions:
                for state_size in state_sizes:
                    for smallestProb in smallestProbs:
                        os.system("python modelGenerator -size %d -numModels %d -numMinActions %d -type %s -smallestProb %f -branchingProb %f" %
                        (state_size, models_per_conf, num_min_action, type, smallestProb, branchingProb))