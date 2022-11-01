import os
import multiprocess as mp
import time

output_dir = "RandomGenBenchmark"
models_per_conf = 1
state_sizes = [10000, 100000, 1000000]
# 2 3 5
# 100 for random
num_min_actions = [1, 10, 100] # In Tree: Branching Factor - In Random: Connectivity Factor (How many Actions lead certainly into state)
branchingProbs = [1] # How probable are multiple transitions per Action? 
backwardsProbs = [1] # How probable are MEC-creating Actions?
smallestProbs = [0.0001] # What's the smallest allowed probability?
fast = [True]

total_iterations = len(state_sizes) * len(num_min_actions) * len(branchingProbs) * len(backwardsProbs) * len(smallestProbs)

def generateBatch(command):
    print(command[1])
    os.system(command[0])

command_list = []
count = 1
for branchingProb in branchingProbs:
    for backwardsProb in backwardsProbs:
        for state_size in state_sizes:
            for num_min_action in num_min_actions:
                for smallestProb in smallestProbs:
                    for fastParam in fast:
                        s = "Iteration %d/%d" % (count, total_iterations)
                        t1 = time.time()
                        command = (f'python modelGenerator.py -size {state_size} -numMinActions {num_min_action} -smallestProb {smallestProb} -branchingProb {branchingProb} -outputDir {output_dir} -backwardsProb {backwardsProb}')
                        if fastParam:
                            command += ' --fastTransitions'
                        os.system(command)
                        t2 = time.time()
                        print(f'[NumStates, NumActions, Fast]: {[state_size, num_min_action, fastParam]}\nSeconds Passed: {t2-t1}')
                        count += 1
