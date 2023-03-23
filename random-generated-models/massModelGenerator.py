import os
import multiprocess as mp
import time

output_dir = "RandomGenBenchmark"
models_per_conf = 2
state_sizes = [10000, 100000, 1000000, 1000000]
# 2 3 5
# 100 for random
num_min_actions = [1, 10, 100] # In Tree: Branching Factor - In Random: Connectivity Factor (How many Actions lead certainly into state)
branchingProbs = [1] # How probable are multiple transitions per Action? 
backwardsProbs = [1] # How probable are MEC-creating Actions?
smallestProbs = [0.0001] # What's the smallest allowed probability?
fast = [True, False]

num_processes = 1
total_iterations = len(state_sizes) * len(num_min_actions) * len(branchingProbs) * len(backwardsProbs) * len(smallestProbs)

def generateBatch(command):
    print(command[1])
    os.system(command[0])

command_list = []
count = 1
process_pool = mp.Pool(num_processes)
for branchingProb in branchingProbs:
    for backwardsProb in backwardsProbs:
        for num_min_action in num_min_actions:
            for state_size in state_sizes:
                for smallestProb in smallestProbs:
                    s = "Iteration %d/%d" % (count, total_iterations)
                    t1 = time.time()
                    command = ("python modelGenerator.py -size %d -numModels %d -numMinActions %d -smallestProb %f -branchingProb %f -outputDir %s -backwardsProb %f -fastTransitions %r" 
                    % (state_size, models_per_conf, num_min_action, smallestProb, branchingProb, output_dir, backwardsProb, fast)
                    )
                    command_list.append([command, s])
                    count += 1

process_pool.map(generateBatch, command_list)