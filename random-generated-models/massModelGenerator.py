import os
import multiprocess as mp

output_dir = "manyAct_steps"
models_per_conf = 20
types = ["tree"]
state_sizes = [100, 1000, 10000]
num_min_actions = [2, 5, 10, 30, 50, 70, 90, 100]
branchingProbs = [0.9]#[0.1, 0.8]
backwardsProbs = [0]
smallestProbs = [0.1]
forceUnknown = True

num_processes = 7
total_iterations = len(types) * len(state_sizes) * len(num_min_actions) * len(branchingProbs) * len(backwardsProbs) * len(smallestProbs)

def generateBatch(command):
    print(command[1])
    os.system(command[0])

command_list = []
count = 1
process_pool = mp.Pool(num_processes)
for type in types:
    for branchingProb in branchingProbs:
        for backwardsProb in backwardsProbs:
            for num_min_action in num_min_actions:
                for state_size in state_sizes:
                    for smallestProb in smallestProbs:
                        s = "Iteration %d/%d" % (count, total_iterations)
                        command = ("python modelGenerator.py -size %d -numModels %d -numMinActions %d -type %s -smallestProb %f -branchingProb %f -outputDir %s -backwardsProb %f -forceUnknown %r" 
                        % (state_size, models_per_conf, num_min_action, type, smallestProb, branchingProb, output_dir, backwardsProb, forceUnknown)
                        )
                        command_list.append([command, s])
                        count += 1

process_pool.map(generateBatch, command_list)