import os
import json

template_path = "../model-extension-configs/config-template/configTemplate.json"
out_dir = "../model-extension-configs/oviBad/"

if (not os.path.isdir(out_dir)):
    os.mkdir(out_dir)

config_json = ""
with open(template_path, 'r') as config_template:
    config_json = json.loads(config_template.read())

mec_model = config_json["MECModel"]
prob_model = config_json["ProbabilisticModel"]
action_model = config_json["ActionTreeModel"]

prob_model["use"] = False
mec_model["use"] = False
action_model["use"] = True

if (mec_model["use"]):
    #3
    for sinkProb in [1, 40, 99]:
        #3
        for initProb in [1, 40, 99]:
            if (sinkProb+initProb > 100):
                continue
            #3
            for c_switchProb in [0, 40, 99]:
                #3
                for numMEC in [40]:#[1, 10, 100]:
                    #3
                    for c_len in [1, 200, 1000]:#[1, 200, 1000]:
                        file_name = "mec_sinkP_%d_initP_%d_cSwitch_%d_numMec_%d_cLen_%d.json" % (
                            sinkProb, initProb, c_switchProb, numMEC, c_len
                        )

                        mec_model["numMECs"] = numMEC
                        mec_model["chainLength"] = c_len
                        mec_model["probabilityLeadingToSink"] = float(sinkProb)/100.0
                        mec_model["probabilityToGoBackToInitialState"] = float(initProb)/100.0
                        mec_model["probabilityForChainSwitch"] = float(c_switchProb)/100.0
                        with open(out_dir+file_name, 'w') as out_file:
                            out_file.write(json.dumps(config_json))

if (prob_model["use"]):
    branch_to_num_comp = dict()
    branch_to_depth = dict()

    branching_factors = [2, 5, 10, 20, 400]

    #Each one gives 400 States
    branch_to_num_comp[2] = 25
    branch_to_num_comp[5] = 16
    branch_to_num_comp[10] = 4
    branch_to_num_comp[20] = 1
    branch_to_num_comp[400] = 1

    branch_to_depth[2] = 4
    branch_to_depth[5] = 2
    branch_to_depth[10] = 2
    branch_to_depth[20] = 2
    branch_to_depth[400] = 1

    #3
    for sinkProb in [1]:
        #3
        for initProb in [0, 1]:
            if (sinkProb+initProb > 100):
                continue
            #3
            for branching_factor in branching_factors:
                #3
                for treeDepth in [branch_to_depth[branching_factor]]:
                    #3
                    for numComponents in [1, 5, 20]: # Gives state sizes 400, 2000, 10000
                        file_name = "probTree_sinkP_%d_initP_%d_branch_%d_depth_%d_numComp_%d.json" % (
                            sinkProb, initProb, branching_factor, treeDepth, numComponents * branch_to_num_comp[branching_factor]
                        )

                        prob_model["numComponents"] = numComponents * branch_to_num_comp[branching_factor]
                        prob_model["componentBranchingFactor"] = branching_factor
                        prob_model["probabilityLeadingToSink"] = float(sinkProb)/100.0
                        prob_model["probabilityToGoBackToInitialState"] = float(initProb)/100.0
                        prob_model["componentTreeDepth"] = treeDepth
                        with open(out_dir+file_name, 'w') as out_file:
                            out_file.write(json.dumps(config_json))

if (action_model["use"]):
    branch_to_num_comp = dict()
    branch_to_depth = dict()

    branching_factors = [1]

    #3
    for sinkProb in [0]:
        #3
        for initProb in [99]:
            if (sinkProb+initProb > 100):
                continue
            #3
            for branching_factor in branching_factors:
                #3
                for numStates in [1]:
                    #3
                    for numComponents in [1, 5, 10, 50, 100, 200]: # Gives state sizes 400, 2000, 10000
                        file_name = "actionTree_sinkP_%d_initP_%d_branch_%d_numStates_%d_numComp_%d.json" % (
                            sinkProb, initProb, branching_factor, numStates, numComponents
                        )

                        action_model["numComponents"] = numComponents
                        action_model["componentBranchingFactor"] = branching_factor
                        action_model["probabilityLeadingToSink"] = float(sinkProb)/100.0
                        action_model["probabilityToGoBackToInitialState"] = float(initProb)/100.0
                        action_model["componentNumberOfStates"] = numStates
                        with open(out_dir+file_name, 'w') as out_file:
                            out_file.write(json.dumps(config_json))