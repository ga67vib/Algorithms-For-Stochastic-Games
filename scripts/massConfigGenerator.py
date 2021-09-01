import os
import json

template_path = "../model-extension-configs/config-template/configTemplate.json"
out_dir = "../model-extension-configs/generated/"

if (not os.path.isdir(out_dir)):
    os.mkdir(out_dir)

config_json = ""
with open(template_path, 'r') as config_template:
    config_json = json.loads(config_template.read())

mec_model = config_json["MECModel"]
prob_model = config_json["ProbabilisticModel"]

prob_model["use"] = False
mec_model["use"] = True

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
