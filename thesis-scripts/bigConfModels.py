import os
from os.path import exists
import sys
import subprocess
import statistics
import subprocess
from shlex import split
from collections import namedtuple
from functools import reduce
import multiprocess as mp
import typing

#TODO: epsilon analysis (but only for best confs)


# If executed in "run" mode,
# this runs prism located at prism_path with the configurations indicated in the dict configurations
# on all the models in the dict models.
# The output is written to the current directory. For every configuration, a folder is created and for every model a file.
# E.g. BVI_0/mdsm1.log
# Also, you can repeat the experiments several times if you want, by setting reps to an int greater than 1.

# If executed in "read" mode,
# this reads the files created by running the benchmarks and creates three csv files: One for the results (values.csv), one for the time taken (times.csv) and one for the iterations (iters.csv)

# If executed in "analyse" mode,
# this analyses all the models enabled and prints a report

# Some general parameters
#prism_path="../../qp/code/prism-games/prism/bin/prism" #Path to PRISM
prism_path="../prism-games-3.0.beta-src/prism/bin/prism"
wp_path="../../CAV20Impl/mycode/WP/bin/prism"
max_processes = 1
TIMEOUT = "30m" # Was 30min! Increased only for T_SI_SI
JAVAMAXMEM = "30g" # "Free: 36g Processes: 6 -> 5g MaxMem and 1g JavaStack"
JAVASTACK = "5g"
reps=1 #Repetitions. If set to 1, it will not appear in filename of log.    
output_dir="bigConfModels"
random_input_dir = "../random-generated-models/moreTrees/"
proc_output = namedtuple('proc_output', 'stdout stderr') #Is needed for the pipeline (i suppose)

#Had to add in these extra commands since subrocess.output seems to make trouble with "|"
#Found on Stackoverflow
#https://stackoverflow.com/questions/24306205/file-not-found-error-when-launching-a-subprocess-containing-piped-commands
def pipeline(starter_command, *commands):
    if not commands:
        try:
            starter_command, *commands = starter_command.split('|')
        except AttributeError:
            pass
    starter_command = _parse(starter_command)
    starter = subprocess.Popen(starter_command, stdout=subprocess.PIPE)
    last_proc = reduce(_create_pipe, map(_parse, commands), starter)
    s = str(proc_output(*last_proc.communicate())[0])
    s = s[2:-3] #usually string has form b'...string-content...\n' => I cut b' and \n'
    return s

def _create_pipe(previous, command):
    proc = subprocess.Popen(command, stdin=previous.stdout, stdout=subprocess.PIPE)
    previous.stdout.close()
    return proc

def _parse(cmd):
    try:
        return split(cmd)
    except Exception:
        return cmd


# Configurations
configurations = dict()

#configurations["BVI_100"] = (prism_path, "-ii -maxiters 100")

#Only normal VI + TOP_VI
#OVI without D
#OVI_opt_1
#WP

#Classis VI
configurations["VI"] = (prism_path, "")
configurations["TOP_VI"] = (prism_path, "-topological")

#BVI
configurations["BVI"] = (prism_path, "-ii -maxiters 1")
configurations["G_BVI"] = (prism_path, "-ii -maxiters 1 -smg_opts 1")
configurations["T_BVI"] = (prism_path, "-ii -maxiters 1 -topological -smg_opts 2")
configurations["TOP_BVI"] = (prism_path, "-ii -maxiters 1 -topological")
configurations["D_BVI"] = (prism_path, "-ii -maxiters 100")

#SVI
configurations["SVI"] = (prism_path, "-svi -maxiters 1")
configurations["G_SVI"] = (prism_path, "-svi -maxiters 1 -smg_opts 1")
configurations["T_SVI"] = (prism_path, "-svi -maxiters 1 -topological -smg_opts 2")
#configurations["TOP_SVI"] = (prism_path, "-svi -maxiters 1 -topological")
configurations["D_SVI"] = (prism_path, "-svi -maxiters 100")

#OVI
configurations["OVI"] = (prism_path, "-ovi -maxiters 1")
configurations["G_OVI"] = (prism_path, "-ovi -maxiters 1 -smg_opts 1")
configurations["T_OVI"] = (prism_path, "-ovi -maxiters 1 -topological -smg_opts 2")
#configurations["TOP_OVI"] = (prism_path, "-ovi -maxiters 1 -topological")
configurations["OPT_OVI"] = (prism_path, "-ovi -maxiters 1 -smg_opts 4")
configurations["WP_OVI"] = (prism_path, "-ovi -maxiters 1 -smg_opts 20")

#WP
configurations["WP"] = (wp_path, "-ex -BVI_A")
configurations["WP_P3"] = (prism_path, "-wp -maxiters 1")

#SI
configurations["SI"] =  (prism_path, "-politer -smg_opts 0 -maxiters 1000000000") # Normal Strategy Iteration - SI Opponent, i.e. MDP Solver is SI
configurations["SI_SI"] =  (prism_path, "-politer -smg_opts 9 -maxiters 1000000000") # Normal Strategy Iteration - SI Opponent, i.e. MDP Solver is SI
configurations["T_SI_SI"] =  (prism_path, "-politer -smg_opts 10 -maxiters 1000000000") # Normal Strategy Iteration - SI Opponent, i.e. MDP Solver is S
configurations["LP_SI"] =  (prism_path, "-politer -smg_opts 13") # Normal Strategy Iteration - SI Opponent, i.e. MDP Solver is S
configurations["T_LP_SI"] =  (prism_path, "-politer -smg_opts 14") # Normal Strategy Iteration - SI Opponent, i.e. MDP Solver is S





#Models
models=dict()

#models["AV18_18_2"]="../case-studies/AV-n.prism ../case-studies/AV.props -const X_MAX=18,Y_MAX=18 -prop 2"
#models["AV20_20_2"]="../case-studies/AV-n.prism ../case-studies/AV.props -const X_MAX=20,Y_MAX=20 -prop 2"
#models["AV25_25_2"]="../case-studies/AV-n.prism ../case-studies/AV.props -const X_MAX=25,Y_MAX=25 -prop 2"
#models["AV30_30_2"]="../case-studies/AV-n.prism ../case-studies/AV.props -const X_MAX=30,Y_MAX=30 -prop 2"
#models["dice200"]="../case-studies/dice-N.prism ../case-studies/dice.props -const N=200 -prop 1"
#models["dice250"]="../case-studies/dice-N.prism ../case-studies/dice.props -const N=250 -prop 1"
#models["dice300"]="../case-studies/dice-N.prism ../case-studies/dice.props -const N=300 -prop 1"
#models["dice400"]="../case-studies/dice-N.prism ../case-studies/dice.props -const N=400 -prop 1"
#models["hallway18_18_2"]="../case-studies/hallway-n.prism ../case-studies/hallway.props -const X_MAX=18,Y_MAX=18 -prop 2"
#models["hallway20_20_2"]="../case-studies/hallway-n.prism ../case-studies/hallway.props -const X_MAX=20,Y_MAX=20 -prop 2"
#models["hallway22_22_2"]="../case-studies/hallway-n.prism ../case-studies/hallway.props -const X_MAX=22,Y_MAX=22 -prop 2"
#models["hallway30_30_2"]="../case-studies/hallway-n.prism ../case-studies/hallway.props -const X_MAX=30,Y_MAX=30 -prop 2"
#models["hallway40_40_2"]="../case-studies/hallway-n.prism ../case-studies/hallway.props -const X_MAX=40,Y_MAX=40 -prop 2"
#models["BigMec_1e6"] = "../case-studies/BigMec.prism ../case-studies/BigMec.props -const N=1000000"
#models["BigMec_1e7"] = "../case-studies/BigMec.prism ../case-studies/BigMec.props -const N=10000000"
#models["BigMec_1e8"] = "../case-studies/BigMec.prism ../case-studies/BigMec.props -const N=100000000"
#models["ManyMECs_1e6"]="../case-studies/ManyMecs.prism ../case-studies/ManyMecs.props -const N=1000000"
#models["ManyMECs_1e7"]="../case-studies/ManyMecs.prism ../case-studies/ManyMecs.props -const N=10000000"
#models["ManyMECs_1e8"]="../case-studies/ManyMecs.prism ../case-studies/ManyMecs.props -const N=100000000"



# Handcrafted Models
###models["dice50MEC"]="../case-studies/dice50MEC.prism ../case-studies/dice.props -prop 1"
###models["cdmsnMEC"]="../case-studies/cdmsnMEC.prism ../case-studies/cdmsn.props"
#models["ManyMECs_1e1"] = "../case-studies/ManyMecs.prism ../case-studies/ManyMecs.props -const N=10"
#models["ManyMECs_1e2"]="../case-studies/ManyMecs.prism ../case-studies/ManyMecs.props -const N=100"
#models["ManyMECs_1e3"]="../case-studies/ManyMecs.prism ../case-studies/ManyMecs.props -const N=1000"
### Fix if you have time models["ManyMECs_1e4"]="../case-studies/ManyMecs.prism ../case-studies/ManyMecs.props -const N=10000"
#models["BigMec_1e1"] = "../case-studies/BigMec.prism ../case-studies/BigMec.props -const N=10"
#models["BigMec_1e2"] = "../case-studies/BigMec.prism ../case-studies/BigMec.props -const N=100"
#models["BigMec_1e3"] = "../case-studies/BigMec.prism ../case-studies/BigMec.props -const N=1000"
#models["BigMec_1e4"] = "../case-studies/BigMec.prism ../case-studies/BigMec.props -const N=10000"
#models["hm_10"]="../case-studies/haddad-monmege-SG.pm ../case-studies/haddad-monmege.prctl -const N=30,p=0.5"
#models["hm_30"]="../case-studies/haddad-monmege-SG.pm ../case-studies/haddad-monmege.prctl -const N=30,p=0.5"
#models["hm_40"]="../case-studies/haddad-monmege-SG.pm ../case-studies/haddad-monmege.prctl -const N=30,p=0.5"


# Random Models
"""
random_model_files_dir = random_input_dir
for random_model_file in os.listdir(random_model_files_dir):
    if random_model_file.startswith("RANDOM") and random_model_file.endswith(".prism"):
        model_name = random_model_file.replace(".prism", "")
        models[model_name] = random_model_files_dir+random_model_file+" ../case-studies/randomModels.props -prop 1"

"""

# Model Extensions
models["simple"]="../case-studies/SimpleModel.prism ../case-studies/randomModels.props -prop 1" # Base Model that is to be extended
# Currently, apply the extension to each Model
extension_config_folder_path = "../model-extension-configs/sccConfs"
extension_config_models = dict()
for config_file in os.listdir(extension_config_folder_path):
    if (not config_file.endswith(".json")):
        continue
    config_name = config_file.replace(".json", "")
    for model_key in models.keys():
        new_key = model_key+"_"+config_name
        extension_config_models[new_key] = models[model_key]+" -smg_extend "+os.path.join(extension_config_folder_path,config_file)
#Unify models with the extensions
models = {**models, **extension_config_models}


# Parse command line to decide whether to run benchmarks or read them
if len(sys.argv) == 0 or str(sys.argv[1]) not in ["run", "read", "analyse", "ovi_iters"]:
    print("This script can only run in the following modes: run, read, analyse or ovi_iters. Call it with one of these as command line parameter")
elif sys.argv[1] == "run":
    command_list = []
    for conf_count, conf in enumerate(sorted(configurations.keys())):
        #print(conf)
        os.system("mkdir -p " + output_dir + "/" + conf)
        for model_count, model in enumerate(sorted(models.keys())):
            counting_str = ("Conf: %s [%d/%d], Model: [%d/%d] - " %(conf, conf_count + 1, len(configurations), model_count + 1, len(models)))
            counting_str = "\t"+counting_str + model
            #print("\t"+counting_str + model)
            for i in range(1, reps+1):
                #print("\t\t"+str(i))
                rep_string = "" if reps == 1 else "_rep" + str(i)
                if exists(output_dir + "/" + conf + "/" + model + rep_string + ".log"):
                    print("\t\tAlready there, skipping")
                    continue
                if (conf != "WP"):
                    prismParams = "-javamaxmem "+JAVAMAXMEM+" -javastack "+JAVASTACK+" " # "-javamaxmem 32g -javastack 16g"  # Change this appropriately
                else:
                    prismParams = "-javamaxmem "+JAVAMAXMEM+" "
                command = "timeout "+TIMEOUT+" " + configurations[conf][0] + " " + \
                    models[model] + " " + configurations[conf][1] + " " + prismParams + \
                    " > " + output_dir + "/" + conf + "/" + model + rep_string + ".log" + \
                    " 2> " + output_dir + "/" + conf + "/" + model + rep_string + ".errors.log"
                command_list.append((command, counting_str))
    process_pool = mp.Pool(processes=max_processes)
    def useCommand(command):
        print(command[1])
        try:
            os.system(command[0])
        except:
            e = sys.exc_info()[0]
            print(e)
    process_pool.map(useCommand, command_list)

elif (sys.argv[1] == "read"):
    # Model, #States, [min/mean/max runtime for each solver]
    with open(output_dir+"/times.csv", "w") as timefile:
        with open(output_dir+"/values.csv", "w") as valuefile:
            with open(output_dir+"/iters.csv", "w") as iterfile:
                timefile.write("Model,#States")
                valuefile.write("Model,#States")
                iterfile.write("Model,#States")
                # Get this here to make sure it has same order for header and rows
                confs = sorted(configurations.keys())
                for conf in confs:
                    timefile.write("," + conf)
                    valuefile.write("," + conf)
                    iterfile.write("," + conf)
                timefile.write("\n")
                valuefile.write("\n")
                iterfile.write("\n")

                for model in sorted(models.keys()):
                    # First print model name and #states. Use first conf and first rep to get #states
                    infile = output_dir + "/" + \
                        confs[0] + "/" + model + \
                        ("" if reps == 1 else "_rep1") + ".log"
                    s1 = "grep 'States' " + infile
                    s2 = "cut -d ' ' -f 7"
                    states = pipeline(s1, s2)
                    timefile.write(str(model) + "," + str(states))
                    valuefile.write(str(model) + "," + str(states))
                    iterfile.write(str(model) + "," + str(states))

                    for conf in confs:
                        if reps == 1:
                            infile = output_dir + "/" + conf + "/" + model + ".log"
                            s1 = "grep 'Time for model checking:' " + infile
                            s2 = "cut -d ' ' -f 5"
                            # s1 = "grep 'Probabilistic reachability took' " + infile
                            # s2 = "cut -d ' ' -f 4"
                            resTime = pipeline(s1,s2)
                            s1 = "grep 'Result:' " + infile
                            s2 = "cut -d ' ' -f 2"
                            resSol = pipeline(s1, s2)
                            s1 = "grep 'Value iteration variant' " + infile
                            s2 = "cut -d ' ' -f 6"
                            resIter = pipeline(s1, s2)


                        else:
                            times = []
                            for i in range(1, reps+1):
                                infile = output_dir + "/" + conf + "/" + \
                                    model + "_rep" + str(i) + ".log"
                                s1 = "grep 'Time for model checking:' " + infile
                                s2 = "cut -d ' ' -f 5"
                                time = pipeline(s1,s2)

                                s1 = "grep 'Result:' " + infile
                                s2 = "cut -d ' ' -f 2"
                                sol = pipeline(s1, s2)

                                s1 = "grep 'Value iteration variant' " + infile
                                s2 = "cut -d ' ' -f 5"
                                iter = pipeline(s1, s2)

                                try:
                                    time = int(time)
                                    sol = float(sol)
                                    iter = int(iter)
                                except (ValueError):
                                    res = "X"
                                    sol = "X"
                                    iter = "X"
                                    break
                                times += [time]
                                sols += [sol]
                                iters += [iter]
                            resTime = min(times) + "/" + statistics.mean(times) + "/" + max(times)
                            resSol = min(sols) + "/" + statistics.mean(sols) + "/" + max(sols)
                            resIter = min(iters) + "/" + statistics.mean(iters) + "/" + max(iters)

                        timefile.write("," + str(resTime))
                        valuefile.write("," + str(resSol))
                        iterfile.write("," + str(resIter))

                    timefile.write("\n")
                    valuefile.write("\n")
                    iterfile.write("\n")

elif (sys.argv[1] == "analyse"):
    conf_params = (prism_path, "-analyse -nopre")
    conf_name = "ANALYSIS"
    print(conf_name)
    os.system("mkdir -p " + output_dir + "/" + conf_name)
    
    #RULES FOR FEATURES: 
    #The key is the name of the column in the resulting .csv file
    #The value contain exactly the string as it is printed in the java-file to avoid errors
    relevantFeatures = dict()

    #Basics
    relevantFeatures["NumStates"] = "Number of States: "
    relevantFeatures["NumActions"] = "Number of Choices: "
    relevantFeatures["NumTargets"] = "Number of Targets (States with trivial value 1): "
    relevantFeatures["NumSinks"] = "Number of Sinks (States with trivial value 0): "
    relevantFeatures["NumUnknown"] = "Number of Unknown States: "
    
    #Actions-related
    relevantFeatures["NumMaxActions"] = "Number of maximal choices per state: "
    relevantFeatures["NumMaxTransitions"] = "Number of maximal transitions per choice: "
    relevantFeatures["SmallestTransProb"] = "Smallest transition probability: "
    relevantFeatures["NumProbActions"] = "Number of Choices with probability: "
    relevantFeatures["AvgNumActionsPerState"] = "Average Number of Choices per state: "
    relevantFeatures["AvgNumTransPerState"] = "Average Number of Transitions per state: "
    relevantFeatures["AvgNumTransPerAction"] = "Average Number of Transitions per choice: "
    relevantFeatures["NumMaxStates"] = "Number of Maximizer States: "
    relevantFeatures["NumMinStates"] = "Number of Minimizer States: "
    relevantFeatures["NumBackwardsTransitions"] = "Percentage of Backwards Transitions: "

    #MEC-related
    relevantFeatures["NumMECs"] = "Number of MECs: "
    relevantFeatures["BiggestMEC"] = "Biggest MEC has size: "
    relevantFeatures["SmallestMEC"] = "Smallest MEC has size: "
    relevantFeatures["AvgMEC"] = "MEC size on average is: "
    relevantFeatures["MedianMEC"] = "MEC size median is: "

    #SCC-related
    relevantFeatures["NumSCCs"] = "Number of SCCs: "
    relevantFeatures["BiggestSCC"] = "Biggest SCC has size: "
    relevantFeatures["SmallestSCC"] = "Smallest SCC has size: "
    relevantFeatures["AvgSCC"] = "Average SCC has size: "
    relevantFeatures["MaxSCCDepth"] = "Longest Chain of SCC has length: "
    #Non-Singleton SCCs
    relevantFeatures["NumNonSingleton"] = "Number of non-Singleton SCCs: "
    relevantFeatures["SmallestSCCNonSing"] = "Smallest non-Singleton SCC has size: "
    relevantFeatures["AvgSccNonSing"] = "Average non-Singleton SCC has size: "

    #Cyclefree-Checks
    relevantFeatures["NearestTarget"] = "Nearest Target from any Initial State: "
    relevantFeatures["FurthestTarget"] = "Furthest Target from any Initial State: "
    relevantFeatures["TargetDistanceAverage"] = "Target-distance Average: "
    relevantFeatures["TargetDistanceMedian"] = "Target-distance Median: "

    #Times
    relevantFeatures["Prob0InSec"] = "Prob0 Time in s: "
    relevantFeatures["Prob1InSec"] = "Prob1 Time in s: "
    relevantFeatures["MECCompInSec"] = "Time to compute MECs (s): "
    relevantFeatures["SCCCompInSec"] = "SCC computation took (s): "


    #Run
    conf_list = []
    def runAnalyse(conf):
        print(conf[1])
        try:
            os.system(conf[0])
        except:
            e = sys.exc_info()[0]
            print(e)

    for model_count, model in enumerate(sorted(models.keys())):
        counting_str = "Model: [%d/%d] - " % (model_count + 1, len(models))
        counting_str = "\t"+counting_str+model
        #print("\t"+counting_str+model)
        for i in range(1, reps+1):
            #print("\t\t"+str(i))
            rep_string = "" if reps == 1 else "_rep" + str(i)
            if exists(output_dir + "/" + conf_name + "/" + model + rep_string + ".log"):
                print("\t\tAlready there, skipping")
                continue
            prismParams = "-javamaxmem "+JAVAMAXMEM+" -javastack "+JAVASTACK+" " # "-javamaxmem 32g -javastack 16g"  # Change this appropriately
            command = f'timeout {TIMEOUT} ' + conf_params[0] + " " + \
                models[model] + " " + conf_params[1] + " " + prismParams + \
                " > " + output_dir + "/" + conf_name + "/" + model + rep_string + ".log" + \
                " 2> " + output_dir + "/" + conf_name + "/" + model + rep_string + ".errors.log"
            conf_list.append([command, counting_str])
    process_pool = mp.Pool(processes=max_processes)
    process_pool.map(runAnalyse, conf_list)

    #Read
    with open(output_dir+"/analysis.csv", "w") as statisticsfile:
        header = "Model"
        for feature in relevantFeatures.keys():
            header += ","+feature
        statisticsfile.write(header+"\n")

        for model in sorted(models.keys()):

            infile = output_dir + "/" + conf_name + "/" + model + ".log"

            #Write name of Model
            statisticsfile.write(str(model))

            for feature in relevantFeatures.keys():
                #Currently don't look for reps
                value = pipeline("grep '"+relevantFeatures[feature]+"' "+infile)
                value = value.replace(relevantFeatures[feature], "")
                statisticsfile.write(", "+value)

            #Print newline for next Model
            statisticsfile.write("\n")

elif (sys.argv[1] == "ovi_iters"):
    #Read
    headers = ["Conf", "Model", "Verification Time", "Number Of Verif. Phases", "Iterations in Verif Phase", "NumAborted", "NumInduceLowerBound", "NumFailedVerif"]

    with open(output_dir+"/ovi_iters.csv", "w") as statisticsfile:
        header = ""
        for h in headers:
            header += h + ","
        header = header[0:-1]
        statisticsfile.write(header+"\n")

        for conf_name in sorted(configurations.keys()):
            if ("OVI" not in conf_name):
                continue

            for model in sorted(models.keys()):
                infile_name = output_dir + "/" + conf_name + "/" + model + ".log"

                num_verifications = 0
                iters_in_verif_phase = "-"
                verif_time = "-"
                num_abort = 0
                num_boost = 0
                num_fails = 0

                with open(infile_name, "r") as infile:
                    for line in infile:
                        if ("Starting a verification phase in iteration" in line):
                            num_verifications+=1
                            modified_line = line.replace("Starting a verification phase in iteration ", "")

                        elif ("Time for model checking:" in line):
                            modified_line = line.replace("Time for model checking: ", "")
                            modified_line = modified_line.replace(" seconds.", "")
                            verif_time = int(float(modified_line))

                        elif (" iterations in the verification phase." in line):
                            modified_line = line.replace("Proved U to be inductive upper bound in iteration", "")
                            modified_line = modified_line.replace("iterations in the verification phase.", "")
                            # There are two numbers left, split them
                            split_iters = modified_line.split(" after ")
                            # Iters in verification phase is the second of those two numbers
                            iters_in_verif_phase = int(split_iters[1])
                            num_fails -= 1

                        elif ("U is inductive lower bound. Using it as new L, stopping verification phase." in line):
                            num_boost+=1

                        elif ("Aborting verification phase after " in line):
                            num_abort+=1

                num_fails += num_verifications
                num_fails -= num_boost
                
                #Write name of Model
                statisticsfile.write(str(conf_name)+","+str(model)+",")
                statisticsfile.write(str(verif_time)+",")
                statisticsfile.write(str(num_verifications)+",")
                statisticsfile.write(str(iters_in_verif_phase)+",")
                statisticsfile.write(str(num_abort)+",")
                statisticsfile.write(str(num_boost)+",")
                statisticsfile.write(str(num_fails)+"\n")