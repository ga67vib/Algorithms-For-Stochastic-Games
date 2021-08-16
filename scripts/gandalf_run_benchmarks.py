import os
from os.path import exists
import sys
import subprocess
import statistics
import subprocess
from shlex import split
from collections import namedtuple
from functools import reduce

# If executed in "run" mode,
# this runs prism located at prism_path with the configurations indicated in the dict configurations
# on all the models in the dict models.
# The output is written to the current directory. For every configuration, a folder is created and for every model a file.
# E.g. BVI_0/mdsm1.log
# Also, you can repeat the experiments several times if you want, by setting reps to an int greater than 1.

# If executed in "read" mode,
# this reads the files created by running the benchmarks and creates a csv file for it.

# Some general parameters
#prism_path="../../qp/code/prism-games/prism/bin/prism" #Path to PRISM
prism_path="../prism-games-3.0.beta-src/prism/bin/prism"
reps=1 #Repetitions. If set to 1, it will not appear in filename of log.
output_dir="output"
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


#Configurations
configurations=dict()

# Configurations
configurations = dict()
configurations["BVI_1"] = "-ii -maxiters 1"
configurations["BVI_100"] = "-ii -maxiters 100"
configurations["TBVI_1"] = "-ii -maxiters 1 -topological"
configurations["TBVI_100"] = "-ii -maxiters 100 -topological"
# QP
configurations["QP_G"] = "-qp -smg_opts 0"    	#QP improved approach
configurations["QP_G_W"] = "-qp -smg_opts 1"    	#Give Gurobi additional start-vector from VI
configurations["QP_C"] = "-qp -smg_opts 3"     	#C = Cplex Model Created (by Gurobi)
configurations["QP_CON"] = "-qp -smg_opts 4"     	#Condon's original QP
# configurations["QP_ECM"] = "-qp -smg_opts 5"     	#Solve ECs with Condon's m-chain
configurations["QP_ECE"] = "-qp -smg_opts 6"     	#Solve ECs with Epsilon-transitions
configurations["QP_G_T"] = "-qp -smg_opts 7"     	#Gurobi, Topological
configurations["QP_G_T_W"] = "-qp -smg_opts 8"    	#Gurobi, Topological, Warm Start with VI
configurations["AMPL_M"] = "-qp -smg_opts 9"  	#AMPL Minos
# configurations["AMPL_W"] = "-qp -smg_opts 8"     	#AMPL Minos, Warm start (not implemented)
configurations["AMPL_T"] = "-qp -smg_opts 11"     	#AMPL With topological sort
# SI
configurations["SI"] = "-politer -smg_opts 0" # Normal Strategy Iteration - SI Opponent, i.e. MDP Solver is SI
configurations["SI_PR"] = "-politer -smg_opts 1" # Value Iteration Precomputation, 20 steps
configurations["SI_T"] = "-politer -smg_opts 2" # MEC decomposition
configurations["SI_PR_T"] = "-politer -smg_opts 3" # Value Iteration Precomputation & MEC decomposition
configurations["SI_BV"] = "-politer -smg_opts 4" # BVI Opponent, i.e. MDP Solver is BVI (Parameters may need adjusting)
configurations["SI_PR_BV"] = "-politer -smg_opts 5" # Value Iteration Precomputation & BVI Opponent (Parameters may need adjusting)
configurations["SI_T_BV"] = "-politer -smg_opts 6" # MEC decomposition & BVI Opponent (Parameters may need adjusting)
configurations["SI_PR_T_BV"] = "-politer -smg_opts 7" # Value Iteration Precomputation & MEC decomposition & BVI Opponent (Parameters may need adjusting)
configurations["SI_VI"] = "-politer -smg_opts 8" # VI Opponent, i.e. MDP Solver is VI
configurations["SI_PR_VI"] = "-politer -smg_opts 9" # Value Iteration Precomputation & VI Opponent
configurations["SI_T_VI"] = "-politer -smg_opts 10" # MEC decomposition & VI Opponent
configurations["SI_PR_T_VI"] = "-politer -smg_opts 11" # Value Iteration Precomputation & MEC decomposition & VI Opponent
configurations["SI_PR_200"] = "-politer -smg_opts 17" # Value Iteration Precomputation, 200 steps
configurations["SI_PR_BV_200"] = "-politer -smg_opts 21" # Value Iteration Precomputation (200 steps) & BVI Opponent (Parameters may need adjusting)
configurations["SI_PR_VI_200"] = "-politer -smg_opts 25" # Value Iteration Precomputation (200 steps) & VI Opponent
configurations["SI_PR_20K"] = "-politer -smg_opts 33" # Value Iteration Precomputation, 20 000 steps
configurations["SI_PR_BV_20K"] = "-politer -smg_opts 37" # Value Iteration Precomputation (20 000 steps) & BVI Opponent (Parameters may need adjusting)
configurations["SI_PR_VI_20K"] = "-politer -smg_opts 41" # Value Iteration Precomputation (20 000 steps) & VI Opponent
configurations["VI"] = ""

#Models
models=dict()
models["cdmsn"]="../case-studies/cdmsn.prism ../case-studies/cdmsn.props"
#models["cloud5"]="../case-studies/cloud_5.prism ../case-studies/cloud.props"
#models["cloud6"]="../case-studies/cloud_6.prism ../case-studies/cloud.props"
#models["mdsm1"]="../case-studies/mdsm.prism ../case-studies/mdsm.props -prop 1"
#models["mdsm2"]="../case-studies/mdsm.prism ../case-studies/mdsm.props -prop 2"
#models["teamform3"]="../case-studies/team-form-3.prism ../case-studies/team-form.props"
#models["teamform4"]="../case-studies/team-form-4.prism ../case-studies/team-form.props"
#models["AV10_10_1"]="../case-studies/AV10_10.prism ../case-studies/AV.props -prop 1"
#models["AV10_10_2"]="../case-studies/AV10_10.prism ../case-studies/AV.props -prop 2"
#models["AV10_10_3"]="../case-studies/AV10_10.prism ../case-studies/AV.props -prop 3"
#models["AV15_15_1"]="../case-studies/AV15_15.prism ../case-studies/AV.props -prop 1"
#models["AV15_15_2"]="../case-studies/AV15_15.prism ../case-studies/AV.props -prop 2"
#models["AV15_15_3"]="../case-studies/AV15_15.prism ../case-studies/AV.props -prop 3"
#models["charlton1"]="../case-studies/charlton.prism ../case-studies/charlton.props -prop 1"
#models["charlton2"]="../case-studies/charlton.prism ../case-studies/charlton.props -prop 2"
#models["dice10"]="../case-studies/dice10.prism ../case-studies/dice.props -prop 1"
#models["dice20"]="../case-studies/dice20.prism ../case-studies/dice.props -prop 1"
#models["dice50"]="../case-studies/dice50.prism ../case-studies/dice.props -prop 1"
#models["hallway5_5_1"]="../case-studies/hallway5_5.prism ../case-studies/hallway.props -prop 1"
#models["hallway5_5_2"]="../case-studies/hallway5_5.prism ../case-studies/hallway.props -prop 2"
#models["hallway8_8_1"]="../case-studies/hallway8_8.prism ../case-studies/hallway.props -prop 1"
#models["hallway8_8_2"]="../case-studies/hallway8_8.prism ../case-studies/hallway.props -prop 2"
#models["hallway10_10_1"]="../case-studies/hallway10_10.prism ../case-studies/hallway.props -prop 1"
#models["hallway10_10_2"]="../case-studies/hallway10_10.prism ../case-studies/hallway.props -prop 2"
#models["dice50MEC"]="../case-studies/dice50MEC.prism ../case-studies/dice.props -prop 1"
#models["cdmsnMEC"]="../case-studies/cdmsnMEC.prism ../case-studies/cdmsn.props"
#models["ManyMECs_1e1"] = "../case-studies/ManyMecs.prism ../case-studies/ManyMecs.props -const N=10"
#models["ManyMECs_1e2"]="../case-studies/ManyMecs.prism ../case-studies/ManyMecs.props -const N=100"
#models["ManyMECs_1e3"]="../case-studies/ManyMecs.prism ../case-studies/ManyMecs.props -const N=1000"
#models["ManyMECs_1e4"]="../case-studies/ManyMecs.prism ../case-studies/ManyMecs.props -const N=10000"
#models["BigMec_1e1"] = "../case-studies/BigMec.prism ../case-studies/BigMec.props -const N=10"
#models["BigMec_1e2"] = "../case-studies/BigMec.prism ../case-studies/BigMec.props -const N=100"
#models["BigMec_1e3"] = "../case-studies/BigMec.prism ../case-studies/BigMec.props -const N=1000"
#models["BigMec_1e4"] = "../case-studies/BigMec.prism ../case-studies/BigMec.props -const N=10000"
#models["hm_30"]="../case-studies/haddad-monmege-SG.pm ../case-studies/haddad-monmege.prctl -const N=30,p=0.5"
#models["hm_100"]="../case-studies/haddad-monmege-SG.pm ../case-studies/haddad-monmege.prctl -const N=100,p=0.5"
#models["hm_200"]="../case-studies/haddad-monmege-SG.pm ../case-studies/haddad-monmege.prctl -const N=200,p=0.5"
#models["adt"]="../case-studies/adt-infect.prism ../case-studies/adt-infect.props -prop 2"
#models["two_investors"]="../case-studies/two_investors.prism ../case-studies/two_investors.props -prop 4"
#models["coins"]="../case-studies/coins.prism ../case-studies/coins.props -prop 1"
#models["prison_dil"]="../case-studies/prisoners_dilemma.prism ../case-studies/prisoners_dilemma.props -prop 9"

# Parse command line to decide whether to run benchmarks or read them
if len(sys.argv) == 0 or str(sys.argv[1]) not in ["run", "read"]:
    print("This script can only run in two modes: run or read. Call it with one of these two as command line parameter")
elif sys.argv[1] == "run":
    for conf in sorted(configurations.keys()):
        print(conf)
        os.system("mkdir -p " + output_dir + "/" + conf)
        for model in sorted(models.keys()):
            print("\t"+model)
            for i in range(1, reps+1):
                print("\t\t"+str(i))
                rep_string = "" if reps == 1 else "_rep" + str(i)
                if exists(output_dir + "/" + conf + "/" + model + rep_string + ".log"):
                    print("\t\tAlready there, skipping")
                    continue
                prismParams = "-javamaxmem 32g -javastack 16g"  # Change this appropriately
                command = "timeout 15m " + prism_path + " " + \
                    models[model] + " " + configurations[conf] + " " + prismParams + \
                    " > " + output_dir + "/" + conf + "/" + model + rep_string + ".log"
                try:
                    os.system(command)
                except:
                    e = sys.exc_info()[0]
                    print(e)

else:  # sys.argv[0] == "read"
    # Model, #States, [min/mean/max runtime for each solver]
    with open(output_dir+"/times.csv", "w") as outfile:
        outfile.write("Model,#States")
        # Get this here to make sure it has same order for header and rows
        confs = sorted(configurations.keys())
        for conf in confs:
            outfile.write("," + conf)
        outfile.write("\n")

        for model in sorted(models.keys()):
            # First print model name and #states. Use first conf and first rep to get #states
            infile = output_dir + "/" + \
                confs[0] + "/" + model + \
                ("" if reps == 1 else "_rep1") + ".log"
            s1 = "grep 'States' " + infile
            s2 = "cut -d ' ' -f 7"
            states = pipeline(s1, s2)
            outfile.write(str(model) + "," + str(states))

            for conf in confs:
                if reps == 1:
                    infile = output_dir + "/" + conf + "/" + model + ".log"
                    s1 = "grep 'Time for model checking:' " + infile
                    s2 = "cut -d ' ' -f 5"
                    # s1 = "grep 'Probabilistic reachability took' " + infile
                    # s2 = "cut -d ' ' -f 4"
                    resTime = pipeline(s1,s2)
                else:
                    times = []
                    for i in range(1, reps+1):
                        infile = output_dir + "/" + conf + "/" + \
                            model + "_rep" + str(i) + ".log"
                        s1 = "grep 'Time for model checking:' " + infile
                        s2 = "cut -d ' ' -f 5"
                        time = pipeline(s1,s2)
                        try:
                            time = int(time)
                        except (ValueError):
                            res = "X"
                            sol = "X"
                            break
                        times += [time]
                    resTime = min(times) + "/" + \
                        statistics.mean(times) + "/" + max(times)
                outfile.write("," + str(resTime))
            outfile.write("\n")

    with open(output_dir+"/values.csv", "w") as outfile:
        outfile.write("Model,#States")
        # Get this here to make sure it has same order for header and rows
        confs = sorted(configurations.keys())
        for conf in confs:
            outfile.write("," + conf)
        outfile.write("\n")

        for model in sorted(models.keys()):
            # First print model name and #states. Use first conf and first rep to get #states
            infile = output_dir + "/" + \
                confs[0] + "/" + model + \
                ("" if reps == 1 else "_rep1") + ".log"
            s1 = "grep 'States' " + infile
            s2 = "cut -d ' ' -f 7"
            states = pipeline(s1, s2)
            outfile.write(str(model) + "," + str(states))

            for conf in confs:
                if reps == 1:
                    infile = output_dir + "/" + conf + "/" + model + ".log"
                    s1 = "grep 'Result:' " + infile
                    s2 = "cut -d ' ' -f 2"
                    resSol = pipeline(s1, s2)
                else:
                    sols = []
                    for i in range(1, reps+1):
                        infile = output_dir + "/" + conf + "/" + \
                            model + "_rep" + str(i) + ".log"
                        s1 = "grep 'Result:' " + infile
                        s2 = "cut -d ' ' -f 2"
                        sol = pipeline(s1, s2)
                        try:
                            sol = float(sol)
                        except (ValueError):
                            res = "X"
                            sol = "X"
                            break
                        sols += [sol]
                    resSol = min(sols) + "/" + \
                        statistics.mean(sols) + "/" + max(sols)
                outfile.write("," + str(resSol))
            outfile.write("\n")
