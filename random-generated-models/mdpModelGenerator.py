import argparse
from genericpath import isdir
import os, os.path
import glob
from graphGenerator import GeneratedGraph
from sccGraphGenerator import SCCGraphGenerator
from treeGraphGenerator import TreeGraphGenerator
from graphGenParams import GraphGenerationParameters
from randomStateGetter import *
import modelGeneratorHelper as helper

model_types = ["tree", "random", "scc"]

def main():
    output_dir = "generatedMDPs"

    #Default parameters
    num_min_actions_default = 2
    model_type_default = "random"
    smallest_transition_probability_default = 0.5
    backwards_probability_default = 0.5
    branching_probability_default = 0.8
    force_unknown_default = False

    parser = argparse.ArgumentParser(description = 'Generate PRISM models randomly')
    parser.add_argument(
        '-outputDir', type=str, default=output_dir, help='Where should the models be saved? Default: '+output_dir
    )
    parser.add_argument(
        '-size', type=int, help='How many states should the model have'
    )
    parser.add_argument(
        '-numModels', type=int, default=1, help='How many models should be created at once'
    )
    parser.add_argument(
        '-numMinActions', type=int, default=num_min_actions_default, help='How many actions should (almost) all states have at least'
    )
    parser.add_argument(
        '-type', type=helper.checkModelType, default=model_type_default, help="Which model template should be used? Available: "+str(model_types)
    )
    parser.add_argument(
        '-smallestProb', type=float, default = smallest_transition_probability_default, help = 'What is the smallest transition probability allowed in the graph?'
    )
    parser.add_argument(
        '-backwardsProb', type = float, default = backwards_probability_default, help=argparse.SUPPRESS#"probability to have backwards actions"
    )
    parser.add_argument(
        '-branchingProb', type = float, default = branching_probability_default, help="probability to add a branch in an action. There can be maximal 10 transitions per action"
    )
    parser.add_argument(
        '-forceUnknown', type = bool, default = force_unknown_default, help='Try to minimize the number of yes- and no-states. This can lead to every state having a small probability of reaching the target.'
    )

    arguments = parser.parse_args()
    output_dir = arguments.outputDir
    num_states = arguments.size
    num_models = arguments.numModels
    num_min_actions = arguments.numMinActions
    model_type_name = arguments.type
    smallest_transition_probability = arguments.smallestProb
    backwards_probability = arguments.backwardsProb
    branching_probability = arguments.branchingProb
    force_unknown = arguments.forceUnknown
    verbose = num_states >= 100000

    choice_permutator=Permutation_AllStatesPossible()
    transition_permutator=Permutation_AllStatesPossible()
    graph_model = GeneratedGraph()
    if (model_type_name == "tree"):
        graph_model = TreeGraphGenerator()
        choice_permutator=Permutation_AllStatesPossible()
        transition_permutator=Permutation_AllStatesPossible()
    elif (model_type_name == "scc"):
        graph_model = SCCGraphGenerator()
        choice_permutator=Permutation_AllStatesPossible()
        transition_permutator=Permutation_AllStatesPossible()

    # Constants
    EMPTY_LINE = "\n\n"

    # Parameters
    max_num_actions_player1 = 3
    max_num_actions_player2 = 3

    num_states_before_iters = num_states #num_states could be overriden to add extra sinks/targets. Need to reset after every modelIteration

    for modelNumber in range(num_models):
        print("Generating Model Number: "+(str(modelNumber+1))+"...")

        num_states = num_states_before_iters

        # Generate Graph
        graph = graph_model
        generation_parameters = GraphGenerationParameters(
                num_states,
                minimum_incoming_edges=num_min_actions,
                maximum_incoming_edges=num_min_actions + 2,
                probability_to_branch=branching_probability,
                probability_for_backwards_action=backwards_probability,
                probability_to_be_maximizer_state=1,
                minimum_outgoing_edges=num_min_actions,
                choice_permutator=choice_permutator,
                transition_permutator=transition_permutator,
                denominator_range=int(1/smallest_transition_probability),
                force_unknown=force_unknown
            )

        graph.generateGraph(
            generation_parameters
        )

        # Could have changed to subvert trivial targets
        num_states = generation_parameters.num_states

        if verbose:
            print("The Graph is generated. Next we need to create a PRISM File from it...")

        max_num_actions_player1 = graph.max_player_1_actions
        max_num_actions_player2 = graph.max_player_2_actions

        output_file_name = "MDP_RANDOM_Size_"+str(num_states)
        if (num_min_actions != num_min_actions_default):
            output_file_name += "_MinAct_"+str(num_min_actions)
        if (model_type_name != model_type_default):
            output_file_name += "_Type-"+model_type_name
        if (smallest_transition_probability != smallest_transition_probability_default):
            output_file_name += "_MinTransProb_"+str(smallest_transition_probability)
        output_file_name += "_Model_"
        output_file_name_index = len(glob.glob1(os.path.join(".",output_dir),output_file_name+"*.prism"))
        output_file_name+=str(output_file_name_index+1)+".prism"

        file_string = ""

        header = "mdp"
        file_string+=header+EMPTY_LINE

        player1_actions_header ="player P1 player1"
        player2_actions_header ="player P2 player2"

        player1_actions_header+=helper.addActionsToActionsHeader(1, max_num_actions_player1)+" endplayer"
        player2_actions_header+=helper.addActionsToActionsHeader(2, max_num_actions_player2)+" endplayer"
        #file_string+=player1_actions_header+"\n"
        #file_string+=player2_actions_header+EMPTY_LINE


        if verbose:
            print("Header is generated. Now we create the player modules...")

        # For efficiency, build a map from state to player
        state_to_player = dict()
        for state in graph.states_of_player1:
            state_to_player[state] = 1
        for state in graph.states_of_player2:
            state_to_player[state] = 2

        # Modules
        player_1_actions = ""
        player_2_actions = ""
        for state in graph.actions_map:
            if verbose and state % 100000 == 0:
                print(f'Processing actions of state {state} / {len(graph.actions_map)}')

            player_identifier = ""
            player_actions = ""

            # We extend either the actions of player 1 or player 2, depending on who the state belongs to
            if (state_to_player[state] == 1):
                player_identifier = "a"
            else:
                player_identifier = "b"

            action_index = 0
            for action in graph.actions_map[state]:
                action_name = "["+player_identifier+str(action_index+1)+"]"
                if (len(graph.actions_map[state]) == 1):
                    action_name = "[]"
                
                action_string = "\t"+action_name+" state="+str(state)+" -> "
                transition_string = ""
                for target_state in graph.actions_map[state][action_index]:
                    transition_string+=graph.actions_map[state][action_index][target_state]+" : (state'="+str(target_state)+")"
                    transition_string+=" + "
                transition_string=transition_string[0:len(transition_string)-3] #Remove last ' + '
                action_string+=transition_string+";\n"

                player_actions+=action_string

                action_index+=1
            
            # Extend the actions of the corresponding player
            if (state_to_player[state] == 1):
                player_1_actions += player_actions
            else:
                player_2_actions += player_actions


        if (verbose):
            print("Modules done! Next create the labels...")

        file_string+="module player1\n"
        states_declaration ="state: [0.."+str(num_states-1)+"];"
        file_string+=states_declaration+EMPTY_LINE
        file_string+=player_1_actions+EMPTY_LINE
        file_string+="endmodule"+EMPTY_LINE

        #file_string+="module player2\n"
        #file_string+=player_2_actions
        #file_string+="endmodule"+EMPTY_LINE

        # Labels
        targetStates = [num_states-1]

        file_string+="//Labels\n"
        label = 'label "p1win" = '
        winCondition = ""
        for state in targetStates:
            winCondition+="state="+str(state)+" & "
        winCondition = winCondition[0:len(winCondition)-3] # Cut out last " & "
        label+=winCondition+";"

        file_string+=label

        if verbose:
            print("PRISM Model is Generated. Saving it...")

        # Write to file
        if (not os.path.isdir(output_dir)):
            os.mkdir(output_dir)
        file = open(os.path.join(output_dir,output_file_name), "w")
        file.write(file_string)
        file.close()

if __name__ == "__main__":
    main()