import argparse
import os, os.path
import glob
from graphGenerator import GeneratedGraph

output_dir = "generatedModels"

parser = argparse.ArgumentParser(description = 'Generate PRISM models randomly')
parser.add_argument(
    '-size', type=int, help='How many states should the model have'
)

arguments = parser.parse_args()

# Constants
EMPTY_LINE = "\n\n"

# Parameters
num_states = arguments.size
max_num_actions_player1 = 3
max_num_actions_player2 = 3

# Generate Graph
graph = GeneratedGraph()
graph.generateGraph(
    num_states,
    minimum_incoming_edges=1,
    maximum_incoming_edges=3,
    probability_to_branch=0.80,
    probability_for_backwards_action=0.2,
    probability_to_be_maximizer_state=0.60
)

max_num_actions_player1 = graph.max_player_1_actions
max_num_actions_player2 = graph.max_player_2_actions

output_file_name = "SIZE_"+str(num_states)+"_MODEL_"
output_file_name_index = len(glob.glob1(os.path.join(".",output_dir),output_file_name+"*.prism"))
output_file_name+=str(output_file_name_index+1)+".prism"

file_string = ""

header = "smg"
file_string+=header+EMPTY_LINE

states_declaration ="global state: [0.."+str(num_states-1)+"];"
file_string+=states_declaration+EMPTY_LINE

player1_actions_header ="player P1 player1"
player2_actions_header ="player P2 player2"

def addActionsToActionsHeader(player, num_actions):
    if (num_actions <= 1):
        return ""

    player_identifier = ""
    if (player==1):
        player_identifier = "a"
    else:
        player_identifier = "b"
    actions = ""

    for action in range(1, num_actions+1):
        actions+=", ["+player_identifier+str(action)+"]"
    return actions

player1_actions_header+=addActionsToActionsHeader(1, max_num_actions_player1)+" endplayer"
player2_actions_header+=addActionsToActionsHeader(2, max_num_actions_player2)+" endplayer"
file_string+=player1_actions_header+"\n"
file_string+=player2_actions_header+EMPTY_LINE

# Modules
player_1_actions = ""
player_2_actions = ""
for state in graph.actions_map:
    player_identifier = ""
    player_actions = ""
    if (state in graph.states_of_player1):
        player_identifier = "a"
        player_actions = player_1_actions
    else:
        player_identifier = "b"
        player_actions = player_2_actions

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
        if (state in graph.states_of_player1):
            player_1_actions = player_actions
        else:
            player_2_actions = player_actions
        action_index+=1

file_string+="module player1\n"
file_string+=player_1_actions
file_string+="endmodule"+EMPTY_LINE

file_string+="module player2\n"
file_string+=player_2_actions
file_string+="endmodule"+EMPTY_LINE

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


# Write to file
file = open(os.path.join(output_dir,output_file_name), "w")
file.write(file_string)
file.close()