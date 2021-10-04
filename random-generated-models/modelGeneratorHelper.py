from argparse import ArgumentError
from modelGenerator import model_types

def checkModelType(value):
    if value not in model_types:
        raise ArgumentError("There is no graph template "+value+" but only "+str(model_types))
    return value

# Helper Functions
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