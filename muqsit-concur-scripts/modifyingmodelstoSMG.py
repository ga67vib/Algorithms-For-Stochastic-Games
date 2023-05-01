import os
import re
import sys

def extract_module_names(file_path):
    pattern = re.compile(r'^module\s+(\w+)')
    result_set = set()

    with open(file_path, 'r') as file:
        for line in file:
            matches = pattern.findall(line.strip())

            for match in matches:
                result_set.add(match)

    return result_set

def extract_square_brackets(file_path):
    pattern = re.compile(r'(^\[.+?\])')
    result_set = set()

    with open(file_path, 'r') as file:
        for line in file:
            matches = pattern.findall(line.strip())
            
            for match in matches:
                result_set.add(match)

    return result_set

def convert_mdp_to_smg(file_path):
    module_names = extract_module_names(file_path)
    action_names = extract_square_brackets(file_path)

    new_file_path = os.path.splitext(file_path)[0] + "-SG.prism"

    with open(file_path, 'r') as old_file, open(new_file_path, 'w') as new_file:
        for line in old_file:
            if line.strip().startswith("mdp"):
                new_file.write("smg\n")
                new_file.write("player P1\n")
                new_file.write(", ".join(module_names))
                if action_names:
                    new_file.write(", " + ", ".join(action_names))
                new_file.write("\nendplayer\n")
            else:
                new_file.write(line)

    print(f"The new file is saved at: {new_file_path}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Please provide a file path as an argument.")
        sys.exit(1)

    file_path = sys.argv[1]
    convert_mdp_to_smg(file_path)

