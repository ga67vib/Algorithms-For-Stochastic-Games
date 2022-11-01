from re import I
from typing import Any, Callable, List, Tuple, Union
import csv

if __name__ != '__main__':
    from .data_entities.file_data_filter import FileDataFilter
    from .data_entities.file_data import FileData
    from .file_loader_utils import is_float, is_property_defined, transpose_array_of_dicts

else:
    from data_entities.file_data_filter import FileDataFilter
    from data_entities.file_data import FileData
    from file_loader_utils import is_float, is_property_defined, transpose_array_of_dicts

def load_file(property_keys : List[str], paths_to_csv_files : Union[List[str], str], primary_key_name : str = "Model", filter_function : Callable[..., bool] = None) -> FileData:
    file_data = None

    if (primary_key_name not in property_keys):
        raise Exception(f"{primary_key_name} is required to be in property_keys {property_keys}.")

    if not isinstance(paths_to_csv_files, list):
        path_to_csv_file = paths_to_csv_files
        file_data = _load_one_file(property_keys=property_keys, path_to_csv_file=path_to_csv_file, primary_key_name=primary_key_name, filter_function=filter_function)
    else:
        loaded_file_data_list = []
        for path_to_csv_file in paths_to_csv_files:
            loaded_file_data = _load_one_file(property_keys=property_keys, path_to_csv_file=path_to_csv_file, primary_key_name=primary_key_name, filter_function=filter_function)
            loaded_file_data_list.append(loaded_file_data)
        
        file_data : FileData = loaded_file_data_list.pop()
        for loaded_file_data in loaded_file_data_list:
            file_data.aggregate_data(loaded_file_data)

    return file_data

def _load_one_file(property_keys : List[str], path_to_csv_file : str, primary_key_name : str = "Model", filter_function : Callable[..., bool] = None) -> FileData:
    table = []
    num_filtered_models = 0
    total_num_of_models = 0
    reason_for_filtering_model = {}
    primary_key_to_index = dict()

    with open(path_to_csv_file, newline='\n') as input_file:
        reader = csv.DictReader(input_file)
        for row in reader:
            # Get the relevant properties and their values from this row as key-value pair 
            conf_name = row["Conf"]
            
            relevant_key_value_pairs = dict()
            for key in property_keys:
                x = key
                if (key != primary_key_name and key != "Model"):
                    # Prepend the conf name to clearly assign the values
                    x = conf_name+"_"+key
                relevant_key_value_pairs[x] = row[key]


            if filter_function != None and not filter_function(relevant_key_value_pairs, reason_for_filtering_model):
                num_filtered_models += 1
                continue

            pk_of_row = relevant_key_value_pairs[primary_key_name]
            if (pk_of_row not in primary_key_to_index.keys()):
                primary_key_to_index[pk_of_row] = len(primary_key_to_index.keys())

            # Prepend the conf to every property
            for key in relevant_key_value_pairs.keys():
                # Every key but the primary_key should be numeric and receives a default value in case it is not
                if (not is_float(relevant_key_value_pairs[key]) and key != primary_key_name and key != "Conf"):
                    relevant_key_value_pairs[key] = -100

            # Check if not already created
            model_already_created = list(filter(lambda x : x[primary_key_name] == relevant_key_value_pairs[primary_key_name], table))
            if (len(model_already_created) > 0):
                index = table.index(model_already_created[0])
                table[index] = {**table[index], **relevant_key_value_pairs}
            else:
                table.append(relevant_key_value_pairs)
                total_num_of_models += 1

    loaded_data = transpose_array_of_dicts(table)

    print(f"{total_num_of_models - num_filtered_models} / {total_num_of_models} Model from Path {path_to_csv_file} remain after filtering")
    if reason_for_filtering_model:
        print("Reasons for filtering:")
        for key, value in reason_for_filtering_model.items():
            print(f'\t{key}: {value}')

    return FileData(primary_key_name=primary_key_name, data=loaded_data)