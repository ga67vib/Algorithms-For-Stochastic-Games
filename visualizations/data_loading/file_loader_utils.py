from typing import Any, Dict, List


def is_float(entry) -> bool:
    try:
        float(entry)
        return True
    except ValueError:
        return False

def is_property_defined(property : str, keys : List[str], row) -> bool:
    return (property in keys and is_float(row[property]))

def transpose_array_of_dicts(list_of_dicts : List[Dict[str, Any]]) -> Dict[str, List]:
    result_dict : Dict[str, List] = dict()
    headers = list_of_dicts[0].keys()

    for header in headers:
        result_dict[header] = []

    for dictionary in list_of_dicts:
        for key in dictionary:
            result_dict[key].append(dictionary[key])
    
    return result_dict