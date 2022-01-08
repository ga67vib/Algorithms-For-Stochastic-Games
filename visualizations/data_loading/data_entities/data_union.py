from dataclasses import dataclass
from typing import Any, Dict, List
import copy

if __name__ != '__main__':
    from .file_data import FileData

else:
    from file_data import FileData

@dataclass
class DataUnionFileInputObject:
    """
    Input Object for DataUnion method "unify"
    """
    file_data_object : FileData
    file_data_prefix: str

class DataUnion:
    """
    Unifies File_Data from different file Types.
    Everything but the primary key values and the key-names is supposed to be 
    """

    _data_union_file_input : List[DataUnionFileInputObject]
    _key_to_input_obj : Dict[str, DataUnionFileInputObject]
    data : Dict[str, List]
    primary_key_values: List[str]
    primary_key : str

    def __init__(self, data_union_file_input : List[DataUnionFileInputObject]) -> None:
        self._data_union_file_input = data_union_file_input
        self._key_to_input_obj = dict()
        self._unify(data_union_file_input=data_union_file_input.copy())

    def _unify(self, data_union_file_input : List[DataUnionFileInputObject]):
        if (len(data_union_file_input) < 1):
            raise Exception("Expect at least one input object")

        self._data_union_file_input = data_union_file_input

        first_input = data_union_file_input.pop(0)
        first_file = first_input.file_data_object
        self.data = dict()
        self.primary_key = first_file.primary_key_name

        # Make deep copies of the loaded data
        prefix = first_input.file_data_prefix
        for key in first_file.data.keys():
            if (key == self.primary_key):
                self.primary_key_values = copy.deepcopy(first_file.data[key])
            else:
                self.data[prefix+key] = copy.deepcopy(first_file[key])
                self._key_to_input_obj[prefix+key] = first_input

        for file_data_input in data_union_file_input:
            # Check if we need to prepend the file_data_prefix
            keys_contained_in_both = [key for key in self.data.keys() if key in file_data_input.file_data_object.data.keys()]
            if (len(keys_contained_in_both) > 0):
                raise Exception(f'When adding {file_data_input.file_data_prefix} to the Data, {keys_contained_in_both} are already provided keys. To not override these, please provide another prefix')
            prefix = file_data_input.file_data_prefix

            # Add keys to the data union
            file_data_keys : List = list(file_data_input.file_data_object.data.keys())

            # Do not include the primary key
            file_data_keys.remove(file_data_input.file_data_object.primary_key_name)

            for key in file_data_keys:
                self.data[prefix+key] = [] 
                self._key_to_input_obj[prefix+key] = file_data_input
            
            # Create an index mapping to join correctly and add file data to data union
            # This is an inner join - only if the key is present in both the data is kept
            index_mapping = dict()
            file_data = file_data_input.file_data_object
            for i, value in enumerate(self.primary_key_values):
                index_mapping[i] = self._indexOf(file_data.data[file_data_input.file_data_object.primary_key_name], value)
                if (index_mapping[i] == -1):
                    continue

                for key in file_data_keys:
                    self.data[prefix+key].append(file_data.data[key][index_mapping[i]])

    def _indexOf(self, list : List, element):
        """Return index of element, else -1"""
        index = -1
        try:
            index = list.index(element)
        except ValueError:
            index = -1
        return index

    def __str__(self) -> str:
        s = 'DataUnion Object'
        s+= f'\nPrimary Key: {self.primary_key}'
        s+= f'\nOther Keys: {self.data.keys()}'
        return(s)

    def replace_values_for_input_object(self, input_object, placeholder, replacement):
        for key in self._key_to_input_obj.keys():
            if (self._key_to_input_obj[key] != input_object):
                continue
            for i in range(len(self.data[key])):
                if (self.data[key][i] == placeholder):
                    self.data[key][i] = replacement

#Test Code:
x_dict = dict()
x_dict['a'] = [1,2]
x_dict['b'] = [2,3]
x = FileData('a', x_dict)
x_input = DataUnionFileInputObject(x, '')

y_dict = dict()
y_dict['c'] = [1,2]
y_dict['d'] = [4,5]
y = FileData('c', y_dict)
y_input = DataUnionFileInputObject(y, '')

z = DataUnion([x_input, y_input])
print(z)
print(z.data)

z.replace_values_for_input_object(y_input, 4, 6)
print(z.data)