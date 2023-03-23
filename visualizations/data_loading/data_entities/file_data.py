from dataclasses import dataclass
from typing import Dict, List


@dataclass
class FileData:
    """
    Data storage object for data derived from one type of file.
    Consists of:
    
    data : A Dict that stores to every key a list. An entry in the list should represent the data of for example the value of the key for a model
    primary_key_name : On the values of this key the data will be joined together in the Data Unification Process. Every entry in the list data[primary_key_name] must be unique
    """
    primary_key_name : str
    data : Dict[str, List]


    def aggregate_data(self, file_data_to_aggregate : 'FileData') -> None:
        self_keys = self.data.keys()
        other_keys = file_data_to_aggregate.data.keys()

        if (self_keys != other_keys):
            print(f"WARNING: Expected that file_data_to_aggregate would have keys {self_keys}, but has keys {other_keys}")

        for key in file_data_to_aggregate:
            if (key in self.data):
                self.data[key] = self.data[key] + file_data_to_aggregate[key]

    # Implement some quality-of-life features 
    def __getitem__(self, k : str):
        ''' Index Data '''
        return self.data.__getitem__(k)

    def __iter__(self):
        ''' Iterate over the Data '''
        return iter(self.data)