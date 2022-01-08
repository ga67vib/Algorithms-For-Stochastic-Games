from typing import Any, Callable, Dict


class FileDataFilter:
    """
    Filter class that is supposed to be handed to the file_loader to decide which data points to filter out.
    The file_loader should then save the FileDataFilter in the returned FileData to propagate which models have been filtered out for what reason

    filter_function: function which receives parameters: 
        row_to_analyse : a row of the csv file
        reason_for_filtering : a dict where you can store how many data points were filtered for which reason

        return: False if model should be filtered out, else True
    """
    def __init__(self, filter_function) -> None:
        self.reason_for_filtering : Dict[str, int] = {}
        self.filtered_models = []
        self.filter_function : Callable[[Any, Dict[str, int]], bool] = filter_function

    def filter(self, row):
        """ 
        Applies the filter_function to the row handed in
        
        return: False is the model should be filtered out according to the filter_function, else True
        """

        value = self.filter_function(row, self.reason_for_filtering)
