from typing import List, Dict
from abc import ABC, abstractmethod


class FileReader(ABC):
    @abstractmethod
    def read(self, file_path: str):
        pass


class FASTAReader(FileReader):
    def __init__(self) -> None:
        self.reader_selector = {
            'dict': self._return_dict,
            'list': self._return_list
        }

    def read(self, file_path: str, output_type: str = 'list'):
        reader = self.reader_selector.get(output_type)
        if reader:
            return reader(file_path)

    @staticmethod
    def _return_dict(file_path: str) -> Dict:
        output = {}
        with open(file_path, 'r') as file:
            for pair in file.read().strip().split(">")[1:]:
                data = pair.split()
                output[data[0]] = "".join(data[1:])

        return output

    @staticmethod
    def _return_list(file_path: str) -> List:
        output = []
        with open(file_path, 'r') as file:
            for pair in file.read().strip().split(">")[1:]:
                data = pair.split()
                output.append((data[0], "".join(data[1:])))

        return output
