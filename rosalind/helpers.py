from typing import Dict, List, Tuple


def read_sequence(file_path: str) -> str:
    with open(file_path, 'r') as file:
        return file.readline().strip()


def read_several_sequences(file_path: str) -> List[str]:
    with open(file_path, 'r') as file:
        return file.read().strip().split("\n")


def fasta_to_dict(file_path: str) -> Dict:
    output = {}
    with open(file_path, 'r') as file:
        for pair in file.read().strip().split(">")[1:]:
            data = pair.split()
            output[data[0]] = "".join(data[1:])

    return output


def fasta_to_list(file_path: str) -> List[Tuple[str, str]]:
    output = []
    with open(file_path, 'r') as file:
        for pair in file.read().strip().split(">")[1:]:
            data = pair.split()
            output.append((data[0], "".join(data[1:])))

    return output
