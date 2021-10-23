from rosalind._cli.analyze_router import Analyze
from rosalind._cli.compare_router import Compare
from rosalind._cli.infer_router import Infer
from rosalind._cli.search_router import Search
from rosalind._cli.simulate_router import Simulate
from rosalind._cli.transform_router import Transform


import fire


class MainRouter:
    """
    The solutions to rosalind as a library and cli tool
    """

    analyze = Analyze
    compare = Compare
    infer = Infer
    search = Search
    simulate = Simulate
    transform = Transform


if __name__ == "__main__":
    fire.Fire(MainRouter, name="rosalind")
