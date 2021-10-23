from _cli.analyze_router import Analyze
from _cli.compare_router import Compare
from _cli.infer_router import Infer
from _cli.search_router import Search
from _cli.simulate_router import Simulate
from _cli.transform_router import Transform


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
