import rosalind.simulate as simulate


class Simulate:
    """
    A collection of solutions to toy problems
    """

    @staticmethod
    def rabbits(month: int, n_offspring: int) -> int:
        """Computes rabbit population at month when rabbits have a given number of offspring"""
        return simulate.number_of_mortal_rabbits(month, n_offspring)

    @staticmethod
    def dominant_allele(homo_dom: int, homo_rec: int, het: int) -> float:
        """
        Computes the probability of an offspring having the domninant allele
        given the number of homogenous and hetergenous people in the population
        """
        return simulate.probability_offspring_has_dominant_allele(homo_dom, het, homo_rec)

    @staticmethod
    def mortal_rabbits(month: int, months_until_death: int) -> int:
        """Computes the population size at month if rabbits have 2 offspring and die after months_until_deaths"""
        return simulate.number_of_mortal_rabbits(month, months_until_death)

    # TODO: write proper docstring
    @staticmethod
    def prob_organisms(generation: int, n_individuals: int) -> float:
        """Computes the probability of at least n organisms being alive after"""
        return simulate.probability_at_least_n_organisms(generation, n_individuals)
