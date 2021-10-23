from math import comb
from typing import Tuple


def number_of_rabbits_with_offspring(number_months: int, number_offspring: int) -> int:
    if 0 < number_months < 3:
        return 1
    i = 3
    r1, r2 = 1, 1
    r3 = r1 + number_offspring * r2

    while i < number_months:
        r1, r2 = r3, r1
        r3 = r1 + number_offspring * r2
        i += 1

    return r3


def probability_offspring_has_dominant_allele(homo_dom: int, hetero: int, homo_rec: int) -> float:
    k, m, n = homo_dom, hetero, homo_rec
    t = homo_dom + hetero + homo_rec
    if t == 0:
        return 0

    output = k * (k - 1) + 2 * k * m + 2 * k * n
    output += m * n + 0.75 * m * (m - 1)
    output /= t * (t - 1)

    return output


def number_of_mortal_rabbits(number_months: int, months_until_death: int) -> int:
    assert months_until_death > 2
    if 0 < number_months < 3:
        return 1

    dp = [0 for _ in range(number_months)]
    dp[0], dp[1] = 1, 1
    dp[months_until_death], dp[months_until_death + 1] = -1, -1

    for i in range(2, number_months):
        dp[i] += dp[i-1] + dp[i-2]
        if i > months_until_death + 1:
            dp[i] -= dp[i - months_until_death - 1]

    return dp[number_months - 1]


def expected_offspring(number_individuals: Tuple) -> float:
    probs = (1.0, 1.0, 1.0, 0.75, 0.5, 0)
    return 2 * sum(map(lambda x: x[0] * x[1],
                       zip(number_individuals, probs)
                       )
                   )


def probability_at_least_n_organisms(generations: int, n_individuals: int) -> float:
    total_population = 2 ** generations
    return 1 - sum(
        map(lambda x: comb(total_population, x) * pow(0.25, x) * pow(0.75, total_population - x),
            range(n_individuals)
            )
    )
