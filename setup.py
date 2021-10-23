from distutils.core import setup

setup(
    name="rosalind",
    version="0.0.1",
    description="The solutions to rosalind as a package and cli tool.",
    author="Jorge Duque",
    install_requires=["fire", "suffix_tree", "requests"],
    scripts=[
        'scripts/rosalind'
    ]
)
