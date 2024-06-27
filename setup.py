from setuptools import setup, find_packages

setup(
    name='dna-sequence-analyzer',
    version='0.1',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'dna-analyzer = dna_sequence_analyzer.sequence:main',
        ],
    },
)

