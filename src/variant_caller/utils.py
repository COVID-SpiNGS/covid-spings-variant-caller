import operator
import functools
import numpy as np

import math

from typing import Dict, List

def from_phred_score(score: float) -> float:
    return math.pow(10, score / -10)

def to_phred_score(probability: float, threshold: int = 99) -> int:
    return min(round(-10 * math.log10(probability)), threshold) if probability > 0.0 else threshold


def genotype_likelihood(hypothesis: str, alleles: Dict[str, List[float]]): 
    hypothesis_value = (1.0 - np.array(alleles[hypothesis])).prod()
    non_hypothesis_value = functools.reduce(operator.mul, {
        allele: np.array(alleles[allele]).prod()
        for allele in alleles
        if allele != hypothesis
    }.values(), 1.0)

    return hypothesis_value * non_hypothesis_value




