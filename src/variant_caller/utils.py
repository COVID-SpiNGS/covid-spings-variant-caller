import operator
import functools
import numpy as np

import math

from typing import Dict, List

def from_phred_score(score: float) -> float:
    return score, math.pow(10, score / -10)

def to_phred_score(probability: float, threshold: int = 99) -> int:
    return min(round(-10 * math.log10(probability)), threshold) if probability > 0.0 else threshold


def genotype_likelihood2(hypothesis: str, alleles: Dict[str, List[float]]): 
    hypothesis_value = (1.0 - np.array(alleles[hypothesis])).prod()
    non_hypothesis_value = functools.reduce(operator.mul, {
        allele: np.array(alleles[allele]).prod()
        for allele in alleles
        if allele != hypothesis
    }.values(), 1.0)

    return round(hypothesis_value * non_hypothesis_value, 2)


## Christians code

index_dict = {'A':0 , 'C': 1, 'G':2, 'T': 3}

def get_likelihood(results):
    num = np.prod(list(map(get_likelihood_for_read,results)),axis=0)
    return num/np.sum(num)

def get_likelihood_for_read(read):
    base,score = read
    prob = 10**(-score)
    likelihood = np.ones((4,))*prob/3
    likelihood[index_dict[base]] = 1-prob
    return likelihood

def extract_base_likelihood(likelihood_array, snvs_tuples, snvs):
    #print("LA", likelihood_array, type(likelihood_array))
    if isinstance(likelihood_array, np.ndarray):
        x = {key: likelihood_array[value] for key, value in index_dict.items()}
        #print("LA", likelihood_array, len(likelihood_array))
        return x
    elif np.isnan(likelihood_array).any():
        print(type(likelihood_array), snvs_tuples, snvs)



