import operator
import functools
from typing import Dict, List
import numpy as np
import math

from structs import Position

def from_phred_scale(score: float) -> float:
    return np.power(10, score / -10)

def to_phred_scale(probability: float, threshold: float = 99):
    if  probability > 0.0:
        score = np.abs(np.multiply(np.log10(probability), -10))

        if score <= threshold:
            return score

    return threshold
# recalibratedAlleleErrorProbabilities , alleleErrorProbabilities
def genotype_likelihood(hypothesis: str, position: Position, key='baseErrorProbabilities'): 
    hypothesisValue = (1.0 - np.array(position[key][hypothesis])).prod()
    alternativeValue = functools.reduce(operator.mul, {
        base: np.array(position[key][base]).prod()
        for base in position[key].keys()
        if base != hypothesis
    }.values())

    return hypothesisValue * alternativeValue

def calculate_genotype_likelihood(hypothesis: str, alleles: Dict[str, List[float]]): 
    hypothesisValue = (1.0 - np.array(alleles[hypothesis])).prod()
    nonHypothesisValue = functools.reduce(operator.mul, {
        allele: np.array(alleles[allele]).prod()
        for allele in alleles
        if allele != hypothesis
    }.values(), 1.0)

    return hypothesisValue * nonHypothesisValue




def error_probability(hypothesis: str, position: Position, key='baseErrorProbabilities'): 
    numHypothesis = len(position[key][hypothesis])
    numAlternative = functools.reduce(operator.add, {
        base: len(position[key][base])
        for base in position[key].keys()
        if base != hypothesis
    }.values())
    numTotal = numHypothesis + numAlternative

    # hypothesisValue = position[key][hypothesis].mean()
    hypothesisValue = (1.0 - position[key][hypothesis]).mean()

    if numAlternative > 0:
        alternativeValue =np.mean(list({
            base: (1.0 - position[key][base]).mean()
            for base in position[key].keys()
            if base != hypothesis and  len(position[key][base]) > 0
        }.values()))

        if position['pos'] == 18873:
            print(numHypothesis, hypothesisValue, numHypothesis / numTotal * hypothesisValue)
            print(numAlternative, alternativeValue, numAlternative / numTotal * alternativeValue )

            weightedHypothesisValue  =  numHypothesis / numTotal * hypothesisValue
            weightedAlternativeValue =  numAlternative / numTotal * alternativeValue
            weigtedTotalValue = weightedHypothesisValue + weightedAlternativeValue
            print(weigtedTotalValue)

            print() 
    else:
        return hypothesisValue


