import operator
import functools
import numpy as np

from structs import Position

def to_error_probability(phredQualityScore: float):
    return np.power(10, phredQualityScore / -10)

def to_phred_quality_score(errorProbability: float, maxThreshold: float = 99):
    if  errorProbability > 0.0:
        score = np.abs(np.multiply(np.log10(errorProbability), -10))

        if score <= maxThreshold:
            return score

    return maxThreshold
# recalibratedAlleleErrorProbabilities , alleleErrorProbabilities
def genotype_likelihood(hypothesis: str, position: Position, key='alleleErrorProbabilities'): 
    hypothesisValue = (1.0 - position[key][hypothesis]).prod()
    alternativeValue = functools.reduce(operator.mul, {
        base: position[key][base].prod()
        for base in position[key].keys()
        if base != hypothesis
    }.values())

    return hypothesisValue * alternativeValue


def error_probability(hypothesis: str, position: Position, key='alleleErrorProbabilities'): 
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



def to_phred_scale(genotypeLikelihood: float):
    return np.round(to_phred_quality_score(1.0 - genotypeLikelihood), 0)
