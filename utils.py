import numpy as np
from config import minDepth

def is_relevant_position(position, minNumReads=minDepth):
    if(position['numReads'] >= minNumReads):
        alt = max(position['calls'], key=position['calls'].get)
        return alt != position['ref']
    else:
        return False

def to_error_probability(phredQualityScore: float):
    return np.power(10, phredQualityScore / -10)

def to_phred_quality_score(errorProbability: float, maxThreshold: float = 99):
    if  errorProbability > 0.0:
        score = np.abs(np.multiply(np.log10(errorProbability), -10))

        if score < maxThreshold:
            return score

    return maxThreshold

def genotype_likelihood(ref: str, position): 
    if ref == 'A':
        return (1.0 - position['errorProbabilities']['A']).prod() \
            * position['errorProbabilities']['T'].prod() \
            * position['errorProbabilities']['C'].prod() \
            * position['errorProbabilities']['G'].prod()
    elif ref == 'T':
        return position['errorProbabilities']['A'].prod() \
            * (1.0 - position['errorProbabilities']['T']).prod() \
            * position['errorProbabilities']['C'].prod() \
            * position['errorProbabilities']['G'].prod()
    elif ref == 'C':
        return position['errorProbabilities']['A'].prod() \
            * position['errorProbabilities']['T'].prod() \
            * (1.0 - position['errorProbabilities']['C']).prod() \
            * position['errorProbabilities']['G'].prod()
    elif ref == 'G':
        return position['errorProbabilities']['A'].prod() \
            * position['errorProbabilities']['T'].prod() \
            * position['errorProbabilities']['C'].prod() \
            * (1.0 - position['errorProbabilities']['G']).prod()

def to_genotype_quality(genotypeLikelihood: float):
    return np.round(to_phred_quality_score(1.0 - genotypeLikelihood), 0)
