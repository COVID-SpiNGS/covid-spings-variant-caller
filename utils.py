import numpy as np

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
        return (1.0 - position['alleleErrorProbabilities']['A']).prod() \
            * position['alleleErrorProbabilities']['T'].prod() \
            * position['alleleErrorProbabilities']['C'].prod() \
            * position['alleleErrorProbabilities']['G'].prod()
    elif ref == 'T':
        return position['alleleErrorProbabilities']['A'].prod() \
            * (1.0 - position['alleleErrorProbabilities']['T']).prod() \
            * position['alleleErrorProbabilities']['C'].prod() \
            * position['alleleErrorProbabilities']['G'].prod()
    elif ref == 'C':
        return position['alleleErrorProbabilities']['A'].prod() \
            * position['alleleErrorProbabilities']['T'].prod() \
            * (1.0 - position['alleleErrorProbabilities']['C']).prod() \
            * position['alleleErrorProbabilities']['G'].prod()
    elif ref == 'G':
        return position['alleleErrorProbabilities']['A'].prod() \
            * position['alleleErrorProbabilities']['T'].prod() \
            * position['alleleErrorProbabilities']['C'].prod() \
            * (1.0 - position['alleleErrorProbabilities']['G']).prod()

def to_genotype_quality(genotypeLikelihood: float):
    return np.round(to_phred_quality_score(1.0 - genotypeLikelihood), 0)
