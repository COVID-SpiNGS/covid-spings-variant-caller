import operator
import functools
import numpy as np

import math

from typing import Dict, List


def from_phred_score(score: float) -> float:
    """
    Converts a Phred quality score to its corresponding error probability.

    @param score: The Phred quality score.
    @type score: float
    @return: A tuple containing the original score and its error probability.
    @rtype: tuple
    """
    return math.pow(10, -score / 10)


def to_phred_score(probability: float, threshold: int = 64) -> int:
    """
    Converts a probability to its corresponding Phred quality score.

    @param probability: The probability of an error.
    @type probability: float
    @param threshold: The maximum Phred score to return. Defaults to 99.
    @type threshold: int, optional
    @return: The Phred score corresponding to the given probability.
    @rtype: int
    """
    if not 0 <= probability <= 1:
        raise ValueError("Probability must be between 0 and 1")

    if probability == 1.0:
        return threshold
    elif probability == 0.0:
        return 0
    else:
        return min(round(-10 * math.log10(1 - probability)), threshold)


def genotype_likelihood_old(hypothesis: str, alleles: Dict[str, List[float]]):
    """
    Calculates the likelihood of a genotype hypothesis given allele frequencies.

    @param hypothesis: The genotype hypothesis to be tested.
    @type hypothesis: str
    @param alleles: A dictionary mapping allele names to their frequencies.
    @type alleles: dict[str, list[float]]
    @return: The likelihood of the given hypothesis.
    @rtype: float

    This function computes the likelihood of a specific genotype hypothesis. It multiplies
    the probability of the hypothesis allele being correct by the product of the
    probabilities of each non-hypothesis allele being incorrect.
    """
    hypothesis_value = (1.0 - np.array(alleles[hypothesis])).prod()
    non_hypothesis_value = functools.reduce(
        operator.mul,
        {
            allele: np.array(alleles[allele]).prod()
            for allele in alleles
            if allele != hypothesis
        }.values(),
        1.0,
    )

    return round(hypothesis_value * non_hypothesis_value, 2)


## Christians code - Bayessian probability

index_dict = {"A": 0, "C": 1, "G": 2, "T": 3}


def get_likelihood(results):
    """
    Calculates the overall likelihood of a set of sequencing reads.

    @param results: A list of tuples representing sequencing reads. Each tuple contains a base (str) and a corresponding quality score (int).
    @type results: list
    @return: An array representing the normalized likelihood of each base (A, C, G, T) across all reads.
    @rtype: np.ndarray

    The function computes the product of the likelihoods for each read (using `get_likelihood_for_read`) and then normalizes
    this product to sum to 1, unless all values in 'results' are 0 or 'nan', in which case it returns the unnormalized product.
    """
    num = np.prod(list(map(get_likelihood_for_read, results)), axis=0)

    if not np.all(results == 0) and not np.all(results == "nan") and np.sum(num) != 0:
        return num / np.sum(num)

    return num


def get_likelihood_for_read(read):
    """
    Calculates the likelihood of a single read.

    @param read: A tuple containing a base (str) and a quality score (int).
    @type read: tuple
    @return: An array of likelihoods for each possible base (A, C, G, T).
    @rtype: np.ndarray

    This function uses the provided quality score to calculate the likelihood of each base
    being the correct one. The 'index_dict' is expected to map base characters to their
    respective indices in the likelihood array.
    """
    base, score = read
    prob = 10 ** (-score)
    likelihood = np.ones((4,)) * prob / 3
    likelihood[index_dict[base]] = 1 - prob
    return likelihood


def extract_base_likelihood(likelihood_array, snvs_tuples, snvs):
    """
    Extracts the likelihood of specific bases from a likelihood array.

    @param likelihood_array: An array of likelihoods for each base.
    @type likelihood_array: np.ndarray
    @param snvs_tuples: Not used in the function but passed as an argument.
    @type snvs_tuples: list
    @param snvs: Not used in the function but passed as an argument.
    @type snvs: list
    @return: A dictionary mapping bases to their respective likelihoods if 'likelihood_array' is an ndarray, otherwise prints the type and content of 'likelihood_array', along with 'snvs_tuples' and 'snvs'.
    @rtype: dict

    This function maps the likelihoods from the array to their corresponding bases using the 'index_dict'. If 'likelihood_array' is not an ndarray, it provides debugging information.
    """
    if isinstance(likelihood_array, np.ndarray):
        x = {key: likelihood_array[value] for key, value in index_dict.items()}
        return x
    else:
        print(type(likelihood_array), likelihood_array, snvs_tuples, snvs)
