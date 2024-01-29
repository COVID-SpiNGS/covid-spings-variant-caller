from typing import Dict, List, Tuple, TypedDict


class Site(TypedDict):
    """
    Represents a data structure for storing information about a specific genomic site.

    Attributes:
        reference (str): A string representing the reference identifier for the genomic site. This could be a gene name, a chromosome number, or any other relevant identifier used to refer to a specific location in the genome.
        total_depth (int): An integer representing the total depth of sequencing at this site. The depth is a measure of how many times a particular site was sequenced, and a higher depth usually indicates higher confidence in the sequencing data.
        snvs (Dict[str, List[int]]): A dictionary mapping single nucleotide variants (SNVs) to their respective positions. The keys are strings representing the nucleotide change (e.g., 'A->T'), and the values are lists of integers representing the positions in the genome where this change occurs.
        indels (Dict[str, List[int]]): A dictionary mapping insertions and deletions (indels) to their respective positions. Similar to 'snvs', the keys are strings representing the indel event (e.g., '+A', '-T'), and the values are lists of integers indicating the positions where these indels occur.
    """
    reference: str
    total_depth: int
    snvs: Dict[str, List[int]]
    indels: Dict[str, List[int]]


class Variant(TypedDict):
    """
    Represents a genomic variant structure, providing key information about a specific genetic variation.

    Attributes:
        start (int): An integer indicating the start position of the variant in the genome. This is usually a base pair number on a specific chromosome.
        stop (int): An integer indicating the stop position of the variant. For single nucleotide variants, this might be the same as the start position. For longer variants, it will be greater than the start position.
        alleles (Tuple[str, str]): A tuple of strings representing the alleles involved in the variant. The first element in the tuple typically represents the reference allele, and the second element represents the variant allele.
        qual (int): An integer representing the quality score of the variant call. A higher score typically indicates greater confidence in the accuracy of the variant call.
        info (Dict): A dictionary containing additional information about the variant. The keys and values in this dictionary can vary depending on the source of the data and the specifics of the variant calling process. This may include information like variant type, effect predictions, population frequencies, etc.
    """

    start: int
    stop: int
    alleles: Tuple[str, str]
    qual: int
    info: Dict
