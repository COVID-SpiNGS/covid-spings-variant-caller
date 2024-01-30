from typing import Dict, List, Tuple, TypedDict


class Site(TypedDict):
    """
    Represents a data structure for storing information about a specific genomic site.

    @ivar reference: A string representing the reference identifier for the genomic site.
    @type reference: str
    @ivar total_depth: An integer representing the total depth of sequencing at this site. Sequencing depth is a critical measure of the data's reliability, with higher values indicating greater confidence.
    @type total_depth: int
    @ivar snvs: A dictionary that maps single nucleotide variants (SNVs) to their respective positions within the genome. The dictionary keys are strings indicating the nucleotide change,
    while the values are lists of integers marking the quality of these changes.
    @type snvs: Dict[str, List[int]]
    @ivar indels: A dictionary that maps insertions and deletions (indels) to their respective positions. -> placeholder for the future!
    @type indels: Dict[str, List[int]]
    """
    reference: str
    total_depth: int
    snvs: Dict[str, List[int]]
    indels: Dict[str, List[int]]


class Variant(TypedDict):
    """
    Represents a genomic variant structure, providing key information about a specific genetic variation.

    @ivar start: An integer indicating the start position of the variant in the genome.
    @type start: int
    @ivar stop: An integer indicating the stop position of the variant.
    @type stop: int
    @ivar alleles: A tuple of strings representing the alleles involved in the variant. The first element is usually the reference allele, and the second is the variant allele.
    @type alleles: Tuple[str, str]
    @ivar qual: An integer representing the quality score of the variant call, where a higher score indicates greater confidence in the variant call's accuracy.
    @type qual: int
    @ivar info: A dictionary containing additional information about the variant, such as variant type, effect predictions, population frequencies, etc.
    @type info: Dict
    """

    start: int
    stop: int
    alleles: Tuple[str, str]
    qual: int
    info: Dict
