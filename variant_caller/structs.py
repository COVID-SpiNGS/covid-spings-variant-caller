from typing import Dict, List, Tuple, TypedDict
class Site(TypedDict):
    reference: str
    totalDepth: int
    snvs: Dict[str, List[int]]
    indels: Dict[str, List[int]]


class Variant(TypedDict):
    start: int
    stop: int
    alleles: Tuple[str, str]
    qual: int
    info: Dict