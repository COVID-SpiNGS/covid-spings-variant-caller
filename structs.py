from re import A, T
import numpy as np
from typing import Dict, List, TypedDict

class Calls(TypedDict):
    A: int
    T: int
    G: int
    C: int

class ErrorProbabilities(TypedDict):
    A: List[float]
    T: List[float]
    G: List[float]
    C: List[float]
class PhredQualityScores(TypedDict):
    A: List[float]
    T: List[float]
    G: List[float]
    C: List[float]

class Position(TypedDict):
    reference: str
    totalDepth: int
    candidatesDepth: int
    baseFrequencies: Calls
    baseQualities: PhredQualityScores
    mappingQualities: PhredQualityScores
    baseErrorProbabilities: ErrorProbabilities


class Alternative(TypedDict):
    isRelevant: bool
    alt: str
    qual: int


class Allele(TypedDict):
    evidenceDepth: int
    baseErrorProbabilities: List[float]
    baseQualities: List[float]
    mappingQualities: List[float]

class MultiPosition(TypedDict):
    reference: str
    totalDepth: int
    alleles: Dict[str, Allele]