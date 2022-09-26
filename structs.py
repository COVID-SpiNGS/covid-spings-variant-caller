import numpy as np
from typing import List, TypedDict

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
