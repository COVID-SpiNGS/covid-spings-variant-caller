import numpy as np
from typing import List, TypedDict

class Calls(TypedDict):
    A: int
    T: int
    G: int
    C: int

class ErrorProbabilities(TypedDict):
    A: np.ndarray
    T: np.ndarray
    G: np.ndarray
    C: np.ndarray

class PhredQualityScores(TypedDict):
    A: np.ndarray
    T: np.ndarray
    G: np.ndarray
    C: np.ndarray

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
