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
    pos: int
    ref: str
    numReads: int
    calls: Calls
    errorProbabilities: ErrorProbabilities
    phredQualityScores: PhredQualityScores

positions: List[Position] = []