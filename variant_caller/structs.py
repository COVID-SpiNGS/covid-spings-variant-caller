from typing import Dict, List, TypedDict
class Site(TypedDict):
    reference: str
    totalDepth: int
    baseQualities: Dict[str, List[int]]