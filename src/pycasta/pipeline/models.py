from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class PocketCandidate:
    pocket_id: int
    tetra_indices: List[int]
    representative_point: Optional[List[float]] = None
    volume: Optional[float] = None
    depth: Optional[float] = None
    mouth_area: Optional[float] = None
    mouth_perimeter: Optional[float] = None
    mouth_center: Optional[List[float]] = None
    rim_atoms: List[int] = field(default_factory=list)
    ligand_distance_to_mouth: Optional[float] = None
    features: Dict[str, Any] = field(default_factory=dict)
    scores: Dict[str, float] = field(default_factory=dict)
    tags: List[str] = field(default_factory=list)