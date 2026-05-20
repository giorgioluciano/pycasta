import math


def compute_basic_score(candidate):
    volume = float(candidate.features.get("volume", 0.0))
    depth = float(candidate.features.get("depth", 0.0))
    mouth_area = float(candidate.features.get("mouth_area", 0.0) or 0.0)
    rim_atom_count = float(candidate.features.get("rim_atom_count", 0.0))

    score = (
        1.0 * volume
        + 0.35 * depth
        + 0.10 * math.log1p(mouth_area)
        + 0.05 * math.log1p(rim_atom_count)
    )
    return float(score)