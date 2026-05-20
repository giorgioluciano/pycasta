from ..scorers.basic_scorer import compute_basic_score


def score_candidates(candidates):
    for candidate in candidates:
        candidate.scores["basic"] = compute_basic_score(candidate)
        candidate.scores["global"] = candidate.scores["basic"]

    return sorted(
        candidates,
        key=lambda c: c.scores.get("global", 0.0),
        reverse=True,
    )