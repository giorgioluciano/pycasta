import numpy as np

from pycasta.pipeline.models import PocketCandidate
from pycasta.pipeline import feature_pipeline as fp


def _make_candidate():
    c = PocketCandidate(
        pocket_id=1,
        tetra_indices=[0, 1, 2],
        representative_point=[0.0, 0.0, 0.0],
    )
    c.features = {"volume": 123.0}
    return c


def _fake_occlusion_features(**kwargs):
    candidate = kwargs["candidate"]
    return {
        "ray_occlusion_fraction": 0.80,
        "ray_escape_fraction": 0.20,
        "ray_mean_hit_distance": 3.0,
        "ray_min_hit_distance": 1.5,
        "hit_residue_sasa_mean": 2.0,
        "low_sasa_hit_fraction": 1.0,
        "hybrid_occlusion_score": 0.75,
        "occlusion_debug": {
            "num_hits": 77,
            "num_rays": 96,
            "hit_atom_indices": [0, 1],
        },
        "pocket_id_echo": candidate.pocket_id,
    }


def test_hybrid_occlusion_baseline_vs_v2_same_output(monkeypatch):
    candidates_base = [_make_candidate()]
    candidates_v2 = [_make_candidate()]

    protein_coords = np.array([
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [0.0, 1.5, 0.0],
        [0.0, 0.0, 1.5],
    ], dtype=float)

    protein_atoms = ["C", "N", "O", "S"]
    atom_radii = np.array([1.7, 1.55, 1.52, 1.8], dtype=float)

    monkeypatch.setattr(fp.config, "USE_HYBRID_OCCLUSION", True, raising=False)
    monkeypatch.setattr(fp.config, "HYBRID_OCCLUSION_METHOD", "baseline", raising=False)
    monkeypatch.setattr(fp, "compute_hybrid_occlusion_features", _fake_occlusion_features)
    monkeypatch.setattr(fp, "compute_hybrid_occlusion_features_v2", _fake_occlusion_features, raising=False)

    out_base = fp.enrich_with_hybrid_occlusion(
        candidates=candidates_base,
        protein_coords=protein_coords,
        protein_atoms=protein_atoms,
        atom_radii=atom_radii,
        sasa_by_residue=None,
    )

    assert "hybrid_occlusion_score" in out_base[0].features, out_base[0].features

    monkeypatch.setattr(fp.config, "HYBRID_OCCLUSION_METHOD", "v2", raising=False)

    out_v2 = fp.enrich_with_hybrid_occlusion(
        candidates=candidates_v2,
        protein_coords=protein_coords,
        protein_atoms=protein_atoms,
        atom_radii=atom_radii,
        sasa_by_residue=None,
    )

    assert "hybrid_occlusion_score" in out_v2[0].features, out_v2[0].features
    assert out_base[0].features["hybrid_occlusion_score"] == 0.75
    assert out_v2[0].features["hybrid_occlusion_score"] == 0.75
    assert out_base[0].features == out_v2[0].features


def test_hybrid_occlusion_disabled_returns_unchanged(monkeypatch):
    candidates = [_make_candidate()]
    before = dict(candidates[0].features)

    protein_coords = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    ], dtype=float)

    protein_atoms = ["C", "N"]
    atom_radii = np.array([1.7, 1.55], dtype=float)

    monkeypatch.setattr(fp.config, "USE_HYBRID_OCCLUSION", False, raising=False)

    out = fp.enrich_with_hybrid_occlusion(
        candidates=candidates,
        protein_coords=protein_coords,
        protein_atoms=protein_atoms,
        atom_radii=atom_radii,
        sasa_by_residue=None,
    )

    assert out[0].features == before