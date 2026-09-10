import numpy as np
import gudhi

from scipy.spatial import Delaunay
from pycasta.utils.geometry_utils import compute_all_circumsphere_radii

rng = np.random.default_rng(20260831)
points = rng.normal(size=(20, 3))

alpha = 1.25
alpha_sq = alpha ** 2

tri = Delaunay(points)
simplices = np.sort(tri.simplices, axis=1)
tetra_positions = points[simplices]

our_radii = compute_all_circumsphere_radii(tetra_positions)
our_selected = {
    tuple(simplex)
    for simplex, radius in zip(simplices, our_radii)
    if np.isfinite(radius) and radius < alpha
}

alpha_complex = gudhi.AlphaComplex(points=points)
simplex_tree = alpha_complex.create_simplex_tree(max_alpha_square=alpha_sq)

gudhi_selected = {
    tuple(sorted(simplex))
    for simplex, filtration in simplex_tree.get_simplices()
    if len(simplex) == 4 and filtration <= alpha_sq
}

only_ours = our_selected - gudhi_selected
only_gudhi = gudhi_selected - our_selected
common = our_selected & gudhi_selected

print(f"alpha pyCASTa = {alpha:.6f} Å")
print(f"alpha GUDHI   = {alpha_sq:.6f} Å²")
print(f"tetra Delaunay totali = {len(simplices)}")
print(f"tetra pyCASTa selezionati = {len(our_selected)}")
print(f"tetra GUDHI selezionati   = {len(gudhi_selected)}")
print(f"tetra comuni              = {len(common)}")
print(f"solo pyCASTa              = {len(only_ours)}")
print(f"solo GUDHI                = {len(only_gudhi)}")

if only_ours:
    print("Esempio solo pyCASTa:", sorted(only_ours)[:5])

if only_gudhi:
    print("Esempio solo GUDHI:", sorted(only_gudhi)[:5])

if not only_ours and not only_gudhi:
    print("PASS: la selezione dei tetraedri coincide esattamente.")
else:
    print("FAIL: la selezione non coincide; confrontare formula, soglia e degenerazioni.")