from __future__ import annotations

import argparse
import csv
from collections import deque
from pathlib import Path

import numpy as np


OFFSETS_26 = [
    (dx, dy, dz)
    for dx in (-1, 0, 1)
    for dy in (-1, 0, 1)
    for dz in (-1, 0, 1)
    if (dx, dy, dz) != (0, 0, 0)
]


def load_hetatm_xyz(pdb_path: Path) -> np.ndarray:
    coords = []

    for line in pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith("HETATM"):
            continue

        resname = line[17:20].strip().upper()
        if resname in {"HOH", "WAT", "DOD"}:
            continue

        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue

        coords.append((x, y, z))

    return np.asarray(coords, dtype=float).reshape((-1, 3))


def connected_components(voxel_keys: set[tuple[int, int, int]]) -> list[list[tuple[int, int, int]]]:
    remaining = set(voxel_keys)
    components = []

    while remaining:
        seed = remaining.pop()
        queue = deque([seed])
        component = [seed]

        while queue:
            x, y, z = queue.popleft()

            for dx, dy, dz in OFFSETS_26:
                neighbour = (x + dx, y + dy, z + dz)

                if neighbour in remaining:
                    remaining.remove(neighbour)
                    queue.append(neighbour)
                    component.append(neighbour)

        components.append(component)

    return sorted(components, key=len, reverse=True)


def tetra_bbox_voxels(tetra: np.ndarray, spacing: float) -> np.ndarray:
    lower = np.floor(tetra.min(axis=0) / spacing).astype(int)
    upper = np.floor(tetra.max(axis=0) / spacing).astype(int)

    xs = np.arange(lower[0], upper[0] + 1)
    ys = np.arange(lower[1], upper[1] + 1)
    zs = np.arange(lower[2], upper[2] + 1)

    grid = np.array(np.meshgrid(xs, ys, zs, indexing="ij"))
    return grid.reshape(3, -1).T


def points_inside_tetra(points: np.ndarray, tetra: np.ndarray, tolerance: float = 1e-7) -> np.ndarray:
    a, b, c, d = tetra
    matrix = np.column_stack((b - a, c - a, d - a))

    try:
        inv_matrix = np.linalg.inv(matrix)
    except np.linalg.LinAlgError:
        return np.zeros(len(points), dtype=bool)

    bary = (points - a) @ inv_matrix.T
    w0 = 1.0 - bary.sum(axis=1)

    return (
        (w0 >= -tolerance)
        & (bary[:, 0] >= -tolerance)
        & (bary[:, 1] >= -tolerance)
        & (bary[:, 2] >= -tolerance)
    )


def voxelize_alpha(alpha_npz: Path, spacing: float) -> np.ndarray:
    data = np.load(alpha_npz)

    if "tetrahedra" not in data or "alpha_mask" not in data:
        raise KeyError("Alpha NPZ must contain 'tetrahedra' and 'alpha_mask'.")

    tetrahedra = np.asarray(data["tetrahedra"], dtype=float)
    alpha_mask = np.asarray(data["alpha_mask"], dtype=bool)

    if len(tetrahedra) != len(alpha_mask):
        raise ValueError("tetrahedra and alpha_mask do not have matching lengths.")

    selected_tetrahedra = tetrahedra[alpha_mask]
    voxel_keys = set()

    for tetra in selected_tetrahedra:
        keys = tetra_bbox_voxels(tetra, spacing)
        centers = (keys.astype(float) + 0.5) * spacing
        inside = points_inside_tetra(centers, tetra)

        for key in keys[inside]:
            voxel_keys.add(tuple(int(v) for v in key))

    if not voxel_keys:
        return np.empty((0, 3), dtype=float)

    keys = np.asarray(sorted(voxel_keys), dtype=float)
    return (keys + 0.5) * spacing


def write_voxel_pdb(points: np.ndarray, output_path: Path, residue_name: str) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as handle:
        for index, (x, y, z) in enumerate(points, start=1):
            atom_name = "C"
            residue_id = min(index, 9999)

            handle.write(
                f"HETATM{index:5d} {atom_name:^4s} {residue_name:>3s} A{residue_id:4d}"
                f"    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
            )

        handle.write("END\n")


def nearest_ligand_distance(points: np.ndarray, ligand_xyz: np.ndarray) -> float:
    if len(points) == 0 or len(ligand_xyz) == 0:
        return float("nan")

    distances = np.linalg.norm(points[:, None, :] - ligand_xyz[None, :, :], axis=2)
    return float(distances.min())


def ligand_coverage(points: np.ndarray, ligand_xyz: np.ndarray, cutoff: float) -> float:
    if len(points) == 0 or len(ligand_xyz) == 0:
        return float("nan")

    distances = np.linalg.norm(points[:, None, :] - ligand_xyz[None, :, :], axis=2)
    covered = np.any(distances <= cutoff, axis=0)

    return float(np.mean(covered))


def write_metrics(
    output_path: Path,
    method: str,
    component_id: int,
    points: np.ndarray,
    ligand_xyz: np.ndarray,
    spacing: float,
) -> dict:
    centroid = points.mean(axis=0) if len(points) else np.array([np.nan, np.nan, np.nan])

    row = {
        "method": method,
        "component_id": component_id,
        "n_voxels": int(len(points)),
        "volume_proxy_a3": float(len(points) * spacing**3),
        "center_x": float(centroid[0]),
        "center_y": float(centroid[1]),
        "center_z": float(centroid[2]),
        "nearest_ligand_atom_distance_a": nearest_ligand_distance(points, ligand_xyz),
        "ligand_coverage_4a": ligand_coverage(points, ligand_xyz, 4.0),
        "ligand_coverage_6a": ligand_coverage(points, ligand_xyz, 6.0),
    }

    return row


def save_pymol_script(
    output_path: Path,
    pdb_path: Path,
    alpha_files: list[Path],
    qs_files: list[Path],
) -> None:
    lines = [
        "reinitialize",
        f'load "{pdb_path}", protein',
        "hide everything, all",
        "show surface, protein",
        "color grey70, protein",
        "set transparency, 0.45, protein",
        "show sticks, organic",
        "color white, organic",
        "set stick_radius, 0.18, organic",
    ]

    for index, path in enumerate(alpha_files):
        name = f"alpha_component_{index}"
        lines.extend(
            [
                f'load "{path}", {name}',
                f"show spheres, {name}",
                f"color green, {name}",
                f"set sphere_scale, 0.45, {name}",
            ]
        )

    for index, path in enumerate(qs_files):
        name = f"quicksilver_component_{index}"
        colour = "cyan" if index == 0 else "magenta"

        lines.extend(
            [
                f'load "{path}", {name}',
                f"show spheres, {name}",
                f"color {colour}, {name}",
                f"set sphere_scale, 0.45, {name}",
            ]
        )

    lines.extend(
        [
            "zoom organic, 16",
            "orient organic",
            "bg_color black",
            "set ray_opaque_background, off",
        ]
    )

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Voxelize Alpha Shape and compare it with QuickSILVER components."
    )
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--alpha-npz", type=Path, required=True)
    parser.add_argument("--quicksilver-npz", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--spacing", type=float, default=1.0)
    parser.add_argument("--min-alpha-component-voxels", type=int, default=10)
    args = parser.parse_args()

    if not args.pdb.is_file():
        raise FileNotFoundError(f"PDB not found: {args.pdb}")

    if not args.alpha_npz.is_file():
        raise FileNotFoundError(f"Alpha NPZ not found: {args.alpha_npz}")

    if not args.quicksilver_npz.is_file():
        raise FileNotFoundError(f"QuickSILVER components NPZ not found: {args.quicksilver_npz}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    molecule_id = args.pdb.stem
    ligand_xyz = load_hetatm_xyz(args.pdb)

    print(f"Molecule: {molecule_id}")
    print(f"Ligand HETATM atoms: {len(ligand_xyz)}")
    print(f"Voxel spacing: {args.spacing:.2f} Å")

    alpha_points = voxelize_alpha(args.alpha_npz, args.spacing)
    alpha_keys = np.floor(alpha_points / args.spacing).astype(int)
    alpha_components = connected_components(
        set(tuple(int(v) for v in key) for key in alpha_keys)
    )

    alpha_files = []
    metric_rows = []

    for component_id, component_keys in enumerate(alpha_components):
        if len(component_keys) < args.min_alpha_component_voxels:
            continue

        points = (np.asarray(component_keys, dtype=float) + 0.5) * args.spacing
        output_pdb = args.output_dir / f"{molecule_id}_alpha_component_{component_id}.pdb"

        write_voxel_pdb(points, output_pdb, "ALP")
        alpha_files.append(output_pdb)

        metric_rows.append(
            write_metrics(
                args.output_dir,
                "alpha_shape",
                component_id,
                points,
                ligand_xyz,
                args.spacing,
            )
        )

        print(
            f"Alpha component {component_id}: "
            f"{len(points)} voxels -> {output_pdb.name}"
        )

    qs_data = np.load(args.quicksilver_npz)
    qs_centers = np.asarray(qs_data["voxel_centers"], dtype=float)
    qs_ids = np.asarray(qs_data["cave_id"], dtype=int)

    qs_files = []

    for cave_id in sorted(np.unique(qs_ids)):
        points = qs_centers[qs_ids == cave_id]
        output_pdb = args.output_dir / f"{molecule_id}_quicksilver_component_{cave_id}.pdb"

        write_voxel_pdb(points, output_pdb, "QSV")
        qs_files.append(output_pdb)

        metric_rows.append(
            write_metrics(
                args.output_dir,
                "quicksilver",
                int(cave_id),
                points,
                ligand_xyz,
                args.spacing,
            )
        )

        print(
            f"QuickSILVER component {cave_id}: "
            f"{len(points)} voxels -> {output_pdb.name}"
        )

    metric_path = args.output_dir / f"{molecule_id}_alpha_vs_quicksilver_metrics.csv"

    fieldnames = [
        "method",
        "component_id",
        "n_voxels",
        "volume_proxy_a3",
        "center_x",
        "center_y",
        "center_z",
        "nearest_ligand_atom_distance_a",
        "ligand_coverage_4a",
        "ligand_coverage_6a",
    ]

    with metric_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(metric_rows)

    pml_path = args.output_dir / f"{molecule_id}_alpha_vs_quicksilver.pml"
    save_pymol_script(pml_path, args.pdb, alpha_files, qs_files)

    print()
    print(f"Wrote metrics: {metric_path}")
    print(f"Wrote PyMOL script: {pml_path}")
    print("Open the PML in PyMOL to view Alpha Shape (green) and QuickSILVER (cyan/magenta).")


if __name__ == "__main__":
    main()