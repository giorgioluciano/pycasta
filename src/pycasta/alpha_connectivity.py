#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections import defaultdict
import numpy as np


def build_alpha_tetra_adjacency(simplices, alpha_mask):
    simplices = np.asarray(simplices, dtype=int)
    alpha_mask = np.asarray(alpha_mask, dtype=bool)

    alpha_indices = np.where(alpha_mask)[0]
    alpha_set = set(alpha_indices.tolist())

    face_to_tetra = defaultdict(list)
    for tetra_idx in alpha_indices:
        tet = simplices[tetra_idx]
        faces = [
            tuple(sorted((tet[0], tet[1], tet[2]))),
            tuple(sorted((tet[0], tet[1], tet[3]))),
            tuple(sorted((tet[0], tet[2], tet[3]))),
            tuple(sorted((tet[1], tet[2], tet[3]))),
        ]
        for face in faces:
            face_to_tetra[face].append(int(tetra_idx))

    graph = {int(i): set() for i in alpha_indices}
    for tetra_list in face_to_tetra.values():
        if len(tetra_list) < 2:
            continue
        for i in range(len(tetra_list)):
            for j in range(i + 1, len(tetra_list)):
                a = tetra_list[i]
                b = tetra_list[j]
                if a in alpha_set and b in alpha_set:
                    graph[a].add(b)
                    graph[b].add(a)

    return graph


def connected_components_from_graph(graph):
    visited = set()
    components = []

    for node in graph:
        if node in visited:
            continue
        stack = [node]
        component = []

        while stack:
            cur = stack.pop()
            if cur in visited:
                continue
            visited.add(cur)
            component.append(cur)
            stack.extend(n for n in graph[cur] if n not in visited)

        components.append(sorted(component))

    components.sort(key=len, reverse=True)
    return components


def compute_alpha_connectivity_metrics(simplices, alpha_mask):
    alpha_mask = np.asarray(alpha_mask, dtype=bool)
    num_alpha_tetrahedra = int(alpha_mask.sum())

    if num_alpha_tetrahedra == 0:
        return {
            "num_alpha_components": 0,
            "largest_component_size": 0,
            "largest_component_ratio": 0.0,
            "small_component_fraction": 0.0,
            "component_size_entropy": 0.0,
        }

    graph = build_alpha_tetra_adjacency(simplices, alpha_mask)
    components = connected_components_from_graph(graph)

    component_sizes = np.array([len(c) for c in components], dtype=float)
    largest_component_size = int(component_sizes.max()) if component_sizes.size else 0
    largest_component_ratio = (
        float(largest_component_size / num_alpha_tetrahedra)
        if num_alpha_tetrahedra > 0
        else 0.0
    )

    if component_sizes.size <= 1:
        small_component_fraction = 0.0
        component_size_entropy = 0.0
    else:
        small_component_fraction = float(
            component_sizes[1:].sum() / num_alpha_tetrahedra
        )
        probs = component_sizes / component_sizes.sum()
        component_size_entropy = float(
            -np.sum(probs * np.log(probs + 1e-12))
        )

    return {
        "num_alpha_components": int(len(components)),
        "largest_component_size": largest_component_size,
        "largest_component_ratio": largest_component_ratio,
        "small_component_fraction": small_component_fraction,
        "component_size_entropy": component_size_entropy,
    }