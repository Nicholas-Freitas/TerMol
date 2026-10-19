# termol/backbone.py
"""
Backbone generation and Catmull-Rom spline interpolation for proteins.
Computes smooth 3D trajectories with secondary structure and sequence annotations.
"""

from dataclasses import dataclass
from typing import List, Tuple
import numpy as np
from .pdb import ProteinStructure, ResidueCA


@dataclass
class BackbonePoint:
    coord: np.ndarray        # 3D coordinate (x, y, z)
    norm_pos: float          # 0.0 to 1.0 along the chain (for rainbow coloring)
    ss: str                  # 'H', 'E', or 'C'
    chain_id: str
    res_num: int


def catmull_rom_spline_segment(coords: np.ndarray, num_points_per_seg: int = 10) -> np.ndarray:
    """
    Interpolate a sequence of 3D points using Catmull-Rom cubic splines.
    Inputs:
        coords: ndarray of shape (N, 3), N >= 3
        num_points_per_seg: number of interpolation samples between consecutive points
    Returns:
        interpolated: ndarray of shape (M, 3)
    """
    n = len(coords)
    if n < 2:
        return coords.copy()
    if n == 2:
        # Linear interpolation
        ts = np.linspace(0, 1, num_points_per_seg)[:, None]
        return coords[0] * (1 - ts) + coords[1] * ts

    # Pad endpoints with virtual control points
    p_start = 2 * coords[0] - coords[1]
    p_end = 2 * coords[-1] - coords[-2]
    padded = np.vstack([p_start, coords, p_end])

    result = []
    # Parametric t values for each segment
    ts = np.linspace(0, 1, num_points_per_seg, endpoint=False)
    t2s = ts * ts
    t3s = t2s * ts

    for i in range(n - 1):
        p0, p1, p2, p3 = padded[i], padded[i+1], padded[i+2], padded[i+3]
        # Catmull-Rom basis:
        # P(t) = 0.5 * ( (2*p1) + (-p0 + p2)*t + (2*p0 - 5*p1 + 4*p2 - p3)*t^2 + (-p0 + 3*p1 - 3*p2 + p3)*t^3 )
        c0 = 2 * p1
        c1 = -p0 + p2
        c2 = 2 * p0 - 5 * p1 + 4 * p2 - p3
        c3 = -p0 + 3 * p1 - 3 * p2 + p3

        pts = 0.5 * (c0 + np.outer(ts, c1) + np.outer(t2s, c2) + np.outer(t3s, c3))
        result.append(pts)

    # Add the final endpoint
    result.append(coords[-1:])
    return np.vstack(result)


def generate_backbone_points(structure: ProteinStructure, samples_per_res: int = 8) -> List[List[BackbonePoint]]:
    """
    Generate smooth backbone trajectory points for all chains in a ProteinStructure.
    Handles chain breaks and returns a list of continuous spline segments.
    """
    all_segments: List[List[BackbonePoint]] = []

    total_residues = max(structure.total_residues, 1)
    global_res_counter = 0

    for chain_id, ca_list in structure.chains.items():
        if not ca_list:
            continue

        # Split chain into continuous fragments at gaps (|CA_{i+1} - CA_i| > 4.5 Å)
        fragments: List[List[ResidueCA]] = []
        current_frag: List[ResidueCA] = [ca_list[0]]

        for i in range(1, len(ca_list)):
            dist = np.linalg.norm(ca_list[i].coord - ca_list[i-1].coord)
            if dist > 4.5:
                fragments.append(current_frag)
                current_frag = [ca_list[i]]
            else:
                current_frag.append(ca_list[i])
        if current_frag:
            fragments.append(current_frag)

        # Interpolate each fragment
        for frag in fragments:
            frag_len = len(frag)
            frag_coords = np.array([res.coord for res in frag])

            if frag_len == 1:
                norm_pos = global_res_counter / total_residues
                pt = BackbonePoint(
                    coord=frag_coords[0],
                    norm_pos=norm_pos,
                    ss=frag[0].ss,
                    chain_id=chain_id,
                    res_num=frag[0].res_num
                )
                all_segments.append([pt])
                global_res_counter += 1
                continue

            # Run Catmull-Rom spline
            spline_coords = catmull_rom_spline_segment(frag_coords, num_points_per_seg=samples_per_res)

            # Map each interpolated point to sequence position and secondary structure
            segment_pts: List[BackbonePoint] = []
            m = len(spline_coords)

            for idx, coord in enumerate(spline_coords):
                # Normalized position along the entire protein (for rainbow)
                frac = idx / max(m - 1, 1)
                res_idx = min(int(frac * frag_len), frag_len - 1)
                current_res = frag[res_idx]

                norm_pos = (global_res_counter + frac * frag_len) / total_residues
                norm_pos = min(max(norm_pos, 0.0), 1.0)

                segment_pts.append(BackbonePoint(
                    coord=coord,
                    norm_pos=norm_pos,
                    ss=current_res.ss,
                    chain_id=chain_id,
                    res_num=current_res.res_num
                ))

            all_segments.append(segment_pts)
            global_res_counter += frag_len

    return all_segments
