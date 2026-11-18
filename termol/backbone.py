# termol/backbone.py
"""
Backbone generation and Catmull-Rom spline interpolation for proteins.
Computes smooth 3D trajectories with secondary structure annotations,
ribbon geometry (with beta-sheet arrowheads and alpha-helix ribbons),
and 3D backbone tube thickness.
"""

from dataclasses import dataclass
from typing import List, Tuple, Optional
import numpy as np
from .pdb import ProteinStructure, ResidueCA


@dataclass
class BackbonePoint:
    coord: np.ndarray        # 3D coordinate (x, y, z)
    norm_pos: float          # 0.0 to 1.0 along the chain (for rainbow coloring)
    ss: str                  # 'H', 'E', or 'C'
    chain_id: str
    res_num: int
    left_coord: Optional[np.ndarray] = None
    right_coord: Optional[np.ndarray] = None
    top_coord: Optional[np.ndarray] = None
    bottom_coord: Optional[np.ndarray] = None


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
        c0 = 2 * p1
        c1 = -p0 + p2
        c2 = 2 * p0 - 5 * p1 + 4 * p2 - p3
        c3 = -p0 + 3 * p1 - 3 * p2 + p3

        pts = 0.5 * (c0 + np.outer(ts, c1) + np.outer(t2s, c2) + np.outer(t3s, c3))
        result.append(pts)

    # Add the final endpoint
    result.append(coords[-1:])
    return np.vstack(result)


def compute_parallel_transport_frame(points: np.ndarray, guide_vector: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute smooth Rotation-Minimizing Frames (RMF / Bishop frame) along a 3D curve.
    Prevents twists and flipping at inflection points.
    Returns: (tangents, normals, binormals)
    """
    N = len(points)
    if N < 2:
        zeros = np.zeros((N, 3))
        return zeros, zeros, zeros

    tangents = np.zeros((N, 3))
    for i in range(N):
        if i == 0:
            t = points[1] - points[0]
        elif i == N - 1:
            t = points[-1] - points[-2]
        else:
            t = points[i+1] - points[i-1]
        norm = np.linalg.norm(t)
        tangents[i] = t / (norm if norm > 1e-6 else 1.0)

    normals = np.zeros((N, 3))
    binormals = np.zeros((N, 3))

    t0 = tangents[0]
    n0 = None
    if guide_vector is not None:
        proj = guide_vector - np.dot(guide_vector, t0) * t0
        if np.linalg.norm(proj) > 1e-4:
            n0 = proj / np.linalg.norm(proj)
    if n0 is None and N > 2:
        curv0 = points[2] - 2 * points[1] + points[0]
        proj = curv0 - np.dot(curv0, t0) * t0
        if np.linalg.norm(proj) > 1e-4:
            n0 = proj / np.linalg.norm(proj)
    if n0 is None:
        up = np.array([0.0, 1.0, 0.0]) if abs(t0[1]) < 0.9 else np.array([1.0, 0.0, 0.0])
        n0 = np.cross(t0, up)
        n0 /= np.linalg.norm(n0)

    b0 = np.cross(t0, n0)
    norm_b0 = np.linalg.norm(b0)
    b0 = b0 / (norm_b0 if norm_b0 > 1e-6 else 1.0)
    normals[0] = n0
    binormals[0] = b0

    for i in range(N - 1):
        t1, t2 = tangents[i], tangents[i+1]
        v = np.cross(t1, t2)
        c = np.dot(t1, t2)
        if np.linalg.norm(v) < 1e-6 or c > 0.9999:
            n = normals[i]
        else:
            axis = v / np.linalg.norm(v)
            angle = np.arccos(np.clip(c, -1.0, 1.0))
            n = normals[i] * np.cos(angle) + np.cross(axis, normals[i]) * np.sin(angle) + axis * np.dot(axis, normals[i]) * (1 - np.cos(angle))
            norm_n = np.linalg.norm(n)
            n = n / (norm_n if norm_n > 1e-6 else 1.0)

        b = np.cross(tangents[i+1], n)
        norm_b = np.linalg.norm(b)
        b = b / (norm_b if norm_b > 1e-6 else 1.0)
        normals[i+1] = n
        binormals[i+1] = b

    return tangents, normals, binormals


def generate_backbone_points(structure: ProteinStructure, samples_per_res: int = 8,
                             representation: str = 'trace') -> List[List[BackbonePoint]]:
    """
    Generate smooth backbone trajectory points for all chains in a ProteinStructure.
    representation options:
        'trace': single 1D Catmull-Rom spline line (default, minimal overhead)
        'cartoon' / 'ribbon': flat ribbons for helices, arrowheads for sheets, thin loops
        'tube': uniform 3D tube thickness across entire backbone
    """
    all_segments: List[List[BackbonePoint]] = []
    total_residues = max(structure.total_residues, 1)
    global_res_counter = 0

    rep = representation.lower()
    is_ribbon = rep in ('cartoon', 'ribbon')
    is_tube = rep == 'tube'

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

            spline_coords = catmull_rom_spline_segment(frag_coords, num_points_per_seg=samples_per_res)
            m = len(spline_coords)

            # Assign SS and residue mapping per spline point
            ss_classes = []
            res_indices = []
            for idx in range(m):
                frac = idx / max(m - 1, 1)
                r_idx = min(int(frac * frag_len), frag_len - 1)
                res_indices.append(r_idx)
                ss_classes.append(frag[r_idx].ss)

            left_coords = None
            right_coords = None
            top_coords = None
            bottom_coords = None

            if is_ribbon or is_tube:
                guide = None
                for res in frag:
                    if res.o_coord is not None:
                        guide = res.o_coord - res.coord
                        break
                tangents, normals, binormals = compute_parallel_transport_frame(spline_coords, guide_vector=guide)

                if is_ribbon:
                    widths = np.zeros(m)
                    for i in range(m):
                        ss = ss_classes[i]
                        if ss == 'H':
                            widths[i] = 1.3
                        elif ss == 'E':
                            widths[i] = 1.8
                        else:
                            widths[i] = 0.35

                    # Arrowheads at the C-terminus of each beta-strand
                    i = 0
                    while i < m:
                        if ss_classes[i] == 'E':
                            start_i = i
                            while i < m and ss_classes[i] == 'E':
                                i += 1
                            end_i = i
                            arrow_len = min(int(samples_per_res * 1.5), end_i - start_i)
                            arrow_start = end_i - arrow_len
                            for k in range(arrow_start, end_i):
                                rem = end_i - 1 - k
                                frac = rem / max(arrow_len, 1)
                                widths[k] = 2.8 * frac
                        else:
                            i += 1

                    left_coords = spline_coords - widths[:, None] * binormals
                    right_coords = spline_coords + widths[:, None] * binormals

                elif is_tube:
                    radius = 0.8
                    left_coords = spline_coords - radius * binormals
                    right_coords = spline_coords + radius * binormals
                    top_coords = spline_coords + radius * normals
                    bottom_coords = spline_coords - radius * normals

            segment_pts: List[BackbonePoint] = []
            for idx in range(m):
                r_idx = res_indices[idx]
                frac = idx / max(m - 1, 1)
                norm_pos = (global_res_counter + frac * frag_len) / total_residues
                norm_pos = min(max(norm_pos, 0.0), 1.0)

                pt = BackbonePoint(
                    coord=spline_coords[idx],
                    norm_pos=norm_pos,
                    ss=ss_classes[idx],
                    chain_id=chain_id,
                    res_num=frag[r_idx].res_num,
                    left_coord=left_coords[idx] if left_coords is not None else None,
                    right_coord=right_coords[idx] if right_coords is not None else None,
                    top_coord=top_coords[idx] if top_coords is not None else None,
                    bottom_coord=bottom_coords[idx] if bottom_coords is not None else None,
                )
                segment_pts.append(pt)

            all_segments.append(segment_pts)
            global_res_counter += frag_len

    return all_segments
