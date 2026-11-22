"""
termol/diagnostics.py
Performance profiling and diagnostics module for TerMol.
Measures per-frame latency, pipeline stage bottlenecks, and rasterization substeps
for small molecules and proteins without affecting core rendering code.
"""

import os
import sys
import time
import copy
import io
import argparse
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Tuple
import numpy as np

from termol.utilities import get_molecule_data, scale_for_canvas, MoleculeData
from termol.canvas import MoleculeCanvas, ProteinCanvas
from termol.render import (
    scale_protein_for_canvas,
    rotation_matrix,
    rotate_points,
    read_available_input,
    generate_backbone_points,
)


@dataclass
class ProfileResult:
    """Diagnostic profile results for a molecular rendering pipeline."""
    input_name: str
    mol_type: str
    representation: str
    style: str
    color: str
    num_frames: int
    total_frame_ms: float
    fps: float
    compute_ms: float
    max_theoretical_fps: float
    stage_times_ms: Dict[str, float] = field(default_factory=dict)
    stage_percentages: Dict[str, float] = field(default_factory=dict)
    compute_percentages: Dict[str, float] = field(default_factory=dict)
    substep_times_ms: Optional[Dict[str, float]] = None


def profile_pipeline(
    input_item: Any,
    num_frames: int = 60,
    representation: str = 'trace',
    style: str = 'high-res',
    color: Optional[str] = None,
    width: int = 80,
    height: int = 40,
    sleep_throttle: float = 0.025,
    warmup: int = 5,
    profile_substeps: bool = True
) -> ProfileResult:
    """
    Profile a molecule or protein through the complete interactive rendering pipeline.
    Measures each stage: rotation, canvas clearing, drawing, ANSI string formatting,
    frame assembly, I/O, sleep throttle, and input reading.
    """
    mol_data = get_molecule_data(input_item, three_d=True, representation=representation)
    is_protein = getattr(mol_data, 'is_protein', False)
    display_name = getattr(mol_data, 'name', None) or (str(input_item)[:20] if isinstance(input_item, str) else "Molecule")

    active_style = style or 'high-res'
    if is_protein:
        active_color = color or 'depth'
        canvas = ProteinCanvas(width, height, float(width), float(height), aspect_ratio=2.0, style=active_style, color_mode=active_color)
        segs = copy.deepcopy(mol_data.backbone_segments)
        ligs = copy.deepcopy(mol_data.ligands)
        scale_protein_for_canvas(segs, ligs, canvas)
        R = rotation_matrix((0, 1, 0), np.radians(1.2))
    else:
        active_color = color or 'cpk'
        canvas = MoleculeCanvas(width, height, float(width), float(height), aspect_ratio=2.0, style=active_style, color_mode=active_color)
        atom_positions, atom_elements, atom_charges, bonds = mol_data
        atom_positions = atom_positions.copy()
        atom_positions = scale_for_canvas(atom_positions, canvas)
        bond_orders = getattr(mol_data, 'bond_orders', None)

    stage_records: Dict[str, List[float]] = {
        'rot': [],
        'clear': [],
        'draw': [],
        'render_ansi': [],
        'frame_assemble': [],
        'out': [],
        'sleep': [],
        'input': [],
        'total_frame': []
    }

    sink = io.StringIO()
    total_iterations = num_frames + warmup

    for idx in range(total_iterations):
        t_frame_start = time.perf_counter_ns()

        # 1. Geometry rotation & coordinate transforms
        t0 = time.perf_counter_ns()
        if is_protein:
            for seg in segs:
                for pt in seg:
                    pt.coord = np.dot(R, pt.coord)
                    if pt.left_coord is not None:
                        pt.left_coord = np.dot(R, pt.left_coord)
                        pt.right_coord = np.dot(R, pt.right_coord)
                    if pt.top_coord is not None:
                        pt.top_coord = np.dot(R, pt.top_coord)
                        pt.bottom_coord = np.dot(R, pt.bottom_coord)
            for lig in ligs:
                lig.coord = np.dot(R, lig.coord)
        else:
            atom_positions = rotate_points(atom_positions, (0, 1, 0), np.radians(1.2))
            stretched = atom_positions.copy()
            stretched[:, 1] /= canvas.aspect_ratio
            pos_3D = [(p[0], p[1], p[2]) for p in stretched]
        stage_records['rot'].append((time.perf_counter_ns() - t0) / 1e6)

        # 2. Canvas buffer clearing
        t0 = time.perf_counter_ns()
        canvas.clear()
        stage_records['clear'].append((time.perf_counter_ns() - t0) / 1e6)

        # 3. Rasterization & drawing
        t0 = time.perf_counter_ns()
        if is_protein:
            canvas.draw_protein(segs, ligs)
        else:
            canvas.draw_molecules(pos_3D, atom_elements, atom_charges, bonds, bond_orders=bond_orders, style=active_style, color_mode=active_color)
        stage_records['draw'].append((time.perf_counter_ns() - t0) / 1e6)

        # 4. Canvas ANSI string serialization
        t0 = time.perf_counter_ns()
        body = canvas.render_ansi()
        stage_records['render_ansi'].append((time.perf_counter_ns() - t0) / 1e6)

        # 5. Frame assembly (header formatting, FPS, string joining)
        t0 = time.perf_counter_ns()
        header_text = f" {display_name} [{active_style.upper()}] "
        header_line = header_text.center(canvas.char_width, "=")
        fps_line = "[30.0 FPS]"
        frame = f"\033[H\033[1;36m{header_line}\033[0m\n{fps_line}\n{body}"
        stage_records['frame_assemble'].append((time.perf_counter_ns() - t0) / 1e6)

        # 6. Terminal I/O simulation (in-memory write & flush)
        t0 = time.perf_counter_ns()
        sink.seek(0)
        sink.write(frame)
        sink.flush()
        stage_records['out'].append((time.perf_counter_ns() - t0) / 1e6)

        # 7. Sleep throttle (replicates real event loop timer)
        t0 = time.perf_counter_ns()
        if sleep_throttle > 0:
            time.sleep(sleep_throttle)
        stage_records['sleep'].append((time.perf_counter_ns() - t0) / 1e6)

        # 8. Non-blocking input polling
        t0 = time.perf_counter_ns()
        read_available_input()
        stage_records['input'].append((time.perf_counter_ns() - t0) / 1e6)

        stage_records['total_frame'].append((time.perf_counter_ns() - t_frame_start) / 1e6)

    # Discard warmup frames
    clean_records = {k: v[warmup:] for k, v in stage_records.items()}
    avg_times = {k: float(np.mean(v)) for k, v in clean_records.items()}
    total_frame_ms = avg_times['total_frame']
    sleep_ms = avg_times['sleep']
    compute_ms = max(0.001, total_frame_ms - sleep_ms)

    # Percentage of total live frametime
    stage_pcts = {k: (avg_times[k] / total_frame_ms) * 100.0 for k in avg_times if k != 'total_frame'}

    # Percentage of computational frametime (excluding sleep)
    compute_keys = ['draw', 'render_ansi', 'rot', 'clear', 'frame_assemble', 'input', 'out']
    compute_pcts = {k: (avg_times[k] / compute_ms) * 100.0 for k in compute_keys}

    fps = 1000.0 / total_frame_ms if total_frame_ms > 0 else 0.0
    max_fps = 1000.0 / compute_ms if compute_ms > 0 else 0.0

    substeps = None
    if profile_substeps:
        substeps = profile_rasterization_substeps(
            input_item,
            representation=representation,
            style=active_style,
            color=active_color,
            width=width,
            height=height,
            num_frames=num_frames,
            warmup=warmup
        )

    return ProfileResult(
        input_name=display_name,
        mol_type="Protein" if is_protein else "Small Molecule",
        representation=representation if is_protein else "N/A",
        style=active_style,
        color=active_color,
        num_frames=num_frames,
        total_frame_ms=total_frame_ms,
        fps=fps,
        compute_ms=compute_ms,
        max_theoretical_fps=max_fps,
        stage_times_ms=avg_times,
        stage_percentages=stage_pcts,
        compute_percentages=compute_pcts,
        substep_times_ms=substeps
    )


def profile_rasterization_substeps(
    input_item: Any,
    representation: str = 'trace',
    style: str = 'high-res',
    color: Optional[str] = None,
    width: int = 80,
    height: int = 40,
    num_frames: int = 50,
    warmup: int = 5
) -> Dict[str, float]:
    """
    Micro-profile the internal sub-steps inside canvas.draw_protein or canvas.draw_molecules.
    Identifies whether time is spent in depth calculations, line Bresenham, bitwise operations,
    or Braille character finalization.
    """
    mol_data = get_molecule_data(input_item, three_d=True, representation=representation)
    is_protein = getattr(mol_data, 'is_protein', False)

    substep_times: Dict[str, List[float]] = {}

    if is_protein:
        canvas = ProteinCanvas(width, height, float(width), float(height), aspect_ratio=2.0, style=style, color_mode=color or 'depth')
        segs = copy.deepcopy(mol_data.backbone_segments)
        ligs = copy.deepcopy(mol_data.ligands)
        scale_protein_for_canvas(segs, ligs, canvas)

        substep_times = {'z_range_scan': [], 'line_rasterization': [], 'finalize_grid': [], 'ligands': []}

        for _ in range(num_frames + warmup):
            # 1. Depth range scan
            t0 = time.perf_counter_ns()
            all_zs = []
            for segment in segs:
                for pt in segment:
                    all_zs.append(pt.coord[2])
                    if pt.left_coord is not None:
                        all_zs.append(pt.left_coord[2])
                        all_zs.append(pt.right_coord[2])
                    if pt.top_coord is not None:
                        all_zs.append(pt.top_coord[2])
                        all_zs.append(pt.bottom_coord[2])
            z_min = float(np.min(all_zs)) if all_zs else -10.0
            z_max = float(np.max(all_zs)) if all_zs else 10.0
            z_range = max(z_max - z_min, 1e-4)
            substep_times['z_range_scan'].append((time.perf_counter_ns() - t0) / 1e6)

            # 2. Line rasterization
            canvas.clear()
            t0 = time.perf_counter_ns()
            for segment in segs:
                if not segment:
                    continue
                has_ribbon = segment[0].left_coord is not None
                if not has_ribbon:
                    for i in range(len(segment)):
                        pt = segment[i]
                        x, y, z = pt.coord
                        prev_pt = segment[i - 1] if i > 0 else pt
                        px, py, pz = prev_pt.coord
                        canvas._draw_line_braille(px, py, pz, x, y, z, pt.norm_pos, pt.ss, z_min, z_range)
                else:
                    for i in range(len(segment)):
                        pt = segment[i]
                        lx, ly, lz = pt.left_coord
                        rx, ry, rz = pt.right_coord
                        canvas._draw_line_braille(lx, ly, lz, rx, ry, rz, pt.norm_pos, pt.ss, z_min, z_range)
                        if i > 0:
                            prev = segment[i - 1]
                            plx, ply, plz = prev.left_coord
                            prx, pry, prz = prev.right_coord
                            canvas._draw_line_braille(plx, ply, plz, lx, ly, lz, pt.norm_pos, pt.ss, z_min, z_range)
                            canvas._draw_line_braille(prx, pry, prz, rx, ry, rz, pt.norm_pos, pt.ss, z_min, z_range)
                        if pt.top_coord is not None and pt.bottom_coord is not None:
                            tx, ty, tz = pt.top_coord
                            bx, by, bz = pt.bottom_coord
                            canvas._draw_line_braille(tx, ty, tz, bx, by, bz, pt.norm_pos, pt.ss, z_min, z_range)
                            if i > 0 and prev.top_coord is not None:
                                ptx, pty, ptz = prev.top_coord
                                pbx, pby, pbz = prev.bottom_coord
                                canvas._draw_line_braille(ptx, pty, ptz, tx, ty, tz, pt.norm_pos, pt.ss, z_min, z_range)
                                canvas._draw_line_braille(pbx, pby, pbz, bx, by, bz, pt.norm_pos, pt.ss, z_min, z_range)
            substep_times['line_rasterization'].append((time.perf_counter_ns() - t0) / 1e6)

            # 3. Finalize grid
            t0 = time.perf_counter_ns()
            for cy in range(canvas.char_height):
                for cx in range(canvas.char_width):
                    bits = canvas.braille_bits[cy, cx]
                    if bits > 0:
                        canvas.data[cy][cx] = chr(0x2800 + bits)
                    else:
                        canvas.data[cy][cx] = ' '
            substep_times['finalize_grid'].append((time.perf_counter_ns() - t0) / 1e6)

            # 4. Ligands
            t0 = time.perf_counter_ns()
            if ligs:
                span_x = canvas.max_x - canvas.min_x
                span_y = canvas.max_y - canvas.min_y
                scale_x = (canvas.char_width - 1) / span_x
                scale_y = (canvas.char_height - 1) / (span_y * canvas.aspect_ratio)
                mid_cx = (canvas.char_width - 1) / 2.0
                mid_cy = (canvas.char_height - 1) / 2.0
                for lig in ligs:
                    lx, ly, lz = lig.coord
                    cx = int(mid_cx + lx * scale_x + 0.5)
                    cy = int(mid_cy + ly * scale_y + 0.5)
                    canvas.set_pixel_cell(cx, cy, lig.element[:1].upper(), z=lz + 0.2)
            substep_times['ligands'].append((time.perf_counter_ns() - t0) / 1e6)

    else:
        # Small molecule
        canvas = MoleculeCanvas(width, height, float(width), float(height), aspect_ratio=2.0, style=style, color_mode=color or 'cpk')
        atom_positions, atom_elements, atom_charges, bonds = mol_data
        atom_positions = scale_for_canvas(atom_positions, canvas)
        bond_orders = getattr(mol_data, 'bond_orders', None)
        stretched = atom_positions.copy()
        stretched[:, 1] /= canvas.aspect_ratio
        pos_3D = [(p[0], p[1], p[2]) for p in stretched]

        substep_times = {'bonds_braille': [], 'finalize_grid': [], 'atoms_braille': []}

        for _ in range(num_frames + warmup):
            canvas.clear()
            # 1. Bonds braille
            t0 = time.perf_counter_ns()
            canvas._draw_bonds_braille(pos_3D, bonds, atom_elements=atom_elements, bond_orders=bond_orders, color_mode=canvas.color_mode)
            substep_times['bonds_braille'].append((time.perf_counter_ns() - t0) / 1e6)

            # 2. Finalize
            t0 = time.perf_counter_ns()
            canvas._finalize_braille_grid()
            substep_times['finalize_grid'].append((time.perf_counter_ns() - t0) / 1e6)

            # 3. Atoms
            t0 = time.perf_counter_ns()
            canvas._draw_atoms_braille(pos_3D, atom_elements, atom_charges, bonds, color_mode=canvas.color_mode)
            substep_times['atoms_braille'].append((time.perf_counter_ns() - t0) / 1e6)

    # Return averaged substep times excluding warmup
    return {k: float(np.mean(v[warmup:])) for k, v in substep_times.items()}


def format_diagnostic_report(results: List[ProfileResult]) -> str:
    """Format a list of ProfileResult objects into a structured diagnostic report."""
    lines = []
    lines.append("=" * 80)
    lines.append(" TerMol Rendering Performance & Frametime Diagnostics ".center(80, "="))
    lines.append("=" * 80)
    lines.append("")

    for res in results:
        lines.append(f"Target: {res.input_name} ({res.mol_type})")
        lines.append(f"  Representation : {res.representation}")
        lines.append(f"  Style / Color  : {res.style} / {res.color}")
        lines.append(f"  Live Frametime : {res.total_frame_ms:6.2f} ms ({res.fps:4.1f} FPS with 25ms sleep throttle)")
        lines.append(f"  Pure Compute   : {res.compute_ms:6.2f} ms ({res.max_theoretical_fps:5.1f} Max Theoretical FPS)")
        lines.append("")

        lines.append("  [End-to-End Pipeline Breakdown]")
        lines.append(f"    {'Stage':<16} {'Time (ms)':>10} {'% Live Frame':>14} {'% Pure Compute':>16}")
        lines.append(f"    {'-'*16} {'-'*10} {'-'*14} {'-'*16}")

        # Order stages logically
        stage_order = ['sleep', 'draw', 'render_ansi', 'rot', 'clear', 'frame_assemble', 'input', 'out']
        for stage in stage_order:
            if stage not in res.stage_times_ms:
                continue
            t_ms = res.stage_times_ms[stage]
            pct_live = res.stage_percentages.get(stage, 0.0)
            pct_comp = res.compute_percentages.get(stage, 0.0) if stage != 'sleep' else 0.0
            comp_str = f"{pct_comp:5.1f}%" if stage != 'sleep' else "N/A (Throttle)"
            lines.append(f"    {stage:<16} {t_ms:9.2f} ms {pct_live:13.1f}% {comp_str:>16}")

        if res.substep_times_ms:
            lines.append("")
            lines.append("  [Rasterization Substep Breakdown]")
            lines.append(f"    {'Substep':<24} {'Time (ms)':>10} {'% Raster Time':>15}")
            lines.append(f"    {'-'*24} {'-'*10} {'-'*15}")
            tot_raster = sum(res.substep_times_ms.values())
            for sub, sub_ms in res.substep_times_ms.items():
                pct_sub = (sub_ms / tot_raster * 100.0) if tot_raster > 0 else 0.0
                lines.append(f"    {sub:<24} {sub_ms:9.2f} ms {pct_sub:14.1f}%")

        lines.append("-" * 80)

    # Diagnostic summary
    lines.append("")
    lines.append("Key Diagnostic Takeaways:")
    for res in results:
        top_compute = max(res.compute_percentages.items(), key=lambda x: x[1])
        sleep_pct = res.stage_percentages.get('sleep', 0.0)
        lines.append(
            f"  * {res.input_name} ({res.representation}): "
            f"Throttle sleep accounts for {sleep_pct:.1f}% of live frame time. "
            f"Within computation ({res.compute_ms:.2f} ms), '{top_compute[0]}' is the largest consumer ({top_compute[1]:.1f}%)."
        )
    lines.append("")
    return "\n".join(lines)


def run_diagnostics(targets: Optional[List[Any]] = None, num_frames: int = 50, verbose: bool = True) -> List[ProfileResult]:
    """Execute complete diagnostics on default or specified targets."""
    if targets is None:
        targets = [
            ("CC(=O)Oc1ccccc1C(=O)O", "trace"),
            ("1UBQ", "trace"),
            ("1UBQ", "cartoon"),
        ]

    results = []
    for item in targets:
        if isinstance(item, tuple):
            inp, rep = item
        else:
            inp, rep = item, "trace"

        res = profile_pipeline(inp, num_frames=num_frames, representation=rep)
        results.append(res)

    if verbose:
        report = format_diagnostic_report(results)
        print(report)

    return results


def main():
    parser = argparse.ArgumentParser(description="TerMol rendering performance & frametime diagnostics.")
    parser.add_argument('inputs', nargs='*', default=[], help="Molecules or PDB IDs to profile (default: Aspirin and 1UBQ).")
    parser.add_argument('--frames', type=int, default=50, help="Number of benchmark frames per input (default: 50).")
    parser.add_argument('--rep', default=None, choices=['trace', 'cartoon', 'tube', 'ribbon'], help="Override protein representation.")
    parser.add_argument('--style', default='high-res', help="Rendering style (default: high-res).")
    parser.add_argument('--no-sleep', action='store_true', help="Disable 25ms frame sleep throttle during benchmark.")
    args = parser.parse_args()

    throttle = 0.0 if args.no_sleep else 0.025

    if not args.inputs:
        targets = [
            ("CC(=O)Oc1ccccc1C(=O)O", "trace"),
            ("1UBQ", "trace"),
            ("1UBQ", "cartoon")
        ]
        if args.rep:
            targets = [("CC(=O)Oc1ccccc1C(=O)O", "trace"), ("1UBQ", args.rep)]
    else:
        targets = [(inp, args.rep or "trace") for inp in args.inputs]

    results = []
    for inp, rep in targets:
        res = profile_pipeline(inp, num_frames=args.frames, representation=rep, style=args.style, sleep_throttle=throttle)
        results.append(res)

    print(format_diagnostic_report(results))


if __name__ == '__main__':
    main()
