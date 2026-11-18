# termol/render.py
import sys
import copy
import time
import numpy as np
import curses
from pathlib import Path
import importlib.resources as pkg_resources
from typing import Optional, List

from .canvas import Canvas, MoleculeCanvas, ProteinCanvas
from .utilities import (
    get_molecule_data, MoleculeData,
    rotate_points, rotate_molecule_for_screen, scale_for_canvas, rotation_matrix
)
from .backbone import BackbonePoint, generate_backbone_points
from .pdb import is_pdb_code, clean_pdb_code


def scale_protein_for_canvas(backbone_segments: List[List[BackbonePoint]], ligands: list, canvas: Canvas):
    """Scale and center protein backbone segments and ligands to fit canvas dimensions."""
    all_pts = []
    for seg in backbone_segments:
        for pt in seg:
            all_pts.append(pt.coord)
            if pt.left_coord is not None:
                all_pts.append(pt.left_coord)
                all_pts.append(pt.right_coord)
            if pt.top_coord is not None:
                all_pts.append(pt.top_coord)
                all_pts.append(pt.bottom_coord)
    for lig in ligands:
        all_pts.append(lig.coord)

    if not all_pts:
        return

    coords = np.array(all_pts)
    min_x, max_x = np.min(coords[:, 0]), np.max(coords[:, 0])
    min_y, max_y = np.min(coords[:, 1]), np.max(coords[:, 1])

    span_x = max(max_x - min_x, 1e-4)
    span_y = max(max_y - min_y, 1e-4)

    sx = (canvas.width - 4) / span_x
    sy = (canvas.height - 4) / span_y
    scale = 0.85 * min(sx, sy * canvas.aspect_ratio)

    for seg in backbone_segments:
        for pt in seg:
            pt.coord *= scale
            if pt.left_coord is not None:
                pt.left_coord *= scale
                pt.right_coord *= scale
            if pt.top_coord is not None:
                pt.top_coord *= scale
                pt.bottom_coord *= scale
    for lig in ligands:
        lig.coord *= scale


def show_molecule_2D(molecule_data, canvas, name=None, style='ascii', color='depth', rep='trace'):
    """Render a static 2D snapshot of a small molecule or protein."""
    # Ensure stdout handles UTF-8 safely (e.g. for Braille)
    if hasattr(sys.stdout, 'reconfigure'):
        try:
            sys.stdout.reconfigure(encoding='utf-8')
        except Exception:
            pass

    # Check if input is a protein
    if isinstance(molecule_data, MoleculeData) and molecule_data.is_protein:
        protein_canvas = ProteinCanvas(
            canvas.char_width, canvas.char_height,
            canvas.width, canvas.height,
            aspect_ratio=canvas.aspect_ratio,
            style=style, color_mode=color
        )

        active_rep = rep or getattr(molecule_data, 'representation', 'trace')
        if active_rep not in ('trace', 'cartoon', 'ribbon', 'tube'):
            active_rep = 'trace'

        # If representation differs from molecule_data's precomputed segments, regenerate
        if active_rep != getattr(molecule_data, 'representation', 'trace') and molecule_data.structure:
            segments = generate_backbone_points(molecule_data.structure, samples_per_res=8, representation=active_rep)
            center = np.mean([pt.coord for seg in segments for pt in seg], axis=0)
            for seg in segments:
                for pt in seg:
                    pt.coord -= center
                    if pt.left_coord is not None:
                        pt.left_coord -= center
                        pt.right_coord -= center
                    if pt.top_coord is not None:
                        pt.top_coord -= center
                        pt.bottom_coord -= center
        else:
            segments = copy.deepcopy(molecule_data.backbone_segments)

        ligands = copy.deepcopy(molecule_data.ligands)

        scale_protein_for_canvas(segments, ligands, protein_canvas)

        protein_canvas.clear()
        protein_canvas.draw_protein(segments, ligands)

        if name or molecule_data.name:
            display_name = name or molecule_data.name
            res_count = molecule_data.structure.total_residues if molecule_data.structure else len(molecule_data.atom_positions)
            header = f" {display_name} ({res_count} res) [{active_rep.upper()} | {style.upper()} | {color.upper()}] "
            if len(header) > protein_canvas.char_width - 6:
                header = header[:protein_canvas.char_width - 10] + "... "
            print(header.center(protein_canvas.char_width, "="))

        if color == 'none':
            print(protein_canvas.render_plain())
        else:
            print(protein_canvas.render_ansi())
        return

    # Small molecule 2D rendering
    atom_positions, atom_elements, atom_charges, bonds = molecule_data
    atom_positions = atom_positions.copy()
    atom_positions = scale_for_canvas(atom_positions, canvas)
    atom_positions[:, 1] /= canvas.aspect_ratio

    atom_positions_2D = [(pos[0], pos[1]) for pos in atom_positions]

    canvas.clear()
    canvas.draw_molecules(atom_positions_2D, atom_elements, atom_charges, bonds)

    if name:
        header = f" {name} "
        if len(header) > canvas.char_width - 10:
            header = header[:canvas.char_width - 10] + "... "
        print(header.center(canvas.char_width, "="))

    print(canvas)


def show_molecule_3D(stdscr, molecule_data, canvas, name=None, timeout=None,
                     auto_resize=True, style='ascii', color='depth', rep='trace'):
    """Interactive 3D curses renderer for small molecules and macromolecules."""
    curses.curs_set(0)
    stdscr.nodelay(1)
    stdscr.timeout(50)

    # Initialize curses color pairs if supported
    has_curses_colors = curses.has_colors()
    if has_curses_colors:
        try:
            curses.start_color()
            curses.use_default_colors()
            # Secondary Structure pairs:
            curses.init_pair(1, curses.COLOR_MAGENTA, -1)  # Helices
            curses.init_pair(2, curses.COLOR_YELLOW, -1)   # Sheets
            curses.init_pair(3, curses.COLOR_CYAN, -1)     # Loops
            # Rainbow spectrum pairs:
            curses.init_pair(4, curses.COLOR_BLUE, -1)
            curses.init_pair(5, curses.COLOR_CYAN, -1)
            curses.init_pair(6, curses.COLOR_GREEN, -1)
            curses.init_pair(7, curses.COLOR_YELLOW, -1)
            curses.init_pair(8, curses.COLOR_RED, -1)
            # Depth cueing pairs:
            curses.init_pair(9, curses.COLOR_WHITE, -1)    # Foreground (bright)
            curses.init_pair(10, curses.COLOR_CYAN, -1)    # Mid-ground
            curses.init_pair(11, curses.COLOR_BLUE, -1)    # Background (dim)
            # Ligand pair:
            curses.init_pair(12, curses.COLOR_YELLOW, -1)
        except Exception:
            has_curses_colors = False

    is_protein = isinstance(molecule_data, MoleculeData) and molecule_data.is_protein

    if auto_resize:
        term_height, term_width = stdscr.getmaxyx()
        canvas.resize_canvas(term_width - 1, term_height - 1)

    if is_protein:
        active_style = style
        active_color = color
        active_rep = rep or getattr(molecule_data, 'representation', 'trace')
        if active_rep not in ('trace', 'cartoon', 'tube'):
            active_rep = 'cartoon' if active_rep == 'ribbon' else 'trace'

        color_modes = ['depth', 'rainbow', 'ss', 'none']
        styles = ['ascii', 'braille']
        rep_modes = ['trace', 'cartoon', 'tube']

        protein_canvas = ProteinCanvas(
            canvas.char_width, canvas.char_height,
            canvas.width, canvas.height,
            aspect_ratio=canvas.aspect_ratio,
            style=active_style, color_mode=active_color
        )

        def build_segments(target_rep):
            if target_rep == getattr(molecule_data, 'representation', 'trace'):
                segs = copy.deepcopy(molecule_data.backbone_segments)
            else:
                segs = generate_backbone_points(molecule_data.structure, samples_per_res=8, representation=target_rep)
                center = np.mean([pt.coord for s in segs for pt in s], axis=0)
                for s in segs:
                    for pt in s:
                        pt.coord -= center
                        if pt.left_coord is not None:
                            pt.left_coord -= center
                            pt.right_coord -= center
                        if pt.top_coord is not None:
                            pt.top_coord -= center
                            pt.bottom_coord -= center
            scale_protein_for_canvas(segs, ligands, protein_canvas)
            return segs

        ligands = copy.deepcopy(molecule_data.ligands)
        segments = build_segments(active_rep)
    else:
        atom_positions, atom_elements, atom_charges, bonds = molecule_data
        atom_positions = atom_positions.copy()
        atom_positions = scale_for_canvas(atom_positions, canvas)

    end_time = time.time() + timeout if timeout else None

    rotation_axis_map = {
        (ord('A'), ord('a'), curses.KEY_LEFT):    (0, 1, 0),
        (ord('D'), ord('d'), curses.KEY_RIGHT):   (0, -1, 0),
        (ord('W'), ord('w'), curses.KEY_UP):      (1, 0, 0),
        (ord('S'), ord('s'), curses.KEY_DOWN):    (-1, 0, 0),
        (ord('Q'), ord('q')):                     (0, 0, -1),
        (ord('E'), ord('e')):                     (0, 0, 1)
    }
    all_rotation_keys = [key for key_list in rotation_axis_map.keys() for key in key_list]
    rotation_axis = (0, 1, 0)
    rotation_paused = False
    theta = 1.2

    stdscr.clear()

    while True:
        stdscr.erase()
        term_height, term_width = stdscr.getmaxyx()

        if auto_resize:
            canvas.resize_canvas(term_width - 1, term_height - 1)
            if is_protein:
                protein_canvas.resize_canvas(term_width - 1, term_height - 1)

        # Header bar
        header_text = ""
        if is_protein:
            display_name = name or molecule_data.name or "Protein"
            res_num = molecule_data.structure.total_residues if molecule_data.structure else len(molecule_data.atom_positions)
            header_text = f" {display_name} ({res_num} res) [Rep: {active_rep.upper()} (M) | Style: {active_style.upper()} (B) | Color: {active_color.upper()} (C)] "
        elif name:
            header_text = f" {name} "

        if header_text:
            if len(header_text) > canvas.char_width - 6:
                header_text = header_text[:canvas.char_width - 10] + "... "
            header_text = header_text.center(canvas.char_width, "=")
            try:
                stdscr.addstr(0, 0, header_text)
            except curses.error:
                pass

        if is_protein:
            # Rotate protein
            if not rotation_paused:
                R = rotation_matrix(rotation_axis, np.radians(theta))
                for seg in segments:
                    for pt in seg:
                        pt.coord = np.dot(R, pt.coord)
                        if pt.left_coord is not None:
                            pt.left_coord = np.dot(R, pt.left_coord)
                            pt.right_coord = np.dot(R, pt.right_coord)
                        if pt.top_coord is not None:
                            pt.top_coord = np.dot(R, pt.top_coord)
                            pt.bottom_coord = np.dot(R, pt.bottom_coord)
                for lig in ligands:
                    lig.coord = np.dot(R, lig.coord)

            protein_canvas.clear()
            protein_canvas.style = active_style
            protein_canvas.color_mode = active_color
            protein_canvas.draw_protein(segments, ligands)

            # Render to curses buffer
            for cy in range(min(protein_canvas.char_height, term_height - 2)):
                for cx in range(min(protein_canvas.char_width, term_width - 1)):
                    char = protein_canvas.data[cy][cx]
                    if char == ' ':
                        continue

                    attr = curses.A_NORMAL
                    if has_curses_colors and active_color != 'none':
                        z = protein_canvas.z_buffer[cy, cx]
                        if active_color == 'depth':
                            if z > 0:
                                attr = curses.color_pair(9) | curses.A_BOLD
                            elif z > -15:
                                attr = curses.color_pair(10)
                            else:
                                attr = curses.color_pair(11) | curses.A_DIM
                        elif active_color == 'rainbow':
                            color_rgb = protein_canvas.color_data[cy][cx]
                            if color_rgb:
                                r, g, b = color_rgb
                                if b > r and b > g:
                                    attr = curses.color_pair(4) | curses.A_BOLD
                                elif g > r:
                                    attr = curses.color_pair(6)
                                else:
                                    attr = curses.color_pair(8) | curses.A_BOLD
                        elif active_color == 'ss':
                            color_rgb = protein_canvas.color_data[cy][cx]
                            if color_rgb:
                                r, g, b = color_rgb
                                if r > 180 and b > 180:
                                    attr = curses.color_pair(1) | curses.A_BOLD  # Helix
                                elif r > 180 and g > 150:
                                    attr = curses.color_pair(2) | curses.A_BOLD  # Sheet
                                else:
                                    attr = curses.color_pair(3)                  # Loop

                    try:
                        stdscr.addstr(cy + 1, cx, char, attr)
                    except curses.error:
                        pass
        else:
            # Small molecule rotation
            if not rotation_paused:
                atom_positions = rotate_points(atom_positions, rotation_axis, np.radians(theta))

            stretched_positions = atom_positions.copy()
            stretched_positions[:, 1] /= canvas.aspect_ratio
            atom_positions_2D = [(pos[0], pos[1]) for pos in stretched_positions]

            canvas.clear()
            canvas.draw_molecules(atom_positions_2D, atom_elements, atom_charges, bonds)
            try:
                stdscr.addstr(1, 0, str(canvas))
            except curses.error:
                pass

        stdscr.refresh()

        key = stdscr.getch()
        if key != -1:
            if key == curses.KEY_RESIZE:
                continue
            elif key in all_rotation_keys:
                for key_list, axis in rotation_axis_map.items():
                    if key in key_list:
                        rotation_axis = axis
                        rotation_paused = False
            elif key == ord(' '):
                rotation_paused = not rotation_paused
            elif key in [ord('R'), ord('r')]:
                if is_protein:
                    for seg in segments:
                        for pt in seg:
                            pt.coord *= 0.96
                            if pt.left_coord is not None:
                                pt.left_coord *= 0.96
                                pt.right_coord *= 0.96
                            if pt.top_coord is not None:
                                pt.top_coord *= 0.96
                                pt.bottom_coord *= 0.96
                    for lig in ligands:
                        lig.coord *= 0.96
                else:
                    atom_positions *= 0.98
            elif key in [ord('F'), ord('f')]:
                if is_protein:
                    for seg in segments:
                        for pt in seg:
                            pt.coord *= 1.04
                            if pt.left_coord is not None:
                                pt.left_coord *= 1.04
                                pt.right_coord *= 1.04
                            if pt.top_coord is not None:
                                pt.top_coord *= 1.04
                                pt.bottom_coord *= 1.04
                    for lig in ligands:
                        lig.coord *= 1.04
                else:
                    atom_positions *= 1.02
            elif key in [ord('M'), ord('m')] and is_protein:
                # Cycle representation mode (trace -> cartoon -> tube)
                idx = rep_modes.index(active_rep) if active_rep in rep_modes else 0
                active_rep = rep_modes[(idx + 1) % len(rep_modes)]
                segments = build_segments(active_rep)
            elif key in [ord('C'), ord('c')] and is_protein:
                # Cycle color mode
                idx = color_modes.index(active_color)
                active_color = color_modes[(idx + 1) % len(color_modes)]
            elif key in [ord('B'), ord('b')] and is_protein:
                # Toggle style (ascii <-> braille)
                idx = styles.index(active_style)
                active_style = styles[(idx + 1) % len(styles)]
            else:
                break

        if timeout and time.time() > end_time and not rotation_paused:
            break


def draw(input_mol, name=None, width=80, height=40, three_d=True,
         add_hydrogens=False, timeout=None, stdscr=None, auto_resize=True,
         style='ascii', color='depth', rep='trace'):
    """
    Main function for TerMol. Renders small molecules or protein structures in 2D or 3D.
    Inputs:
        input_mol: SMILES string, .sdf/.mol/.pdb file path, or 4-letter RCSB PDB ID (e.g. '1ubq').
        name: Optional title to display.
        width: Canvas width in characters (default 80).
        height: Canvas height in characters (default 40).
        three_d: Whether to display interactive 3D viewer (default True).
        style: For proteins: 'ascii' (default density ramp) or 'braille' (subpixel dots).
        color: For proteins: 'depth' (default luminance ramp), 'rainbow' (N->C), 'ss' (secondary structure), or 'none'.
        rep: For proteins: 'trace' (default 1D line), 'cartoon'/'ribbon' (flat ribbons and arrows), or 'tube'.
    """
    if three_d:
        curses.wrapper(draw_persistent, input_mol, name=name, width=width, height=height,
                       three_d=three_d, add_hydrogens=add_hydrogens, timeout=timeout,
                       auto_resize=auto_resize, style=style, color=color, rep=rep)
    else:
        draw_persistent(None, input_mol, name=name, width=width, height=height,
                        three_d=three_d, add_hydrogens=add_hydrogens, timeout=timeout,
                        auto_resize=auto_resize, style=style, color=color, rep=rep)


def draw_persistent(stdscr, input_mol, name=None, width=80, height=40,
                    three_d=True, add_hydrogens=False, timeout=None,
                    auto_resize=True, style='ascii', color='depth', rep='trace'):
    """Persistent curses wrapper for draw."""
    molecule_data = get_molecule_data(input_mol, three_d=three_d, add_hydrogens=add_hydrogens, representation=rep)
    canvas = MoleculeCanvas(width, height, float(width), float(height), aspect_ratio=2.0)

    if three_d:
        if stdscr:
            show_molecule_3D(stdscr, molecule_data, canvas, name=name, timeout=timeout,
                             auto_resize=auto_resize, style=style, color=color, rep=rep)
        else:
            curses.wrapper(show_molecule_3D, molecule_data, canvas, name=name,
                           timeout=timeout, auto_resize=auto_resize, style=style, color=color, rep=rep)
    else:
        show_molecule_2D(molecule_data, canvas, name=name, style=style, color=color, rep=rep)


def draw_multi(input_mols, names=None, width=80, height=40, three_d=True,
               add_hydrogens=False, timeout=None, auto_resize=True,
               style='ascii', color='depth', rep='trace'):
    """Render multiple molecules or proteins in sequence."""
    if names and len(names) != len(input_mols):
        raise ValueError("The number of names must match the number of inputs!")

    if names is None:
        names = []
        for i, m in enumerate(input_mols):
            if isinstance(m, str) and is_pdb_code(m):
                names.append(clean_pdb_code(m))
            elif isinstance(m, (str, Path)) and Path(m).exists():
                names.append(Path(m).stem)
            else:
                names.append(f"Molecule {i + 1}")

    if not three_d:
        for i in range(len(input_mols)):
            draw(input_mols[i], name=names[i], width=width, height=height,
                 three_d=three_d, add_hydrogens=add_hydrogens, style=style, color=color, rep=rep)
        return

    def loop_molecules(stdscr, molecules, mol_names, loop_timeout=None):
        for i in range(len(molecules)):
            try:
                draw_persistent(stdscr, molecules[i], name=mol_names[i], timeout=loop_timeout,
                                width=width, height=height, three_d=three_d,
                                add_hydrogens=add_hydrogens, auto_resize=auto_resize,
                                style=style, color=color, rep=rep)
            except Exception:
                continue

    curses.wrapper(loop_molecules, input_mols, names, loop_timeout=timeout)


def showcase(timeout=5, width=80, height=40):
    """Showcase TerMol's capabilities with a random gallery."""
    smiles_dict = {}
    with pkg_resources.open_text(__package__, 'smiles_1000.csv') as file:
        lines = file.readlines()
        for line in lines:
            split_line = line.split('\t')
            if len(split_line) != 2:
                continue
            name, smiles = split_line
            smiles_dict[name] = smiles

    molecules = list(smiles_dict.values())
    names = list(smiles_dict.keys())
    draw_multi(molecules, names, width=width, height=height, timeout=timeout, auto_resize=True)
