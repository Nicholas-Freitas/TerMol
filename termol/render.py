# termol/render.py
import sys
import os
import re
import copy
import time
import ctypes
import traceback
import numpy as np
import curses
from pathlib import Path
import importlib.resources as pkg_resources
from typing import Optional, List, Tuple

from .canvas import Canvas, MoleculeCanvas, ProteinCanvas
from .utilities import (
    get_molecule_data, MoleculeData,
    rotate_points, rotate_molecule_for_screen, scale_for_canvas, rotation_matrix
)
from .backbone import BackbonePoint, generate_backbone_points
from .pdb import is_pdb_code, clean_pdb_code


def scale_protein_for_canvas(backbone_segments: List[List[BackbonePoint]], ligands: list, canvas: Canvas) -> float:
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
        return 1.0

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

    return float(scale)


def show_molecule_2D(molecule_data, canvas, name=None, style=None, color=None, rep='trace'):
    """Render a static 2D snapshot of a small molecule or protein."""
    # Ensure stdout handles UTF-8 safely (e.g. for Braille/Unicode glyphs)
    if hasattr(sys.stdout, 'reconfigure'):
        try:
            sys.stdout.reconfigure(encoding='utf-8')
        except Exception:
            pass

    # Check if input is a protein
    if isinstance(molecule_data, MoleculeData) and molecule_data.is_protein:
        active_style = style or 'braille'
        active_color = color or 'depth'
        protein_canvas = ProteinCanvas(
            canvas.char_width, canvas.char_height,
            canvas.width, canvas.height,
            aspect_ratio=canvas.aspect_ratio,
            style=active_style, color_mode=active_color
        )

        active_rep = rep or getattr(molecule_data, 'representation', 'trace')
        if active_rep not in ('trace', 'cartoon', 'ribbon', 'tube'):
            active_rep = 'trace'

        # If representation differs from molecule_data's precomputed segments, regenerate
        if active_rep != getattr(molecule_data, 'representation', 'trace') and molecule_data.structure:
            segments = generate_backbone_points(molecule_data.structure, samples_per_res=8, representation=active_rep)
        else:
            segments = copy.deepcopy(molecule_data.backbone_segments)

        ligands = copy.deepcopy(molecule_data.ligands)

        scale_protein_for_canvas(segments, ligands, protein_canvas)

        protein_canvas.clear()
        protein_canvas.draw_protein(segments, ligands)

        if name or molecule_data.name:
            display_name = name or molecule_data.name
            res_count = molecule_data.structure.total_residues if molecule_data.structure else len(molecule_data.atom_positions)
            header = f" {display_name} ({res_count} res) [{active_rep.upper()} | {active_style.upper()} | {active_color.upper()}] "
            if len(header) > protein_canvas.char_width - 6:
                header = header[:protein_canvas.char_width - 10] + "... "
            print(header.center(protein_canvas.char_width, "="))

        if active_color == 'none':
            print(protein_canvas.render_plain())
        else:
            print(protein_canvas.render_ansi())
        return

    # Small molecule 2D rendering
    active_style = style or 'symbols'
    active_color = color or 'cpk'

    atom_positions, atom_elements, atom_charges, bonds = molecule_data
    atom_positions = atom_positions.copy()
    atom_positions = scale_for_canvas(atom_positions, canvas)
    atom_positions[:, 1] /= canvas.aspect_ratio

    atom_positions_2D = [(pos[0], pos[1], 0.0) for pos in atom_positions]
    bond_orders = getattr(molecule_data, 'bond_orders', None)

    canvas.clear()
    canvas.draw_molecules(atom_positions_2D, atom_elements, atom_charges, bonds,
                          bond_orders=bond_orders, style=active_style, color_mode=active_color)

    display_name = name or getattr(molecule_data, 'name', None)
    if display_name:
        header = f" {display_name} [{active_style.upper()}] "
        if len(header) > canvas.char_width - 10:
            header = header[:canvas.char_width - 10] + "... "
        if active_color == 'none':
            print(header.center(canvas.char_width, "="))
        else:
            print(f"\033[1;36m{header.center(canvas.char_width, '=')}\033[0m")

    if active_color == 'none':
        print(canvas.render_plain())
    else:
        print(canvas.render_ansi())


_RGB_256_CACHE = {}


def rgb_to_256(r: int, g: int, b: int) -> int:
    """Map 24-bit RGB (0-255) to the closest xterm-256 color index (16-255)."""
    r_c = max(0, min(255, int(r)))
    g_c = max(0, min(255, int(g)))
    b_c = max(0, min(255, int(b)))
    key = (r_c, g_c, b_c)
    cached = _RGB_256_CACHE.get(key)
    if cached is not None:
        return cached

    # Grayscale check: if all channels are within 8 units, map to grayscale ramp
    if abs(r_c - g_c) <= 8 and abs(g_c - b_c) <= 8 and abs(r_c - b_c) <= 8:
        gray = (r_c + g_c + b_c) // 3
        if gray < 4:
            res = 16  # Pure black in 6x6x6 cube
        elif gray > 246:
            res = 231  # Pure white in 6x6x6 cube
        else:
            # Grayscale ramp 232..255 spans 8..238 in steps of 10
            idx = 232 + int(round((gray - 8) / 10.0))
            res = max(232, min(255, idx))
        if len(_RGB_256_CACHE) > 4096:
            _RGB_256_CACHE.clear()
        _RGB_256_CACHE[key] = res
        return res

    # 6x6x6 color cube: indices 16..231
    def to_6cube(v: int) -> int:
        if v < 48:
            return 0
        elif v < 115:
            return 1
        elif v < 155:
            return 2
        elif v < 195:
            return 3
        elif v < 235:
            return 4
        return 5

    res = 16 + 36 * to_6cube(r_c) + 6 * to_6cube(g_c) + to_6cube(b_c)
    if len(_RGB_256_CACHE) > 4096:
        _RGB_256_CACHE.clear()
    _RGB_256_CACHE[key] = res
    return res


def show_molecule_3D(stdscr, molecule_data, canvas, name=None, timeout=None,
                     auto_resize=True, style='braille', color='depth', rep='trace'):
    """Interactive 3D curses renderer for small molecules and macromolecules."""
    curses.curs_set(0)
    stdscr.nodelay(1)
    stdscr.timeout(50)

    # Initialize curses color pairs if supported
    has_curses_colors = curses.has_colors()
    use_256 = False
    num_colors = 0
    num_pairs = 0
    pair_map = {}
    next_pair_idx = 1

    if has_curses_colors:
        try:
            curses.start_color()
            try:
                curses.use_default_colors()
                bg_color = -1
            except Exception:
                bg_color = curses.COLOR_BLACK

            num_colors = getattr(curses, 'COLORS', 8)
            num_pairs = getattr(curses, 'COLOR_PAIRS', 64)
            # On Windows, PDCurses wincon uses a 4-bit foreground / 4-bit background console attribute.
            # Passing color indices >= 16 overflows into the background attribute nibble,
            # causing character backgrounds to turn gray/colored.
            # True 256-color curses pairs are only supported on non-Windows ncurses platforms.
            use_256 = (sys.platform != 'win32') and (num_colors >= 256) and (num_pairs >= 256)

            # Always initialize the standard 16-color palette (pairs 1..14 with indices 0..7)
            curses.init_pair(1, curses.COLOR_MAGENTA, bg_color)  # Helices
            curses.init_pair(2, curses.COLOR_YELLOW, bg_color)   # Sheets
            curses.init_pair(3, curses.COLOR_CYAN, bg_color)     # Loops
            curses.init_pair(4, curses.COLOR_BLUE, bg_color)     # Rainbow blue
            curses.init_pair(5, curses.COLOR_CYAN, bg_color)     # Rainbow cyan
            curses.init_pair(6, curses.COLOR_GREEN, bg_color)    # Rainbow green
            curses.init_pair(7, curses.COLOR_YELLOW, bg_color)   # Rainbow yellow
            curses.init_pair(8, curses.COLOR_RED, bg_color)      # Rainbow red
            curses.init_pair(9, curses.COLOR_WHITE, bg_color)    # Depth white
            curses.init_pair(10, curses.COLOR_CYAN, bg_color)    # Depth cyan
            curses.init_pair(11, curses.COLOR_BLUE, bg_color)    # Depth blue
            curses.init_pair(12, curses.COLOR_YELLOW, bg_color)  # Ligand / Twilight yellow
            curses.init_pair(13, curses.COLOR_MAGENTA, bg_color) # Twilight magenta
            curses.init_pair(14, curses.COLOR_BLACK, bg_color)   # Monochrome charcoal / dark

            if use_256:
                max_init = min(256, num_pairs)
                for c in range(1, max_init):
                    try:
                        curses.init_pair(c, c, bg_color)
                        pair_map[c] = c
                    except Exception:
                        break
                next_pair_idx = len(pair_map) + 1
        except Exception:
            has_curses_colors = False
            use_256 = False

    def get_color_pair(color_idx: int) -> int:
        nonlocal next_pair_idx
        if color_idx in pair_map:
            return pair_map[color_idx]
        if next_pair_idx < num_pairs:
            try:
                curses.init_pair(next_pair_idx, color_idx, bg_color)
                pair_map[color_idx] = next_pair_idx
                pid = next_pair_idx
                next_pair_idx += 1
                return pid
            except Exception:
                pass
        return pair_map.get(color_idx % max(1, len(pair_map)), 0)

    is_protein = isinstance(molecule_data, MoleculeData) and molecule_data.is_protein

    if auto_resize:
        term_height, term_width = stdscr.getmaxyx()
        canvas.resize_canvas(term_width - 1, term_height - 1)

    if is_protein:
        active_style = style
        active_color = color
        if active_color == 'depth':
            active_color = 'depth1'
        active_rep = rep or getattr(molecule_data, 'representation', 'trace')
        if active_rep not in ('trace', 'cartoon', 'tube'):
            active_rep = 'cartoon' if active_rep == 'ribbon' else 'trace'

        color_modes = ['depth1', 'depth2', 'depth3', 'rainbow', 'ss', 'none']
        styles = ['ascii', 'braille']
        rep_modes = ['trace', 'cartoon', 'tube']

        protein_canvas = ProteinCanvas(
            canvas.char_width, canvas.char_height,
            canvas.width, canvas.height,
            aspect_ratio=canvas.aspect_ratio,
            style=active_style, color_mode=active_color
        )

        initial_rep = active_rep
        ligands = copy.deepcopy(molecule_data.ligands)
        canonical_cache = {}

        # Scale initial representation and ligands once to fit canvas
        initial_segs = copy.deepcopy(molecule_data.backbone_segments)
        protein_scale = scale_protein_for_canvas(initial_segs, ligands, protein_canvas)
        canonical_cache[initial_rep] = initial_segs

        # Track cumulative interactive rotation matrix and zoom factor
        cum_R = np.eye(3)
        zoom_scale = 1.0

        def get_segments(target_rep):
            if target_rep not in canonical_cache:
                raw_segs = generate_backbone_points(molecule_data.structure, samples_per_res=8, representation=target_rep)
                for seg in raw_segs:
                    for pt in seg:
                        pt.coord *= protein_scale
                        if pt.left_coord is not None:
                            pt.left_coord *= protein_scale
                            pt.right_coord *= protein_scale
                        if pt.top_coord is not None:
                            pt.top_coord *= protein_scale
                            pt.bottom_coord *= protein_scale
                canonical_cache[target_rep] = raw_segs

            target_segs = copy.deepcopy(canonical_cache[target_rep])
            for seg in target_segs:
                for pt in seg:
                    pt.coord = np.dot(cum_R, pt.coord * zoom_scale)
                    if pt.left_coord is not None:
                        pt.left_coord = np.dot(cum_R, pt.left_coord * zoom_scale)
                        pt.right_coord = np.dot(cum_R, pt.right_coord * zoom_scale)
                    if pt.top_coord is not None:
                        pt.top_coord = np.dot(cum_R, pt.top_coord * zoom_scale)
                        pt.bottom_coord = np.dot(cum_R, pt.bottom_coord * zoom_scale)
            return target_segs

        segments = get_segments(active_rep)
    else:
        active_style = style or 'symbols'
        active_color = color or 'cpk'
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
            color_cap = "256c" if use_256 else ("16c" if sys.platform == 'win32' else f"{num_colors}c")
            header_text = f" {display_name} ({res_num} res) [Rep: {active_rep.upper()} (M) | Style: {active_style.upper()} (B) | Color: {active_color.upper()} (C) | {color_cap}] "
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
                cum_R = np.dot(R, cum_R)
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
                        color_rgb = protein_canvas.color_data[cy][cx]
                        if use_256 and color_rgb is not None:
                            c256 = rgb_to_256(color_rgb[0], color_rgb[1], color_rgb[2])
                            attr = curses.color_pair(get_color_pair(c256))
                        else:
                            # 8/16-color fallback with dynamic depth normalization
                            z = protein_canvas.z_buffer[cy, cx]
                            z_range = getattr(protein_canvas, 'last_z_range', 1.0)
                            z_min = getattr(protein_canvas, 'last_z_min', 0.0)
                            z_norm = (z - z_min) / z_range if z_range > 0 else 0.5
                            if z_norm < 0.0:
                                z_norm = 0.0
                            elif z_norm > 1.0:
                                z_norm = 1.0

                            # Ligand color check
                            if color_rgb == (255, 230, 80):
                                attr = curses.color_pair(12) | curses.A_BOLD
                            elif active_color in ('depth', 'depth1'):
                                if z_norm >= 0.85:
                                    attr = curses.color_pair(9) | curses.A_BOLD
                                elif z_norm >= 0.70:
                                    attr = curses.color_pair(9)
                                elif z_norm >= 0.55:
                                    attr = curses.color_pair(10) | curses.A_BOLD
                                elif z_norm >= 0.40:
                                    attr = curses.color_pair(10)
                                elif z_norm >= 0.25:
                                    attr = curses.color_pair(10) | curses.A_DIM
                                elif z_norm >= 0.12:
                                    attr = curses.color_pair(11)
                                else:
                                    attr = curses.color_pair(11) | curses.A_DIM
                            elif active_color == 'depth2':
                                if z_norm >= 0.75:
                                    attr = curses.color_pair(9) | curses.A_BOLD   # Bold White
                                elif z_norm >= 0.50:
                                    attr = curses.color_pair(9)                  # Normal White
                                elif z_norm >= 0.25:
                                    attr = curses.color_pair(9) | curses.A_DIM   # Dim White
                                else:
                                    attr = curses.color_pair(14) | curses.A_BOLD # Bold Black (Charcoal)
                            elif active_color == 'depth3':
                                if z_norm >= 0.85:
                                    attr = curses.color_pair(12) | curses.A_BOLD
                                elif z_norm >= 0.70:
                                    attr = curses.color_pair(12)
                                elif z_norm >= 0.55:
                                    attr = curses.color_pair(13) | curses.A_BOLD
                                elif z_norm >= 0.40:
                                    attr = curses.color_pair(13)
                                elif z_norm >= 0.25:
                                    attr = curses.color_pair(13) | curses.A_DIM
                                elif z_norm >= 0.12:
                                    attr = curses.color_pair(11)
                                else:
                                    attr = curses.color_pair(11) | curses.A_DIM
                            elif active_color == 'rainbow':
                                if color_rgb:
                                    r, g, b = color_rgb
                                    depth_bold = curses.A_BOLD if z_norm > 0.45 else 0
                                    if b > r and b > g:
                                        attr = curses.color_pair(4) | depth_bold
                                    elif g > r and g > b:
                                        attr = curses.color_pair(6) | depth_bold
                                    elif g > 150 and r > 150:
                                        attr = curses.color_pair(7) | depth_bold
                                    elif r > g and r > b:
                                        attr = curses.color_pair(8) | depth_bold
                                    else:
                                        attr = curses.color_pair(5) | depth_bold
                            elif active_color == 'ss':
                                if color_rgb:
                                    r, g, b = color_rgb
                                    depth_bold = curses.A_BOLD if z_norm > 0.45 else 0
                                    if r > 180 and b > 180:
                                        attr = curses.color_pair(1) | depth_bold  # Helix
                                    elif r > 180 and g > 150:
                                        attr = curses.color_pair(2) | depth_bold  # Sheet
                                    else:
                                        attr = curses.color_pair(3) | depth_bold  # Loop

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
                    zoom_scale *= 0.96
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
                    zoom_scale *= 1.04
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
                segments = get_segments(active_rep)
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


MOUSE_SGR_PATTERN = re.compile(r'\x1b\[<(\d+);(\d+);(\d+)([Mm])')
CSI_PATTERN = re.compile(r'\x1b\[[0-9;?<=]*[a-zA-Z~]')
INCOMPLETE_CSI_PATTERN = re.compile(r'\x1b\[[0-9;?<=]*$')
INCOMPLETE_ESC_PATTERN = re.compile(r'\x1b(\[[0-9;?<=]*)?$')


class AnsiTerminalSession:
    """Context manager for full-screen ANSI terminal session with raw input."""

    def __init__(self):
        self.h_out = None
        self.h_in = None
        self.orig_out_mode = None
        self.orig_in_mode = None
        self.orig_termios = None

    def __enter__(self):
        if hasattr(sys.stdout, 'reconfigure'):
            try:
                sys.stdout.reconfigure(encoding='utf-8')
            except Exception:
                pass

        if sys.platform == 'win32':
            try:
                kernel32 = ctypes.windll.kernel32
                h_out = kernel32.GetStdHandle(-11)  # STD_OUTPUT_HANDLE
                h_in = kernel32.GetStdHandle(-10)   # STD_INPUT_HANDLE
                out_m = ctypes.c_ulong()
                in_m = ctypes.c_ulong()
                if kernel32.GetConsoleMode(h_out, ctypes.byref(out_m)):
                    self.h_out = h_out
                    self.orig_out_mode = out_m.value
                    # ENABLE_VIRTUAL_TERMINAL_PROCESSING = 0x0004
                    kernel32.SetConsoleMode(h_out, self.orig_out_mode | 0x0004 | 0x0008)
                if kernel32.GetConsoleMode(h_in, ctypes.byref(in_m)):
                    self.h_in = h_in
                    self.orig_in_mode = in_m.value
                    # ENABLE_VIRTUAL_TERMINAL_INPUT = 0x0200, ENABLE_WINDOW_INPUT = 0x0008
                    kernel32.SetConsoleMode(h_in, self.orig_in_mode | 0x0200 | 0x0008)
            except Exception:
                pass
        else:
            try:
                import termios
                import tty
                self.orig_termios = termios.tcgetattr(sys.stdin)
                tty.setcbreak(sys.stdin.fileno())
            except Exception:
                pass

        # Switch to alternate screen, hide cursor, clear screen
        # \033[?1049h: Alternate screen buffer
        # \033[?25l: Hide cursor
        # \033[2J\033[H: Clear screen & move cursor to top-left
        try:
            sys.stdout.write("\033[?1049h\033[?25l\033[2J\033[H")
            sys.stdout.flush()
        except Exception:
            pass
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Disable mouse tracking, show cursor, switch back to main screen buffer
        try:
            sys.stdout.write("\033[?1006l\033[?1003l\033[?1002l\033[?1000l\033[?25h\033[?1049l")
            sys.stdout.flush()
        except Exception:
            pass

        if sys.platform == 'win32':
            try:
                if self.h_out is not None and self.orig_out_mode is not None:
                    ctypes.windll.kernel32.SetConsoleMode(self.h_out, self.orig_out_mode)
                if self.h_in is not None and self.orig_in_mode is not None:
                    ctypes.windll.kernel32.SetConsoleMode(self.h_in, self.orig_in_mode)
            except Exception:
                pass
        else:
            try:
                if self.orig_termios is not None:
                    import termios
                    termios.tcsetattr(sys.stdin, termios.TCSADRAIN, self.orig_termios)
            except Exception:
                pass


def read_available_input() -> str:
    """Read any pending characters from stdin non-blockingly, assembling complete escape sequences."""
    chars = []
    if sys.platform == 'win32':
        import msvcrt
        try:
            while msvcrt.kbhit():
                c = msvcrt.getwch()
                if c in ('\x00', '\xe0') and msvcrt.kbhit():
                    code = msvcrt.getwch()
                    if code == 'H':
                        chars.append('\x1b[A')
                    elif code == 'P':
                        chars.append('\x1b[B')
                    elif code == 'M':
                        chars.append('\x1b[C')
                    elif code == 'K':
                        chars.append('\x1b[D')
                else:
                    chars.append(c)

            # Assemble any in-flight escape sequence so mouse packets are never split
            for _ in range(5):
                raw = "".join(chars)
                if not INCOMPLETE_ESC_PATTERN.search(raw):
                    break
                time.sleep(0.002)
                if not msvcrt.kbhit():
                    break
                while msvcrt.kbhit():
                    chars.append(msvcrt.getwch())
        except Exception:
            pass
    else:
        import select
        try:
            while select.select([sys.stdin], [], [], 0)[0]:
                c = sys.stdin.read(1)
                chars.append(c)
            for _ in range(5):
                raw = "".join(chars)
                if not INCOMPLETE_ESC_PATTERN.search(raw):
                    break
                time.sleep(0.002)
                if select.select([sys.stdin], [], [], 0)[0]:
                    chars.append(sys.stdin.read(1))
                else:
                    break
        except Exception:
            pass
    return "".join(chars)


def show_molecule_3D_ansi(molecule_data, canvas, name=None, timeout=None,
                          auto_resize=True, style=None, color=None, rep='trace'):
    """
    Interactive 3D TrueColor ANSI renderer with responsive keyboard navigation.
    Directly streams 24-bit RGB ANSI escape sequences, completely bypassing curses.
    """
    is_protein = isinstance(molecule_data, MoleculeData) and molecule_data.is_protein

    try:
        ts = os.get_terminal_size()
        term_width, term_height = ts.columns, ts.lines
    except Exception:
        term_width, term_height = canvas.char_width, canvas.char_height

    if auto_resize:
        canvas.resize_canvas(term_width, max(10, term_height - 3))

    if is_protein:
        active_style = style or 'braille'
        active_color = color or 'depth'
        if active_color == 'depth':
            active_color = 'depth1'
        active_rep = rep or getattr(molecule_data, 'representation', 'trace')
        if active_rep not in ('trace', 'cartoon', 'tube'):
            active_rep = 'cartoon' if active_rep == 'ribbon' else 'trace'

        color_modes = ['depth1', 'depth2', 'depth3', 'rainbow', 'ss', 'none']
        styles = ['ascii', 'braille']
        rep_modes = ['trace', 'cartoon', 'tube']

        protein_canvas = ProteinCanvas(
            canvas.char_width, canvas.char_height,
            canvas.width, canvas.height,
            aspect_ratio=canvas.aspect_ratio,
            style=active_style, color_mode=active_color
        )

        initial_rep = active_rep
        ligands = copy.deepcopy(molecule_data.ligands)
        canonical_cache = {}

        # Scale initial representation and ligands once to fit canvas
        initial_segs = copy.deepcopy(molecule_data.backbone_segments)
        protein_scale = scale_protein_for_canvas(initial_segs, ligands, protein_canvas)
        canonical_cache[initial_rep] = initial_segs

        # Track cumulative interactive rotation matrix and zoom factor
        cum_R = np.eye(3)
        zoom_scale = 1.0

        def get_segments(target_rep):
            if target_rep not in canonical_cache:
                raw_segs = generate_backbone_points(molecule_data.structure, samples_per_res=8, representation=target_rep)
                for seg in raw_segs:
                    for pt in seg:
                        pt.coord *= protein_scale
                        if pt.left_coord is not None:
                            pt.left_coord *= protein_scale
                            pt.right_coord *= protein_scale
                        if pt.top_coord is not None:
                            pt.top_coord *= protein_scale
                            pt.bottom_coord *= protein_scale
                canonical_cache[target_rep] = raw_segs

            target_segs = copy.deepcopy(canonical_cache[target_rep])
            for seg in target_segs:
                for pt in seg:
                    pt.coord = np.dot(cum_R, pt.coord * zoom_scale)
                    if pt.left_coord is not None:
                        pt.left_coord = np.dot(cum_R, pt.left_coord * zoom_scale)
                        pt.right_coord = np.dot(cum_R, pt.right_coord * zoom_scale)
                    if pt.top_coord is not None:
                        pt.top_coord = np.dot(cum_R, pt.top_coord * zoom_scale)
                        pt.bottom_coord = np.dot(cum_R, pt.bottom_coord * zoom_scale)
            return target_segs

        def apply_protein_zoom(factor: float):
            nonlocal zoom_scale
            new_scale = zoom_scale * factor
            if new_scale < 0.05 or new_scale > 50.0:
                return
            zoom_scale = new_scale
            for seg in segments:
                for pt in seg:
                    pt.coord *= factor
                    if pt.left_coord is not None:
                        pt.left_coord *= factor
                        pt.right_coord *= factor
                    if pt.top_coord is not None:
                        pt.top_coord *= factor
                        pt.bottom_coord *= factor
            for lig in ligands:
                lig.coord *= factor

        def apply_protein_rotation_step(axis, angle_deg=5.0):
            nonlocal cum_R
            R = rotation_matrix(axis, np.radians(angle_deg))
            cum_R = np.dot(R, cum_R)
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

        segments = get_segments(active_rep)
    else:
        active_style = style or 'symbols'
        active_color = color or 'cpk'
        atom_positions, atom_elements, atom_charges, bonds = molecule_data
        atom_positions = atom_positions.copy()
        atom_positions = scale_for_canvas(atom_positions, canvas)
        mol_zoom_scale = 1.0

        def apply_mol_zoom(factor: float):
            nonlocal atom_positions, mol_zoom_scale
            new_scale = mol_zoom_scale * factor
            if new_scale < 0.05 or new_scale > 50.0:
                return
            mol_zoom_scale = new_scale
            atom_positions *= factor

        def apply_mol_rotation_step(axis, angle_deg=5.0):
            nonlocal atom_positions
            atom_positions = rotate_points(atom_positions, axis, np.radians(angle_deg))

    end_time = time.time() + timeout if timeout else None
    rotation_axis = (0, 1, 0)
    rotation_paused = False
    theta = 1.2

    fps_history = []
    last_frame_time = time.perf_counter()

    try:
        with AnsiTerminalSession():
            while True:
                # Rolling FPS calculation (15-frame rolling average)
                now = time.perf_counter()
                dt = now - last_frame_time
                last_frame_time = now
                if dt > 0:
                    fps_history.append(1.0 / dt)
                    if len(fps_history) > 15:
                        fps_history.pop(0)

                fps_val = sum(fps_history) / len(fps_history) if fps_history else 0.0
                fps_str = f"[{fps_val:4.1f} FPS]"

                # Check terminal size
                if auto_resize:
                    try:
                        curr_ts = os.get_terminal_size()
                        cw, ch = curr_ts.columns, curr_ts.lines
                        if cw != term_width or ch != term_height:
                            term_width, term_height = cw, ch
                            canvas.resize_canvas(term_width, max(10, term_height - 3))
                            if is_protein:
                                protein_canvas.resize_canvas(term_width, max(10, term_height - 3))
                            sys.stdout.write("\033[2J")
                    except Exception:
                        pass

                if is_protein:
                    # Rotate protein if auto-spinning
                    if not rotation_paused:
                        R = rotation_matrix(rotation_axis, np.radians(theta))
                        cum_R = np.dot(R, cum_R)
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

                    display_name = name or molecule_data.name or "Protein"
                    res_num = molecule_data.structure.total_residues if molecule_data.structure else len(molecule_data.atom_positions)
                    header_text = f" {display_name} ({res_num} res) [{active_rep.upper()} (M) | {active_style.upper()} (B) | {active_color.upper()} (C) | Arrows/WASD: steer | Space: pause | +/-: zoom] "
                    if len(header_text) > protein_canvas.char_width - 6:
                        header_text = f" {display_name} [{active_rep.upper()} | {active_style.upper()} | {active_color.upper()} | Arrows: steer] "
                    if len(header_text) > protein_canvas.char_width - 6:
                        header_text = header_text[:protein_canvas.char_width - 10] + "... "
                    header_line = header_text.center(protein_canvas.char_width, "=")

                    if active_color == 'none':
                        body = protein_canvas.render_plain()
                    else:
                        body = protein_canvas.render_ansi()

                    fps_line = f"\033[1;32m{fps_str}\033[0m"
                    frame = f"\033[H\033[1;36m{header_line}\033[0m\n{fps_line}\n{body}"
                    sys.stdout.write(frame)
                    sys.stdout.flush()
                else:
                    # Small molecule
                    if not rotation_paused:
                        atom_positions = rotate_points(atom_positions, rotation_axis, np.radians(theta))

                    stretched = atom_positions.copy()
                    stretched[:, 1] /= canvas.aspect_ratio
                    positions_3D = [(p[0], p[1], p[2]) for p in stretched]
                    bond_orders = getattr(molecule_data, 'bond_orders', None)

                    canvas.clear()
                    canvas.draw_molecules(positions_3D, atom_elements, atom_charges, bonds,
                                          bond_orders=bond_orders, style=active_style, color_mode=active_color)

                    display_name = name or getattr(molecule_data, 'name', None) or "Molecule"
                    header_text = f" {display_name} [Small Mol | {active_style.upper()} (B) | {active_color.upper()} (C) | Arrows/WASD: steer | Space: pause | +/-: zoom | Q: quit] "
                    if len(header_text) > canvas.char_width - 6:
                        header_text = f" {display_name} [{active_style.upper()} | {active_color.upper()} | Arrows: steer] "
                    if len(header_text) > canvas.char_width - 6:
                        header_text = header_text[:canvas.char_width - 10] + "... "
                    header_line = header_text.center(canvas.char_width, "=")
                    if active_color == 'none':
                        body = canvas.render_plain()
                    else:
                        body = canvas.render_ansi()

                    fps_line = f"\033[1;32m{fps_str}\033[0m"
                    frame = f"\033[H\033[1;36m{header_line}\033[0m\n{fps_line}\n{body}"
                    sys.stdout.write(frame)
                    sys.stdout.flush()

                # Non-blocking input loop
                time.sleep(0.025)
                raw_input = read_available_input()
                if raw_input:
                    # 1. Parse arrow keys (CSI sequences)
                    if '\x1b[A' in raw_input:
                        if rotation_paused:
                            if is_protein:
                                apply_protein_rotation_step((1, 0, 0), 5.0)
                            else:
                                apply_mol_rotation_step((1, 0, 0), 5.0)
                        else:
                            rotation_axis = (1, 0, 0)
                    elif '\x1b[B' in raw_input:
                        if rotation_paused:
                            if is_protein:
                                apply_protein_rotation_step((-1, 0, 0), 5.0)
                            else:
                                apply_mol_rotation_step((-1, 0, 0), 5.0)
                        else:
                            rotation_axis = (-1, 0, 0)
                    elif '\x1b[C' in raw_input:
                        if rotation_paused:
                            if is_protein:
                                apply_protein_rotation_step((0, -1, 0), 5.0)
                            else:
                                apply_mol_rotation_step((0, -1, 0), 5.0)
                        else:
                            rotation_axis = (0, -1, 0)
                    elif '\x1b[D' in raw_input:
                        if rotation_paused:
                            if is_protein:
                                apply_protein_rotation_step((0, 1, 0), 5.0)
                            else:
                                apply_mol_rotation_step((0, 1, 0), 5.0)
                        else:
                            rotation_axis = (0, 1, 0)

                    # Strip all CSI escape sequences and trailing incomplete sequences
                    key_chars = INCOMPLETE_CSI_PATTERN.sub('', CSI_PATTERN.sub('', raw_input))

                    # Standalone Escape key exits
                    if raw_input == '\x1b':
                        return

                    for ch in key_chars:
                        if ch in ('q', 'Q', 'x', 'X'):
                            return
                        elif ch == ' ':
                            rotation_paused = not rotation_paused
                        elif ch in ('a', 'A'):
                            if rotation_paused:
                                if is_protein:
                                    apply_protein_rotation_step((0, 1, 0), 5.0)
                                else:
                                    apply_mol_rotation_step((0, 1, 0), 5.0)
                            else:
                                rotation_axis = (0, 1, 0)
                        elif ch in ('d', 'D'):
                            if rotation_paused:
                                if is_protein:
                                    apply_protein_rotation_step((0, -1, 0), 5.0)
                                else:
                                    apply_mol_rotation_step((0, -1, 0), 5.0)
                            else:
                                rotation_axis = (0, -1, 0)
                        elif ch in ('w', 'W'):
                            if rotation_paused:
                                if is_protein:
                                    apply_protein_rotation_step((1, 0, 0), 5.0)
                                else:
                                    apply_mol_rotation_step((1, 0, 0), 5.0)
                            else:
                                rotation_axis = (1, 0, 0)
                        elif ch in ('s', 'S'):
                            if rotation_paused:
                                if is_protein:
                                    apply_protein_rotation_step((-1, 0, 0), 5.0)
                                else:
                                    apply_mol_rotation_step((-1, 0, 0), 5.0)
                            else:
                                rotation_axis = (-1, 0, 0)
                        elif ch in ('e', 'E'):
                            rotation_axis = (0, 0, 1)
                        elif ch in ('+', '=', 'r', 'R'):
                            if is_protein:
                                apply_protein_zoom(1.05)
                            else:
                                apply_mol_zoom(1.05)
                        elif ch in ('-', '_', 'f', 'F'):
                            if is_protein:
                                apply_protein_zoom(0.95)
                            else:
                                apply_mol_zoom(0.95)
                        elif ch in ('[', '{'):
                            theta = max(0.2, theta - 0.3)
                        elif ch in (']', '}'):
                            theta = min(8.0, theta + 0.3)
                        elif ch in ('m', 'M') and is_protein:
                            idx = rep_modes.index(active_rep) if active_rep in rep_modes else 0
                            active_rep = rep_modes[(idx + 1) % len(rep_modes)]
                            segments = get_segments(active_rep)
                        elif ch in ('c', 'C'):
                            if is_protein:
                                idx = color_modes.index(active_color) if active_color in color_modes else 0
                                active_color = color_modes[(idx + 1) % len(color_modes)]
                            else:
                                small_mol_color_modes = ['cpk', 'none']
                                idx = small_mol_color_modes.index(active_color) if active_color in small_mol_color_modes else 0
                                active_color = small_mol_color_modes[(idx + 1) % len(small_mol_color_modes)]
                        elif ch in ('b', 'B'):
                            if is_protein:
                                idx = styles.index(active_style) if active_style in styles else 0
                                active_style = styles[(idx + 1) % len(styles)]
                            else:
                                small_mol_styles = ['symbols', 'balls', 'braille']
                                idx = small_mol_styles.index(active_style) if active_style in small_mol_styles else 0
                                active_style = small_mol_styles[(idx + 1) % len(small_mol_styles)]

                if timeout and time.time() > end_time and not rotation_paused:
                    break
    except Exception as e:
        traceback.print_exc(file=sys.stderr)
        raise e


def draw(input_mol, name=None, width=80, height=40, three_d=True,
         add_hydrogens=False, timeout=None, stdscr=None, auto_resize=True,
         style=None, color=None, rep='trace', backend='ansi'):
    """
    Main function for TerMol. Renders small molecules or protein structures in 2D or 3D.
    Inputs:
        input_mol: SMILES string, .sdf/.mol/.pdb file path, or 4-letter RCSB PDB ID (e.g. '1ubq').
        name: Optional title to display.
        width: Canvas width in characters (default 80).
        height: Canvas height in characters (default 40).
        three_d: Whether to display interactive 3D viewer (default True).
        style: For proteins: 'braille' (default) or 'ascii'; for small molecules: 'symbols' (default), 'balls', or 'braille'.
        color: For proteins: 'depth' (default), 'rainbow', 'ss', or 'none'; for small molecules: 'cpk' (default) or 'none'.
        rep: For proteins: 'trace' (default 1D line), 'cartoon'/'ribbon' (flat ribbons and arrows), or 'tube'.
        backend: Interactive 3D backend: 'ansi' (default, 24-bit TrueColor) or 'curses' (legacy curses).
    """
    if three_d:
        if backend == 'curses':
            curses.wrapper(draw_persistent, input_mol, name=name, width=width, height=height,
                           three_d=three_d, add_hydrogens=add_hydrogens, timeout=timeout,
                           auto_resize=auto_resize, style=style, color=color, rep=rep, backend='curses')
        else:
            draw_persistent(None, input_mol, name=name, width=width, height=height,
                            three_d=three_d, add_hydrogens=add_hydrogens, timeout=timeout,
                            auto_resize=auto_resize, style=style, color=color, rep=rep, backend='ansi')
    else:
        draw_persistent(None, input_mol, name=name, width=width, height=height,
                        three_d=three_d, add_hydrogens=add_hydrogens, timeout=timeout,
                        auto_resize=auto_resize, style=style, color=color, rep=rep, backend=backend)


def draw_persistent(stdscr, input_mol, name=None, width=80, height=40,
                    three_d=True, add_hydrogens=False, timeout=None,
                    auto_resize=True, style=None, color=None, rep='trace',
                    backend='ansi'):
    """Persistent wrapper for draw supporting both ANSI and curses backends."""
    molecule_data = get_molecule_data(input_mol, three_d=three_d, add_hydrogens=add_hydrogens, representation=rep)
    canvas = MoleculeCanvas(width, height, float(width), float(height), aspect_ratio=2.0)

    if three_d:
        if backend == 'curses':
            if stdscr:
                show_molecule_3D(stdscr, molecule_data, canvas, name=name, timeout=timeout,
                                 auto_resize=auto_resize, style=style, color=color, rep=rep)
            else:
                curses.wrapper(show_molecule_3D, molecule_data, canvas, name=name,
                               timeout=timeout, auto_resize=auto_resize, style=style, color=color, rep=rep)
        else:
            show_molecule_3D_ansi(molecule_data, canvas, name=name, timeout=timeout,
                                  auto_resize=auto_resize, style=style, color=color, rep=rep)
    else:
        show_molecule_2D(molecule_data, canvas, name=name, style=style, color=color, rep=rep)


def draw_multi(input_mols, names=None, width=80, height=40, three_d=True,
               add_hydrogens=False, timeout=None, auto_resize=True,
               style=None, color=None, rep='trace', backend='ansi'):
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
                 three_d=three_d, add_hydrogens=add_hydrogens, style=style, color=color, rep=rep, backend=backend)
        return

    if backend == 'curses':
        def loop_molecules(stdscr, molecules, mol_names, loop_timeout=None):
            for i in range(len(molecules)):
                draw_persistent(stdscr, molecules[i], name=mol_names[i], timeout=loop_timeout,
                                width=width, height=height, three_d=three_d,
                                add_hydrogens=add_hydrogens, auto_resize=auto_resize,
                                style=style, color=color, rep=rep, backend='curses')

        curses.wrapper(loop_molecules, input_mols, names, loop_timeout=timeout)
    else:
        for i in range(len(input_mols)):
            draw_persistent(None, input_mols[i], name=names[i], timeout=timeout,
                            width=width, height=height, three_d=three_d,
                            add_hydrogens=add_hydrogens, auto_resize=auto_resize,
                            style=style, color=color, rep=rep, backend='ansi')


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
