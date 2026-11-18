# termol/canvas.py
"""
Terminal canvas implementation supporting:
1. 2D ASCII character grid with Aspect Ratio correction.
2. Z-buffering (depth buffer) for realistic 3D occlusion.
3. Depth cueing via character density ramps (" .:-=+*#%@").
4. High-resolution Braille subpixel rendering (2x4 dots per character cell).
5. Depth, Rainbow, and Secondary Structure color shading via TrueColor ANSI escape codes.
"""

import math
import colorsys
import numpy as np
from typing import List, Tuple, Optional, Dict, Any

ASCII_RAMP = " .:-=+*#%@"

# Braille sub-dot bit masks:
# 1 4
# 2 5
# 3 6
# 7 8
BRAILLE_DOT_MAP = {
    (0, 0): 0x01, (0, 1): 0x02, (0, 2): 0x04, (0, 3): 0x40,
    (1, 0): 0x08, (1, 1): 0x10, (1, 2): 0x20, (1, 3): 0x80
}


def rgb_to_ansi(r: int, g: int, b: int, text: str) -> str:
    """Format text with 24-bit TrueColor ANSI escape sequence."""
    return f"\033[38;2;{r};{g};{b}m{text}\033[0m"


class Canvas:
    def __init__(self, char_width: int, char_height: int, width: float, height: float, aspect_ratio: float = 2.0):
        """
        An ASCII/Unicode canvas mapping real-world coordinates to terminal character cells.
        aspect_ratio: accounts for terminal fonts being approximately 2x taller than wide.
        """
        self.char_width = char_width
        self.char_height = char_height

        self.min_x = -width / 2.0
        self.max_x = width / 2.0
        self.min_y = -height / 2.0
        self.max_y = height / 2.0

        self.width = width
        self.height = height
        self.aspect_ratio = aspect_ratio

        self.data: List[List[str]] = [[' '] * char_width for _ in range(char_height)]
        self.z_buffer: np.ndarray = np.full((char_height, char_width), -np.inf, dtype=float)
        self.color_data: List[List[Optional[Tuple[int, int, int]]]] = [[None] * char_width for _ in range(char_height)]

    def resize_canvas(self, char_width: int, char_height: int):
        """Resize the canvas to match new character dimensions and clear content."""
        if char_width % 2 == 1:
            char_width -= 1
        if char_height % 2 == 1:
            char_height -= 1

        scale_x = char_width / max(self.char_width, 1)
        scale_y = char_height / max(self.char_height, 1)

        self.char_width = char_width
        self.char_height = char_height

        self.width = self.width * scale_x
        self.height = self.height * scale_y

        self.min_x = -self.width / 2.0
        self.max_x = self.width / 2.0
        self.min_y = -self.height / 2.0
        self.max_y = self.height / 2.0

        self.clear()

    def map_coords(self, x: float, y: float) -> Tuple[int, int]:
        """Maps real-world coordinates to discrete character coordinates."""
        span_x = self.max_x - self.min_x
        span_y = self.max_y - self.min_y
        cx = int(((x - self.min_x) / span_x) * (self.char_width - 1) + 0.5)
        cy = int(((y - self.min_y) / span_y) * (self.char_height - 1) + 0.5)
        return cx, cy

    def clear(self):
        """Reset character, depth, and color buffers."""
        self.data = [[' '] * self.char_width for _ in range(self.char_height)]
        self.z_buffer = np.full((self.char_height, self.char_width), -np.inf, dtype=float)
        self.color_data = [[None] * self.char_width for _ in range(self.char_height)]

    def set_pixel_cell(self, cx: int, cy: int, char: str, z: float = 0.0, color: Optional[Tuple[int, int, int]] = None) -> bool:
        """Place a character directly at grid cell (cx, cy) with Z-buffering."""
        if 0 <= cx < self.char_width and 0 <= cy < self.char_height:
            if z >= self.z_buffer[cy, cx]:
                self.z_buffer[cy, cx] = z
                self.data[cy][cx] = char
                self.color_data[cy][cx] = color
                return True
        return False

    def set_pixel(self, x: float, y: float, char: str, z: float = 0.0, color: Optional[Tuple[int, int, int]] = None) -> bool:
        """
        Place a character on the canvas with Z-buffering.
        Returns True if the pixel was placed (closer than existing depth).
        """
        cx, cy = self.map_coords(x, y)
        return self.set_pixel_cell(cx, cy, char, z=z, color=color)

    def __getitem__(self, coords: tuple):
        x, y = coords
        cx, cy = self.map_coords(x, y)
        if 0 <= cx < self.char_width and 0 <= cy < self.char_height:
            return self.data[cy][cx]
        return None

    def __setitem__(self, coords: tuple, value):
        x, y = coords
        cx, cy = self.map_coords(x, y)
        if 0 <= cx < self.char_width and 0 <= cy < self.char_height:
            self.data[cy][cx] = value

    def render_ansi(self) -> str:
        """Render canvas with TrueColor ANSI escape sequences."""
        lines = []
        for cy in range(self.char_height):
            line_parts = []
            for cx in range(self.char_width):
                char = self.data[cy][cx]
                color = self.color_data[cy][cx]
                if color is not None and char != ' ':
                    r, g, b = color
                    line_parts.append(f"\033[38;2;{r};{g};{b}m{char}\033[0m")
                else:
                    line_parts.append(char)
            lines.append("".join(line_parts))
        return "\n".join(lines)

    def render_plain(self) -> str:
        """Render canvas in pure text without any ANSI sequences."""
        return "\n".join("".join(row) for row in self.data)

    def __str__(self) -> str:
        # Default to plain text for backwards-compatibility with existing code and curses
        return self.render_plain()

    def __repr__(self) -> str:
        return self.render_plain()


CPK_COLORS = {
    'H': (240, 240, 240),
    'C': (195, 195, 200),
    'N': (55, 125, 255),
    'O': (255, 50, 50),
    'F': (110, 240, 110),
    'Cl': (40, 220, 40),
    'Br': (180, 50, 30),
    'I': (155, 35, 185),
    'P': (255, 140, 20),
    'S': (245, 225, 30),
    'B': (255, 180, 180),
    'Si': (200, 180, 140),
    'Fe': (220, 120, 0),
    'Zn': (125, 128, 168),
    'Mg': (140, 175, 110),
    'Ca': (130, 130, 145),
    'Na': (170, 95, 240),
    'K': (145, 65, 215),
}
DEFAULT_ATOM_COLOR = (200, 200, 200)


class MoleculeCanvas(Canvas):
    def __init__(self, char_width: int, char_height: int, width: float, height: float,
                 aspect_ratio: float = 2.0, style: str = 'symbols', color_mode: str = 'cpk'):
        super().__init__(char_width, char_height, width, height, aspect_ratio)
        self.style = style.lower() if style else 'symbols'
        self.color_mode = color_mode.lower() if color_mode else 'cpk'

    def draw_molecules(self, atom_positions, atom_elements, atom_charges, bonds,
                       bond_orders=None, style=None, color_mode=None):
        active_style = (style or self.style).lower()
        active_color = (color_mode or self.color_mode).lower()

        # Determine depth range across all atom positions for depth-cueing
        zs = [p[2] for p in atom_positions if len(p) > 2]
        if zs and (max(zs) - min(zs)) > 1e-4:
            min_z = min(zs)
            max_z = max(zs)
            z_range = max_z - min_z
        else:
            min_z = 0.0
            max_z = 0.0
            z_range = 0.0

        def get_depth_scale(z: float) -> float:
            if z_range > 1e-4:
                t = (z - min_z) / z_range
                return 0.40 + 0.60 * t
            return 1.0

        self.draw_bonds(atom_positions, bonds, atom_elements=atom_elements,
                        bond_orders=bond_orders, style=active_style, color_mode=active_color,
                        get_depth_scale=get_depth_scale)
        self.draw_atoms(atom_positions, atom_elements, atom_charges,
                        style=active_style, color_mode=active_color,
                        get_depth_scale=get_depth_scale)

    def draw_bonds(self, atom_positions, bonds, atom_elements=None,
                   bond_orders=None, style='symbols', color_mode='cpk', get_depth_scale=None):
        if get_depth_scale is None:
            get_depth_scale = lambda z: 1.0

        def bresenham(x1, y1, x2, y2):
            points = []
            dx = abs(x2 - x1)
            dy = abs(y2 - y1)
            sx = 1 if x1 < x2 else -1
            sy = 1 if y1 < y2 else -1
            err = dx - dy
            while True:
                points.append((x1, y1))
                if x1 == x2 and y1 == y2:
                    break
                e2 = err * 2
                if e2 > -dy:
                    err -= dy
                    x1 += sx
                if e2 < dx:
                    err += dx
                    y1 += sy
            return points

        def get_bond_character(x1, y1, x2, y2, order=1.0):
            dx = x2 - x1
            dy = y2 - y1
            angle = math.atan2(dy, dx) * 180 / math.pi
            if angle < 0:
                angle += 180

            is_aromatic = (abs(order - 1.5) < 1e-3)
            is_triple = (order >= 2.8)
            is_double = (not is_aromatic and order >= 1.8 and order < 2.8)

            # Triple bond
            if is_triple:
                if angle < 30 or angle > 150:
                    return '≡'
                elif 60 <= angle < 120:
                    return '⦀'
                else:
                    return '⫽'

            # Double bond
            if is_double:
                if angle < 30 or angle > 150:
                    return '='
                elif 60 <= angle < 120:
                    return '‖'
                elif 30 <= angle < 60:
                    return '⫻'
                else:
                    return '⫽'

            # Aromatic bond
            if is_aromatic:
                if angle < 30 or angle > 150:
                    return '╌'
                elif 60 <= angle < 120:
                    return '┆'
                elif 30 <= angle < 60:
                    return '⋱'
                else:
                    return '⋰'

            # Standard single bond
            if angle < 30 or angle > 150:
                return '-'
            elif 30 <= angle < 60:
                return '\\'
            elif 60 <= angle < 120:
                return '|'
            else:
                return '/'

        for b_idx, bond in enumerate(bonds):
            begin_idx = bond[0]
            end_idx = bond[1]
            order = bond[2] if len(bond) > 2 else (bond_orders[b_idx] if bond_orders and b_idx < len(bond_orders) else 1.0)
            p1 = atom_positions[begin_idx]
            p2 = atom_positions[end_idx]
            x1, y1 = int(p1[0]), int(p1[1])
            x2, y2 = int(p2[0]), int(p2[1])
            z1 = p1[2] if len(p1) > 2 else 0.0
            z2 = p2[2] if len(p2) > 2 else 0.0

            elem1 = atom_elements[begin_idx] if atom_elements and begin_idx < len(atom_elements) else 'C'
            elem2 = atom_elements[end_idx] if atom_elements and end_idx < len(atom_elements) else 'C'
            color1 = CPK_COLORS.get(elem1, DEFAULT_ATOM_COLOR)
            color2 = CPK_COLORS.get(elem2, DEFAULT_ATOM_COLOR)

            points = bresenham(x1, y1, x2, y2)
            bond_char = get_bond_character(x1, y1, x2, y2, order=order)
            num_pts = len(points)

            for step_idx, (x, y) in enumerate(points):
                alpha = step_idx / (num_pts - 1) if num_pts > 1 else 0.5
                curr_z = z1 + alpha * (z2 - z1)
                if color_mode == 'none':
                    bond_color = None
                else:
                    scale = get_depth_scale(curr_z)
                    base_color = color1 if alpha < 0.5 else color2
                    bond_color = (
                        int(base_color[0] * scale * 0.8),
                        int(base_color[1] * scale * 0.8),
                        int(base_color[2] * scale * 0.8)
                    )
                self.set_pixel(x, y, bond_char, z=curr_z - 0.1, color=bond_color)

    def draw_atoms(self, atom_positions, atom_elements, atom_charges,
                   style='symbols', color_mode='cpk', get_depth_scale=None):
        if get_depth_scale is None:
            get_depth_scale = lambda z: 1.0

        for i in range(len(atom_positions)):
            pos = atom_positions[i]
            element = atom_elements[i]
            x, y = int(pos[0]), int(pos[1])
            z = pos[2] if len(pos) > 2 else 0.0

            if color_mode == 'none':
                atom_color = None
            else:
                base_color = CPK_COLORS.get(element, DEFAULT_ATOM_COLOR)
                scale = get_depth_scale(z)
                atom_color = (
                    int(base_color[0] * scale),
                    int(base_color[1] * scale),
                    int(base_color[2] * scale)
                )

            cx, cy = self.map_coords(x, y)

            # Choose atom glyph according to style
            if style in ('balls', 'spheres', 'ball'):
                if element == 'H':
                    glyph = '•'
                elif element in ('P', 'S', 'Cl', 'Br', 'I', 'Fe', 'Zn', 'Mg', 'Ca'):
                    glyph = '◉'
                else:
                    glyph = '●'
                self.set_pixel_cell(cx, cy, glyph, z=z, color=atom_color)
            elif style in ('dots', 'dot'):
                glyph = '·' if element == 'H' else '•'
                self.set_pixel_cell(cx, cy, glyph, z=z, color=atom_color)
            else:
                # Classic chemical symbol (C, N, O, Cl, etc.)
                for ch_offset, ch in enumerate(element):
                    self.set_pixel_cell(cx + ch_offset, cy, ch, z=z, color=atom_color)

            if atom_charges[i] != 0:
                charge = atom_charges[i]
                ch_x = cx - 1
                ch_y = cy - 1
                sign_char = '+' if charge > 0 else '-'
                if color_mode == 'none':
                    self.set_pixel_cell(ch_x, ch_y, '[', z=z)
                    self.set_pixel_cell(ch_x + 1, ch_y, sign_char, z=z)
                    self.set_pixel_cell(ch_x + 2, ch_y, str(abs(charge)), z=z)
                    self.set_pixel_cell(ch_x + 3, ch_y, ']', z=z)
                else:
                    scale = get_depth_scale(z)
                    sign_color = (255, 80, 80) if charge > 0 else (80, 160, 255)
                    bracket_color = (int(180 * scale), int(180 * scale), int(180 * scale))
                    num_color = (int(240 * scale), int(240 * scale), int(240 * scale))

                    self.set_pixel_cell(ch_x, ch_y, '[', z=z, color=bracket_color)
                    self.set_pixel_cell(ch_x + 1, ch_y, sign_char, z=z, color=sign_color)
                    self.set_pixel_cell(ch_x + 2, ch_y, str(abs(charge)), z=z, color=num_color)
                    self.set_pixel_cell(ch_x + 3, ch_y, ']', z=z, color=bracket_color)


class ProteinCanvas(Canvas):
    def __init__(self, char_width: int, char_height: int, width: float, height: float,
                 aspect_ratio: float = 2.0, style: str = 'braille', color_mode: str = 'depth'):
        """
        Interactive Protein Canvas with Z-buffering, depth cueing, and Braille support.
        style: 'braille' (default subpixel dots) or 'ascii' (density ramp)
        color_mode: 'depth' (default), 'rainbow', 'ss', or 'none'
        """
        super().__init__(char_width, char_height, width, height, aspect_ratio)
        self.style = style.lower()
        self.color_mode = color_mode.lower()

        # Braille subpixel grid (2 horizontal dots x 4 vertical dots per character cell)
        self.sub_width = self.char_width * 2
        self.sub_height = self.char_height * 4
        self.braille_bits = np.zeros((self.char_height, self.char_width), dtype=int)
        self.sub_z_buffer = np.full((self.sub_height, self.sub_width), -np.inf, dtype=float)
        self.last_z_min: float = -10.0
        self.last_z_range: float = 20.0

    def clear(self):
        super().clear()
        self.sub_width = self.char_width * 2
        self.sub_height = self.char_height * 4
        self.braille_bits = np.zeros((self.char_height, self.char_width), dtype=int)
        self.sub_z_buffer = np.full((self.sub_height, self.sub_width), -np.inf, dtype=float)
        self.last_z_min: float = -10.0
        self.last_z_range: float = 20.0

    def draw_protein(self, backbone_segments: List[Any], ligands: Optional[List[Any]] = None):
        """
        Render 3D backbone segments and optional ligands onto the canvas.
        Each point in backbone_segments contains:
            coord: (x, y, z)
            norm_pos: 0.0 to 1.0 (N to C terminus)
            ss: 'H' (helix), 'E' (sheet), 'C' (loop)
            chain_id: str
            res_num: int
        """
        # Collect all points to determine global depth range for depth cueing
        all_zs = []
        for segment in backbone_segments:
            for pt in segment:
                all_zs.append(pt.coord[2])
                if pt.left_coord is not None:
                    all_zs.append(pt.left_coord[2])
                    all_zs.append(pt.right_coord[2])
                if pt.top_coord is not None:
                    all_zs.append(pt.top_coord[2])
                    all_zs.append(pt.bottom_coord[2])
        if ligands:
            for lig in ligands:
                all_zs.append(lig.coord[2])

        if not all_zs:
            return

        z_min = float(np.min(all_zs))
        z_max = float(np.max(all_zs))
        z_range = max(z_max - z_min, 1e-4)

        self.last_z_min = z_min
        self.last_z_range = z_range

        if self.style == 'braille':
            self._draw_protein_braille(backbone_segments, ligands, z_min, z_range)
        else:
            self._draw_protein_ascii(backbone_segments, ligands, z_min, z_range)

    def _get_color(self, z_norm: float, norm_pos: float, ss: str) -> Optional[Tuple[int, int, int]]:
        """Compute RGB color based on active color_mode and normalized depth."""
        if self.color_mode == 'none':
            return None

        # Clamp z_norm strictly to [0.0, 1.0] to prevent negative roots or float overflow
        if z_norm <= 0.0:
            z_norm = 0.0
        elif z_norm >= 1.0:
            z_norm = 1.0
        else:
            z_norm = float(z_norm)

        # Depth luminance scaling: closer = brighter
        val_depth = 0.35 + 0.65 * z_norm

        if self.color_mode == 'rainbow':
            # Blue (0.66) at N-terminus down to Red (0.0) at C-terminus
            hue = 0.66 * (1.0 - norm_pos)
            sat = 0.85
            val = val_depth
            r, g, b = colorsys.hsv_to_rgb(hue, sat, val)
            return int(255 * r), int(255 * g), int(255 * b)

        elif self.color_mode == 'ss':
            # Secondary structure color coding:
            # Helices: Magenta/Pink (#E040FB)
            # Sheets: Golden Yellow (#FFD54F)
            # Loops: Cyan (#40C4FF)
            if ss == 'H':
                base_rgb = (235, 60, 240)
            elif ss == 'E':
                base_rgb = (255, 215, 60)
            else:
                base_rgb = (50, 210, 255)
            r = int(base_rgb[0] * val_depth)
            g = int(base_rgb[1] * val_depth)
            b = int(base_rgb[2] * val_depth)
            return r, g, b

        elif self.color_mode == 'depth2':
            # Option 2: High-contrast monochrome / pure grayscale (charcoal -> silver -> bright white)
            t = z_norm ** 1.3
            val = int(28 + 227 * t)
            return val, val, val

        elif self.color_mode == 'depth3':
            # Option 3: Twilight / neon plasma (midnight violet -> vibrant magenta -> golden peach/white)
            k0 = (20, 10, 45)       # dark violet
            k1 = (185, 45, 145)     # vibrant magenta
            k2 = (255, 238, 200)    # warm golden white
            if z_norm < 0.55:
                t = max(0.0, z_norm / 0.55) ** 1.2
                return (
                    int(k0[0] + (k1[0] - k0[0]) * t),
                    int(k0[1] + (k1[1] - k0[1]) * t),
                    int(k0[2] + (k1[2] - k0[2]) * t)
                )
            else:
                t = max(0.0, min(1.0, (z_norm - 0.55) / 0.45)) ** 0.9
                return (
                    int(k1[0] + (k2[0] - k1[0]) * t),
                    int(k1[1] + (k2[1] - k1[1]) * t),
                    int(k1[2] + (k2[2] - k1[2]) * t)
                )

        else:
            # Default 'depth' / 'depth1': Deep cyan / oceanic ice fog (midnight navy -> electric cyan -> cyan-white)
            k0 = (15, 30, 60)       # deep midnight navy
            k1 = (35, 185, 230)     # vibrant electric cyan
            k2 = (235, 252, 255)    # glowing cyan-white
            if z_norm < 0.60:
                t = max(0.0, z_norm / 0.60) ** 1.3
                return (
                    int(k0[0] + (k1[0] - k0[0]) * t),
                    int(k0[1] + (k1[1] - k0[1]) * t),
                    int(k0[2] + (k1[2] - k0[2]) * t)
                )
            else:
                t = max(0.0, min(1.0, (z_norm - 0.60) / 0.40)) ** 0.85
                return (
                    int(k1[0] + (k2[0] - k1[0]) * t),
                    int(k1[1] + (k2[1] - k1[1]) * t),
                    int(k1[2] + (k2[2] - k1[2]) * t)
                )

    def _draw_line_ascii(self, x1: float, y1: float, z1: float,
                          x2: float, y2: float, z2: float,
                          norm_pos: float, ss: str, z_min: float, z_range: float):
        """Draw a 3D line segment on the ASCII canvas with depth testing and color (fast incremental)."""
        span_x = self.max_x - self.min_x
        span_y = self.max_y - self.min_y
        scale_x = (self.char_width - 1) / span_x
        scale_y = (self.char_height - 1) / (span_y * self.aspect_ratio)
        mid_cx = (self.char_width - 1) / 2.0
        mid_cy = (self.char_height - 1) / 2.0

        cx1 = int(mid_cx + x1 * scale_x + 0.5)
        cy1 = int(mid_cy + y1 * scale_y + 0.5)
        cx2 = int(mid_cx + x2 * scale_x + 0.5)
        cy2 = int(mid_cy + y2 * scale_y + 0.5)

        dx = cx2 - cx1
        dy = cy2 - cy1
        steps = abs(dx) if abs(dx) > abs(dy) else abs(dy)
        if steps == 0:
            steps = 1

        inv_steps = 1.0 / steps
        x_inc = dx * inv_steps
        y_inc = dy * inv_steps
        z_inc = (z2 - z1) * inv_steps

        curr_x = float(cx1)
        curr_y = float(cy1)
        curr_z = z1

        char_w = self.char_width
        char_h = self.char_height
        z_buf = self.z_buffer
        data = self.data
        color_data = self.color_data
        inv_z_range = 1.0 / z_range if z_range > 0 else 1.0
        ramp_len_m1 = len(ASCII_RAMP) - 1

        for _ in range(steps + 1):
            ix = int(curr_x + 0.5)
            iy = int(curr_y + 0.5)
            if 0 <= ix < char_w and 0 <= iy < char_h:
                if curr_z >= z_buf[iy, ix]:
                    z_buf[iy, ix] = curr_z
                    iz_norm = (curr_z - z_min) * inv_z_range
                    char_idx = int(iz_norm * ramp_len_m1)
                    if char_idx > ramp_len_m1:
                        char_idx = ramp_len_m1
                    elif char_idx < 0:
                        char_idx = 0
                    data[iy][ix] = ASCII_RAMP[char_idx]
                    color_data[iy][ix] = self._get_color(iz_norm, norm_pos, ss)
            curr_x += x_inc
            curr_y += y_inc
            curr_z += z_inc

    def _draw_line_braille(self, x1: float, y1: float, z1: float,
                            x2: float, y2: float, z2: float,
                            norm_pos: float, ss: str, z_min: float, z_range: float):
        """Draw a 3D line segment on the Braille subpixel grid with depth testing and color (fast incremental)."""
        span_x = self.max_x - self.min_x
        span_y = self.max_y - self.min_y
        scale_x = (self.sub_width - 1) / span_x
        scale_y = (self.sub_height - 1) / (span_y * self.aspect_ratio)
        mid_sub_x = (self.sub_width - 1) / 2.0
        mid_sub_y = (self.sub_height - 1) / 2.0

        sub_x1 = int(mid_sub_x + x1 * scale_x + 0.5)
        sub_y1 = int(mid_sub_y + y1 * scale_y + 0.5)
        sub_x2 = int(mid_sub_x + x2 * scale_x + 0.5)
        sub_y2 = int(mid_sub_y + y2 * scale_y + 0.5)

        dx = sub_x2 - sub_x1
        dy = sub_y2 - sub_y1
        steps = abs(dx) if abs(dx) > abs(dy) else abs(dy)
        if steps == 0:
            steps = 1

        inv_steps = 1.0 / steps
        x_inc = dx * inv_steps
        y_inc = dy * inv_steps
        z_inc = (z2 - z1) * inv_steps

        curr_x = float(sub_x1)
        curr_y = float(sub_y1)
        curr_z = z1

        sub_w = self.sub_width
        sub_h = self.sub_height
        sub_z = self.sub_z_buffer
        braille_bits = self.braille_bits
        z_buf = self.z_buffer
        color_data = self.color_data
        inv_z_range = 1.0 / z_range if z_range > 0 else 1.0

        for _ in range(steps + 1):
            ix = int(curr_x + 0.5)
            iy = int(curr_y + 0.5)
            if 0 <= ix < sub_w and 0 <= iy < sub_h:
                if curr_z >= sub_z[iy, ix]:
                    sub_z[iy, ix] = curr_z
                    cx = ix >> 1
                    cy = iy >> 2
                    dx_bit = ix & 1
                    dy_bit = iy & 3
                    braille_bits[cy, cx] |= BRAILLE_DOT_MAP.get((dx_bit, dy_bit), 0)

                    iz_norm = (curr_z - z_min) * inv_z_range
                    if curr_z >= z_buf[cy, cx]:
                        z_buf[cy, cx] = curr_z
                        color_data[cy][cx] = self._get_color(iz_norm, norm_pos, ss)
            curr_x += x_inc
            curr_y += y_inc
            curr_z += z_inc

    def _draw_protein_ascii(self, backbone_segments, ligands, z_min: float, z_range: float):
        """Rasterize backbone segments onto ASCII grid using density ramp and Z-buffer."""
        for segment in backbone_segments:
            if not segment:
                continue

            has_ribbon = segment[0].left_coord is not None
            if not has_ribbon:
                # Standard trace mode
                for i in range(len(segment)):
                    pt = segment[i]
                    x, y, z = pt.coord
                    if i > 0:
                        prev_pt = segment[i - 1]
                        px, py, pz = prev_pt.coord
                        self._draw_line_ascii(px, py, pz, x, y, z, pt.norm_pos, pt.ss, z_min, z_range)
                    else:
                        span_x = self.max_x - self.min_x
                        span_y = self.max_y - self.min_y
                        scale_x = (self.char_width - 1) / span_x
                        scale_y = (self.char_height - 1) / (span_y * self.aspect_ratio)
                        mid_cx = (self.char_width - 1) / 2.0
                        mid_cy = (self.char_height - 1) / 2.0
                        cx = int(mid_cx + x * scale_x + 0.5)
                        cy = int(mid_cy + y * scale_y + 0.5)
                        z_norm = (z - z_min) / z_range
                        char_idx = int(z_norm * (len(ASCII_RAMP) - 1))
                        char = ASCII_RAMP[max(0, min(char_idx, len(ASCII_RAMP) - 1))]
                        color = self._get_color(z_norm, pt.norm_pos, pt.ss)
                        self.set_pixel_cell(cx, cy, char, z=z, color=color)
            else:
                # Ribbon / Tube / Cartoon mode
                for i in range(len(segment)):
                    pt = segment[i]
                    lx, ly, lz = pt.left_coord
                    rx, ry, rz = pt.right_coord
                    self._draw_line_ascii(lx, ly, lz, rx, ry, rz, pt.norm_pos, pt.ss, z_min, z_range)

                    if i > 0:
                        prev = segment[i - 1]
                        plx, ply, plz = prev.left_coord
                        prx, pry, prz = prev.right_coord
                        self._draw_line_ascii(plx, ply, plz, lx, ly, lz, pt.norm_pos, pt.ss, z_min, z_range)
                        self._draw_line_ascii(prx, pry, prz, rx, ry, rz, pt.norm_pos, pt.ss, z_min, z_range)

                    if pt.top_coord is not None and pt.bottom_coord is not None:
                        tx, ty, tz = pt.top_coord
                        bx, by, bz = pt.bottom_coord
                        self._draw_line_ascii(tx, ty, tz, bx, by, bz, pt.norm_pos, pt.ss, z_min, z_range)
                        if i > 0 and prev.top_coord is not None:
                            ptx, pty, ptz = prev.top_coord
                            pbx, pby, pbz = prev.bottom_coord
                            self._draw_line_ascii(ptx, pty, ptz, tx, ty, tz, pt.norm_pos, pt.ss, z_min, z_range)
                            self._draw_line_ascii(pbx, pby, pbz, bx, by, bz, pt.norm_pos, pt.ss, z_min, z_range)

        # Draw any bound ligands as prominent element markers
        if ligands:
            span_x = self.max_x - self.min_x
            span_y = self.max_y - self.min_y
            scale_x = (self.char_width - 1) / span_x
            scale_y = (self.char_height - 1) / (span_y * self.aspect_ratio)
            mid_cx = (self.char_width - 1) / 2.0
            mid_cy = (self.char_height - 1) / 2.0
            for lig in ligands:
                lx, ly, lz = lig.coord
                cx = int(mid_cx + lx * scale_x + 0.5)
                cy = int(mid_cy + ly * scale_y + 0.5)
                lig_color = (255, 230, 80) if self.color_mode != 'none' else None
                self.set_pixel_cell(cx, cy, lig.element[:1].upper(), z=lz + 0.2, color=lig_color)

    def _draw_protein_braille(self, backbone_segments, ligands, z_min: float, z_range: float):
        """Rasterize backbone onto high-resolution 2x4 Braille subpixel grid."""
        for segment in backbone_segments:
            if not segment:
                continue

            has_ribbon = segment[0].left_coord is not None
            if not has_ribbon:
                # Standard trace mode
                for i in range(len(segment)):
                    pt = segment[i]
                    x, y, z = pt.coord
                    prev_pt = segment[i - 1] if i > 0 else pt
                    px, py, pz = prev_pt.coord
                    self._draw_line_braille(px, py, pz, x, y, z, pt.norm_pos, pt.ss, z_min, z_range)
            else:
                # Ribbon / Tube / Cartoon mode
                for i in range(len(segment)):
                    pt = segment[i]
                    lx, ly, lz = pt.left_coord
                    rx, ry, rz = pt.right_coord
                    self._draw_line_braille(lx, ly, lz, rx, ry, rz, pt.norm_pos, pt.ss, z_min, z_range)

                    if i > 0:
                        prev = segment[i - 1]
                        plx, ply, plz = prev.left_coord
                        prx, pry, prz = prev.right_coord
                        self._draw_line_braille(plx, ply, plz, lx, ly, lz, pt.norm_pos, pt.ss, z_min, z_range)
                        self._draw_line_braille(prx, pry, prz, rx, ry, rz, pt.norm_pos, pt.ss, z_min, z_range)

                    if pt.top_coord is not None and pt.bottom_coord is not None:
                        tx, ty, tz = pt.top_coord
                        bx, by, bz = pt.bottom_coord
                        self._draw_line_braille(tx, ty, tz, bx, by, bz, pt.norm_pos, pt.ss, z_min, z_range)
                        if i > 0 and prev.top_coord is not None:
                            ptx, pty, ptz = prev.top_coord
                            pbx, pby, pbz = prev.bottom_coord
                            self._draw_line_braille(ptx, pty, ptz, tx, ty, tz, pt.norm_pos, pt.ss, z_min, z_range)
                            self._draw_line_braille(pbx, pby, pbz, bx, by, bz, pt.norm_pos, pt.ss, z_min, z_range)

        # Convert accumulated Braille bits to Unicode characters
        for cy in range(self.char_height):
            for cx in range(self.char_width):
                bits = self.braille_bits[cy, cx]
                if bits > 0:
                    self.data[cy][cx] = chr(0x2800 + bits)
                else:
                    self.data[cy][cx] = ' '

        # Overlay ligands
        if ligands:
            span_x = self.max_x - self.min_x
            span_y = self.max_y - self.min_y
            scale_x = (self.char_width - 1) / span_x
            scale_y = (self.char_height - 1) / (span_y * self.aspect_ratio)
            mid_cx = (self.char_width - 1) / 2.0
            mid_cy = (self.char_height - 1) / 2.0
            for lig in ligands:
                lx, ly, lz = lig.coord
                cx = int(mid_cx + lx * scale_x + 0.5)
                cy = int(mid_cy + ly * scale_y + 0.5)
                lig_color = (255, 230, 80) if self.color_mode != 'none' else None
                self.set_pixel_cell(cx, cy, lig.element[:1].upper(), z=lz + 0.2, color=lig_color)
