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
        cx = int(np.round(((x - self.min_x) / span_x) * (self.char_width - 1)))
        cy = int(np.round(((y - self.min_y) / span_y) * (self.char_height - 1)))
        return cx, cy

    def clear(self):
        """Reset character, depth, and color buffers."""
        self.data = [[' '] * self.char_width for _ in range(self.char_height)]
        self.z_buffer = np.full((self.char_height, self.char_width), -np.inf, dtype=float)
        self.color_data = [[None] * self.char_width for _ in range(self.char_height)]

    def set_pixel(self, x: float, y: float, char: str, z: float = 0.0, color: Optional[Tuple[int, int, int]] = None) -> bool:
        """
        Place a character on the canvas with Z-buffering.
        Returns True if the pixel was placed (closer than existing depth).
        """
        cx, cy = self.map_coords(x, y)
        if 0 <= cx < self.char_width and 0 <= cy < self.char_height:
            if z >= self.z_buffer[cy, cx]:
                self.z_buffer[cy, cx] = z
                self.data[cy][cx] = char
                self.color_data[cy][cx] = color
                return True
        return False

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


class MoleculeCanvas(Canvas):
    def __init__(self, char_width: int, char_height: int, width: float, height: float, aspect_ratio: float = 2.0):
        super().__init__(char_width, char_height, width, height, aspect_ratio)

    def draw_molecules(self, atom_positions, atom_elements, atom_charges, bonds):
        self.draw_bonds(atom_positions, bonds)
        self.draw_atoms(atom_positions, atom_elements, atom_charges)

    def draw_bonds(self, atom_positions, bonds):
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

        def get_bond_character(x1, y1, x2, y2):
            dx = x2 - x1
            dy = y2 - y1
            angle = math.atan2(dy, dx) * 180 / math.pi
            if angle < 0:
                angle += 180
            if angle < 30 or angle > 150:
                return '-'
            elif 30 <= angle < 60:
                return '\\'
            elif 60 <= angle < 120:
                return '|'
            else:
                return '/'

        for bond in bonds:
            begin_idx, end_idx = bond
            p1 = atom_positions[begin_idx]
            p2 = atom_positions[end_idx]
            x1, y1 = int(p1[0]), int(p1[1])
            x2, y2 = int(p2[0]), int(p2[1])
            z1 = p1[2] if len(p1) > 2 else 0.0
            z2 = p2[2] if len(p2) > 2 else 0.0
            avg_z = 0.5 * (z1 + z2)

            points = bresenham(x1, y1, x2, y2)
            bond_char = get_bond_character(x1, y1, x2, y2)
            for x, y in points:
                if self[x, y] == ' ':
                    self.set_pixel(x, y, bond_char, z=avg_z - 0.1)

    def draw_atoms(self, atom_positions, atom_elements, atom_charges):
        for i in range(len(atom_positions)):
            pos = atom_positions[i]
            element = atom_elements[i]
            x, y = int(pos[0]), int(pos[1])
            z = pos[2] if len(pos) > 2 else 0.0

            self.set_pixel(x, y, element, z=z)

            if atom_charges[i] != 0:
                charge = atom_charges[i]
                cx, cy = self.map_coords(x, y)
                ch_x = cx - 1
                ch_y = cy - 1
                if 0 <= ch_x and ch_x + 3 < self.char_width and 0 <= ch_y < self.char_height:
                    self.data[ch_y][ch_x] = '['
                    self.data[ch_y][ch_x + 1] = '+' if charge > 0 else '-'
                    self.data[ch_y][ch_x + 2] = str(abs(charge))
                    self.data[ch_y][ch_x + 3] = ']'


class ProteinCanvas(Canvas):
    def __init__(self, char_width: int, char_height: int, width: float, height: float,
                 aspect_ratio: float = 2.0, style: str = 'ascii', color_mode: str = 'depth'):
        """
        Interactive Protein Canvas with Z-buffering, depth cueing, and Braille support.
        style: 'ascii' (default density ramp) or 'braille' (subpixel dots)
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

    def clear(self):
        super().clear()
        self.sub_width = self.char_width * 2
        self.sub_height = self.char_height * 4
        self.braille_bits = np.zeros((self.char_height, self.char_width), dtype=int)
        self.sub_z_buffer = np.full((self.sub_height, self.sub_width), -np.inf, dtype=float)

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
        if ligands:
            for lig in ligands:
                all_zs.append(lig.coord[2])

        if not all_zs:
            return

        z_min = float(np.min(all_zs))
        z_max = float(np.max(all_zs))
        z_range = max(z_max - z_min, 1e-4)

        if self.style == 'braille':
            self._draw_protein_braille(backbone_segments, ligands, z_min, z_range)
        else:
            self._draw_protein_ascii(backbone_segments, ligands, z_min, z_range)

    def _get_color(self, z_norm: float, norm_pos: float, ss: str) -> Optional[Tuple[int, int, int]]:
        """Compute RGB color based on active color_mode and normalized depth."""
        if self.color_mode == 'none':
            return None

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

        else:
            # Default 'depth' mode: Ice-blue to bright white luminance ramp
            # Farthest: dark slate blue (30, 60, 110)
            # Nearest: bright glowing cyan-white (220, 245, 255)
            r = int(30 + 190 * val_depth)
            g = int(60 + 185 * val_depth)
            b = int(110 + 145 * val_depth)
            return r, g, b

    def _draw_protein_ascii(self, backbone_segments, ligands, z_min: float, z_range: float):
        """Rasterize backbone segments onto ASCII grid using density ramp and Z-buffer."""
        for segment in backbone_segments:
            if not segment:
                continue

            for i in range(len(segment)):
                pt = segment[i]
                x, y, z = pt.coord
                z_norm = (z - z_min) / z_range
                color = self._get_color(z_norm, pt.norm_pos, pt.ss)

                char_idx = int(z_norm * (len(ASCII_RAMP) - 1))
                char = ASCII_RAMP[max(0, min(char_idx, len(ASCII_RAMP) - 1))]

                # Connect consecutive points with linear steps if needed
                if i > 0:
                    prev_pt = segment[i - 1]
                    px, py, pz = prev_pt.coord
                    cx1, cy1 = self.map_coords(px, py)
                    cx2, cy2 = self.map_coords(x, y)
                    steps = max(abs(cx2 - cx1), abs(cy2 - cy1), 1)
                    for s in range(steps + 1):
                        frac = s / steps
                        ix = px + frac * (x - px)
                        iy = py + frac * (y - py)
                        iz = pz + frac * (z - pz)
                        iz_norm = (iz - z_min) / z_range
                        ic_idx = int(iz_norm * (len(ASCII_RAMP) - 1))
                        ic_char = ASCII_RAMP[max(0, min(ic_idx, len(ASCII_RAMP) - 1))]
                        ic_color = self._get_color(iz_norm, pt.norm_pos, pt.ss)
                        self.set_pixel(ix, iy, ic_char, z=iz, color=ic_color)
                else:
                    self.set_pixel(x, y, char, z=z, color=color)

        # Draw any bound ligands as prominent element markers
        if ligands:
            for lig in ligands:
                lx, ly, lz = lig.coord
                lz_norm = (lz - z_min) / z_range
                # Golden/amber highlight for ligand atoms
                lig_color = (255, 230, 80) if self.color_mode != 'none' else None
                self.set_pixel(lx, ly, lig.element[:1].upper(), z=lz + 0.2, color=lig_color)

    def _draw_protein_braille(self, backbone_segments, ligands, z_min: float, z_range: float):
        """Rasterize backbone onto high-resolution 2x4 Braille subpixel grid."""
        span_x = self.max_x - self.min_x
        span_y = self.max_y - self.min_y

        for segment in backbone_segments:
            if not segment:
                continue

            for i in range(len(segment)):
                pt = segment[i]
                x, y, z = pt.coord
                prev_pt = segment[i - 1] if i > 0 else pt
                px, py, pz = prev_pt.coord

                # Map real coords to subpixel coords
                sub_x1 = int(np.round(((px - self.min_x) / span_x) * (self.sub_width - 1)))
                sub_y1 = int(np.round(((py - self.min_y) / span_y) * (self.sub_height - 1)))
                sub_x2 = int(np.round(((x - self.min_x) / span_x) * (self.sub_width - 1)))
                sub_y2 = int(np.round(((y - self.min_y) / span_y) * (self.sub_height - 1)))

                steps = max(abs(sub_x2 - sub_x1), abs(sub_y2 - sub_y1), 1)

                for s in range(steps + 1):
                    frac = s / steps
                    sub_x = int(np.round(sub_x1 + frac * (sub_x2 - sub_x1)))
                    sub_y = int(np.round(sub_y1 + frac * (sub_y2 - sub_y1)))
                    iz = pz + frac * (z - pz)

                    if 0 <= sub_x < self.sub_width and 0 <= sub_y < self.sub_height:
                        if iz >= self.sub_z_buffer[sub_y, sub_x]:
                            self.sub_z_buffer[sub_y, sub_x] = iz
                            cx = sub_x // 2
                            cy = sub_y // 4
                            dx = sub_x % 2
                            dy = sub_y % 4
                            self.braille_bits[cy, cx] |= BRAILLE_DOT_MAP.get((dx, dy), 0)

                            iz_norm = (iz - z_min) / z_range
                            if iz >= self.z_buffer[cy, cx]:
                                self.z_buffer[cy, cx] = iz
                                self.color_data[cy][cx] = self._get_color(iz_norm, pt.norm_pos, pt.ss)

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
            for lig in ligands:
                lx, ly, lz = lig.coord
                lz_norm = (lz - z_min) / z_range
                lig_color = (255, 230, 80) if self.color_mode != 'none' else None
                self.set_pixel(lx, ly, lig.element[:1].upper(), z=lz + 0.2, color=lig_color)
