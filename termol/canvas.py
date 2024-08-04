import math
import numpy as np

class Canvas():
    def __init__(self, char_width, char_height, width, height, aspect_ratio):
        '''
        This object is an ascii-canvas of char_width x char_height characters.
        The canvas maps real-world coordinates to the ascii-canvas coordinates, so the center of the canvas has coords (0,0).
        The aspect ratio is used to scale the coordinates to accound for the fact that characters are not square.        
        '''
        self.char_width = char_width
        self.char_height = char_height

        self.min_x = -width/2
        self.max_x = width/2
        self.min_y = -height/2
        self.max_y = height/2

        self.width = width
        self.height = height

        self.aspect_ratio = aspect_ratio

        self.data = [[' '] * char_width for i in range(char_height)]

    def resize_canvas(self, char_width, char_height):
        '''
        Resize the canvas. This changes the character width and height,
        and also scales the "real-world" width in height.
        Also clears canvas data.
        '''
        # If odd, make it even
        if char_width % 2 == 1:
            char_width -= 1
        if char_height % 2 == 1:
            char_height -= 1

        # Save old values
        char_width_old = self.char_width
        char_height_old = self.char_height

        # Update values
        self.char_width = char_width
        self.char_height = char_height

        # Scale width and height:
        scale_x = char_width / char_width_old
        scale_y = char_height / char_height_old

        # Update the real-world width and height
        self.width  = self.width * scale_x
        self.height = self.height * scale_y

        self.min_x = -self.width/2
        self.max_x = self.width/2
        self.min_y = -self.height/2
        self.max_y = self.height/2

        self.data = [[' '] * char_width for i in range(char_height)]

    def map_coords(self, x, y):
        '''
        Maps real-world coordinates to ascii-canvas coordinates.
        '''
        x = int(np.round( (x - self.min_x) / (self.max_x - self.min_x) * self.char_width) )
        y = int(np.round( (y - self.min_y) / (self.max_y - self.min_y) * self.char_height))
        return x, y
    
    def clear(self):
        '''
        Clears the canvas.
        '''
        self.data = [[' '] * self.char_width for i in range(self.char_height)]
    
    def __getitem__(self, coords: tuple):
        x, y = coords
        x, y = self.map_coords(x, y)
        if x < 0 or x >= self.char_width or y < 0 or y >= self.char_height:
            return None
        return self.data[y][x]
    
    def __setitem__(self, coords: tuple, value):
        x, y = coords
        x, y = self.map_coords(x, y)
        if x < 0 or x >= self.char_width or y < 0 or y >= self.char_height:
            return 
        self.data[y][x] = value

    def __repr__(self):
        return "\n".join(["".join(row) for row in self.data])
    

class MoleculeCanvas(Canvas):
    def __init__(self, char_width, char_height, width, height, aspect_ratio):
        super().__init__(char_width, char_height, width, height, aspect_ratio)
        self.bonds = []
        self.atoms = []

    def draw_molecules(self, atom_positions, atom_elements, atom_charges, bonds):
        '''
        Inputs:
            atom_positions: list of tuples, each tuple is the x,y position of an atom
            atom_elements: list of strings, each string is the element of an atom
            atom_charges: list of integers, each integer is the charge of an atom
            bonds: list of tuples, each tuple is the index of two atoms that are bonded
        Returns:
            None
            Updates the internal canvas to draw the atoms and bonds.
        '''
        self.draw_atoms(atom_positions, atom_elements, atom_charges)
        self.draw_bonds(atom_positions, bonds)

    def draw_bonds(self, atom_positions, bonds):
        '''
        Inputs:
            atom_positions: list of tuples, each tuple is the x,y position of an atom
            bonds: list of tuples, each tuple is the index of two atoms that are bonded
        '''

        def bresenham(x1, y1, x2, y2):
            """Bresenham's Line Algorithm to find the points on a line between two points.
            Currently, this takes a step of length 1 along the path and assumes integer inputs."""
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
            """Determine the character to use for the bond based on the angle"""
            dx = x2 - x1
            dy = y2 - y1
            angle = math.atan2(dy, dx) * 180 / math.pi
            
            # this correction flips the angles on the lower two quadrants to the upper two, for simplicity.
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
            x1, y1 = atom_positions[begin_idx]
            x2, y2 = atom_positions[end_idx]
            x1, y1, x2, y2 = int(x1), int(y1), int(x2), int(y2)
            
            points = bresenham(x1, y1, x2, y2)
            bond_char = get_bond_character(x1, y1, x2, y2)
            for x, y in points:
                #if 0 <= y < height and 0 <= x < width:
                if self[x,y] == ' ':
                    self[x,y] = bond_char

    def draw_atoms(self, atom_positions, atom_elements, atom_charges):
        '''
        Inputs:
            atom_positions: list of tuples, each tuple is the x,y position of an atom
            atom_elements: list of strings, each string is the element of an atom
            atom_charges: list of integers, each integer is the charge of an atom
        Returns:
            None
            Updats the internal canvas to draw the atoms with elements and charges.
        '''
        # Draw atoms
        for i in range(len(atom_positions)):
            pos = atom_positions[i]
            element = atom_elements[i]
            x, y = int(pos[0]), int(pos[1])
            self[x,y] = element

            if atom_charges[i] != 0:
                charge = atom_charges[i]
                canvas_coords_x, canvas_coords_y = self.map_coords(x, y)
                charge_x = canvas_coords_x - 1
                charge_y = canvas_coords_y - 1

                # Show the charge. In this case, we manually retrieve the canvas x,y position:
                if charge_x >= 0 and charge_x+3 < self.char_width and charge_y >= 0 and charge_y < self.char_height:
                    self.data[charge_y][charge_x] = '['
                    self.data[charge_y][charge_x+1] = '+' if charge > 0 else '-'
                    self.data[charge_y][charge_x+2] = str(abs(charge))
                    self.data[charge_y][charge_x+3] = ']'
