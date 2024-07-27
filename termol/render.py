from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger 
import numpy as np
import curses
from pathlib import Path
from .canvas import MoleculeCanvas
import time

def get_molecule_data(input_mol, three_d=True, add_hydrogens=False):

    # Disable those pesky RDKit warnings
    RDLogger.DisableLog('rdApp.*')

    # Is the input_mol a filepath or a SMILES string?
    if Path(input_mol).suffix in ['.sdf', '.mol']:
        mol = Chem.MolFromMolFile(input_mol)
    else:
        mol = Chem.MolFromSmiles(input_mol)
    
    if not mol:
        raise ValueError("Invalid SMILES string or file path!")
    
    if add_hydrogens:
        mol = Chem.AddHs(mol)

    if three_d:
        # Generate 3D coordinates
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
    else:
        AllChem.Compute2DCoords(mol)
    
    
    # Get bond and atom positions
    conf = mol.GetConformer()
    atom_positions = [conf.GetAtomPosition(i) for i in range(len(mol.GetAtoms()))]
    # Get positions in the form of [(x,y,z), ...]
    atom_positions = np.array([[pos.x, pos.y, pos.z] for pos in atom_positions])
    atom_elements  = [atom.GetSymbol() for atom in mol.GetAtoms()]
    atom_charges   = [atom.GetFormalCharge() for atom in mol.GetAtoms()]
    bonds          = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()]

    return atom_positions, atom_elements, atom_charges, bonds

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation
    about the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([
        [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
        [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
        [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]
    ])

def rotate_points(points, axis, theta):
    R = rotation_matrix(axis, theta)
    return np.array([np.dot(R, point) for point in points])

def show_molecule_2D(molecule_data, canvas, name=None):
    # Get molecule data:
    atom_positions, atom_elements, atom_charges, bonds = molecule_data

    # Scale to fit the canvas:
    # What's the maximum distance between any two atoms?
    max_distance = 0
    for i in range(len(atom_positions)):
        for j in range(i+1, len(atom_positions)):
            distance = np.linalg.norm(atom_positions[i] - atom_positions[j])
            max_distance = max(max_distance, distance)
    
    scaling_factor = 1.8 * min(canvas.width, canvas.height) / max_distance

    # Scale all positions:
    atom_positions *= scaling_factor

    # Stretch to aspect ratio:
    atom_positions[:, 1] /= canvas.aspect_ratio

    # Get 2D positions:
    atom_positions_2D = [(pos[0], pos[1]) for pos in atom_positions]

    canvas.clear()
    canvas.draw_molecules(atom_positions_2D, atom_elements, atom_charges, bonds)

    if name:
        # Create a stylish header:
        header = f" {name} "
        if len(header) > canvas.char_width-10:
            header = header[:canvas.char_width-10] + "... "
        header = header.center(canvas.char_width, "=")
        print(header)

    print(canvas)

def show_molecule_3D(stdscr, molecule_data, canvas, name=None, timeout=None):
    curses.curs_set(0)  # Hide cursor
    stdscr.nodelay(1)  # Make getch non-blocking
    stdscr.timeout(50)  # Refresh every 50ms
    
    #Get molecule data:
    atom_positions, atom_elements, atom_charges, bonds = molecule_data

    #Do an initial scaling of the molecule:
    # Scale to fit the canvas:
    # What's the maximum distance between any two atoms?
    max_distance = 0
    for i in range(len(atom_positions)):
        for j in range(i+1, len(atom_positions)):
            distance = np.linalg.norm(atom_positions[i] - atom_positions[j])
            max_distance = max(max_distance, distance)
    
    scaling_factor = 2* min(canvas.width, canvas.height) / max_distance

    # Scale all positions:
    atom_positions *= scaling_factor

    # When we'll quit, if we have a timeout:
    timeout = time.time() + timeout if timeout else None

    theta = 0
    while True:
        stdscr.clear()
        
        # Get terminal size
        term_height, term_width = stdscr.getmaxyx()
        
        # Ensure the canvas fits within the terminal size
        if canvas.char_width > term_width or canvas.char_height > term_height+1:
            error_message = "Please increase the size of your terminal window."
            stdscr.addstr(0, 0, error_message)
            stdscr.refresh()
            key = stdscr.getch()
            if key != -1:
                if key == curses.KEY_RESIZE:
                    continue  # Ignore resize keypress
                break  # Exit on any other key press
            continue
        
        if name:
            # Create a stylish header:
            header = f" {name} "
            if len(header) > canvas.char_width-10:
                header = header[:canvas.char_width-10] + "... "
            header = header.center(canvas.char_width, "=")
            stdscr.addstr(0, 0, header)
        
        # Rotate the 3D coordinates:
        rotated_positions = rotate_points(atom_positions, (0, 1, 0), np.radians(theta))

        # Stretch to aspect ratio:
        rotated_positions[:, 1] /= canvas.aspect_ratio

        # Get 2D positions:
        atom_positions_2D = [(pos[0], pos[1]) for pos in rotated_positions]

        canvas.clear()
        canvas.draw_molecules(atom_positions_2D, atom_elements, atom_charges, bonds)
        try:
            stdscr.addstr(1, 1, str(canvas))
        except curses.error:
            pass  # Handle the error gracefully

        stdscr.refresh()
        theta += np.radians(100)  # For some reason this isn't working as expected...
        
        key = stdscr.getch()
        if key != -1:
            if key == curses.KEY_RESIZE:
                continue  # Ignore resize keypress
            break  # Exit on any other key press

        if timeout and time.time() > timeout:
            break

def draw(input_mol, name=None, width=80, height=40, three_d=True, add_hydrogens=False, timeout=None):
    '''
    Main function for TerMol:
    Inputs:
        input_mol: Either a SMILES string or a file path to a .sdf or .mol file.
        name: Optional name for the molecule.
        width: Width of the canvas in characters.'
        height: Height of the canvas in characters.
        three_d: Whether to show the molecule in 3D.
        add_hydrogens: Whether to add hydrogens to the molecule.
        timeout: Time in seconds to show the molecule. If None, the molecule will be shown indefinitely Only applies for 3D viewing.
    Returns:
        None
        Renders 2D or 3D ASCII art of the molecule.
    '''
    # Get the molecule data:
    molecule_data = get_molecule_data(input_mol, three_d=three_d, add_hydrogens=add_hydrogens)

    # Create a canvas:
    canvas = MoleculeCanvas(width, height, width, height, aspect_ratio=2.0)

    # Show the molecule:
    if three_d:
        curses.wrapper(show_molecule_3D, molecule_data, canvas, name=name, timeout=timeout) 
    else:
        show_molecule_2D(molecule_data, canvas, name=name)


def main():
    ### Run a showcase of the program ###
    # Load CSV as dictionary:
    smiles_dict = {}
    with open('smiles.csv', 'r') as file:
        lines = file.readlines()
        for line in lines:
            split_line = line.split('\t')
            if len(split_line) != 2:
                continue
            name, smiles = split_line
            smiles_dict[name] = smiles
    
    while True:

        # choose a random molecule:
        name = np.random.choice(list(smiles_dict.keys()))
        smiles = smiles_dict[name]
        
        try:
            draw(smiles, name=name, three_d=True, add_hydrogens=False, timeout=10)
        except Exception as e:
            if Exception == KeyboardInterrupt:
                break
            print(f"Failed to render {name}")
            continue

if __name__ == "__main__":
    ### Example Usage ###
    # smiles = "c1cc2c(cc1[N+](=O)[O-])[nH]nn2" # Nitrobenzotriazole
    # name = "Nitrobenzotriazole"
    # draw(smiles, name=name, three_d=True, add_hydrogens=True)

    ### Run a showcase of the program ###
    main()
