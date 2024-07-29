from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger 
import numpy as np
import curses
from pathlib import Path
from .canvas import MoleculeCanvas
import time
import importlib.resources as pkg_resources

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

    # Rotate the molecule so it's flat on the camera plane:
    if three_d:
       atom_positions = rotate_molecule_for_screen(atom_positions)

    return atom_positions, atom_elements, atom_charges, bonds

def rotate_molecule_for_screen(atom_positions):
    '''
    Here, we rotate the molecule so it's longest axis is in the x-direction, it's second longest in the y-direction, and it's shortest in the z-direction.
    This avoids the molecule lying flat in the camera view.
    Thank you to ChatGPT lol.
    Inputs:
        atom_positions: The 3D coordinates of the molecule.
    Returns:
        rotated_positions: The rotated 3D coordinates of the molecule.
    '''

    # Step 1: Center the data
    mean_centered = atom_positions - np.mean(atom_positions, axis=0)

    # Step 2: Compute the covariance matrix
    cov_matrix = np.cov(mean_centered, rowvar=False)

    # Step 3: Perform eigen decomposition
    eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)

    # Step 4: Sort the eigenvectors by eigenvalues in descending order
    sorted_indices = np.argsort(eigenvalues)[::-1]
    sorted_eigenvectors = eigenvectors[:, sorted_indices]

    # Step 5: Rotate the data
    rotated_positions = np.dot(mean_centered, sorted_eigenvectors)

    return rotated_positions



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
    # What's the minimum and maximum x and y values?
    min_x, max_x = np.min(atom_positions[:, 0]), np.max(atom_positions[:, 0])
    min_y, max_y = np.min(atom_positions[:, 1]), np.max(atom_positions[:, 1])

    # How much would we need to scale on the x and y axes, to fit the width and height?
    x_scaling_factor = canvas.width / (max_x - min_x)
    y_scaling_factor = canvas.height / (max_y - min_y)

    # Scale all positions:
    atom_positions *= 0.9 * min(x_scaling_factor, canvas.aspect_ratio*y_scaling_factor)

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

    # Which way will we rotate?
    rotation_axis_map = {
        (ord('A'), ord('a'), curses.KEY_LEFT):    (0, 1, 0),
        (ord('D'), ord('d'), curses.KEY_RIGHT):   (0, -1, 0),
        (ord('W'), ord('w'), curses.KEY_UP):      (1, 0, 0),
        (ord('S'), ord('s'), curses.KEY_DOWN):    (-1, 0, 0),
        (ord('Q'), ord('q')):                     (0, 0, -1),
        (ord('E'), ord('e')):                     (0, 0, 1)
    }
    all_rotation_keys = [key for key_list in rotation_axis_map.keys() for key in key_list]
    rotation_axis = (0,1,0)
    rotation_paused = False

    theta = 1
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
        
        if not rotation_paused:
            atom_positions = rotate_points(atom_positions, rotation_axis, np.radians(theta))

        # Stretch to aspect ratio:
        stretched_positions = atom_positions.copy()
        stretched_positions[:, 1] /= canvas.aspect_ratio

        # Get 2D positions:
        atom_positions_2D = [(pos[0], pos[1]) for pos in stretched_positions]

        canvas.clear()
        canvas.draw_molecules(atom_positions_2D, atom_elements, atom_charges, bonds)
        try:
            stdscr.addstr(1, 1, str(canvas))
        except curses.error:
            pass  # Handle the error gracefully

        stdscr.refresh()
        
        key = stdscr.getch()
        if key != -1:
            if key == curses.KEY_RESIZE:
                continue  # Ignore resize keypress
            elif key in all_rotation_keys:
                for key_list, axis in rotation_axis_map.items():
                    if key in key_list:
                        rotation_axis = axis
                        rotation_paused = False
            elif key == ord(' '):
                rotation_paused = not rotation_paused
            else:
                break  # Exit on any other key press

        if timeout and time.time() > timeout:
            # If we've paused using the spacebar, we don't want to exit on timeout
            if not rotation_paused:
                break

def draw(input_mol, name=None, width=80, height=40, three_d=True, add_hydrogens=False, timeout=None, stdscr=None):
    '''
    Main function for TerMol. This wraps the draw_persistent() function in a curses wrapper.
    This is so the user can utilize the draw_persistent() function directly if they want to keep the curses window open between renders.
    Inputs:
        input_mol: Either a SMILES string or a file path to a .sdf or .mol file.
        name: Optional name for the molecule.
        width: Width of the canvas in characters.'
        height: Height of the canvas in characters.
        three_d: Whether to show the molecule in 3D.
        add_hydrogens: Whether to add hydrogens to the molecule.
        timeout: Time in seconds to show the molecule. If None, the molecule will be shown indefinitely. Only applies for 3D viewing.
        stdscr: The curses stdscr object. If None, the function will create a new curses window. This only needs to be used if you're keeping the curses Window between renders.
    Returns:
        None
        Renders 2D or 3D ASCII art of the molecule.
    '''
    if three_d:
        curses.wrapper(draw_persistent, input_mol, name=name, width=width, height=height, three_d=three_d, add_hydrogens=add_hydrogens, timeout=timeout)
    else:
        draw_persistent(None, input_mol, name=name, width=width, height=height, three_d=three_d, add_hydrogens=add_hydrogens, timeout=timeout)

def draw_persistent(stdscr, input_mol, name=None, width=80, height=40, three_d=True, add_hydrogens=False, timeout=None):
    '''
    Main function for TerMol:
    Inputs:
        stdscr: The curses stdscr object. Allows for a persistent window between molecules.
        input_mol: Either a SMILES string or a file path to a .sdf or .mol file.
        name: Optional name for the molecule.
        width: Width of the canvas in characters.'
        height: Height of the canvas in characters.
        three_d: Whether to show the molecule in 3D.
        add_hydrogens: Whether to add hydrogens to the molecule.
        timeout: Time in seconds to show the molecule. If None, the molecule will be shown indefinitely. Only applies for 3D viewing.
        stdscr: The curses stdscr object. If None, the function will create a new curses window. This only needs to be used if you're keeping the curses Window between renders.
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
        if stdscr:
            show_molecule_3D(stdscr, molecule_data, canvas, name=name, timeout=timeout)
        else:
            curses.wrapper(show_molecule_3D, molecule_data, canvas, name=name, timeout=timeout) 
    else:
        show_molecule_2D(molecule_data, canvas, name=name)


def showcase(timeout=5):
    ### Run a showcase of the program ###

    def loop_molecules(stdscr, smiles_dict, name=None, timeout=None):
        while True:
            # choose a random molecule:
            name = np.random.choice(list(smiles_dict.keys()))
            smiles = smiles_dict[name]
            
            try:
                draw_persistent(stdscr, smiles, name=name, timeout=timeout)
            except Exception as e:
                if Exception == KeyboardInterrupt:
                    break
                print(f"Failed to render {name}")
                continue

    # Load CSV as dictionary:
    smiles_dict = {}
    with pkg_resources.open_text(__package__, 'smiles_1000.csv') as file:
        lines = file.readlines()
        for line in lines:
            split_line = line.split('\t')
            if len(split_line) != 2:
                continue
            name, smiles = split_line
            smiles_dict[name] = smiles
    
    curses.wrapper(loop_molecules, smiles_dict, timeout=timeout)
