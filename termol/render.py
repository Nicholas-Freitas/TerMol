import numpy as np
import curses
from pathlib import Path
from .canvas import MoleculeCanvas
from .utilities import get_molecule_data, rotate_points, rotate_molecule_for_screen, scale_for_canvas, rotation_matrix
import time
import importlib.resources as pkg_resources

def show_molecule_2D(molecule_data, canvas, name=None):
    # Get molecule data:
    atom_positions, atom_elements, atom_charges, bonds = molecule_data

    # Scale to fit the canvas:
    atom_positions = scale_for_canvas(atom_positions, canvas)

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

def show_molecule_3D(stdscr, molecule_data, canvas, name=None, timeout=None, auto_resize=True):
    curses.curs_set(0)  # Hide cursor
    stdscr.nodelay(1)  # Make getch non-blocking
    stdscr.timeout(50)  # Refresh every 50ms
    
    #Get molecule data:
    atom_positions, atom_elements, atom_charges, bonds = molecule_data

    # Scale to fit the canvas:
    if auto_resize:
        term_height, term_width = stdscr.getmaxyx()
        canvas.resize_canvas(term_width-1, term_height-1)
    atom_positions = scale_for_canvas(atom_positions, canvas)

    # When we'll quit, if we have a timeout:
    end_time = time.time() + timeout if timeout else None

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

    stdscr.clear()
    
    theta = 1
    while True:
        stdscr.erase()

        # Get terminal size
        term_height, term_width = stdscr.getmaxyx()

        # Resize canvas:
        if auto_resize:
            canvas.resize_canvas(term_width-1, term_height-1)
        
        # Ensure the canvas fits within the terminal size
        if canvas.char_width > term_width-1 or canvas.char_height > term_height-1:
            error_message = "Please increase the size of your terminal window."
            stdscr.addstr(0, 0, error_message)
            stdscr.refresh()
            key = stdscr.getch()
            if key != -1:
                if key == curses.KEY_RESIZE:
                    continue  # Ignore resize keypress
                break  # Exit on any other key press

            # Pause timeout:
            end_time = time.time() + timeout if timeout else None
            continue
        
        if name:
            # Create a stylish header:
            header = f" {name} "
            #header = f" {canvas.char_width, canvas.char_height, canvas.width, canvas.height} "
            if len(header) > canvas.char_width-5:
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
            stdscr.addstr(1, 0, str(canvas))
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
            elif key in [ord('R'), ord('r')]:
                # Scale the molecule down (further away):
                atom_positions *= 0.99
            elif key in [ord('F'), ord('f')]:
                # Scale the molecule up (closer):
                atom_positions *= 1.01
            else:
                break  # Exit on any other key press

        if timeout and time.time() > end_time:
            # If we've paused using the spacebar, we don't want to exit on timeout
            if not rotation_paused:
                break

def draw(input_mol, name=None, width=80, height=40, three_d=True, add_hydrogens=False, timeout=None, stdscr=None, auto_resize=True):
    '''
    Main function for TerMol. This wraps the draw_persistent() function in a curses wrapper.
    This is so the user can utilize the draw_persistent() function directly if they want to keep the curses window open between renders.
    Inputs:
        input_mol: Either a SMILES string, file path to a .sdf or .mol file, or rdkit Mol.
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
        curses.wrapper(draw_persistent, input_mol, name=name, width=width, height=height, three_d=three_d, add_hydrogens=add_hydrogens, timeout=timeout, auto_resize=auto_resize)
    else:
        draw_persistent(None, input_mol, name=name, width=width, height=height, three_d=three_d, add_hydrogens=add_hydrogens, timeout=timeout, auto_resize=auto_resize)

def draw_persistent(stdscr, input_mol, name=None, width=80, height=40, three_d=True, add_hydrogens=False, timeout=None, auto_resize=True):
    '''
    Main function for TerMol:
    Inputs:
        stdscr: The curses stdscr object. Allows for a persistent window between molecules.
        input_mol: Either a SMILES string,a file path to a .sdf or .mol file, or rdkit Mol.
        name: Optional name for the molecule.
        width: Width of the canvas in characters.'
        height: Height of the canvas in characters.
        three_d: Whether to show the molecule in 3D.
        add_hydrogens: Whether to add hydrogens to the molecule.
        timeout: Time in seconds to show the molecule. If None, the molecule will be shown indefinitely. Only applies for 3D viewing.
        stdscr: The curses stdscr object. If None, the function will create a new curses window. This only needs to be used if you're keeping the curses Window between renders.
        auto_resize: Whether to automatically resize the canvas to fit the terminal window. Only applies for 3D viewing.
    Returns:
        None
        Renders 2D or 3D ASCII art of the molecule.
    '''

    # Get the molecule data:
    molecule_data = get_molecule_data(input_mol, three_d=three_d, add_hydrogens=add_hydrogens)

    # Create a canvas:
    # The canvas has a width and height in characters, and a width and height in Angstroms. Here, we're making them the same.
    canvas = MoleculeCanvas(width, height, width*1, height*1, aspect_ratio=2.0)

    # Show the molecule:
    if three_d:
        if stdscr:
            show_molecule_3D(stdscr, molecule_data, canvas, name=name, timeout=timeout, auto_resize=auto_resize)
        else:
            curses.wrapper(show_molecule_3D, molecule_data, canvas, name=name, timeout=timeout, auto_resize=auto_resize) 
    else:
        show_molecule_2D(molecule_data, canvas, name=name)

def draw_multi(input_mols, names=None, width=80, height=40, three_d=True, add_hydrogens=False, timeout=None, auto_resize=True):
    '''
    Wraps termol.draw_persistent() to allow for multiple molecules to be rendered in sequence.
    Inputs:
        input_mols: A list of SMILES strings or file paths to .sdf or .mol files.
        names: Optional names for the molecules. If None, the names will be the file names or "Molecule X" where X is the index.
                If a list is provided, it must be the same length as the number of molecules.
        width: Width of the canvas in characters.'
        height: Height of the canvas in characters.
        three_d: Whether to show the molecule in 3D.
        add_hydrogens: Whether to add hydrogens to the molecule.
        timeout: Time in seconds to show the molecule. If None, the molecule will be shown indefinitely. Only applies for 3D viewing.
    Returns:
        None
        Renders 2D or 3D ASCII art of the molecules in sequence.
    '''

    if names and len(names) != len(input_mols):
        raise ValueError("The number of names must be the same as the number of molecules!")
    
    if names is None:
        names = [f"Molecule {i+1}" for i in range(len(input_mols))] # 1-indexing for the biologists. Forgive me.

    ### For 2D rendering, just loop:
    if not three_d:
        for i in range(len(input_mols)):
            # choose a random molecule:
            name = names[i]
            molecule = input_mols[i]
            draw(molecule, name=name, width=width, height=height, three_d=three_d, add_hydrogens=add_hydrogens)
        return
    
    ### For 3D rendering, create a persistent window, then loop using draw_persistent():
    
    def loop_molecules(stdscr, molecules, names, timeout=None):
        for i in range(len(molecules)):
            name = names[i]
            molecule = molecules[i]
            
            try:
                draw_persistent(stdscr, molecule, name=name, timeout=timeout, width=width, height=height, three_d=three_d, add_hydrogens=add_hydrogens, auto_resize=auto_resize)
            except Exception as e:
                if Exception == KeyboardInterrupt:
                    break
                print(f"Failed to render {name}")
                continue
    
    curses.wrapper(loop_molecules, input_mols, names, timeout=timeout)

def showcase(timeout=5, width=80, height=40):
    ### Run a showcase of the program ###

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
    
    molecules = list(smiles_dict.values())
    names = list(smiles_dict.keys())
    draw_multi(molecules, names, width=width, height=height, timeout=timeout, auto_resize=True)
