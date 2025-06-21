from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger 
import numpy as np
from pathlib import Path

def get_molecule_data(input_mol, three_d=True, add_hydrogens=False):

    # Disable those pesky RDKit warnings
    RDLogger.DisableLog('rdApp.*')

    # Is input_mol an RDKit mol object?
    if isinstance(input_mol, Chem.Mol):
        mol = input_mol
    # else, assume it's a string or file path
    elif Path(input_mol).suffix in ['.sdf', '.mol']:
        mol = Chem.MolFromMolFile(input_mol)
    elif Path(input_mol).suffix in ['.pdb']:
        mol = Chem.MolFromPDBFile(input_mol)
    else:
        mol = Chem.MolFromSmiles(input_mol)
    
    if not mol:
        raise ValueError("Invalid SMILES string or file path!")
    
    if add_hydrogens:
        mol = Chem.AddHs(mol)

    if isinstance(input_mol, Chem.Mol) or not Path(input_mol).suffix in ['.pdb']:
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

def scale_for_canvas(atom_positions, canvas):
    '''
    Given a set of atom positions and canvas object, increase the size of the molecule to fill up the screen space.
    Assumes molecules is centered at (0,0) and sitting roughly flat on the x-y plane.
    '''
    # Scale to fit the canvas:
    # What's the minimum and maximum x and y values?
    min_x, max_x = np.min(atom_positions[:, 0]), np.max(atom_positions[:, 0])
    min_y, max_y = np.min(atom_positions[:, 1]), np.max(atom_positions[:, 1])

    # How much would we need to scale on the x and y axes, to fit the width and height?
    x_scaling_factor = canvas.width / (max_x - min_x)
    y_scaling_factor = canvas.height / (max_y - min_y)

    # Scale all positions:
    atom_positions *= 0.8 * min(x_scaling_factor, canvas.aspect_ratio*y_scaling_factor)

    return atom_positions

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