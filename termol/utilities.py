# termol/utilities.py
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger 
import numpy as np
from pathlib import Path
from typing import List, Tuple, Optional, Any

from .pdb import is_pdb_code, load_pdb, ProteinStructure
from .backbone import generate_backbone_points, BackbonePoint


class MoleculeData:
    """
    Container for molecule or protein representation.
    Supports unpacking as (atom_positions, atom_elements, atom_charges, bonds)
    for 100% backwards compatibility with existing TerMol code.
    """
    def __init__(self,
                 is_protein: bool = False,
                 atom_positions: Optional[np.ndarray] = None,
                 atom_elements: Optional[List[str]] = None,
                 atom_charges: Optional[List[int]] = None,
                 bonds: Optional[List[Tuple[int, int]]] = None,
                 structure: Optional[ProteinStructure] = None,
                 backbone_segments: Optional[List[List[BackbonePoint]]] = None,
                 ligands: Optional[List[Any]] = None,
                 name: Optional[str] = None):
        self.is_protein = is_protein
        self.atom_positions = atom_positions if atom_positions is not None else np.empty((0, 3))
        self.atom_elements = atom_elements if atom_elements is not None else []
        self.atom_charges = atom_charges if atom_charges is not None else []
        self.bonds = bonds if bonds is not None else []
        self.structure = structure
        self.backbone_segments = backbone_segments if backbone_segments is not None else []
        self.ligands = ligands if ligands is not None else []
        self.name = name

    def __iter__(self):
        # Allows legacy unpacking: atom_positions, atom_elements, atom_charges, bonds = data
        return iter([self.atom_positions, self.atom_elements, self.atom_charges, self.bonds])

    def __getitem__(self, idx):
        return [self.atom_positions, self.atom_elements, self.atom_charges, self.bonds][idx]

    def __len__(self):
        return 4


def get_molecule_data(input_mol, three_d=True, add_hydrogens=False):
    """
    Load molecule data from SMILES string, chemical file (.sdf, .mol),
    PDB file (.pdb), or RCSB PDB 4-letter accession code (e.g. '1ubq').
    Automatically detects proteins and generates backbone spline data.
    """
    RDLogger.DisableLog('rdApp.*')

    # 1. Check if input is a 4-character RCSB PDB code or .pdb file
    is_code = is_pdb_code(str(input_mol))
    is_pdb_file = False
    if not is_code and isinstance(input_mol, (str, Path)):
        p = Path(input_mol)
        if p.suffix.lower() in ['.pdb', '.ent']:
            is_pdb_file = True

    if is_code or is_pdb_file:
        try:
            structure = load_pdb(str(input_mol))
            if structure.is_protein:
                # Generate Catmull-Rom backbone splines
                backbone_segments = generate_backbone_points(structure, samples_per_res=8)
                ligands = structure.ligands

                # Extract all coordinates to center and PCA-align
                all_coords_list = []
                for seg in backbone_segments:
                    for pt in seg:
                        all_coords_list.append(pt.coord)
                if ligands:
                    for lig in ligands:
                        all_coords_list.append(lig.coord)

                if all_coords_list:
                    all_coords = np.array(all_coords_list)
                    center = np.mean(all_coords, axis=0)

                    # Center all points around (0, 0, 0)
                    for seg in backbone_segments:
                        for pt in seg:
                            pt.coord = pt.coord - center
                    for lig in ligands:
                        lig.coord = lig.coord - center

                    if three_d:
                        # PCA rotation so longest inertia axis is along X, second along Y
                        all_centered = all_coords - center
                        cov = np.cov(all_centered, rowvar=False)
                        eigvals, eigvecs = np.linalg.eigh(cov)
                        rot = eigvecs[:, np.argsort(eigvals)[::-1]]

                        for seg in backbone_segments:
                            for pt in seg:
                                pt.coord = np.dot(pt.coord, rot)
                        for lig in ligands:
                            lig.coord = np.dot(lig.coord, rot)

                ca_coords = structure.all_ca_coords
                if len(ca_coords) > 0 and all_coords_list:
                    ca_coords = ca_coords - center
                    if three_d:
                        ca_coords = np.dot(ca_coords, rot)

                return MoleculeData(
                    is_protein=True,
                    atom_positions=ca_coords,
                    atom_elements=['CA'] * len(ca_coords),
                    atom_charges=[0] * len(ca_coords),
                    bonds=[],
                    structure=structure,
                    backbone_segments=backbone_segments,
                    ligands=ligands,
                    name=structure.title or str(input_mol)
                )
        except Exception as e:
            # Fall through to RDKit if it's a small molecule PDB or load_pdb failed
            pass

    # 2. Standard RDKit pathway for small molecules
    if isinstance(input_mol, Chem.Mol):
        mol = input_mol
    elif Path(input_mol).suffix in ['.sdf', '.mol']:
        mol = Chem.MolFromMolFile(str(input_mol))
    elif Path(input_mol).suffix in ['.pdb']:
        mol = Chem.MolFromPDBFile(str(input_mol))
    else:
        mol = Chem.MolFromSmiles(str(input_mol))
    
    if not mol:
        raise ValueError(f"Invalid SMILES string, PDB code, or file path: '{input_mol}'")
    
    if add_hydrogens:
        mol = Chem.AddHs(mol)

    if isinstance(input_mol, Chem.Mol) or not Path(input_mol).suffix in ['.pdb']:
        if three_d:
            AllChem.EmbedMolecule(mol)
            AllChem.MMFFOptimizeMolecule(mol)
        else:
            AllChem.Compute2DCoords(mol)
    
    # Get bond and atom positions
    conf = mol.GetConformer()
    atom_positions = [conf.GetAtomPosition(i) for i in range(len(mol.GetAtoms()))]
    atom_positions = np.array([[pos.x, pos.y, pos.z] for pos in atom_positions])
    atom_elements  = [atom.GetSymbol() for atom in mol.GetAtoms()]
    atom_charges   = [atom.GetFormalCharge() for atom in mol.GetAtoms()]
    bonds          = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()]

    # Center coordinates around (0, 0, 0)
    if len(atom_positions) > 0:
        center = np.mean(atom_positions, axis=0)
        atom_positions = atom_positions - center

    # Rotate the molecule so it's flat on the camera plane:
    if three_d and len(atom_positions) > 0:
        atom_positions = rotate_molecule_for_screen(atom_positions)
    
    return MoleculeData(
        is_protein=False,
        atom_positions=atom_positions,
        atom_elements=atom_elements,
        atom_charges=atom_charges,
        bonds=bonds,
        name=Path(input_mol).stem if isinstance(input_mol, (str, Path)) and Path(input_mol).exists() else None
    )


def rotate_molecule_for_screen(atom_positions):
    '''
    Rotate coordinates using PCA so longest axis is X, second is Y, shortest is Z.
    Avoids molecules lying flat in the camera view.
    '''
    if len(atom_positions) < 3:
        return atom_positions

    mean_centered = atom_positions - np.mean(atom_positions, axis=0)
    cov_matrix = np.cov(mean_centered, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
    sorted_indices = np.argsort(eigenvalues)[::-1]
    sorted_eigenvectors = eigenvectors[:, sorted_indices]
    rotated_positions = np.dot(mean_centered, sorted_eigenvectors)
    return rotated_positions


def scale_for_canvas(atom_positions, canvas):
    '''
    Given atom positions and canvas object, center and scale coordinates to fill screen space.
    '''
    if len(atom_positions) == 0:
        return atom_positions

    # Always ensure centroid is at (0, 0)
    center = np.mean(atom_positions[:, :2], axis=0)
    atom_positions[:, 0] -= center[0]
    atom_positions[:, 1] -= center[1]

    min_x, max_x = np.min(atom_positions[:, 0]), np.max(atom_positions[:, 0])
    min_y, max_y = np.min(atom_positions[:, 1]), np.max(atom_positions[:, 1])

    span_x = max(max_x - min_x, 1e-4)
    span_y = max(max_y - min_y, 1e-4)

    x_scaling_factor = canvas.width / span_x
    y_scaling_factor = canvas.height / span_y

    atom_positions *= 0.85 * min(x_scaling_factor, canvas.aspect_ratio * y_scaling_factor)
    return atom_positions


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation
    about the given axis by theta radians.
    """
    axis = np.asarray(axis, dtype=float)
    norm = np.sqrt(np.dot(axis, axis))
    if norm == 0:
        return np.eye(3)
    axis = axis / norm
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