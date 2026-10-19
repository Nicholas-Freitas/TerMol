# termol/pdb.py
"""
Zero-dependency PDB downloader and parser for macromolecules and proteins.
Extracts backbone CA coordinates, chain boundaries, secondary structure,
and non-water ligands using only Python standard library and NumPy.
"""

import os
import re
import urllib.request
import urllib.error
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
import numpy as np

# Standard amino acid 3-letter codes
STANDARD_AMINO_ACIDS = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    # Common modified amino acids
    'MSE', 'SEC', 'PYL', 'HYP', 'SEP', 'TPO', 'PTR'
}

COMMON_NUCLEIC_ACIDS = {
    'A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT', 'DI', 'I'
}

WATER_RESIDUES = {'HOH', 'WAT', 'H2O', 'DOD', 'TIP', 'TIP3'}


@dataclass
class ResidueCA:
    chain_id: str
    res_num: int
    res_name: str
    coord: np.ndarray  # shape (3,)
    ss: str = 'C'       # 'H' (helix), 'E' (sheet/strand), 'C' (coil/loop)
    b_factor: float = 0.0


@dataclass
class LigandAtom:
    element: str
    res_name: str
    res_num: int
    chain_id: str
    coord: np.ndarray  # shape (3,)


@dataclass
class ProteinStructure:
    chains: Dict[str, List[ResidueCA]] = field(default_factory=dict)
    ligands: List[LigandAtom] = field(default_factory=list)
    helices: List[Tuple[str, int, int]] = field(default_factory=list)  # (chain, start, end)
    sheets: List[Tuple[str, int, int]] = field(default_factory=list)   # (chain, start, end)
    title: str = ""
    is_protein: bool = False

    @property
    def total_residues(self) -> int:
        return sum(len(ca_list) for ca_list in self.chains.values())

    @property
    def all_ca_coords(self) -> np.ndarray:
        coords = []
        for ca_list in self.chains.values():
            for res in ca_list:
                coords.append(res.coord)
        if not coords:
            return np.empty((0, 3), dtype=float)
        return np.array(coords)


def is_pdb_code(query: str) -> bool:
    """Check if query is a 4-character RCSB PDB code (e.g. '1ubq', '4HHB')."""
    if not isinstance(query, str):
        return False
    query = query.strip()
    return bool(re.match(r'^[1-9][0-9a-zA-Z]{3}$', query)) or bool(re.match(r'^(pdb|rcsb):[1-9][0-9a-zA-Z]{3}$', query, re.IGNORECASE))


def clean_pdb_code(query: str) -> str:
    """Normalize PDB code to uppercase 4 characters."""
    query = query.strip()
    if ':' in query:
        query = query.split(':', 1)[1]
    return query.upper()


def fetch_pdb_from_rcsb(pdb_id: str, cache: bool = True) -> str:
    """
    Fetch PDB file text from RCSB PDB database.
    Caches downloaded files to ~/.cache/termol/{PDB_ID}.pdb
    """
    clean_id = clean_pdb_code(pdb_id)
    cache_dir = Path.home() / ".cache" / "termol"
    cache_path = cache_dir / f"{clean_id}.pdb"

    if cache and cache_path.exists():
        try:
            return cache_path.read_text(encoding='utf-8')
        except Exception:
            pass

    url = f"https://files.rcsb.org/download/{clean_id}.pdb"
    req = urllib.request.Request(url, headers={'User-Agent': 'TerMol/0.2 (Terminal Molecular Viewer)'})
    try:
        with urllib.request.urlopen(req, timeout=12) as response:
            pdb_text = response.read().decode('utf-8', errors='replace')
    except urllib.error.HTTPError as e:
        raise ValueError(f"Failed to fetch PDB '{clean_id}' from RCSB (HTTP {e.code})") from e
    except urllib.error.URLError as e:
        raise ValueError(f"Network error fetching PDB '{clean_id}': {e.reason}") from e

    if cache:
        try:
            cache_dir.mkdir(parents=True, exist_ok=True)
            cache_path.write_text(pdb_text, encoding='utf-8')
        except Exception:
            pass

    return pdb_text


def parse_pdb_text(pdb_text: str, title: str = "") -> ProteinStructure:
    """
    Parse PDB text and extract protein chains, CA atoms, secondary structures,
    and ligands without external dependencies.
    """
    structure = ProteinStructure(title=title)
    lines = pdb_text.splitlines()

    chains: Dict[str, List[ResidueCA]] = {}
    ligands: List[LigandAtom] = []
    helices: List[Tuple[str, int, int]] = []
    sheets: List[Tuple[str, int, int]] = []
    header_title = ""

    seen_ca: set = set()  # (chain, res_num, alt_loc)

    for line in lines:
        if line.startswith('TITLE') and not header_title:
            header_title = line[10:].strip()

        elif line.startswith('HELIX'):
            parsed = False
            try:
                ch = line[19].strip() or 'A'
                start_res = int(line[21:25].strip())
                end_res = int(line[33:37].strip())
                helices.append((ch, start_res, end_res))
                parsed = True
            except (ValueError, IndexError):
                pass
            if not parsed:
                # Regex fallback for non-standard spaced PDBs
                m = re.search(r'HELIX\s+\S+\s+\S+\s+\S+\s+(\S)\s+(\d+)\s+\S+\s+(\S)\s+(\d+)', line)
                if m:
                    helices.append((m.group(1), int(m.group(2)), int(m.group(4))))

        elif line.startswith('SHEET'):
            parsed = False
            try:
                ch = line[21].strip() or 'A'
                start_res = int(line[22:26].strip())
                end_res = int(line[33:37].strip())
                sheets.append((ch, start_res, end_res))
                parsed = True
            except (ValueError, IndexError):
                pass
            if not parsed:
                # Regex fallback for non-standard spaced PDBs
                m = re.search(r'SHEET\s+\S+\s+\S+\s+\d+\s+\S+\s+(\S)\s+(\d+)\s+\S+\s+(\S)\s+(\d+)', line)
                if m:
                    sheets.append((m.group(1), int(m.group(2)), int(m.group(4))))

        elif line.startswith('ATOM') or line.startswith('HETATM'):
            atom_name = line[12:16].strip()
            alt_loc = line[16]
            res_name = line[17:20].strip()
            chain_id = line[21].strip() or 'A'

            # Ignore non-default alternate locations
            if alt_loc not in (' ', 'A', '1'):
                continue

            try:
                res_num = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                b_factor = float(line[60:66].strip()) if len(line) >= 66 else 0.0
                element = line[76:78].strip() if len(line) >= 78 else atom_name[:1]
            except (ValueError, IndexError):
                continue

            coord = np.array([x, y, z], dtype=float)

            # Is it a C-alpha atom of an amino acid (or C4'/P of a nucleic acid)?
            if atom_name in ('CA', 'C4*', "C4'", 'P'):
                key = (chain_id, res_num)
                if key not in seen_ca:
                    seen_ca.add(key)
                    if chain_id not in chains:
                        chains[chain_id] = []
                    chains[chain_id].append(ResidueCA(
                        chain_id=chain_id,
                        res_num=res_num,
                        res_name=res_name,
                        coord=coord,
                        b_factor=b_factor
                    ))

            # Non-water ligand atom
            elif line.startswith('HETATM') and res_name not in WATER_RESIDUES:
                ligands.append(LigandAtom(
                    element=element,
                    res_name=res_name,
                    res_num=res_num,
                    chain_id=chain_id,
                    coord=coord
                ))

    structure.chains = chains
    structure.ligands = ligands
    structure.helices = helices
    structure.sheets = sheets
    if not structure.title and header_title:
        structure.title = header_title

    # Determine if this structure qualifies as a protein/macromolecule
    total_ca = structure.total_residues
    structure.is_protein = total_ca >= 10

    # Assign secondary structure to each residue
    if structure.is_protein:
        _assign_secondary_structures(structure)

    return structure


def _assign_secondary_structures(structure: ProteinStructure):
    """
    Annotate each CA residue with secondary structure ('H', 'E', or 'C').
    Prefers PDB header HELIX/SHEET records; falls back to CA geometry.
    """
    has_header_ss = bool(structure.helices or structure.sheets)

    for ch, ca_list in structure.chains.items():
        if has_header_ss:
            # Map by residue number
            helix_set = set()
            for h_ch, h_start, h_end in structure.helices:
                if h_ch == ch:
                    for r in range(h_start, h_end + 1):
                        helix_set.add(r)

            sheet_set = set()
            for s_ch, s_start, s_end in structure.sheets:
                if s_ch == ch:
                    for r in range(s_start, s_end + 1):
                        sheet_set.add(r)

            for res in ca_list:
                if res.res_num in helix_set:
                    res.ss = 'H'
                elif res.res_num in sheet_set:
                    res.ss = 'E'
                else:
                    res.ss = 'C'
        else:
            # Geometric assignment from CA-CA distance profile
            coords = np.array([r.coord for r in ca_list])
            n = len(coords)
            ss_classes = ['C'] * n
            for i in range(n):
                if i >= 1 and i + 2 < n:
                    d_i_i2 = np.linalg.norm(coords[i+2] - coords[i])
                    if i + 3 < n:
                        d_i_i3 = np.linalg.norm(coords[i+3] - coords[i])
                        if 4.7 < d_i_i3 < 6.2 and 4.9 < d_i_i2 < 6.0:
                            ss_classes[i] = 'H'
                            continue
                    if d_i_i2 > 6.4:
                        ss_classes[i] = 'E'

            # Smooth single-residue noise
            for i in range(1, n - 1):
                if ss_classes[i] == 'H' and ss_classes[i-1] != 'H' and ss_classes[i+1] != 'H':
                    ss_classes[i] = 'C'
                elif ss_classes[i] == 'E' and ss_classes[i-1] != 'E' and ss_classes[i+1] != 'E':
                    ss_classes[i] = 'C'

            for res, ss_code in zip(ca_list, ss_classes):
                res.ss = ss_code


def load_pdb(input_target: str) -> ProteinStructure:
    """
    Load a PDB structure from a file path or 4-character RCSB PDB ID.
    """
    # 1. Is it a 4-character PDB code?
    if is_pdb_code(input_target):
        pdb_id = clean_pdb_code(input_target)
        pdb_text = fetch_pdb_from_rcsb(pdb_id)
        return parse_pdb_text(pdb_text, title=pdb_id)

    # 2. Is it a local file?
    path = Path(input_target)
    if path.exists() and path.is_file():
        pdb_text = path.read_text(encoding='utf-8', errors='replace')
        return parse_pdb_text(pdb_text, title=path.stem)

    raise ValueError(f"Input '{input_target}' is neither an existing file nor a valid PDB ID.")
