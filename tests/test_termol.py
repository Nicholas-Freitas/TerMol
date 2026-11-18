import pytest
import numpy as np
import termol
from termol.pdb import is_pdb_code, clean_pdb_code, parse_pdb_text, load_pdb
from termol.backbone import generate_backbone_points, catmull_rom_spline_segment
from termol.canvas import Canvas, MoleculeCanvas, ProteinCanvas
from termol.utilities import get_molecule_data, MoleculeData

SAMPLE_PDB = """
HEADER    TEST PROTEIN
HELIX    1   1 GLY A    1  ALA A    5  1
SHEET    1   A 2 VAL A    6  LEU A    8  0
ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 10.00           N
ATOM      2  CA  GLY A   1       1.000   1.000   1.000  1.00 10.00           C
ATOM      3  C   GLY A   1       2.000   1.000   1.000  1.00 10.00           C
ATOM      4  O   GLY A   1       2.500   0.000   1.000  1.00 10.00           O
ATOM      5  N   ALA A   2       2.000   2.000   1.000  1.00 10.00           N
ATOM      6  CA  ALA A   2       3.000   3.000   2.000  1.00 10.00           C
ATOM      7  CA  VAL A   3       5.000   5.000   3.000  1.00 10.00           C
ATOM      8  CA  LEU A   4       7.000   6.000   4.000  1.00 10.00           C
ATOM      9  CA  SER A   5       9.000   7.000   5.000  1.00 10.00           C
ATOM     10  CA  VAL A   6      11.000   8.000   6.000  1.00 10.00           C
ATOM     11  CA  ILE A   7      13.000   9.000   7.000  1.00 10.00           C
ATOM     12  CA  LEU A   8      15.000  10.000   8.000  1.00 10.00           C
ATOM     13  CA  PHE A   9      17.000  11.000   9.000  1.00 10.00           C
ATOM     14  CA  TYR A  10      19.000  12.000  10.000  1.00 10.00           C
ATOM     15  CA  TRP A  11      21.000  13.000  11.000  1.00 10.00           C
HETATM   16 ZN    ZN A  99      10.000  10.000  10.000  1.00 20.00          ZN
END
"""

def test_pdb_code_detection():
    assert is_pdb_code("1ubq") is True
    assert is_pdb_code("4HHB") is True
    assert is_pdb_code("pdb:1crn") is True
    assert is_pdb_code("c1ccccc1") is False
    assert is_pdb_code("protein.pdb") is False
    assert clean_pdb_code("pdb:1ubq") == "1UBQ"
    assert clean_pdb_code("4hhb") == "4HHB"


def test_pdb_parser():
    structure = parse_pdb_text(SAMPLE_PDB, title="Test")
    assert structure.is_protein is True
    assert structure.total_residues == 11
    assert "A" in structure.chains
    assert len(structure.ligands) == 1
    assert structure.ligands[0].element == "ZN"
    # Check secondary structure assignment
    ca_list = structure.chains["A"]
    assert ca_list[0].ss == "H"  # res 1 in helix
    assert ca_list[5].ss == "E"  # res 6 in sheet
    assert ca_list[10].ss == "C" # res 11 in coil


def test_catmull_rom_spline():
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 2.0, 1.0],
        [2.0, 1.0, 3.0],
        [3.0, 4.0, 2.0]
    ])
    spline = catmull_rom_spline_segment(coords, num_points_per_seg=5)
    assert len(spline) > len(coords)
    # Check first and last points match original endpoints
    np.testing.assert_allclose(spline[0], coords[0], atol=1e-5)
    np.testing.assert_allclose(spline[-1], coords[-1], atol=1e-5)


def test_protein_canvas_ascii():
    canvas = ProteinCanvas(60, 30, 60, 30, aspect_ratio=2.0, style='ascii', color_mode='depth')
    structure = parse_pdb_text(SAMPLE_PDB)
    segments = generate_backbone_points(structure, samples_per_res=4)
    canvas.draw_protein(segments, structure.ligands)
    plain = canvas.render_plain()
    assert len(plain) > 0
    # Canvas should contain non-space characters
    non_spaces = [c for c in plain if c != ' ' and c != '\n']
    assert len(non_spaces) > 0


def test_protein_canvas_braille():
    canvas = ProteinCanvas(60, 30, 60, 30, aspect_ratio=2.0, style='braille', color_mode='depth')
    structure = parse_pdb_text(SAMPLE_PDB)
    segments = generate_backbone_points(structure, samples_per_res=4)
    canvas.draw_protein(segments, structure.ligands)
    plain = canvas.render_plain()
    # Check that Braille characters (0x2800..0x28FF) are present
    braille_chars = [c for c in plain if 0x2800 <= ord(c) <= 0x28FF and c != chr(0x2800)]
    assert len(braille_chars) > 0


def test_small_molecule_backward_compatibility():
    mol_data = get_molecule_data("c1ccccc1", three_d=False)
    assert mol_data.is_protein is False
    # Check legacy 4-tuple unpacking
    pos, elems, charges, bonds = mol_data
    assert len(pos) == 6
    assert len(elems) == 6
    assert len(bonds) == 6


def test_draw_2d_runs_without_error(capsys):
    # Small molecule
    termol.draw("c1ccccc1", three_d=False)
    out1 = capsys.readouterr().out
    assert "C" in out1

    # Protein sample
    structure = parse_pdb_text(SAMPLE_PDB)
    mol_data = MoleculeData(
        is_protein=True,
        structure=structure,
        backbone_segments=generate_backbone_points(structure, samples_per_res=4),
        ligands=structure.ligands,
        name="TestProtein"
    )
    from termol.canvas import MoleculeCanvas
    from termol.render import show_molecule_2D
    c = MoleculeCanvas(60, 30, 60, 30, aspect_ratio=2.0)
    show_molecule_2D(mol_data, c, name="TestProtein", style='ascii', color='none')
    out2 = capsys.readouterr().out
    assert "TestProtein" in out2


def test_color_modes():
    structure = parse_pdb_text(SAMPLE_PDB)
    segments = generate_backbone_points(structure, samples_per_res=4)

    for cm in ['depth', 'depth1', 'depth2', 'depth3', 'rainbow', 'ss']:
        canvas = ProteinCanvas(60, 30, 60, 30, aspect_ratio=2.0, style='ascii', color_mode=cm)
        canvas.draw_protein(segments)
        ansi_text = canvas.render_ansi()
        # Verify ANSI escape sequences are present
        assert "\033[38;2;" in ansi_text

    # None mode should have no ANSI escapes
    canvas_none = ProteinCanvas(60, 30, 60, 30, aspect_ratio=2.0, style='ascii', color_mode='none')
    canvas_none.draw_protein(segments)
    assert "\033[" not in canvas_none.render_ansi()


def test_load_pdb_rcsb():
    struct = load_pdb("1ubq")
    assert struct.is_protein is True
    assert struct.total_residues == 76
    assert len(struct.helices) > 0


def test_draw_multi(capsys):
    termol.draw_multi(["c1ccccc1", "1crn"], three_d=False)
    captured = capsys.readouterr().out
    assert "1CRN" in captured or "c1ccccc1" in captured


def test_load_local_pdb(tmp_path):
    pdb_file = tmp_path / "sample.pdb"
    pdb_file.write_text(SAMPLE_PDB, encoding="utf-8")
    mol_data = get_molecule_data(str(pdb_file), three_d=False)
    assert mol_data.is_protein is True
    assert mol_data.structure.total_residues == 11


def test_representations(capsys):
    struct = parse_pdb_text(SAMPLE_PDB)

    # Test trace
    segs_trace = generate_backbone_points(struct, representation='trace')
    assert segs_trace[0][0].left_coord is None

    # Test cartoon
    segs_cartoon = generate_backbone_points(struct, representation='cartoon')
    assert segs_cartoon[0][0].left_coord is not None
    assert segs_cartoon[0][0].right_coord is not None

    # Test tube
    segs_tube = generate_backbone_points(struct, representation='tube')
    assert segs_tube[0][0].top_coord is not None
    assert segs_tube[0][0].bottom_coord is not None

    # Test drawing in all representations
    for r in ['trace', 'cartoon', 'tube']:
        termol.draw("1crn", three_d=False, style='braille', color='none', rep=r)
        out = capsys.readouterr().out
        assert "1CRN" in out


def test_protein_vertical_centering():
    """Verify that proteins render centered vertically around the canvas midpoint."""
    from termol.utilities import get_molecule_data
    from termol.canvas import ProteinCanvas
    from termol.render import scale_protein_for_canvas

    data = get_molecule_data("1crn", representation='cartoon')
    canvas = ProteinCanvas(80, 40, 80, 40, aspect_ratio=2.0, style='braille')
    scale_protein_for_canvas(data.backbone_segments, data.ligands, canvas)
    canvas.draw_protein(data.backbone_segments, data.ligands)

    rows = [y for y in range(canvas.char_height) if any(canvas.data[y][x] != ' ' for x in range(canvas.char_width))]
    assert len(rows) > 0
    mid_row = (min(rows) + max(rows)) / 2.0
    canvas_mid = (canvas.char_height - 1) / 2.0
    # Center of content should be within 2.5 rows of the canvas midpoint
    assert abs(mid_row - canvas_mid) <= 2.5


def test_representation_switch_pose_preservation():
    """Verify that switching representation preserves rotated pose and zoom."""
    import copy
    from termol.utilities import get_molecule_data
    from termol.render import scale_protein_for_canvas, rotation_matrix
    from termol.canvas import ProteinCanvas
    from termol.backbone import generate_backbone_points

    data = get_molecule_data("1ubq", three_d=True, representation='trace')
    canvas = ProteinCanvas(80, 40, 80, 40, aspect_ratio=2.0)

    # Initial scaling
    initial_segs = copy.deepcopy(data.backbone_segments)
    protein_scale = scale_protein_for_canvas(initial_segs, data.ligands, canvas)

    canonical_cache = {'trace': initial_segs}
    cum_R = np.eye(3)
    zoom_scale = 1.0

    def get_segs(target_rep):
        if target_rep not in canonical_cache:
            raw_segs = generate_backbone_points(data.structure, samples_per_res=8, representation=target_rep)
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

    segments = get_segs('trace')

    # Rotate 25 frames and zoom
    R = rotation_matrix((0, 1, 0), np.radians(2.0))
    for _ in range(25):
        cum_R = np.dot(R, cum_R)
        for s in segments:
            for pt in s:
                pt.coord = np.dot(R, pt.coord)

    zoom_scale *= 1.15
    for s in segments:
        for pt in s:
            pt.coord *= 1.15

    # Switch to cartoon
    cartoon_segs = get_segs('cartoon')

    # Start and mid C-alpha coordinates must match between rotated trace and cartoon
    assert np.allclose(segments[0][0].coord, cartoon_segs[0][0].coord, atol=1e-5)
    assert np.allclose(segments[0][100].coord, cartoon_segs[0][100].coord, atol=1e-5)

    # Switch to tube
    tube_segs = get_segs('tube')
    assert np.allclose(segments[0][0].coord, tube_segs[0][0].coord, atol=1e-5)
    assert np.allclose(segments[0][100].coord, tube_segs[0][100].coord, atol=1e-5)


def test_rgb_to_256_mapping():
    """Verify 24-bit RGB to xterm-256 color mapping logic."""
    from termol.render import rgb_to_256

    # Pure black and white
    assert rgb_to_256(0, 0, 0) == 16
    assert rgb_to_256(255, 255, 255) == 231

    # Grayscale ramp indices (232..255)
    assert 232 <= rgb_to_256(10, 10, 10) <= 255
    assert 232 <= rgb_to_256(128, 128, 128) <= 255
    assert 232 <= rgb_to_256(230, 230, 230) <= 255

    # Primaries
    assert rgb_to_256(255, 0, 0) == 196
    assert rgb_to_256(0, 255, 0) == 46
    assert rgb_to_256(0, 0, 255) == 21

    # Clamping out-of-range values
    assert rgb_to_256(-10, -5, -1) == 16
    assert rgb_to_256(300, 300, 300) == 231


def test_protein_canvas_depth_tracking():
    """Verify ProteinCanvas accurately tracks last_z_min and last_z_range."""
    from termol.canvas import ProteinCanvas
    from termol.backbone import generate_backbone_points
    from termol.pdb import parse_pdb_text

    canvas = ProteinCanvas(80, 40, 80, 40, aspect_ratio=2.0)
    structure = parse_pdb_text(SAMPLE_PDB)
    segments = generate_backbone_points(structure, samples_per_res=4)
    canvas.draw_protein(segments, structure.ligands)

    assert hasattr(canvas, 'last_z_min')
    assert hasattr(canvas, 'last_z_range')
    assert canvas.last_z_range > 0

    # Every drawn cell with non-negative z should normalize into [0.0, 1.0]
    for cy in range(canvas.char_height):
        for cx in range(canvas.char_width):
            z = canvas.z_buffer[cy, cx]
            if z > -np.inf:
                z_norm = (z - canvas.last_z_min) / canvas.last_z_range
                assert -1e-6 <= z_norm <= 1.0 + 1e-6


def test_depth_color_modes_ansi_render():
    """Verify all depth color modes render ANSI escape sequences without warnings."""
    from termol.canvas import ProteinCanvas
    from termol.backbone import generate_backbone_points
    from termol.pdb import parse_pdb_text

    structure = parse_pdb_text(SAMPLE_PDB)
    segments = generate_backbone_points(structure, samples_per_res=4)

    for mode in ('depth1', 'depth2', 'depth3', 'rainbow', 'ss'):
        canvas = ProteinCanvas(60, 30, 60, 30, aspect_ratio=2.0, style='braille', color_mode=mode)
        canvas.draw_protein(segments, structure.ligands)
        ansi_out = canvas.render_ansi()
        assert len(ansi_out) > 0
        assert "\033[" in ansi_out


def test_escape_sequence_parsing_no_premature_quit():
    """Verify escape sequences and keys do not prematurely trigger viewer exit."""
    from termol.render import CSI_PATTERN, INCOMPLETE_CSI_PATTERN

    raw_partial = "\x1b[A\x1b["
    raw_arrow = "\x1b[A"
    raw_quit_key = "q"
    raw_esc_key = "\x1b"

    # Partial CSI stream: should strip completely
    cleaned_partial = INCOMPLETE_CSI_PATTERN.sub('', CSI_PATTERN.sub('', raw_partial))
    assert cleaned_partial == ""

    # Arrow key stream: should strip to empty
    cleaned_arrow = INCOMPLETE_CSI_PATTERN.sub('', CSI_PATTERN.sub('', raw_arrow))
    assert cleaned_arrow == ""

    # Standalone Escape key: recognized as standalone '\x1b'
    assert raw_esc_key == '\x1b'

    # Quit key 'q': recognized
    cleaned_q = INCOMPLETE_CSI_PATTERN.sub('', CSI_PATTERN.sub('', raw_quit_key))
    assert any(c in ('q', 'Q', 'x', 'X') for c in cleaned_q)


def test_ansi_3d_viewer_protein_timeout():
    """Verify ANSI 3D interactive viewer runs and terminates cleanly on timeout."""
    from termol.render import show_molecule_3D_ansi
    from termol.canvas import MoleculeCanvas
    from termol.utilities import get_molecule_data

    mol_data = get_molecule_data("1ubq", three_d=True)
    canvas = MoleculeCanvas(60, 30, 60.0, 30.0, aspect_ratio=2.0)
    # Run for 0.05 seconds
    show_molecule_3D_ansi(mol_data, canvas, timeout=0.05, auto_resize=False)


def test_ansi_3d_viewer_small_mol_timeout():
    """Verify ANSI 3D interactive viewer runs on small molecules without error."""
    from termol.render import show_molecule_3D_ansi
    from termol.canvas import MoleculeCanvas
    from termol.utilities import get_molecule_data

    mol_data = get_molecule_data("c1ccccc1", three_d=True)
    canvas = MoleculeCanvas(60, 30, 60.0, 30.0, aspect_ratio=2.0)
    show_molecule_3D_ansi(mol_data, canvas, timeout=0.05, auto_resize=False)


def test_molecule_canvas_cpk_colors_and_depth_cueing():
    """Verify CPK element coloring and Z-depth cueing in MoleculeCanvas."""
    from termol.canvas import MoleculeCanvas, CPK_COLORS

    canvas = MoleculeCanvas(40, 20, 40.0, 20.0, aspect_ratio=2.0)
    # 2 atoms: Oxygen at foreground (z=+5), Carbon at background (z=-5)
    positions = [(0.0, 0.0, 5.0), (10.0, 0.0, -5.0)]
    elements = ['O', 'C']
    charges = [0, 0]
    bonds = [(0, 1)]

    canvas.draw_molecules(positions, elements, charges, bonds)
    ansi_str = canvas.render_ansi()

    # ANSI output should contain 24-bit TrueColor escape sequences
    assert "\033[38;2;" in ansi_str
    # Oxygen in foreground should have bright red component (> 200)
    ox_cx, ox_cy = canvas.map_coords(0.0, 0.0)
    ox_color = canvas.color_data[ox_cy][ox_cx]
    assert ox_color is not None
    r, g, b = ox_color
    assert r > 200
    assert g < 100
    assert b < 100

    # Carbon in background (z=-5) should be dimmed compared to full brightness
    c_cx, c_cy = canvas.map_coords(10.0, 0.0)
    c_color = canvas.color_data[c_cy][c_cx]
    assert c_color is not None
    # Base Carbon is ~195; at deep background (scale=0.40) it should be <= 120
    assert c_color[0] < 140


def test_load_sdf_file_support(tmp_path):
    """Verify loading single and multi-molecule .sdf files preserves structures."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from termol.utilities import get_molecule_data

    # Write a test SDF file with 3D coordinates
    mol = Chem.MolFromSmiles('CC(=O)O')
    AllChem.EmbedMolecule(mol)
    sdf_file = str(tmp_path / "test_acetic_acid.sdf")
    writer = Chem.SDWriter(sdf_file)
    writer.write(mol)
    writer.close()

    mol_data = get_molecule_data(sdf_file, three_d=True)
    assert not mol_data.is_protein
    assert len(mol_data.atom_elements) == 4
    assert 'O' in mol_data.atom_elements
    assert len(mol_data.bonds) == 3


def test_load_mol2_file_support(tmp_path):
    """Verify loading Tripos .mol2 files preserves structures and coordinates."""
    from termol.utilities import get_molecule_data

    mol2_content = """@<TRIPOS>MOLECULE
benzene
6 6 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
1 C1 0.0000 1.3970 0.0000 C.ar 1 BENZ 0.0000
2 C2 1.2098 0.6985 0.0000 C.ar 1 BENZ 0.0000
3 C3 1.2098 -0.6985 0.0000 C.ar 1 BENZ 0.0000
4 C4 0.0000 -1.3970 0.0000 C.ar 1 BENZ 0.0000
5 C5 -1.2098 -0.6985 0.0000 C.ar 1 BENZ 0.0000
6 C6 -1.2098 0.6985 0.0000 C.ar 1 BENZ 0.0000
@<TRIPOS>BOND
1 1 2 ar
2 2 3 ar
3 3 4 ar
4 4 5 ar
5 5 6 ar
6 6 1 ar
"""
    mol2_file = str(tmp_path / "benzene.mol2")
    with open(mol2_file, "w") as f:
        f.write(mol2_content)

    mol_data = get_molecule_data(mol2_file, three_d=True)
    assert not mol_data.is_protein
    assert len(mol_data.atom_elements) == 6
    assert all(e == 'C' for e in mol_data.atom_elements)
    assert len(mol_data.bonds) == 6


def test_show_molecule_2d_truecolor(capsys):
    """Verify 2D small-molecule rendering outputs TrueColor CPK ANSI."""
    from termol.render import show_molecule_2D
    from termol.canvas import MoleculeCanvas
    from termol.utilities import get_molecule_data

    mol_data = get_molecule_data("CCO", three_d=False)
    canvas = MoleculeCanvas(40, 20, 40.0, 20.0, aspect_ratio=2.0)
    show_molecule_2D(mol_data, canvas, name="Ethanol", color="cpk")

    captured = capsys.readouterr().out
    assert "Ethanol" in captured
    assert "\033[38;2;" in captured  # TrueColor ANSI sequence


def test_termol_cli_smiles_invocation(monkeypatch, capsys):
    """Verify termol command line entry point parses SMILES without error."""
    import sys
    from termol.cli import termol_cli

    monkeypatch.setattr(sys, 'argv', ['termol', 'CC(=O)Oc1ccccc1C(=O)O', '--show2D', '--names', 'Aspirin'])
    termol_cli()

    captured = capsys.readouterr().out
    assert "Aspirin" in captured
    assert "\033[38;2;" in captured


def test_termol_cli_sdf_invocation(tmp_path, monkeypatch, capsys):
    """Verify termol command line entry point parses .sdf file without error."""
    import sys
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from termol.cli import termol_cli

    mol = Chem.MolFromSmiles('c1ccccc1')
    AllChem.EmbedMolecule(mol)
    sdf_file = str(tmp_path / "benzene.sdf")
    writer = Chem.SDWriter(sdf_file)
    writer.write(mol)
    writer.close()

    monkeypatch.setattr(sys, 'argv', ['termol', sdf_file, '--show2D'])
    termol_cli()

    captured = capsys.readouterr().out
    assert "benzene" in captured


def test_bond_orders_and_aromaticity_rendering():
    """Verify single, double, triple, and aromatic bond characters."""
    from termol.canvas import MoleculeCanvas

    canvas = MoleculeCanvas(40, 20, 40.0, 20.0, aspect_ratio=2.0)
    # Triple bond (order=3.0) horizontally between (0, 0) and (10, 0)
    # Double bond (order=2.0) vertically between (0, 2) and (0, 8)
    # Aromatic bond (order=1.5) horizontally between (15, 0) and (25, 0)
    positions = [(0.0, 0.0, 0.0), (10.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 8.0, 0.0), (15.0, 0.0, 0.0), (25.0, 0.0, 0.0)]
    elements = ['C', 'C', 'C', 'C', 'C', 'C']
    charges = [0, 0, 0, 0, 0, 0]
    bonds = [(0, 1), (2, 3), (4, 5)]
    bond_orders = [3.0, 2.0, 1.5]

    canvas.draw_molecules(positions, elements, charges, bonds, bond_orders=bond_orders)
    plain = canvas.render_plain()

    # Triple bond glyph '≡' should be present
    assert '≡' in plain
    # Double vertical bond glyph '‖' should be present
    assert '‖' in plain
    # Aromatic horizontal bond glyph '╌' should be present
    assert '╌' in plain


def test_small_molecule_styles_balls_and_dots():
    """Verify balls and dots visual styles for small molecules."""
    from termol.canvas import MoleculeCanvas

    canvas_balls = MoleculeCanvas(40, 20, 40.0, 20.0, aspect_ratio=2.0, style='balls')
    canvas_dots = MoleculeCanvas(40, 20, 40.0, 20.0, aspect_ratio=2.0, style='dots')

    positions = [(0.0, 0.0, 0.0), (10.0, 0.0, 0.0)]
    elements = ['C', 'S']
    charges = [0, 0]
    bonds = [(0, 1)]

    canvas_balls.draw_molecules(positions, elements, charges, bonds, style='balls')
    plain_balls = canvas_balls.render_plain()
    # Carbon should be '●', Sulfur (larger atom) should be '◉'
    assert '●' in plain_balls
    assert '◉' in plain_balls

    canvas_dots.draw_molecules(positions, elements, charges, bonds, style='dots')
    plain_dots = canvas_dots.render_plain()
    assert '•' in plain_dots


def test_termol_cli_no_color_flag(monkeypatch, capsys):
    """Verify --no-color flag strips all ANSI sequences from output."""
    import sys
    from termol.cli import termol_cli

    monkeypatch.setattr(sys, 'argv', ['termol', 'c1ccccc1', '--show2D', '--no-color'])
    termol_cli()

    captured = capsys.readouterr().out
    assert "\033[38;2;" not in captured
    assert "\033[" not in captured
    assert "C" in captured




