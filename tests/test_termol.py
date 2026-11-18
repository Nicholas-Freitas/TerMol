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

    for cm in ['depth', 'rainbow', 'ss']:
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
