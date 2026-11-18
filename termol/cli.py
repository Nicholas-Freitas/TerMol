import argparse
import csv
from pathlib import Path
from .render import draw_multi, showcase

def termol_cli():
    parser = argparse.ArgumentParser(
        description='TerMol: A terminal molecular renderer for small molecules and proteins.'
    )
    parser.add_argument(
        'molecule_inputs', nargs='+',
        help='One or more SMILES strings, PDB accession codes (e.g. 1ubq), or molecule files (.sdf, .mol, .mol2, .pdb). '
             'Separate with spaces. Can also be a text file (.txt) or CSV file (name, smiles/code/file).'
    )
    parser.add_argument('--names', nargs='+', help='Names for each molecule. Separate with spaces. Must be same length as inputs.')
    parser.add_argument('--show2D', action='store_true', help='Render static 2D view (default interactive 3D).')
    parser.add_argument('--width', type=int, help='Width of the output display in characters (default 80, or auto-resized in 3D).')
    parser.add_argument('--height', type=int, help='Height of the output display in characters (default 40, or auto-resized in 3D).')
    parser.add_argument('--add_hydrogens', action='store_true', help='Add hydrogens to small molecules via RDKit (default False).')
    parser.add_argument('--timeout', type=int, default=None, help='3D window will close after this number of seconds (default indefinitely).')
    parser.add_argument(
        '--style', default=None,
        help='Rendering style: for small molecules: "symbols" (chemical symbols, default), "balls", "dots", or "braille" (high-res subpixel line drawing); '
             'for proteins: "braille" (subpixel dots, default) or "ascii".'
    )
    parser.add_argument(
        '--bond-spacing', '--spacing', dest='bond_spacing', type=float, default=3.0,
        help='Spacing between parallel lines for double/triple/aromatic bonds in braille mode (default 3.0 subpixels).'
    )
    parser.add_argument(
        '--color', default=None,
        help='Color mode: for small molecules: "cpk" (default) or "none"; '
             'for proteins: "depth"/"depth1" (default), "depth2", "depth3", "rainbow", "ss", or "none".'
    )
    parser.add_argument(
        '--no-color', action='store_true',
        help='Disable ANSI TrueColor escape sequences and render in monochrome text.'
    )
    parser.add_argument(
        '--rep', '--representation', dest='rep', choices=['trace', 'cartoon', 'ribbon', 'tube'], default='trace',
        help='Representation for proteins: "trace" (1D spline line, default), "cartoon"/"ribbon" (secondary structure ribbons & arrowheads), or "tube" (3D uniform thickness).'
    )
    parser.add_argument(
        '--backend', choices=['ansi', 'curses'], default='ansi',
        help='Interactive 3D rendering backend: "ansi" (24-bit TrueColor with responsive keyboard navigation, default) or "curses" (legacy curses).'
    )
    args = parser.parse_args()

    if args.no_color:
        args.color = 'none'

    width_default = 80
    height_default = 40
    if args.width is None and args.height is None:
        auto_resize = True
        args.width = width_default
        args.height = height_default
    else:
        if args.width is None:
            args.width = width_default
        if args.height is None:
            args.height = height_default
        auto_resize = False

    names_list = None
    if args.names:
        names_list = args.names

    molecules_list = []
    for input_item in args.molecule_inputs:
        if input_item.endswith('.txt'):
            with open(input_item, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    val = line.strip().removeprefix('"').removesuffix('"')
                    if val:
                        molecules_list.append(val)
        elif input_item.endswith('.csv'):
            with open(input_item, 'r') as file:
                reader = csv.reader(file)
                for row in reader:
                    if len(row) >= 2:
                        molecules_list.append(row[1].strip().removeprefix('"').removesuffix('"'))
                        if names_list is None:
                            names_list = [row[0]]
                        else:
                            names_list.append(row[0])
        elif Path(input_item).suffix.lower() in ('.sdf', '.sd') and Path(input_item).exists():
            # Expand multi-molecule SDF file into individual structures
            try:
                from rdkit import Chem
                suppl = Chem.SDMolSupplier(input_item)
                valid_mols = [m for m in suppl if m is not None]
                if len(valid_mols) > 1:
                    for idx, m in enumerate(valid_mols):
                        molecules_list.append(m)
                        m_name = m.GetProp('_Name') if m.HasProp('_Name') and m.GetProp('_Name').strip() else f"{Path(input_item).stem}_{idx + 1}"
                        if names_list is None:
                            names_list = []
                        names_list.append(m_name)
                    continue
            except Exception:
                pass
            molecules_list.append(input_item)
        else:
            molecules_list.append(input_item)

    if names_list and len(names_list) != len(molecules_list):
        raise ValueError('Number of names must match number of molecules.')

    three_d = not args.show2D

    draw_multi(
        molecules_list,
        names=names_list,
        width=args.width,
        height=args.height,
        three_d=three_d,
        add_hydrogens=args.add_hydrogens,
        timeout=args.timeout,
        auto_resize=auto_resize,
        style=args.style,
        color=args.color,
        rep=args.rep,
        backend=args.backend,
        bond_spacing=args.bond_spacing
    )


def showcase_cli():
    parser = argparse.ArgumentParser(description='Showcase of 3D molecular rendering in the terminal using TerMol.')
    parser.add_argument('--width', type=int, default=80, help='Width of the output display in characters (default 80).')
    parser.add_argument('--height', type=int, default=40, help='Height of the output display in characters (default 40).')
    parser.add_argument('--timeout', type=int, default=5, help='3D window will cycle after this number of seconds (default 5).')

    args = parser.parse_args()
    showcase(timeout=args.timeout, width=args.width, height=args.height)


if __name__ == '__main__':
    termol_cli()

