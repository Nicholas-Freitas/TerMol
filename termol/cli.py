import argparse
from .render import draw_multi, showcase
import csv

def termol_cli():
    parser = argparse.ArgumentParser(
        description='TerMol: A terminal molecular renderer for small molecules and proteins.'
    )
    parser.add_argument(
        'molecule_inputs', nargs='+',
        help='One or more SMILES strings, PDB accession codes (e.g. 1ubq), or molecule files (.sdf, .mol, .pdb). '
             'Separate with spaces. Can also be a text file (.txt) or CSV file (name, smiles/code).'
    )
    parser.add_argument('--names', nargs='+', help='Names for each molecule. Separate with spaces. Must be same length as inputs.')
    parser.add_argument('--show2D', action='store_true', help='Render static 2D view (default interactive 3D).')
    parser.add_argument('--width', type=int, help='Width of the output display in characters (default 80, or auto-resized in 3D).')
    parser.add_argument('--height', type=int, help='Height of the output display in characters (default 40, or auto-resized in 3D).')
    parser.add_argument('--add_hydrogens', action='store_true', help='Add hydrogens to small molecules via RDKit (default False).')
    parser.add_argument('--timeout', type=int, default=None, help='3D window will close after this number of seconds (default indefinitely).')
    parser.add_argument(
        '--style', choices=['ascii', 'braille'], default='ascii',
        help='Rendering style for proteins: "ascii" (depth ramp, default) or "braille" (high-resolution subpixel dots).'
    )
    parser.add_argument(
        '--color', choices=['depth', 'rainbow', 'ss', 'none'], default='depth',
        help='Color mode for proteins: "depth" (default luminance ramp), "rainbow" (N->C spectrum), "ss" (secondary structure), or "none".'
    )
    parser.add_argument(
        '--rep', '--representation', dest='rep', choices=['trace', 'cartoon', 'ribbon', 'tube'], default='trace',
        help='Representation for proteins: "trace" (1D spline line, default), "cartoon"/"ribbon" (secondary structure ribbons & arrowheads), or "tube" (3D uniform thickness).'
    )
    args = parser.parse_args()

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
        rep=args.rep
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

