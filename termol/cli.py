import argparse
from .render import draw_multi, showcase
import csv

def termol_cli():
    parser = argparse.ArgumentParser(description='A simple molecular renderer for the terminal using RDKit.')
    parser.add_argument('molecule_inputs', nargs='+', help='One or more SMILES strings, RDKit compatible molecule files (.sdf, .mol2). Separate with spaces.\n \
                        Can also be a text file with one input per line (file or smiles), or a CSV file in the format "name, smiles".')    
    parser.add_argument('--names', nargs='+', help='Names for each molecule. Separate with spaces. Must be same length as inputs.')
    parser.add_argument('--show2D', action='store_true', help='Render molecules in 2D (default 3D).')
    parser.add_argument('--width', type=int, default=80, help='Width of the output display, in characters (default 80).')
    parser.add_argument('--height', type=int, default=40, help='Height of the output display in characters (default 40).')
    parser.add_argument('--add_hydrogens', action='store_true', help='Add hydrogens to the molecule (default no hydrogens).')
    parser.add_argument('--timeout', type=int, default=None, help='3D window will close after this number of seconds (default no timeout).')

    args = parser.parse_args()

    print('Names', args.names)

    # Process inputs
    names_list = None
    if args.names:
        names_list = args.names
    
    molecules_list = []
    for input_item in args.molecule_inputs:
        if input_item.endswith('.txt'):
            with open(input_item, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    molecules_list.append(line.strip().removeprefix('"').removesuffix('"'))
        elif input_item.endswith('.csv'):
            with open(input_item, 'r') as file:
                reader = csv.reader(file)
                for row in reader:
                    molecules_list.append(row[1].strip().removeprefix('"').removesuffix('"'))
                    if names_list is None:
                        names_list = [row[0]]
                    else:
                        names_list.append(row[0])
        else:
            molecules_list.append(input_item)
    
    print(molecules_list)

    # Process names
    if names_list and len(names_list) != len(molecules_list):
        raise ValueError('Number of names must match number of molecules.')
    
    three_d = not args.show2D

    draw_multi(molecules_list, names=names_list, width=args.width, height=args.height, three_d=three_d, add_hydrogens=args.add_hydrogens, timeout=args.timeout)


def showcase_cli(timeout=None, width=80, height=40):
    # Call the showcase function with the arguments
    parser = argparse.ArgumentParser(description='A showcase of 3D molecular rendering in the terminal using TerMol.')
    parser.add_argument('--width', type=int, default=80, help='Width of the output display, in characters (default 80).')
    parser.add_argument('--height', type=int, default=40, help='Height of the output display in characters (default 40).')
    parser.add_argument('--timeout', type=int, default=5, help='3D window will close after this number of seconds (default 5 seconds).')

    args = parser.parse_args()

    showcase(timeout=args.timeout, width=args.width, height=args.height)
