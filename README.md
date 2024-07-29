# TerMol
A simple molecular renderer for the terminal using RDKit.

![molecules](https://github.com/user-attachments/assets/a80a2287-8aeb-4055-a71c-88a9c4474998)

## Table of Contents
1. [Overview](#overview)
3. [Installation](#installation)
4. [Usage](#usage)
8. [License](#license)
9. [Contact](#contact)

## Overview
This python package allows for the 2D or 3D rendering of molecules in the terminal, using RDKit and Curses.

## Installation
Install TerMol with:
```bash
pip install termol
```

On Windows only, the Curses python package must be installed manually.

## Command-Line Usage:
Once installed, you can use the `termol` command to render molecules in 2D or 3D.

It takes as input one or more SMILES strings or RDKit compatible molecule files (.sdf, .mol2, separated with spaces. Inputs can also be a text file with one input per line (file or smiles), or a CSV file in the format "name, smiles".

Example usage: `termol ESR.sdf "O=C4N(Cc1ccccc1)CC(C(=O)NC(Cc2ccccc2)C(O)CNCc3cccc(N(C)C)c3)(N4)Cc5ccccc5" "O=C(O)CCCCC(=O)O"`

Additional arguments include:
- `--names`: The molecule names to be displayed (separated by spaces).
- `--width, --height`: The size of the 'screen' (in number of characters) Default 80x40.
- `--show2D`: Display in 2D instead of 3D.
- `--add_hydrogens`: Have RDKit attempt to add hydrogens automatically. Default False.
- `--timeout`: In the 3D viewer, this will automatically close after this number of seconds. Default None, which allows the viewer to stay open indefinitely.

### Controlling the viewer
In the 3D viewer, the molecule will initially rotate about the Y-axis. To change the direction of the rotation, use the arrow keys, or WASD+QE. To pause rotation, hit the spacebar.

To make the molecule smaller or larger, you can use the R and F keys, respectively.

To exit the 3D viewer, hit any other key.

## Python Usage:

### 3D Rendering:
Import the package and call the draw function as shown:
```python
import termol

smiles = "c1cc2c(cc1[N+](=O)[O-])[nH]nn2"
name   = "Nitrobenzotriazole"

termol.draw(smiles, name=name)
```

The molecule input can be a smiles string or an RDKit-compatible filepath (.mol2, .sdf, etc). 

### 2D Rendering:
Sometimes, a simple 2D graphic is sufficient. To render in 2D, use the flag `three_d=False`:
```python
termol.draw(smiles, name=name, three_d=False)
```
![image](https://github.com/user-attachments/assets/63694895-b34c-4166-8815-0da9afc6bc62)

### Showcase:
To display a showcase of termol's capabilities, you may run:
```python
import termol

termol.showcase()
```

Want a fun screensaver? Use the `timeout=60` argument to cycle through a random molecule every 60 seconds.

### Other Options
The draw function only requires the molecule SMILES/file as input. Other options include:
- `name`: A molecule name to be displayed
- `width, height`: The size of the 'screen' (in number of characters) Default 80x40.
- `three_d`: defaults to True to display in 3D. Set to False to print a 2D view.
- `add_hydrogens`: Have RDKit attempt to add hydrogens automatically. Default False.
- `timeout`: In the 3D viewer, this will automatically close after this number of seconds. Default None, which allows the viewer to stay open indefinitely.

```python
termol.draw(input_mol, name=None, width=80, height=40, three_d=True, add_hydrogens=False, timeout=None)
```

### Drawing Many Molecules
The `termol.draw_multi()` function takes as input a list of SMILES or molecule files as input, and displays them in succession (In 2D or 3D). If the optional argument `names` is provided, it must be the same length as the number of inputs. All other arguments are the same.

## License:
This software is provided under the MIT License.

## Contact:
[Nicholas Freitas](https://github.com/Nicholas-Freitas)

