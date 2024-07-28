# TerMol
A simple molecular renderer for the terminal using RDKit.

https://github.com/user-attachments/assets/6fa67f32-f310-4851-aabf-9476feff8eca

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

## Usage:

Import the package and call the draw function as shown:
```python
import termol

smiles = "c1cc2c(cc1[N+](=O)[O-])[nH]nn2"
name   = "Nitrobenzotriazole"

termol.draw(smiles, name=name)
```

The molecule input can be a smiles string or an RDKit-compatible filepath. Optional arguments for the draw function include a molecule name to be displayed, the width and height (in characters) of your preferred canvas, whether to render in animated 3D or static 2D, and whether to have RDKit add hydrogens to the molecule:

```python
termol.draw(input_mol, name=None, width=80, height=40, three_d=True, add_hydrogens=False, timeout=None)
```
## License:
This software is provided under the MIT License.

## Contact:
[Nicholas Freitas](https://github.com/Nicholas-Freitas)

