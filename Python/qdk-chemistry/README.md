# QDK chemistry

Prototype for chemistry library's Python application layer, contains tools for creating 2D molecular diagrams and calculating their 3D geometry using RDKit.

## How to install (development mode)

```
jupyter nbextension install jupyter_jsme
jupyter nbextension enable jupyter_jsme/extension
```

## How to use

To create a new JSME widget, run the following in a Jupyter notebook:

```python
from qdk_chemistry.widets import JsmeWidget
w = JsmeWidget()
w
```

This will display an interactive JSME widget in which you can design a molecule. When finished, click "save". We can then load the data into an RDKit Mol object

```python
mol = w.to_mol(add_hs=True)
```

where `add_hs` is used to optionally add Hydrogen atoms. Note that setting this value to `True` triggers a recalculation to find the optimal conformer for the molecule geometry.

To convert the 2D format to 3D and visualize it, run

```python
from qdk_chemistry.widgets import JsmolWidget
JsmolWidget.from_mol(mol)
```

This will display an interactive JSMol widget.
