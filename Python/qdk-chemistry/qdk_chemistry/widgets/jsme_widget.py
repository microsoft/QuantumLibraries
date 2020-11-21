# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.

"""Module for Jupyter Widget that displays the JSME editor.
"""
import uuid
import ipywidgets

from collections import namedtuple
from typing import TYPE_CHECKING
from notebook.nbextensions import check_nbextension
from varname import varname, VarnameRetrievingError
from IPython.core.display import Javascript, display, HTML

from rdkit.Chem import AllChem as Chem

from .jsmol_widget import JsmolWidget

if TYPE_CHECKING:
    from rdkit.Chem import Mol, Conformer

# Named tuple for storing serialized JSME widget value in several popular formats
JsmeValue = namedtuple("JsmeValue", ["jme", "molblock", "smiles"])

# Rel file path for Javascript source
JS_SOURCE = "jsme/jsme.nocache.js"

# Base URL for widget source javascript
BASE_URL = "/nbextensions/jupyter_jsme"

if check_nbextension("jupyter_jsme") is False:
    # Fallback option in case user didn't install Jupyter extension
    BASE_URL = "https://peter-ertl.com/jsme/JSME_2020-06-11"

_HTML_STR_FORMAT = '''
<script type="text/javascript" src="{base_url}/{js_source}"></script>
<script type="text/javascript">
function jsmeOnLoad() {{
    jsmeApp_{uid} = new JSApplet.JSME("JSApp_{uid}", "{width}px", "{height}px", {{
        "options" : "{options}",
        "{data_type}": "{data}"
    }});
}}

function saveAs_{uid}(var_name) {{
    var molblock = jsmeApp_{uid}.molFile();
    var smiles = jsmeApp_{uid}.smiles();
    var jme = jsmeApp_{uid}.jmeFile();
    var command = var_name + '.set_value(jme="' + jme + '", smiles="' + smiles + '", molblock="""' + molblock + '""")';
    console.log(command);
    IPython.notebook.kernel.execute(command);
}}

</script>
<div id="JSApp_{uid}"></div>
'''

class JsmeWidgetButton(ipywidgets.Button):
    def __init__(self, uid: str, *args, **kwargs):
        """Create "save" button for JsmeWidget

        :param uid: Unique ID of widget
        :type uid: str
        """
        self.uid = uid
        return super(JsmeWidgetButton, self).__init__(description="Save", *args, **kwargs)


class JsmeWidget:
    """Jupyter widget for displaying JSME molecular editor.
    Allows to either build molecule from scratch and save the JME, SMILES and
    MolBlock values, or set either of the aforementioned values manually.
    """
    n = 0

    def __init__(self, width=400, height=350, jme="", smiles="", molblock="", options="query,hydrogens", parent_varname: str = ""):
        """Create JsmeWidget instance

        :param width: Widget width in pixels, defaults to 400
        :type width: int, optional
        :param height: Widget height in pixels, defaults to 350
        :type height: int, optional
        :param jme: JME value string, defaults to ""
        :type jme: str, optional
        :param smiles: SMILES value string, defaults to ""
        :type smiles: str, optional
        :param molblock: MolBlock value string, defaults to ""
        :type molblock: str, optional
        :param options: Options to pass to widget, defaults to "query,hydrogens"
        :type options: str, optional
        """
        try:
            self.name = f"{parent_varname}.{varname()}"
        except VarnameRetrievingError:
            self.name = "_"

        self.width = width
        self.height = height
        self.value = JsmeValue(jme=jme, molblock=molblock, smiles=smiles)
        self.options = options
        self._uids = []
        self._updated = False

    def set_value(self, jme: str, smiles: str, molblock: str):
        """Set JSME value (JME, SMILES and MolBlock)

        :param jme: JME value string
        :type jme: str
        :param smiles: SMILES value string
        :type smiles: str
        :param molblock: MolBlock value string
        :type molblock: str
        """
        self._updated = True
        self.value = JsmeValue(jme=jme, molblock=molblock, smiles=smiles)

    @property
    def was_updated(self):
        return self._updated

    def reset_updated(self):
        self._updated = False

    def _gen_uid(self):
        """Generate unique identifier for javascript applet"""
        uid = str(uuid.uuid1()).replace("-", "")
        # Keep track of all UIDs
        self._uids.append(uid)
        return uid

    def html_str(self, uid: str) -> str:
        """Returns an HTML string that contains the widget.

        :param uid: Unique identifier of widget
        :type uid: str
        :return: HTML string for displaying widget
        :rtype: str
        """
        JsmeWidget.n += 1
        data_type = "jme" if self.value.jme else "smiles"
        data = self.value.jme if self.value.jme else self.value.smiles

        return _HTML_STR_FORMAT.format(
            base_url=BASE_URL,
            js_source=JS_SOURCE,
            var_name=self.name,
            uid=uid,
            width=self.width,
            height=self.height,
            options=self.options,
            data_type=data_type,
            data=data
        )

    def _ipython_display_(self):
        """Display the widget"""
        uid = self._gen_uid()
        jsme = HTML(self.html_str(uid=uid))
        button = JsmeWidgetButton(uid=uid)
        
        def _save(button):
            self.save(button.uid)
        
        button.on_click(_save)
        display(jsme, button)

    def save(self, uid):
        """Save the value of the molecular diagram drawn in the widget to the Python namespace
        """
        display(Javascript(f"saveAs_{uid}('{self.name}')"))

    def to_mol(self, add_hs: bool=False, num_confs: int=10) -> "Mol":
        """Convert widget value to RDKit molecule.
        If Hydrogen atoms are added, calculate the optimal conformer to get H coordinates.

        :param add_hs: Add Hydrogen atoms
        :type add_hs: bool
        :param num_confs: Number of conformers to generate
        :type num_confs: int
        :return: RDKit molecule object
        :rtype: Mol
        """
        if self.value.smiles != "":
            mol = Chem.MolFromSmiles(self.value.smiles)
        elif self.value.molblock != "":
            mol = Chem.MolFromMolBlock(self.value.molblock)
        else:
            raise ValueError("Cannot create Mol object: JSME value is empty")

        if add_hs is True and mol is not None:
            mol = Chem.AddHs(mol)
            # Calculate conformers after adding hydrogens
            Chem.EmbedMultipleConfs(mol, numConfs=num_confs)

        if mol is None:
            raise ValueError(f"Cannot convert JSME widget value to Mol object:\
                \n{self.value}\
                \nPlease try clicking 'Save' button on widget or remove H atoms from SMILES string and try again.")

        return mol
