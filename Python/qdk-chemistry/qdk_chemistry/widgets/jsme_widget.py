from collections import namedtuple
from typing import TYPE_CHECKING
from notebook.nbextensions import check_nbextension
from varname import varname, VarnameRetrievingError

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
    jsmeApp_ = new JSApplet.JSME("JSApp_{uid}", "{width}px", "{height}px", {{
        "options" : "{options}",
        "{data_type}": "{data}"
    }});
}}

function saveAs(var_name) {{
    var molblock = jsmeApp_.molFile();
    var smiles = jsmeApp_.smiles();
    var jme = jsmeApp_.jmeFile();
    var command = var_name + '.set_value(jme="' + jme + '", smiles="' + smiles + '", molblock="""' + molblock + '""")';
    console.log(command);
    IPython.notebook.kernel.execute(command);
}}
</script>
<div id="JSApp_{uid}"></div>
<button onclick="saveAs('{var_name}')">save</button>
'''

class JsmeWidget:
    """Jupyter widget for displaying JSME molecular editor.
    Allows to either build molecule from scratch and save the JME, SMILES and
    MolBlock values, or set either of the aforementioned values manually.
    """
    n = 0

    def __init__(self, width=400, height=350, jme="", smiles="", molblock="", options="query,hydrogens"):
        """Create JsmeWidget instance

        Args:
            width (int, optional): Widget width in pixels. Defaults to 400.
            height (int, optional): Widget height in pixels. Defaults to 350.
            jme (str, optional): JME value string. Defaults to "".
            smiles (str, optional): SMILES value string. Defaults to "".
            molblock (str, optional): MolBlock value string. Defaults to "".
            options (str, optional): Options to pass to widget. Defaults to "query,hydrogens".
        """
        try:
            self.name = varname()
        except VarnameRetrievingError:
            self.name = "_"

        self.width = width
        self.height = height
        self.value = JsmeValue(jme=jme, molblock=molblock, smiles=smiles)
        self.options = options

    def set_value(self, jme: str, smiles: str, molblock: str):
        """Set JSME value (JME, SMILES and MolBlock)

        Args:
            jme (str, optional): JME value string. Defaults to "".
            smiles (str, optional): SMILES value string. Defaults to "".
            molblock (str, optional): MolBlock value string. Defaults to "".
        """
        self.value = JsmeValue(jme=jme, molblock=molblock, smiles=smiles)

    @property
    def html_str(self) -> str:
        """Returns an HTML string that contains the widget.

        Returns:
            str: HTML string for displaying widget
        """
        JsmeWidget.n += 1
        data_type = "jme" if self.value.jme else "smiles"
        data = self.value.jme if self.value.jme else self.value.smiles

        return _HTML_STR_FORMAT.format(
            base_url=BASE_URL,
            js_source=JS_SOURCE,
            var_name=self.name,
            uid=str(JsmeWidget.n),
            width=self.width,
            height=self.height,
            options=self.options,
            data_type=data_type,
            data=data
        )

    def _repr_html_(self):
        return self.html_str

    def to_mol(self, add_hs: bool=False, num_confs: int=10) -> "Mol":
        """Convert widget value to RDKit molecule.
        If Hydrogen atoms are added, calculate the optimal conformer to get H coordinates.

        Args:
            add_hs (bool, optional): Add Hydrogen atoms
            num_confs (int, optional): Number of conformers to generate

        Returns:
            Mol: RDKit Mol object
        """
        if self.value.smiles is not "":
            mol = Chem.MolFromSmiles(self.value.smiles)
        elif self.value.molblock is not "":
            mol = Chem.MolFromMolBlock(self.value.molblock)
        else:
            raise ValueError("Cannot create Mol object: JSME value is empty")

        if add_hs is True:
            mol = Chem.AddHs(mol)
            # Calculate conformers after adding hydrogens
            Chem.EmbedMultipleConfs(mol, numConfs=num_confs)

        return mol
