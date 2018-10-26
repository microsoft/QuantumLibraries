: This script is called automatically by conda-build to test
: a newly built package.
: We use this script as opposed to specifying tests in meta.yaml directly
: so as to allow a chance to install pip requirements.
python -m pip install qinfer

python -c "import qsharp"
nosetests qsharp
