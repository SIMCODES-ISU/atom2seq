# atom2seq

Detects primary structure of proteins given the Cartesian coordinates of the 
atoms.

# Installation

Installation instructions for users are coming soon!

## Development
All commands assume a Linux-like command line.

1. Obtain the repo from github
```
git clone https://github.com/SIMCODES-ISU/atom2seq.git
```

2. Create a Python virtual environment
```
cd atom2seq
python3 -m venv .venv
```

3. Activate the virtual environment
```
source .venv/bin/activate
```

4. Install the development dependencies
```
pip install ".[dev]"
```

5. Editable install the repo
```
pip install -e .
```

You can then test the repo by running
```
pytest
```

# Acknowledgements

This material is based upon work supported by the National Science Foundation under Grant No. 2348724.