from argparse import ArgumentParser
import subprocess, os

# note: the github actions will insert a
# __pypi_dochash__ on top of this file
try:
    _docHash = __pypi_dochash__
except NameError:
    # fallback for default
    _docHash = "master"

parser = ArgumentParser(description='Execute DUNE-FEM commands', prog='dune.fem')
subparsers = parser.add_subparsers(dest='command')

# Reader
parserConfigure = subparsers.add_parser('reader',
                  help='Environment variable export needed for paraview to find reader')
parserConfigure = subparsers.add_parser('readerpath',
                  help='Path to paraview file reader')
parserConfigure = subparsers.add_parser('tutorial',
                  help='Path to paraview file reader')

args = parser.parse_args()

if args.command == 'reader':
    file_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "reader")
    print(f"Paste the following or run (in your .bashrc for example):",end="")
    print(f"export PV_PLUGIN_PATH=`python -m dune.fem readerpath`")
    print(f"export PV_PLUGIN_PATH={file_path}")
elif args.command == 'readerpath':
    file_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "reader")
    print(file_path)
else: # download tutorial
    print("Downloading dune-fem tutorial examples, this will take a short time.")
    commands=f'''
    mkdir fem_tutorial
    cd fem_tutorial

    TMPNAME=`mktemp -d ./tmptutorial.XXXXXX`

    # clone repo without history
    git clone --quiet https://gitlab.dune-project.org/dune-fem/dune-fempy.git $TMPNAME
    cd $TMPNAME
    git checkout --quiet {_docHash}

    cp doc/*.py doc/*.ipynb doc/*.hh doc/*.dgf doc/*.msh ..
    cd ../

    rm -rf $TMPNAME
    rm pandoc-formatting.py gitlab-formatting.py interpolation.py svg2pdf.py
    '''
    subprocess.check_output(commands, shell=True)

    print("############################################################")
    print("## The tutorial is now located in the 'fem_tutorial' folder.")
    try:
        import matplotlib
    except ImportError:
        print("## Note: some of the examples require the installation of 'matplotlib'.")
    try:
        import scipy
    except ImportError:
        print("## Note: some of the examples require the installation of 'scipy'.")
    print("############################################################")
