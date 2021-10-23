import subprocess
commands='''
mkdir fem_tutorial
cd fem_tutorial

TMPNAME=`mktemp -d ./tmptutorial.XXXXXX`

# clone repo without history
git clone --quiet --depth 1 https://gitlab.dune-project.org/dune-fem/dune-fempy.git $TMPNAME
cd $TMPNAME

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
