try:
    from dune.packagemetadata import metaData
except ImportError:
    from packagemetadata import metaData
from skbuild import setup
setup(**metaData('2.8.0.1')[1])
