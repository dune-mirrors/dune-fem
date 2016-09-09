from .grid             import leafGrid as leafGrid
from .gridpart         import create as gridpart
from .space            import create as space
from .discretefunction import create as discretefunction
from .scheme           import create as scheme

from .grid             import myGenerator as gridGenerator
from .gridpart         import myGenerator as gridpartGenerator
from .space            import myGenerator as spaceGenerator
from .discretefunction import myGenerator as discretefunctionGenerator
from .scheme           import myGenerator as schemeGenerator

def ellipticModel(*args,**kwargs):
    from ..models.elliptic import create
    return create(*args,**kwargs)

def rebuild():
    gridGenerator.force=gridpartGenerator.force=spaceGenerator.force=\
            discretefunctionGenerator.force=schemeGenerator.force = True
