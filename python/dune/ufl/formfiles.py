import ufl
import dune.ufl

from ufl.algorithms.formfiles import read_ufl_file as readUFLFile

def executeUFLCode(code, predefined=None):
    namespace = {}
    namespace.update(vars(ufl))
    namespace.update(vars(dune.ufl))
    if predefined is not None:
        namespace.update(predefined)
    exec(code, namespace)
    return namespace
