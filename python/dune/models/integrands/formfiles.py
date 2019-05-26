from ufl import Form

import dune.ufl.formfiles

from dune.ufl import DirichletBC
from dune.ufl.formfiles import readUFLFile

from .ufl import compileUFL

class Model:
    def __init__(self, name, form):
        self.name = name
        if not isinstance(form, Form):
            raise Exception("Expecting 'form' to be a Form instance")
        self.form = form


def executeUFLCode(code, predefined=None):
    predefined = {} if predefined is None else predefined
    predefined["Model"] = Model
    return dune.ufl.formfiles.executeUFLCode(code, predefined)


def interpretUFLNamespace(namespace):
    models = namespace.get("models")
    if models is None:
        form = namespace.get("F")
        if isinstance(form, Form):
            models = [Model("Integrands", form)]
        elif form is None:
            models = []
        else:
            raise Exception("'F' is not a form")

    if not isinstance(models, (list, tuple)):
        raise Exception("Expecting 'models' to be a list of a tuple, not '%s'" % type(models))
    if not all(isinstance(m, Model) for m in models):
        raise Exception("Expecting 'models' to be a list of Model instances")

    return models


def loadUFLFile(filename, predefined=None):
    code = readUFLFile(filename)
    namespace = executeUFLCode(code, predefined=predefined)
    return interpretUFLNamespace(namespace)


def compileUFLFile(filename, predefined=None, tempVars=True):
    models = loadUFLFile(filename, predefined=predefined)
    return [compileUFL(model.form, tempVars=tempVars).code(model.name) for model in models]
