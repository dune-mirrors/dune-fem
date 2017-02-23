from ufl import Form

import dune.ufl.formfiles

from dune.ufl import DirichletBC
from dune.ufl.formfiles import readUFLFile

from .ufl import compileUFL

class Model:
    def __init__(self, name, form, *args):
        self.name = name
        if not isinstance(form, Form):
            raise Exception("Expecting 'form' to be a Form instance")
        self.form = form
        if not all(isinstance(arg, DirichletBC) for arg in args):
            raise Exception("Expecting each constraint to be a DirichletBC instances")
        self.constraints = args


def executeUFLCode(code):
    return dune.ufl.formfiles.executeUFLCode(code, {"Model": Model})


def interpretUFLNamespace(namespace):
    models = namespace.get("models")
    if models is None:
        form = namespace.get("F")
        if isinstance(form, Form):
            constraints = namespace.get("constraints", [])
            if not isinstance(constraints, (list, tuple)):
                raise Exception("Expecting 'constraints' to be a list of a tuple, not '%s'" % type(constraints))
            if not all(isinstance(c, DirichletBC) for c in constraints):
                raise Exception("Expecting 'constraints' to be a list of DirichletBC instances")
            models = [Model("Model", form, *list(constraints))]
        elif form is None:
            model = []
        else:
            raise Exception("'F' is not a form")

    if not isinstance(models, (list, tuple)):
        raise Exception("Expecting 'models' to be a list of a tuple, not '%s'" % type(models))
    if not all(isinstance(m, Model) for m in models):
        raise Exception("Expecting 'models' to be a list of Model instances")

    return models


def loadUFLFile(filename):
    code = readUFLFile(filename)
    namespace = executeUFLCode(code)
    return interpretUFLNamespace(namespace)


def compileUFLFile(filename, tempVars=True):
    models = loadUFLFile(filename)
    return [compileUFL(model.form, *model.constraints, tempVars=tempVars).code(model.name) for model in models]
