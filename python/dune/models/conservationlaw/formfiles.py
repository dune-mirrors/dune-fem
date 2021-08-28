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

def to_alphabetic(i,base=26):
    if base < 0 or 62 < base:
        raise ValueError("Invalid base")

    if i < 0:
        return '-'+to_alphabetic(-i-1)

    quot = int(i/base)
    rem = i%base
    if rem < 26:
        letter = chr( ord("A") + rem)
    elif rem < 36:
        letter = str( rem-26)
    else:
        letter = chr( ord("a") + rem - 36)
    if quot == 0:
        return letter
    else:
        return to_alphabetic(quot-1,base) + lett

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
            models = []
        else:
            raise Exception("'F' is not a form")
    else:
        if isinstance(models,dict):
            models = [
                Model(name, models[name]) if isinstance(models[name],Form) else
                Model(name, models[name][0], *models[name][1:])
                for name in models ]
        elif all(isinstance(m, Form) for m in models):
            constraints = namespace.get("constraints", [])
            if not isinstance(constraints, (list, tuple)):
                raise Exception("Expecting 'constraints' to be a list of a tuple, not '%s'" % type(constraints))
            if not all(isinstance(c, DirichletBC) for c in constraints):
                raise Exception("Expecting 'constraints' to be a list of DirichletBC instances")
            models = [Model("Model"+to_alphabetic(i), m, *list(constraints)) for i,m in enumerate(models)]

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
    return [compileUFL(model.form, *model.constraints, tempVars=tempVars)[0].code(model.name) for model in models]

def loadModels(view,filename, *args, **kwargs):
    from dune.fem.model import conservationlaw
    models = loadUFLFile(filename)
    return [ conservationlaw(view, model.form,
             *model.constraints, *args, **kwargs) for model in models]
