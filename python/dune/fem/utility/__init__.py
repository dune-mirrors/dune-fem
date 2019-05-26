from __future__ import print_function

# just an idea

def evaluate(expression, coordinate):
    assert isinstance(expression, ufl.Expr)
    return value
def integrate(form,grid=None):
    assert check_form_arity(expression, arguments) == 0
    return value
def assemble(form,space=None):
    arity = check_form_arity(form, arguments)
    assert arity == 1 or arity == 2
    assert space or hasattr(form.testFunction,"space")
    if arity == 1:
        return functional
    else:
        return matrix
def solve(equation,solutionGF):
    arity = check_form_arity(equation, arguments)
    assert arity == 1 or arity == 2
    if arity == 2:
        return info
    else:
        assert solutionGF
        # assert that solutionGF is a coefficient in equation and replace
        # it by trialfunction
        return info

# or all in one?
def evaluate(expression, grid=None, space=None, target=None, coordinate=None)
    if isinstance(expression, ufl.Expr):
        assert coordinate
        return expression(coordinate)
    elif isinstance(expression, ufl.Form):
        if check_form_arity(expression, arguments) == 0:
            assert grid or space
            if not grid: grid = space.grid
            pass # integrate function
        if check_form_arity(expression, arguments) == 1:
            assert space
            pass # return a df
        elif check_form_arity(expression, arguments) == 2:
            assert space
            pass # return a matrix
    elif isinstance(expression, ufl.Equation):
        assert target
        pass # return solver info
