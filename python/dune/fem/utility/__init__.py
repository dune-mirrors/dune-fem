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
def evaluate(expression, grid=None, space=None, target=None, coordinate=None):
    if isinstance(expression, ufl.Expr):
        assert coordinate
        return expression(coordinate)
    elif isinstance(expression, ufl.Form):
        if check_form_arity(expression, arguments) == 0:
            assert grid or space
            if not grid: grid = space.gridView
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

class Sampler:
    def __init__(self, gridFunction ):
        self.gridFunction = gridFunction
        self.lineSampler = None
        self.pointSampler = None

    def lineSample(self,x0,x1,N):
        from dune.generator import algorithm, path
        from dune.common import FieldVector
        import numpy
        x0, x1 = FieldVector(x0), FieldVector(x1)
        if self.lineSampler is None:
            import io
            _lineSampleCode = \
            """
            #ifndef FEMPY_UTILITY_HH
            #define FEMPY_UTILITY_HH
            #include <vector>
            #include <utility>
            #include <dune/fem/misc/linesegmentsampler.hh>
            #include <dune/fem/gridpart/common/entitysearch.hh>
            #include <dune/fem/function/localfunction/const.hh>

            template <class GF, class DT>
            std::pair<std::vector<DT>, std::vector<typename GF::RangeType>>
            sample(const GF &gf, DT &start, DT &end, int n)
            {
              Dune::Fem::LineSegmentSampler<typename GF::GridPartType> sampler(gf.gridPart(),start,end);
              std::vector<DT> coords(n);
              std::vector<typename GF::RangeType> values(n);
              sampler(gf,values);
              sampler.samplePoints(coords);
              return std::make_pair(coords,values);
            }
            #endif

            """
            self.lineSampler = algorithm.load('sample', io.StringIO(_lineSampleCode), self.gridFunction, x0, x1, N)

        p,v = self.lineSampler( self.gridFunction, x0, x1, N )
        if self.gridFunction.scalar:
            x,y = numpy.zeros(len(p)), numpy.zeros(len(p))
        else:
            x,y = numpy.zeros(len(p)), numpy.zeros(( len(p),self.gridFunction.dimRange) )
        length = (x1-x0).two_norm
        for i in range(len(x)):
            x[i] = (p[i]-x0).two_norm / length
            if self.gridFunction.scalar:
                y[i] = v[i][0]
            else:
                y[i] = v[i]
        return x,y

    def pointSample(self, x0):
        from dune.generator import algorithm, path
        from dune.common import FieldVector
        import numpy
        x0 = FieldVector(x0)
        if self.pointSampler is None:
            import io
            _pointSampleCode = \
            """
            #ifndef FEMPY_UTILITY_HH
            #define FEMPY_UTILITY_HH
            #include <dune/fem/misc/linesegmentsampler.hh>
            #include <dune/fem/gridpart/common/entitysearch.hh>
            #include <dune/fem/function/localfunction/const.hh>

            template <class GF, class DT>
            typename GF::RangeType sample(const GF &gf, DT &point)
            {
              typedef typename GF::DiscreteFunctionSpaceType::GridPartType GridPartType;
              Dune::Fem::EntitySearch<GridPartType> search(gf.space().gridPart());
              const auto &entity = search(point);
              const auto localPoint = entity.geometry().local(point);
              return constLocalFunction(gf,entity).evaluate(localPoint);
            }
            #endif // FEMPY_UTILITY_HH

            """
            self.pointSampler = algorithm.load('sample', io.StringIO(_pointSampleCode), self.gridFunction, x0)

        v = self.pointSampler(self.gridFunction, x0 )
        if self.gridFunction.scalar:
            return v[0]
        else:
            return v

def lineSample(gridFunction,x0,x1,N):
    return Sampler(gridFunction).lineSample(x0, x1, N)

def pointSample(gridFunction,x0):
    return Sampler(gridFunction).pointSample( x0 )
