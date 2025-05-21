import io
from dune.generator import algorithm
def interpolate(fromGF, toDF):
    try:
        toFilter = fromGF.gridView.cppTypeName in toDF.gridView.hostCppTypeName
    except AttributeError:
        try:
            if not toDF.gridView.cppTypeName in fromGF.gridView.hostCppTypeName:
                raise ValueError("provided functions to not have identical underlying grids or both are not over filtered grid views.")
            toFilter = False
        except AttributeError:
            raise ValueError("provided functions to not have identical underlying grids or both are not over filtered grid views.")
    if not hasattr(toDF,"interpolate"):
        raise ValueError("second argument has to be a discrete function")

    if toFilter: # the df is over the filter:
        code = """
#include <dune/fem/space/common/localinterpolation.hh>
#include <dune/fem/common/localcontribution.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/common/bindguard.hh>
template <class GF, class DF>
void filteredInterpolate( const GF &gf, DF& df )
{
  Dune::Fem::LocalInterpolation< typename DF::DiscreteFunctionSpaceType > interpolation( df.space() );
  Dune::Fem::ConstLocalFunction< GF > gfLocal( gf );
  Dune::Fem::SetSelectedLocalContribution< DF > dfLocal( df );
  for( const auto &entity : df.space() ) // iterate over filter coming from df
  {
    auto gfGuard = Dune::Fem::bindGuard( gfLocal, entity ); // note: entities on filtergp are not wrapped
    auto dfGuard = Dune::Fem::bindGuard( dfLocal, entity );
    auto iGuard  = Dune::Fem::bindGuard( interpolation, entity );
    interpolation(gfLocal,dfLocal);
  }
}
        """
    else:
        code = """
#include <dune/fem/space/common/localinterpolation.hh>
#include <dune/fem/common/localcontribution.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/common/bindguard.hh>
template <class GF, class DF>
void filteredInterpolate( const GF &gf, DF& df )
{
  Dune::Fem::LocalInterpolation< typename DF::DiscreteFunctionSpaceType > interpolation( df.space() );
  Dune::Fem::ConstLocalFunction< GF > gfLocal( gf );
  Dune::Fem::SetLocalContribution< DF > dfLocal( df );
  for( const auto &entity : gf.space() ) // iterate over filter coming from gf
  {
    auto gfGuard = Dune::Fem::bindGuard( gfLocal, entity );
    auto dfGuard = Dune::Fem::bindGuard( dfLocal, entity ); // note: entities on filtergp are not wrapped
    auto iGuard  = Dune::Fem::bindGuard( interpolation, entity );
    interpolation(gfLocal,dfLocal);
  }
}
        """
    algorithm.run("filteredInterpolate", io.StringIO(code), fromGF, toDF )
