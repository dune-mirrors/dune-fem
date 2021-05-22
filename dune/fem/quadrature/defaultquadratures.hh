#ifndef DUNE_FEM_DEFAULTQUADRATURES_HH
#define DUNE_FEM_DEFAULTQUADRATURES_HH

//#include <vector>
#include <cassert>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/grid.hh>

#include <dune/fem/quadrature/idprovider.hh>

#include "femquadratures.hh"

namespace Dune
{

  namespace Fem
  {

    // default defines for used quadratures
    template< typename FieldType, int dim >
    struct DefaultQuadratureTraits
    {
      typedef CubeQuadrature< FieldType, dim > CubeQuadratureType;
      typedef QuadratureImp< FieldType, dim > IntegrationPointListType;

      // dummy types for d > 3
      typedef CubeQuadratureType SimplexQuadratureType;

      typedef int QuadratureKeyType ;
    };



    // quadratures for points
    template< typename FieldType >
    struct DefaultQuadratureTraits< FieldType, 0 >
    {
      typedef CubeQuadrature< FieldType, 0 > PointQuadratureType;
      typedef QuadratureImp< FieldType, 0 > IntegrationPointListType;
      typedef int QuadratureKeyType ;
    };



    // quadratures for lines
    template< typename FieldType >
    struct DefaultQuadratureTraits< FieldType, 1 >
    {
      typedef CubeQuadrature< FieldType, 1 > LineQuadratureType;
      typedef QuadratureImp< FieldType, 1 > IntegrationPointListType;

      typedef int QuadratureKeyType ;
    };



    // quadratures for simplex and cubes
    template< typename FieldType >
    struct DefaultQuadratureTraits< FieldType, 2 >
    {
      typedef CubeQuadrature< FieldType, 2 >       CubeQuadratureType;
      typedef SimplexQuadrature< FieldType, 2 >    SimplexQuadratureType;
      typedef PolyhedronQuadrature< FieldType, 2 > PolyhedronQuadratureType;
      typedef QuadratureImp< FieldType, 2 > IntegrationPointListType;

      typedef int QuadratureKeyType ;
    };



    // quadratures for simplex, cubes, prisms, and pyramids
    template< typename FieldType >
    struct DefaultQuadratureTraits< FieldType , 3 >
    {
      typedef CubeQuadrature< FieldType, 3 >    CubeQuadratureType;
      typedef SimplexQuadrature< FieldType, 3 > SimplexQuadratureType;

      typedef PrismQuadrature< FieldType >   PrismQuadratureType;
      typedef PyramidQuadrature< FieldType > PyramidQuadratureType;

      typedef PolyhedronQuadrature< FieldType, 3 > PolyhedronQuadratureType;

      typedef QuadratureImp< FieldType, 3 > IntegrationPointListType;

      typedef int QuadratureKeyType ;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DEFAULTQUADRATURES_HH
