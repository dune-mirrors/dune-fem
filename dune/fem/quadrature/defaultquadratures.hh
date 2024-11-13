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
    struct DefaultQuadratureImplementationTraits
    {
      //! type of integration point list / quadrature implementation
      typedef QuadratureImp< FieldType, dim > QuadratureImplementationType;

      // deprecated typedef
      typedef QuadratureImplementationType  IntegrationPointListType;
    };

    // default defines for used quadratures
    template< typename FieldType, int dim >
    struct DefaultQuadratureTraits
      : public DefaultQuadratureImplementationTraits< FieldType, dim >
    {
      typedef CubeQuadrature< FieldType, dim > CubeQuadratureType;

      // dummy types for d > 3
      typedef CubeQuadratureType SimplexQuadratureType;

      typedef int QuadratureKeyType ;
    };



    // quadratures for points
    template< typename FieldType >
    struct DefaultQuadratureTraits< FieldType, 0 >
      : public DefaultQuadratureImplementationTraits< FieldType, 0 >
    {
      typedef CubeQuadrature< FieldType, 0 > PointQuadratureType;
      typedef int QuadratureKeyType ;
    };



    // quadratures for lines
    template< typename FieldType >
    struct DefaultQuadratureTraits< FieldType, 1 >
      : public DefaultQuadratureImplementationTraits< FieldType, 1 >
    {
      typedef CubeQuadrature< FieldType, 1 > LineQuadratureType;
      typedef QuadratureImp< FieldType, 1 > IntegrationPointListType;

      typedef int QuadratureKeyType ;
    };



    // quadratures for simplex and cubes
    template< typename FieldType >
    struct DefaultQuadratureTraits< FieldType, 2 >
      : public DefaultQuadratureImplementationTraits< FieldType, 2 >
    {
      typedef CubeQuadrature< FieldType, 2 >       CubeQuadratureType;
      typedef SimplexQuadrature< FieldType, 2 >    SimplexQuadratureType;
      typedef PolyhedronQuadrature< FieldType, 2 > PolyhedronQuadratureType;

      typedef int QuadratureKeyType ;
    };



    // quadratures for simplex, cubes, prisms, and pyramids
    template< typename FieldType >
    struct DefaultQuadratureTraits< FieldType , 3 >
      : public DefaultQuadratureImplementationTraits< FieldType, 3 >
    {
      typedef CubeQuadrature< FieldType, 3 >    CubeQuadratureType;
      typedef SimplexQuadrature< FieldType, 3 > SimplexQuadratureType;

      typedef PrismQuadrature< FieldType >   PrismQuadratureType;
      typedef PyramidQuadrature< FieldType > PyramidQuadratureType;

      typedef PolyhedronQuadrature< FieldType, 3 > PolyhedronQuadratureType;

      typedef int QuadratureKeyType ;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DEFAULTQUADRATURES_HH
