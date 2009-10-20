#ifndef DUNE_FEM_DEFAULTQUADRATURES_HH
#define DUNE_FEM_DEFAULTQUADRATURES_HH

//#include <vector>
#include <cassert>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/grid.hh>

#include <dune/fem/storage/array.hh>
#include <dune/fem/quadrature/idprovider.hh>

// don't use quadratures from dune-grid 
//#define USE_DUNE_QUADRATURES

// include quadrature points 
#ifdef USE_DUNE_QUADRATURES
#warning "Don't use DUNE Quadratures!!!" 
#include "dunequadratures.hh"
#else
#include "femquadratures.hh"
#endif

namespace Dune
{

  // default defines for used quadratures 
  template< typename FieldType, int dim >
  struct DefaultQuadratureTraits
  {
#ifdef USE_DUNE_QUADRATURES
    typedef QuadratureRulesFactory< FieldType, dim > CubeQuadratureType;
#else 
    typedef CubeQuadrature< FieldType, dim > CubeQuadratureType; 
#endif
    typedef QuadratureImp< FieldType, dim > IntegrationPointListType;
  }; 



  // quadratures for points 
  template< typename FieldType >
  struct DefaultQuadratureTraits< FieldType, 0 >  
  {
#ifdef USE_DUNE_QUADRATURES
    typedef QuadratureRulesFactory< FieldType, 0 > PointQuadratureType;
#else 
    typedef CubeQuadrature< FieldType, 0 > PointQuadratureType;     
#endif
    typedef QuadratureImp< FieldType, 0 > IntegrationPointListType;
  };
 


  // quadratures for lines 
  template< typename FieldType >
  struct DefaultQuadratureTraits< FieldType, 1 >  
  {
#ifdef USE_DUNE_QUADRATURES
    typedef QuadratureRulesFactory< FieldType, 1 > LineQuadratureType;
#else 
    typedef CubeQuadrature< FieldType, 1 > LineQuadratureType;     
#endif
    typedef QuadratureImp< FieldType, 1 > IntegrationPointListType;
  };
 


  // quadratures for simplex and cubes 
  template< typename FieldType >
  struct DefaultQuadratureTraits< FieldType, 2 >  
  {
#ifdef USE_DUNE_QUADRATURES
    typedef QuadratureRulesFactory< FieldType, 2 > SimplexQuadratureType;
    typedef QuadratureRulesFactory< FieldType, 2 > CubeQuadratureType;
#else
    typedef CubeQuadrature< FieldType, 2 > CubeQuadratureType; 
    typedef SimplexQuadrature< FieldType, 2 > SimplexQuadratureType;     
#endif
    typedef QuadratureImp< FieldType, 2 > IntegrationPointListType;
  };


  
  // quadratures for simplex, cubes, prisms, and pyramids
  template< typename FieldType >
  struct DefaultQuadratureTraits< FieldType , 3 >  
  {
#ifdef USE_DUNE_QUADRATURES
    typedef QuadratureRulesFactory< FieldType, 3 > SimplexQuadratureType;
    typedef QuadratureRulesFactory< FieldType, 3 > CubeQuadratureType;

    typedef QuadratureRulesFactory< FieldType, 3 > PrismQuadratureType;
    typedef QuadratureRulesFactory< FieldType, 3 > PyramidQuadratureType;
#else 
    typedef CubeQuadrature< FieldType, 3 > CubeQuadratureType;     
    typedef SimplexQuadrature< FieldType, 3 > SimplexQuadratureType;     

    typedef PrismQuadrature< FieldType > PrismQuadratureType;
    typedef PyramidQuadrature< FieldType > PyramidQuadratureType;
#endif

    typedef QuadratureImp< FieldType, 3 > IntegrationPointListType;
  };

} // end namespace Dune

#undef USE_DUNE_QUADRATURES
#endif
