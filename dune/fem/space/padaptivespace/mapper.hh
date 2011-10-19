#ifndef DUNE_PADAPTIVELAGRANGESPACE_MAPPER_HH
#define DUNE_PADAPTIVELAGRANGESPACE_MAPPER_HH

//- Dune includes 
#include <dune/common/geometrytype.hh>
#include <dune/common/exceptions.hh>

//- Dune-Fem includes 
#include <dune/fem/misc/capabilities.hh>
#include <dune/fem/misc/codimmap.hh>
#include <dune/fem/misc/metaprogramming.hh>
#include <dune/fem/space/common/dofmanager.hh>

#include <dune/fem/space/mapper/dofmapper.hh>
#include <dune/fem/space/mapper/codimensionmapper.hh>

//- local includes 
#include <dune/fem/space/lagrangespace/lagrangepoints.hh>

#include <dune/grid/utility/persistentcontainer.hh>

// include generic adaptive dof mapper 
#include <dune/fem/space/mapper/genericadaptivedofmapper.hh>

namespace Dune
{
  
  template< class GridPart, int polOrder >
  class PAdaptiveLagrangeMapper;

  template< class GridPart, int polOrder >
  struct PAdaptiveLagrangeMapperTraits
  {
    typedef GridPart GridPartType;
    
    static const int polynomialOrder = polOrder;

    typedef typename GridPartType::template Codim< 0 >::IteratorType::Entity EntityType;
    typedef PAdaptiveLagrangeMapper< GridPartType, polynomialOrder > DofMapperType;
    typedef DefaultDofMapIterator< EntityType, DofMapperType > DofMapIteratorType;

    //! type of the compiled local key 
    typedef LagrangePointSet< GridPartType, polynomialOrder >  CompiledLocalKeyType;
  };



  // First Order Lagrange Mapper
  // ---------------------------

  template< class GridPart >
  class PAdaptiveLagrangeMapper< GridPart, 1 >
  : public CodimensionMapper< GridPart, GridPart::GridType::dimension >
  {
    typedef PAdaptiveLagrangeMapper< GridPart, 1 > ThisType;
    typedef CodimensionMapper< GridPart, GridPart::GridType::dimension > BaseType;

  public:
    //! type of the grid part
    typedef typename BaseType::GridPartType GridPartType;

    //! type of entities (codim 0)
    typedef typename BaseType::EntityType EntityType;

    //! type of the underlying grid
    typedef typename GridPartType::GridType GridType;

    //! type of coordinates within the grid
    typedef typename GridType::ctype FieldType;

    //! dimension of the grid
    static const int dimension = GridType::dimension;

    //! order of the Lagrange polynoms
    static const int polynomialOrder = 1; 

    //! type of the Lagrange point set
    typedef LagrangePointSet< GridPartType, polynomialOrder >
      CompiledLocalKeyType;
    //! type of the map for the Lagrange point sets
    typedef std::map< const GeometryType, const CompiledLocalKeyType* >
      LagrangePointSetMapType;

    typedef std::vector< LagrangePointSetMapType > LagrangePointSetMapVectorType;

  public:
    //! constructor
    PAdaptiveLagrangeMapper ( const GridPartType &gridPart, LagrangePointSetMapVectorType &lagrangePointSet )
    : BaseType( gridPart )
    {}

    bool fixedDataSize ( const int codim ) const
    {
      return true;
    }

    int polynomOrder( const EntityType& entity ) const 
    {
      return 1;
    }

    void setPolynomOrder( const EntityType& entity, const int polOrd ) 
    {
    }
  };



  // Higher Order Lagrange Mapper
  // ----------------------------

  template< class GridPart, int polOrder >
  class PAdaptiveLagrangeMapper 
    : public GenericAdaptiveDofMapper< PAdaptiveLagrangeMapperTraits< GridPart, polOrder > >
  {
  public:   
    // my traits class 
    typedef PAdaptiveLagrangeMapperTraits< GridPart, polOrder > Traits;

  private:   
    typedef PAdaptiveLagrangeMapper< GridPart, polOrder > ThisType;
    typedef GenericAdaptiveDofMapper< Traits > BaseType;

  public:
    //! type of the grid part
    typedef typename Traits::GridPartType GridPartType;

    //! type of the compiled local key 
    typedef typename Traits :: CompiledLocalKeyType CompiledLocalKeyType;

    //! type of the map for the Lagrange point sets
    typedef std::map< const GeometryType, const CompiledLocalKeyType* >
      CompiledLocalKeyMapType;

    //! to be revised 
    typedef std::vector< CompiledLocalKeyMapType > CompiledLocalKeyVectorType;

  public:
    //! constructor
    PAdaptiveLagrangeMapper ( const GridPartType &gridPart,
                              CompiledLocalKeyVectorType &compiledLocalKeys )
      : BaseType( gridPart, compiledLocalKeys )
    {
    }

    //! sort of copy constructor
    PAdaptiveLagrangeMapper ( const PAdaptiveLagrangeMapper& other,
                              CompiledLocalKeyVectorType &compiledLocalKeys )
      : BaseType( other, compiledLocalKeys ) 
    {} 
  };

} // end namespace Dune 

#endif // #ifndef DUNE_LAGRANGESPACE_MAPPER_HH
