#ifndef DUNE_FEM_COMBINEDADAPTMANAGER_HH
#define DUNE_FEM_COMBINEDADAPTMANAGER_HH

#include <dune/common/exceptions.hh>
//- local includes  
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>




namespace Dune{

//***********************************************************************
/** \brief This is a restriction/prolongation operator for DG data. 
 */

  namespace Fem{
    template< class ,
              class >
    class CombinedDiscreteFunctionSpace; 
  }

  
  /** \brief specialization of RestrictProlongDefault for
      LegendreDiscontinuousGalerkinSpace.
  */
  template < class DF, class DFS1, class DFS2 > 
  class RestrictProlongDefaultImplementation< DF, Fem::CombinedDiscreteFunctionSpace< DFS1, DFS2 > >
  : public RestrictProlongInterfaceDefault<   
     RestrictProlongTraits<  RestrictProlongDefaultImplementation < DF, Fem::CombinedDiscreteFunctionSpace< DFS1, DFS2 > >,
                             typename DFS1::GridType::ctype  > >
  {
     typedef  RestrictProlongInterfaceDefault< 
     RestrictProlongTraits<  RestrictProlongDefaultImplementation < DF, Fem::CombinedDiscreteFunctionSpace< DFS1, DFS2 > >,
                             typename DFS1::GridType::ctype > > BaseType;

    typedef typename BaseType::DomainFieldType DomainFieldType;

  public:
    typedef DF DiscreteFunctionType;
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    typedef CachingQuadrature<GridPartType,0> QuadratureType;
    typedef typename GridType::template Codim<0>::Entity::LocalGeometry LocalGeometry;

  protected:
    using BaseType :: calcWeight;
    using BaseType :: entitiesAreCopies;
    
  public:  
    //! Constructor
    explicit RestrictProlongDefaultImplementation( DiscreteFunctionType &df )
    {}

    /** \brief explicit set volume ratio of son and father
     *
     *  \param[in]  weight  volume of son / volume of father
     *
     *  \note If this ratio is set, it is assume to be constant.
     */
    void setFatherChildWeight ( const DomainFieldType &weight ) const
    {
      DUNE_THROW(NotImplemented, "Not yet impl."); 
    }

    //! restrict data to father 
    template< class EntityType >
    void restrictLocal ( const EntityType &father, const EntityType &son, bool initialize ) const
    {
      DUNE_THROW(NotImplemented, "Not yet impl."); 
    }

    //! prolong data to children 
    template< class EntityType >
    void prolongLocal ( const EntityType &father, const EntityType &son, bool initialize ) const
    {
      DUNE_THROW(NotImplemented, "Not yet impl."); 
    }

    //! add discrete function to communicator 
    template <class CommunicatorImp>
    void addToList(CommunicatorImp& comm)
    {
      DUNE_THROW(NotImplemented, "Not yet impl."); 
    }

  };
} // end namespace Dune 
#endif
