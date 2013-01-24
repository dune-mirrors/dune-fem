#ifndef DUNE_FEM_SPACE_PADAPTIVE_LAGRANGE_HH
#define DUNE_FEM_SPACE_PADAPTIVE_LAGRANGE_HH

#include <dune/fem/operator/lagrangeinterpolation.hh>
#include <dune/fem/space/basefunctions/basefunctionproxy.hh>
#include <dune/fem/space/basefunctions/basefunctionsets.hh>
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>

#include "adaptmanager.hh"
#include "declaration.hh"
#include "generic.hh"
#include "mapper.hh"
#include "restrictprolong.hh"


namespace Dune
{

  namespace Fem 
  {

    /** \addtogroup PAdaptiveLagrangeSpace
     *
     *  Provides access to base function sets for different element types in
     *  one grid and size of function space and maps from local to global dof
     *  number.
     *
     *  \note This space can only be used with special index sets. If you want
     *  to use the PAdaptiveLagrangeSpace with an index set only
     *  supporting the index set interface you will have to use the
     *  IndexSetWrapper class to provide the required functionality.
     *
     *  \note For adaptive calculations one has to use index sets that are
     *  capable of adaption (i.e. the method adaptive returns true). See also
     *  AdaptiveLeafIndexSet.
     */

    // PAdaptiveLagrangeSpaceTraits
    // ----------------------------

    template< class FunctionSpace, class GridPart, unsigned int polOrder,
              template< class > class BaseFunctionStorage >
    struct PAdaptiveLagrangeSpaceTraits
    {
      dune_static_assert((polOrder > 0), "LagrangeSpace only defined for polOrder > 0" );
      
      static const int codimension = 0;

      typedef FunctionSpace FunctionSpaceType;
      typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType :: DomainType DomainType;
      typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
      typedef typename FunctionSpaceType :: RangeType RangeType;
      typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

      enum { dimRange = FunctionSpaceType :: dimRange };
      
      typedef GridPart GridPartType;
      typedef typename GridPartType :: GridType GridType;
      typedef typename GridPartType :: IndexSetType IndexSetType;
      typedef typename GridPartType :: template Codim< codimension > :: IteratorType
        IteratorType;

      // get dimension of local coordinate 
      enum { dimLocal = GridType :: dimension };

      typedef typename ToLocalFunctionSpace< FunctionSpaceType, dimLocal > :: Type 
        BaseFunctionSpaceType;

      enum { polynomialOrder = polOrder };
      
      typedef PAdaptiveLagrangeSpace
        < FunctionSpaceType, GridPartType, polynomialOrder, BaseFunctionStorage >
        DiscreteFunctionSpaceType;

      enum { localBlockSize = dimRange };

      //! this is a continuous space 
      static const bool continuousSpace = true ;

      // mapper for block
      typedef PAdaptiveLagrangeMapper< GridPartType, polynomialOrder > BlockMapperType;
      typedef NonBlockMapper< BlockMapperType, localBlockSize > MapperType;
      
      // implementation of shapefunction set 
      typedef VectorialBaseFunctionSet< BaseFunctionSpaceType, BaseFunctionStorage >
        ShapeFunctionSetType ;

      // exported type of base function set 
      typedef SimpleBaseFunctionProxy< ShapeFunctionSetType > BaseFunctionSetType;

      //! type of a compiled local key 
      typedef LagrangePointSet< GridPartType, polynomialOrder >
        CompiledLocalKeyType;

      /** \brief defines type of communication data handle for this type of space
       */
      template< class DiscreteFunction,
                class Operation = DFCommunicationOperation :: Add >
      struct CommDataHandle
      {
        //! type of data handle 
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        //! type of operatation to perform on scatter 
        typedef Operation OperationType;
      };
    };



    // PAdaptiveLagrangeSpace
    // ----------------------

    /** \class   PAdaptiveLagrangeSpace
     *  \ingroup PAdaptiveLagrangeSpace
     *  \brief   Lagrange discrete function space
     */
    template< class FunctionSpaceImp,
              class GridPartImp,
              int polOrder,
              template< class > class BaseFunctionStorageImp = CachingStorage >
    class PAdaptiveLagrangeSpace
    : public GenericDiscreteFunctionSpace
             < PAdaptiveLagrangeSpaceTraits< FunctionSpaceImp,
                                             GridPartImp,
                                             polOrder,
                                             BaseFunctionStorageImp > >
    {
    public:
      //! traits for the discrete function space
      typedef PAdaptiveLagrangeSpaceTraits< FunctionSpaceImp,
                                            GridPartImp,
                                            polOrder,
                                            BaseFunctionStorageImp >
        Traits;

      typedef GenericDiscreteFunctionSpace< Traits > BaseType ;

      //! type of the discrete function space
      typedef PAdaptiveLagrangeSpace< FunctionSpaceImp,
                                      GridPartImp,
                                      polOrder,
                                      BaseFunctionStorageImp >
              PAdaptiveLagrangeSpaceType;

      typedef typename Traits :: GridPartType GridPartType;
      typedef typename Traits :: GridType GridType;
      typedef typename Traits :: IndexSetType IndexSetType;
      typedef typename Traits :: IteratorType IteratorType;

      //! maximum polynomial order of functions in this space
      enum { polynomialOrder = Traits :: polynomialOrder };
      
      //! type of compiled local key 
      typedef typename Traits :: CompiledLocalKeyType  CompiledLocalKeyType;

      // deprecated name 
      typedef CompiledLocalKeyType LagrangePointSetType;

      //! mapper used to implement mapToGlobal
      typedef typename Traits :: MapperType MapperType;

      //! mapper used to for block vector function 
      typedef typename Traits :: BlockMapperType BlockMapperType;

      //! size of local blocks
      enum { localBlockSize = Traits :: localBlockSize };

      //! dimension of a value
      enum { dimVal = 1 };

      //! type of DoF manager
      typedef DofManager< GridType > DofManagerType;

      //! type of intersections
      typedef typename BaseType ::IntersectionType IntersectionType;

    public:
      using BaseType :: gridPart;
      using BaseType :: blockMapper;
      using BaseType :: compiledLocalKey;
      using BaseType :: order;

      //! default communication interface 
      static const InterfaceType defaultInterface = InteriorBorder_InteriorBorder_Interface;

      //! default communication direction 
      static const CommunicationDirection defaultDirection = ForwardCommunication;

      /** \brief constructor
       *
       *  \param[in]  gridPart       grid part for the Lagrange space
       *  \param[in]  commInterface  communication interface to use (optional)
       *  \param[in]  commDirection  communication direction to use (optional)
       */
      explicit PAdaptiveLagrangeSpace
        ( GridPartType &gridPart,
          const InterfaceType commInterface = defaultInterface,
          const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, commInterface, commDirection )
      {
      }

      //! copy constructor needed for p-adaption 
      PAdaptiveLagrangeSpace( const PAdaptiveLagrangeSpace& other ) 
      : BaseType( other )
      {
      }

    protected:
      using BaseType :: dfList_ ;
      using BaseType :: searchFunction ;

    public:
      /*! \brief add function to discrete function space for p-adaptation 
          (currently only supported by AdaptiveDiscreteFunction )
       */
      template <class DiscreteFunction> 
      void addFunction( DiscreteFunction& df ) const
      {
        assert( searchFunction( df ) == dfList_.end() );
        // select LagrangeInterpolation to be the LocalInterpolation 
        typedef typename BaseType :: template PAdaptiveDiscreteFunctionEntry< 
            DiscreteFunction, LagrangeInterpolation< DiscreteFunction, DiscreteFunction > > RealEntryType ;
        typedef typename BaseType :: PAdaptiveDiscreteFunctionEntryInterface
          EntryInterface;
        EntryInterface* entry = new RealEntryType( df );

        assert( entry );
        dfList_.push_front( entry );
      }

      //! deprecated method 
      template< class EntityType >
      inline const CompiledLocalKeyType &lagrangePointSet( const EntityType &entity ) const
      {
        return compiledLocalKey( entity.type(),
                                 blockMapper().polynomOrder( entity ) );
      }

      //! deprecated method 
      inline const CompiledLocalKeyType &lagrangePointSet( const GeometryType type ) const
      {
        return compiledLocalKey( type, polynomialOrder );
      }

      //! deprecated method 
      inline const CompiledLocalKeyType &lagrangePointSet( const GeometryType type, const int order ) const
      {
        return compiledLocalKey( type, order );
      }

      using BaseType::continuous;
      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      inline bool continuous (const IntersectionType &intersection) const
      { 
        if ( order() > 0 && intersection.conforming())
        {
          return true;
          if (intersection.neighbor())
            return (order(*(intersection.inside())) == order(*(intersection.outside())));
          else
            return true;
        }
        return false;
      }
    };

  } // namespace Fem
    
} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_PADAPTIVE_LAGRANGE_HH
