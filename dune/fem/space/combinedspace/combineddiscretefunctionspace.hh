#ifndef DUNE_COMBINEDDISCRETFUNCTIONSPACE_HH
#define DUNE_COMBINEDDISCRETFUNCTIONSPACE_HH

//- system includes 
#include <algorithm>

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/grid/common/grid.hh>


//- Dune-Fem includes 
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/basefunctions/basefunctionsets.hh>
#include <dune/fem/space/basefunctions/basefunctionproxy.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>

//- local includes 
#include "combinedbasefunctions.hh"
#include "combinedmapper.hh"
#include "adaptmanager.hh"

namespace Dune
{
  namespace Fem
  {

    template< class DFunctionSpaceImp1,
              class DFunctionSpaceImp2>
    class CombinedDiscreteFunctionSpace;



    template< class DFunctionSpace1, class DFunctionSpace2> 
    struct CombinedDiscreteFunctionSpaceTraits
    {

      typedef DFunctionSpace1   DiscreteFunctionSpace1;
      typedef DFunctionSpace2   DiscreteFunctionSpace2;   

      typedef typename DiscreteFunctionSpace1 :: FunctionSpaceType     FunctionSpaceType1; 
      typedef typename DiscreteFunctionSpace2 :: FunctionSpaceType     FunctionSpaceType2; 
      
      typedef typename FunctionSpaceType1 :: DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType1 :: DomainType DomainType;
      typedef typename FunctionSpaceType1 :: RangeFieldType RangeFieldType;
      typedef typename FunctionSpaceType1 :: RangeType RangeType1;
      typedef typename FunctionSpaceType1 :: JacobianRangeType JacobianRangeType1;
      typedef typename FunctionSpaceType2 :: RangeType RangeType2;
      typedef typename FunctionSpaceType2 :: JacobianRangeType JacobianRangeType2;

      enum { dimRange = FunctionSpaceType1 :: dimRange + FunctionSpaceType2 :: dimRange };

      typedef Dune :: FunctionSpace< DomainFieldType, RangeFieldType, DomainType :: dimension, dimRange> FunctionSpaceType;

      typedef typename FunctionSpaceType :: RangeType RangeType;
      typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

      // assume for the first that DFSpace1 and DFSpace2 have the same gridpart
      typedef typename DiscreteFunctionSpace1 :: GridPartType GridPartType;
      typedef typename GridPartType :: GridType GridType;
      typedef typename GridPartType :: IndexSetType IndexSetType;
      typedef typename GridPartType :: template Codim< 0 > :: IteratorType
        IteratorType;

      typedef typename IteratorType :: Entity EntityType;    

      // get dimension of local coordinate 
      enum { dimLocal = GridType :: dimension };

      typedef typename ToLocalFunctionSpace< FunctionSpaceType, dimLocal > :: Type
        BaseFunctionSpaceType;

      enum { polynomialOrder1 =  DiscreteFunctionSpace1 :: polynomialOrder, 
             polynomialOrder2 =  DiscreteFunctionSpace2 :: polynomialOrder,
             polynomialOrder =  ( polynomialOrder1 > polynomialOrder2 ? polynomialOrder1 : polynomialOrder2 )   };
      
      typedef CombinedDiscreteFunctionSpace < DiscreteFunctionSpace1, DiscreteFunctionSpace2 >
        DiscreteFunctionSpaceType;

      enum { localBlockSize = 1 };

      private:
      // mapper for block
      typedef typename DiscreteFunctionSpace1 :: MapperType     MapperType1;
      typedef typename DiscreteFunctionSpace2 :: MapperType     MapperType2;

      // get base function sets of the two spaces
      typedef typename DiscreteFunctionSpace1 :: BaseFunctionSetType   BaseFunctionSetType1;
      typedef typename DiscreteFunctionSpace2 :: BaseFunctionSetType   BaseFunctionSetType2;

      public:
      // atm no other possibility
      typedef CombinedMapper< GridType, MapperType1, MapperType2 > BlockMapperType;
      typedef BlockMapperType MapperType;

      // implementation of basefunction set 
      typedef CombinedBaseFunctionSet< FunctionSpaceType, BaseFunctionSetType1, BaseFunctionSetType2 >
          BaseFunctionSetType;   

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



    /** \addtogroup LagrangeDiscreteFunctionSpace
     *
     *  Provides access to bse function sets for different element types in
     *  one grid and size of function space and maps from local to global dof
     *  number.
     *
     *  \note This space can only be used with special index sets. If you want
     *  to use the LagrangeDiscreteFunctionSpace with an index set only
     *  supporting the index set interface you will have to use the
     *  IndexSetWrapper class to provide the required functionality.
     *
     *  \note For adaptive calculations one has to use index sets that are
     *  capable of adaption (i.e. the method adaptive returns true). See also
     *  AdaptiveLeafIndexSet.
     */



    /** \class   LagrangeDiscreteFunctionSpace
     *  \ingroup LagrangeDiscreteFunctionSpace
     *  \brief   Lagrange discrete function space
     */
    template< class DFSpace1, class DFSpace2 >
    class CombinedDiscreteFunctionSpace
    : public DiscreteFunctionSpaceDefault
             < CombinedDiscreteFunctionSpaceTraits< DFSpace1,
                                                    DFSpace2 > >
    {
    public:
      typedef DFSpace1    DiscreteFunctionSpaceType1;
      typedef DFSpace2    DiscreteFunctionSpaceType2;

      //! traits for the discrete function space
      typedef CombinedDiscreteFunctionSpaceTraits< DFSpace1,
                                                   DFSpace2 > 
        Traits;

      //! type of the discrete function space
      typedef CombinedDiscreteFunctionSpace< DFSpace1,
                                             DFSpace2 >
       CombinedDiscreteFunctionSpaceType;

      typedef typename Traits :: GridPartType GridPartType;
      typedef typename Traits :: GridType GridType;
      typedef typename Traits :: IndexSetType IndexSetType;
      typedef typename Traits :: IteratorType IteratorType;
      //! dimension of the grid (not the world)
      enum { dimension = GridType :: dimension };

      typedef typename Traits :: FunctionSpaceType FunctionSpaceType;
      //! field type for function space's domain
      typedef typename Traits :: DomainFieldType DomainFieldType;
      //! type for function space's domain
      typedef typename Traits :: DomainType DomainType;
      //! field type for function space's range
      typedef typename Traits :: RangeFieldType RangeFieldType;
      //! type for function space's range
      typedef typename Traits :: RangeType RangeType;
      //! dimension of function space's range
      enum { dimRange = FunctionSpaceType :: dimRange };
      //! type of scalar function space
      typedef typename Traits :: BaseFunctionSpaceType BaseFunctionSpaceType;

      //! maximum polynomial order of functions in this space
      enum { polynomialOrder1 = Traits :: polynomialOrder1,
             polynomialOrder2 = Traits :: polynomialOrder2,
             polynomialOrder = Traits :: polynomialOrder
      };
      
      //! type of the base function set(s)
      typedef typename Traits :: BaseFunctionSetType BaseFunctionSetType;

      typedef typename DiscreteFunctionSpaceType1 :: MapperType MapperType1;
      typedef typename DiscreteFunctionSpaceType2 :: MapperType MapperType2;

#if 0
      //! type of base function factory
      typedef CombinedBaseFunctionFactory
        < typename BaseFunctionSpaceType :: ScalarFunctionSpaceType, dimension, polynomialOrder >
        ScalarFactoryType;
      //! type of singleton base function factory
      typedef BaseFunctionSetSingletonFactory
        < GeometryType, BaseFunctionSetImp, ScalarFactoryType >
        BaseFunctionSetSingletonFactoryType;
      //! type of singleton list (singleton provider) for base functions
      typedef SingletonList
        < GeometryType, BaseFunctionSetImp, BaseFunctionSetSingletonFactoryType >
        BaseFunctionSetSingletonProviderType;


      //! type of a Lagrange point set
      typedef LagrangePointSet< GridPartType, polynomialOrder >
        LagrangePointSetType;
      //! type of Lagrange point set map
      typedef std :: map< const GeometryType, const LagrangePointSetType* >
        LagrangePointSetMapType;
#endif

      //! mapper used to implement mapToGlobal
      typedef typename Traits :: MapperType MapperType;

      //! mapper used to for block vector function 
      typedef typename Traits :: BlockMapperType BlockMapperType;

      //! size of local blocks
      enum { localBlockSize = Traits :: localBlockSize };

      //! type for DoF
      typedef RangeFieldType DofType;
      //! dimension of a value
      enum { dimVal = 1 };

#if 0
      //! mapper factory 
      typedef CombinedMapperSingletonFactory
        < BlockMapperType1, BlockMapperType2 >
        BlockMapperSingletonFactoryType;

      typedef typename BlockMapperSingletonFactoryType :: Key KeyType;

      //! singleton list of mappers 
      typedef SingletonList
        < KeyType, BlockMapperType, BlockMapperSingletonFactoryType >
        BlockMapperProviderType;
#endif

    public:
      //! type of identifier for this discrete function space
      typedef int IdentifierType;
      //! identifier of this discrete function space
      static const IdentifierType id = 667;
      
    private:
      typedef CombinedDiscreteFunctionSpaceType ThisType;
      typedef DiscreteFunctionSpaceDefault< Traits > BaseType;


    public:
      using BaseType :: gridPart;

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
      explicit CombinedDiscreteFunctionSpace
        ( GridPartType &gridPart, 
          const InterfaceType commInterface = defaultInterface,
          const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        space1_( gridPart, commInterface, commDirection ),
        space2_( gridPart, commInterface, commDirection ),
        mapper_( 0 )
      {
#if 0
        const IndexSetType &indexSet = gridPart.indexSet();

        AllGeomTypes< IndexSetType, GridType > allGeometryTypes( indexSet );
        const std :: vector< GeometryType >& geometryTypes
          = allGeometryTypes.geomTypes( 0 );
        for( unsigned int i = 0; i < geometryTypes.size(); ++i )
        {
          const GeometryType &geometryType = geometryTypes[ i ];
          
          if( baseFunctionSet_.find( geometryType ) == baseFunctionSet_.end() )
          {
            const BaseFunctionSetImp *baseFunctionSet; 
            baseFunctionSet = new BaseFunctionSetImp( 
                space1_.baseFunctionSet( geometryType ), 
                space2_.baseFunctionSet( geometryType ) );
            /*
              = &(BaseFunctionSetSingletonProviderType
                  :: getObject( geometryType ));
            assert( baseFunctionSet != NULL );
            */
            baseFunctionSet_[ geometryType ] = baseFunctionSet;
          }
        }
#endif

        /*
        MapperSingletonKeyType key( gridPart, lagrangePointSet_, polynomialOrder );
        blockMapper_ = &BlockMapperProviderType :: getObject( key );
        */
        mapper_ = new MapperType( gridPart.grid(), space1_.mapper(), space2_.mapper());
        assert( mapper_ != 0 );
      }

    private:
      // forbid the copy constructor
      CombinedDiscreteFunctionSpace ( const ThisType& );

    public:
      /** \brief Destructor (freeing base functions and mapper)
          \return 
      **/
      ~CombinedDiscreteFunctionSpace ()
      {
        delete mapper_;
  /*
        BlockMapperProviderType::removeObject( *blockMapper_ );
        */

#if 0
        typedef typename BaseFunctionMapType :: iterator BFIteratorType;
        BFIteratorType bfend = baseFunctionSet_.end();
        for( BFIteratorType it = baseFunctionSet_.begin(); it != bfend; ++it ) 
        {
          const BaseFunctionSetImp *baseFunctionSet = (*it).second;        
          if( baseFunctionSet != NULL )
            delete baseFunctionSet;
          /*
            BaseFunctionSetSingletonProviderType
            :: removeObject( *baseFunctionSet );
            */
        }
#endif
      }

      /** \copydoc Dune::DiscreteFunctionSpaceInterface::contains */
      inline bool contains ( const int codim ) const
      {
        // forward to mapper since this information is held there 
        return blockMapper().contains( codim );
      }

      /** \copydoc Dune::DiscreteFunctionSpaceInterface::continuous */
      inline bool continuous () const
      {
        return space1_.continuous() && space2_.continuous();
      }

      /** \brief get the type of this discrete function space 
          \return DFSpaceIdentifier
      **/
      inline DFSpaceIdentifier type () const
      {
        return CombinedSpace_id;
      }

      /** \copydoc Dune::DiscreteFunctionSpaceInterface::order */
      inline int order () const
      {
        return polynomialOrder;
      }

      /** \copydoc Dune::DiscreteFunctionSpaceInterface::baseFunctionSet(const EntityType &entity) const */
      template< class EntityType >
      BaseFunctionSetType baseFunctionSet ( const EntityType &entity ) const
      {
        return BaseFunctionSetType( space1_.baseFunctionSet( entity ), space2_.baseFunctionSet( entity ) );
      }

      /** \brief get dimension of value
          \return int
      **/
      inline int dimensionOfValue () const
      {
        return dimVal;
      }

      /** \copydoc Dune::DiscreteFunctionSpaceInterface::mapper */
      MapperType &mapper () const
      {
        assert( mapper_ != 0 );
        return *mapper_;
      }

      /** \brief obtain the DoF block mapper of this space
          \return BlockMapperType
      **/
      BlockMapperType &blockMapper () const
      {
        return mapper();
      }

      const DiscreteFunctionSpaceType1 &space1 () const
      {
        return space1_;
      }

      const DiscreteFunctionSpaceType2 &space2 () const
      {
        return space2_;
      }

    protected:
      //! space 1
      DiscreteFunctionSpaceType1  space1_;
      //! space 2
      DiscreteFunctionSpaceType2  space2_;

#if 0
      //! map for the base function sets
      mutable BaseFunctionMapType baseFunctionSet_;
#endif

      //! corresponding mapper
      MapperType *mapper_;
    };

  } // end namespace Fem
  template< class DFunctionSpaceImp1,
            class DFunctionSpaceImp2,
            class NewFunctionSpace >
  class DifferentDiscreteFunctionSpace< Fem::CombinedDiscreteFunctionSpace<
        DFunctionSpaceImp1, DFunctionSpaceImp2 >, NewFunctionSpace >
  {
    static const int dimRange1 = DFunctionSpaceImp1 :: dimRange;
    static const int dimRange2 = DFunctionSpaceImp2 :: dimRange;
    static const int newDimRange = NewFunctionSpace :: dimRange;

    static const int newDimRange1 = (newDimRange * dimRange1)/( dimRange1 + dimRange2 );
    static const int newDimRange2 = (newDimRange * dimRange2)/( dimRange1 + dimRange2 );

    //////////  Gurke in gruen /////////
    typedef typename DFunctionSpaceImp1 :: template ToNewDimRange< newDimRange1 > :: Type Type1;
    typedef typename DFunctionSpaceImp2 :: template ToNewDimRange< newDimRange2 > :: Type Type2;

    public:
    typedef Fem::CombinedDiscreteFunctionSpace< Type1, Type2 > Type;
  };
} // end namespace Dune 


#endif // #ifndef DUNE_COMBINDEDDISCRETFUNCTIONSPACE_HH
