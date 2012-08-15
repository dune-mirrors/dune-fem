#ifndef DUNE_FEM_COMBINEDDISCRETFUNCTIONSPACE_HH
#define DUNE_FEM_COMBINEDDISCRETFUNCTIONSPACE_HH

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

      typedef FunctionSpace< DomainFieldType, RangeFieldType, DomainType :: dimension, dimRange> FunctionSpaceType;

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
      //! mapper for block
      typedef typename DiscreteFunctionSpace1 :: MapperType     MapperType1;
      typedef typename DiscreteFunctionSpace2 :: MapperType     MapperType2;

      //! get base function sets of the two spaces
      typedef typename DiscreteFunctionSpace1 :: BaseFunctionSetType   BaseFunctionSetType1;
      typedef typename DiscreteFunctionSpace2 :: BaseFunctionSetType   BaseFunctionSetType2;

      public:
      //! define a combined DofMapper and the block mapper
      typedef CombinedMapper< GridType, MapperType1, MapperType2 > BlockMapperType;
      typedef BlockMapperType MapperType;

      //! implementation of basefunction set 
      typedef CombinedBaseFunctionSet< FunctionSpaceType, BaseFunctionSetType1, BaseFunctionSetType2 >
          BaseFunctionSetType;   

      template< class DiscreteFunction,
                class Operation = DFCommunicationOperation :: Copy >
      struct CommDataHandle
      {
        //! type of data handle 
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        //! type of operatation to perform on scatter 
        typedef Operation OperationType;
      };
    };



    /** \addtogroup CombinedDiscreteFunctionSpace
     *
     *  Provides a DiscreteFunctionSpace combined from two arbitrary 
     *  DiscreteFunctionSpaces U_h and V_h into a single \ref Dune::Fem::DiscreteFunctionSpaceInterface ( U_h times V_h ). 
     *
     *  \note It is assumed that the spaces U_h and V_h  are constructed on the same
     *        gridpart!!!
     *
     *  \note adaptivity is not yet implemented!!!
     */

    /** \class   CombinedDiscreteFunctionSpace
     *  \ingroup CombinedDiscreteFunctionSpace
     *  \brief   Combined discrete function space
     */
    template< class DFSpace1, class DFSpace2 >
    class CombinedDiscreteFunctionSpace
    : public DiscreteFunctionSpaceDefault
             < CombinedDiscreteFunctionSpaceTraits< DFSpace1,
                                                    DFSpace2 > >
    {
    public:
      //! types of Discrete Subspace 
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

      //! extract grid informations, it is assumed the both spaces are living on the
      //! same gridPart
      typedef typename Traits :: GridPartType GridPartType;
      typedef typename Traits :: GridType GridType;
      //! extract informations about IndexSet and Iterators
      typedef typename Traits :: IndexSetType IndexSetType;
      typedef typename Traits :: IteratorType IteratorType;
      //! dimension of the grid (not the world)
      enum { dimension = GridType :: dimension };

      //! the underlaying Analytical function space
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

      //! types of the dof mapper for the two subspaces
      typedef typename DiscreteFunctionSpaceType1 :: MapperType MapperType1;
      typedef typename DiscreteFunctionSpaceType2 :: MapperType MapperType2;

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

    public:
      //! type of identifier for this discrete function space
      typedef int IdentifierType;
      //! identifier of this discrete function space
      static const IdentifierType id = 669;
      
    private:
      typedef CombinedDiscreteFunctionSpaceType ThisType;
      typedef DiscreteFunctionSpaceDefault< Traits > BaseType;


    public:
      using BaseType :: gridPart;

      //! default communication interface 
      static const InterfaceType defaultInterface = InteriorBorder_All_Interface;

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
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::contains */
      inline bool contains ( const int codim ) const
      {
        // forward to mapper since this information is held there 
        return blockMapper().contains( codim );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      inline bool continuous () const
      {
        // forward to the subsapces
        return space1_.continuous() && space2_.continuous();
      }

      /** \brief get the type of this discrete function space 
          \return DFSpaceIdentifier
      **/
      inline DFSpaceIdentifier type () const
      {
        return CombinedSpace_id;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      inline int order () const
      {
        return polynomialOrder;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      template<class Entity>
      inline int order ( const Entity &entity ) const
      {
        return polynomialOrder;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::baseFunctionSet(const EntityType &entity) const */
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

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::mapper */
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

      /** \brief obtain the first subspace
       *  \return DiscreteFunctionSpaceType1
       **/
      const DiscreteFunctionSpaceType1 &space1 () const
      {
        return space1_;
      }

      /** \brief obtain the second subspace
       *  \return DiscreteFunctionSpaceType2
       **/
      const DiscreteFunctionSpaceType2 &space2 () const
      {
        return space2_;
      }

    protected:
      //! space 1
      DiscreteFunctionSpaceType1  space1_;
      //! space 2
      DiscreteFunctionSpaceType2  space2_;

      //! corresponding mapper
      MapperType *mapper_;
    };

    //! specialization of DifferentDiscreteFunctionSpace for CombinedDiscreteFunctionSpace
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

  } // namespace Fem

} // namespace Dune 

#endif // #ifndef DUNE_FEM_COMBINDEDDISCRETFUNCTIONSPACE_HH
