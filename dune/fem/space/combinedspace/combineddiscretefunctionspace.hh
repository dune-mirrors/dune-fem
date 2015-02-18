#ifndef DUNE_FEM_COMBINEDDISCRETFUNCTIONSPACE_HH
#define DUNE_FEM_COMBINEDDISCRETFUNCTIONSPACE_HH

#include <algorithm>

#include <dune/common/math.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/grid.hh>

#include <dune/fem/space/combinedspace/tuplespace.hh>


namespace Dune
{

  namespace Fem
  {

    /** \addtogroup CombinedDiscreteFunctionSpace
     *
     *  Provides a DiscreteFunctionSpace combined from two arbitrary
     *  DiscreteFunctionSpaces U_h and V_h into a single \ref Dune::Fem::DiscreteFunctionSpaceInterface ( U_h times V_h ).
     *
     *  \note It is assumed that the spaces U_h and V_h  are constructed on the same
     *        gridpart!!!
     */

    /** \class   CombinedDiscreteFunctionSpace
     *  \ingroup CombinedDiscreteFunctionSpace
     *  \brief   Combined discrete function space
     */
    template< class DFSpace1, class DFSpace2 >
    class CombinedDiscreteFunctionSpace
      : public TupleDiscreteFunctionSpace< DFSpace1, DFSpace2 >
    {
      typedef CombinedDiscreteFunctionSpace< DFSpace1, DFSpace2 > ThisType;
      typedef TupleDiscreteFunctionSpace< DFSpace1, DFSpace2 > BaseType;

    public:
      //! types of Discrete Subspace
      typedef DFSpace1 DiscreteFunctionSpaceType1;
      typedef DFSpace2 DiscreteFunctionSpaceType2;

      //! maximum polynomial order of functions in this space
      enum { polynomialOrder1 =  DiscreteFunctionSpaceType1::polynomialOrder,
             polynomialOrder2 =  DiscreteFunctionSpaceType2::polynomialOrder,
             polynomialOrder = BaseType::polynomialOrder };

      template< int newDimRange >
      struct ToNewDimRange
      {
        typedef typename conditional< (newDimRange == 1),
                                      typename DiscreteFunctionSpaceType1::template ToNewDimRange< 1 >::Type,
                                      typename BaseType::template ToNewDimRange< newDimRange >::Type
                                      >::type Type;
      };

      // type of gridPart
      typedef typename BaseType::GridPartType GridPartType;

      /** \brief constructor
       *
       *  \param[in]  gridPart       grid part for the Lagrange space
       *  \param[in]  commInterface  communication interface to use (optional)
       *  \param[in]  commDirection  communication direction to use (optional)
       */
      explicit CombinedDiscreteFunctionSpace ( GridPartType &gridPart,
                                               const InterfaceType commInterface = InteriorBorder_All_Interface,
                                               const CommunicationDirection commDirection = ForwardCommunication )
      DUNE_DEPRECATED_MSG( "CommunicationDirection is Deprecated, us the more general TupleDiscreteFunctionSpace instead." )
        : BaseType( gridPart, commInterface, commDirection )
      {}

      // forbid the copy constructor
      CombinedDiscreteFunctionSpace ( const ThisType & ) = delete;

      /** \brief obtain the first subspace
       *  \return DiscreteFunctionSpaceType1
       **/
      const DiscreteFunctionSpaceType1 &space1 () const
      {
        return BaseType::template subDiscreteFunctionSpace< 0 >();
      }

      /** \brief obtain the second subspace
       *  \return DiscreteFunctionSpaceType2
       **/
      const DiscreteFunctionSpaceType2 &space2 () const
      {
        return BaseType::template subDiscreteFunctionSpace< 1 >();
      }
    };


    //! specialization of DifferentDiscreteFunctionSpace for CombinedDiscreteFunctionSpace
    template< class DFunctionSpaceImp1,
              class DFunctionSpaceImp2,
              class NewFunctionSpace >
    struct DifferentDiscreteFunctionSpace< Fem::CombinedDiscreteFunctionSpace<
                                             DFunctionSpaceImp1, DFunctionSpaceImp2 >, NewFunctionSpace >
    {
    private:
      static const int dimRange1 = DFunctionSpaceImp1::dimRange;
      static const int dimRange2 = DFunctionSpaceImp2::dimRange;
      static const int newDimRange = NewFunctionSpace::dimRange;

      static const int newDimRange1 = (newDimRange * dimRange1)/( dimRange1 + dimRange2 );
      static const int newDimRange2 = (newDimRange * dimRange2)/( dimRange1 + dimRange2 );

      //////////  Gurke in gruen /////////
      typedef typename DFunctionSpaceImp1::template ToNewDimRange< newDimRange1 >::Type Type1;
      typedef typename DFunctionSpaceImp2::template ToNewDimRange< newDimRange2 >::Type Type2;

    public:
      typedef CombinedDiscreteFunctionSpace< Type1, Type2 > Type;
    };


    // DefaultLocalRestrictProlong for CombinedDiscreteFunctionSpace

    template< class SP1, class SP2 >
    class DefaultLocalRestrictProlong< CombinedDiscreteFunctionSpace< SP1, SP2 > >
      : public TupleLocalRestrictProlong< SP1, SP2 >
    {
      typedef DefaultLocalRestrictProlong< CombinedDiscreteFunctionSpace< SP1, SP2 > > ThisType;
      typedef TupleLocalRestrictProlong< SP1, SP2 > BaseType;

    public:
      DefaultLocalRestrictProlong ( const CombinedDiscreteFunctionSpace< SP1, SP2 > &space )
        : BaseType( space.space1(), space.space2() )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_COMBINDEDDISCRETFUNCTIONSPACE_HH
