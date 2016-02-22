#ifndef DUNE_FEM_SPACE_COMBINEDSPACE_COMBINEDSPACE_HH
#define DUNE_FEM_SPACE_COMBINEDSPACE_COMBINEDSPACE_HH

#include <map>
#include <type_traits>
#include <vector>

#include <dune/fem/version.hh>

#include <dune/fem/space/combinedspace/powerspace.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/dofstorage.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>


namespace Dune
{

  namespace Fem
  {

    /** @addtogroup CombinedSpace
        Class to combine N scalar discrete function spaces.
        Policies PointBased and VariableBased decide, how dof are stored in
        vectors. PointBased stores all local dofs consecutive,
        VectorBased stores all dofs for one component consecutive.
       @{

     */

    template< class DiscreteFunctionSpace, int N, DofStoragePolicy policy >
    class CombinedSpace;

    /** @brief
        Combined Space Function Space
     **/

    template< class DiscreteFunctionSpace, int N >
    class CombinedSpace< DiscreteFunctionSpace, N, PointBased >
      : public DiscreteFunctionSpace::template ToNewDimRange< DiscreteFunctionSpace::dimRange *N >::Type
    {
      typedef CombinedSpace< DiscreteFunctionSpace, N, PointBased > ThisType;
      typedef typename DiscreteFunctionSpace::template ToNewDimRange< DiscreteFunctionSpace::dimRange *N >::Type BaseType;

      static const DofStoragePolicy policy = PointBased;

    public:
      typedef typename BaseType::GridPartType GridPartType;
      typedef DiscreteFunctionSpace ContainedDiscreteFunctionSpaceType;

      CombinedSpace ( GridPartType &gridPart,
                      const InterfaceType commInterface = InteriorBorder_All_Interface,
                      const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, commInterface, commDirection ),
          containedSpace_( gridPart, commInterface, commDirection )
      {}

      CombinedSpace ( const ThisType& ) = delete;
      ThisType& operator= ( const ThisType& ) = delete;

      //- Additional methods
      //! number of components
      int numComponents () const
      {
        return N;
      }

      //! policy of this space
      DofStoragePolicy myPolicy () const
      {
        return policy;
      }

      //! contained space
      const ContainedDiscreteFunctionSpaceType &containedSpace () const
      {
        return containedSpace_;
      }

    private:
      ContainedDiscreteFunctionSpaceType containedSpace_;
    };


    template< class DiscreteFunctionSpace, int N >
    class CombinedSpace< DiscreteFunctionSpace, N, VariableBased >
      : public PowerDiscreteFunctionSpace< DiscreteFunctionSpace, N >
    {
      typedef CombinedSpace< DiscreteFunctionSpace, N, VariableBased > ThisType;
      typedef PowerDiscreteFunctionSpace< DiscreteFunctionSpace, N > BaseType;

      static const DofStoragePolicy policy = VariableBased;

    public:
      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::ContainedDiscreteFunctionSpaceType ContainedDiscreteFunctionSpaceType;

      CombinedSpace ( GridPartType &gridPart,
                      const InterfaceType commInterface = InteriorBorder_All_Interface,
                      const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, commInterface, commDirection )
      {}

      CombinedSpace ( const ThisType& ) = delete;
      ThisType& operator= ( const ThisType& ) = delete;

      //- Additional methods
      //! number of components
      int numComponents () const
      {
        return N;
      }

      //! policy of this space
      DofStoragePolicy myPolicy () const
      {
        return policy;
      }

      //! contained space
      const ContainedDiscreteFunctionSpaceType &containedSpace () const
      {
        return BaseType::containedSpace();
      }
    };



    // specialization of DifferentDiscreteFunctionSpace for this CombinedSapce
    template< class ContainedSpace, int N, DofStoragePolicy policy, class NewFunctionSpace >
    struct DifferentDiscreteFunctionSpace< CombinedSpace< ContainedSpace, N, policy >, NewFunctionSpace >
    {
      typedef CombinedSpace< ContainedSpace, NewFunctionSpace::dimRange, policy > Type;
    };


    // DefaultLocalRestrictProlong ( specialization for CombinedSpace< DiscreteFunctionSpace, N, VariableBased > )
    template< class DiscreteFunctionSpace, int N >
    class DefaultLocalRestrictProlong< CombinedSpace< DiscreteFunctionSpace, N, PointBased > >
      : public DefaultLocalRestrictProlong< typename DiscreteFunctionSpace::template ToNewDimRange< DiscreteFunctionSpace::dimRange *N >::Type >
    {
      typedef DefaultLocalRestrictProlong< CombinedSpace< DiscreteFunctionSpace, N, PointBased > > ThisType;
      typedef DefaultLocalRestrictProlong< typename DiscreteFunctionSpace::template ToNewDimRange< DiscreteFunctionSpace::dimRange *N >::Type > BaseType;

    public:
      DefaultLocalRestrictProlong ( const CombinedSpace< DiscreteFunctionSpace, N, PointBased > &space )
        : BaseType( space )
      {}
    };


    // DefaultLocalRestrictProlong ( specialization for CombinedSpace< DiscreteFunctionSpace, N, VariableBased > )
    template< class DiscreteFunctionSpace, int N >
    class DefaultLocalRestrictProlong< CombinedSpace< DiscreteFunctionSpace, N, VariableBased > >
      : public PowerLocalRestrictProlong< DiscreteFunctionSpace, N >
    {
      typedef DefaultLocalRestrictProlong< CombinedSpace< DiscreteFunctionSpace, N, VariableBased > > ThisType;
      typedef PowerLocalRestrictProlong< DiscreteFunctionSpace, N > BaseType;

    public:
      DefaultLocalRestrictProlong ( const CombinedSpace< DiscreteFunctionSpace, N, VariableBased > &space )
        : BaseType( space.containedSpace() )
      {}
    };

    /** @} **/

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMBINEDSPACE_COMBINEDSPACE_HH
