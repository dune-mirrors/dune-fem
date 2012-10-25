#ifndef DUNE_FEM_SPACE_LAGRANGE_SPACE_HH
#define DUNE_FEM_SPACE_LAGRANGE_SPACE_HH

// dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/nullptr.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/common/basesetlocalkeystorage.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/lagrangespace/lagrangepoints.hh>

// local includes
#include "adaptmanager.hh"


namespace Dune
{

  namespace Fem
  {

    // LagrangeDiscreteFunctionSpaceTraits
    // -----------------------------------

    template< class FunctionSpace, class GridPart, unsigned int polOrder, template< class > class Storage = CachingStorage >
    struct LagrangeDiscreteFunctionSpaceTraits;



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

    // LagrangeDiscreteFunctionSpace
    // -----------------------------

    /** \class   LagrangeDiscreteFunctionSpace
     *  \ingroup LagrangeDiscreteFunctionSpace
     *  \brief   Lagrange discrete function space
     */

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    class LagrangeDiscreteFunctionSpace
    : public DiscreteFunctionSpaceDefault< LagrangeDiscreteFunctionSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
    {
      typedef LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType;
      typedef DiscreteFunctionSpaceDefault< LagrangeDiscreteFunctionSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > > BaseType;

    public:
      static const int polynomialOrder = polOrder;

      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;
      typedef typename BaseType::IntersectionType IntersectionType;
      typedef typename BaseType::MapperType MapperType;
      typedef typename BaseType::BlockMapperType BlockMapperType;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::EntityType EntityType;

      typedef LagrangePointSet< GridPartType, polynomialOrder > LagrangePointSetType;

    private:
      typedef CompiledLocalKeyContainer< LagrangePointSetType, polynomialOrder, polynomialOrder > LagrangePointSetContainerType;

    public:
      ///////////////////////
      // Interface methods //
      ///////////////////////

      using BaseType::order;

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::type */
      DFSpaceIdentifier type () const
      {
        return LagrangeSpace_id;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::basisFunctionSet */
      const BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        DUNE_THROW( NotImplemented, "Method basisFunctionSet() not implemented yet." );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous () const
      {
        return true;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous (const IntersectionType &intersection) const
      {
        return intersection.conforming();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order () const
      {
        return polOrder;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::mapper */
      MapperType &mapper () const
      {
        assert( mapper_ );
        return *mapper_;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::blockMapper */
      BlockMapperType &blockMapper () const
      {
        assert( blockMapper_ );
        return *blockMapper_;
      }

      ///////////////////////////
      // Non-interface methods //
      ///////////////////////////
      
      /** \brief provide access to the base function set for a geometry type
       *
       *  \param[in]  type  type of geometry the base function set is requested for
       *
       *  \returns base function set for the specified geometry
       */
      const BasisFunctionSetType basisFunctionSet ( const GeometryType type ) const
      {
        DUNE_THROW( NotImplemented, "Method basisFunctionSet() not implemented yet." );
      }

      /** \brief provide access to the Lagrange point set for an entity
       *
       *  \note This method is not part of the DiscreteFunctionSpaceInterface. It
       *        is unique to the LagrangeDiscreteFunctionSpace.
       *
       *  \param[in]  entity  entity the Lagrange point set is requested for
       *  
       *  \returns LagrangePointSet
       */
      template< class EntityType >
      const LagrangePointSetType &lagrangePointSet ( const EntityType &entity ) const
      {
        return lagrangePointSet( entity.type() );
      }

      /** \brief provide access to the Lagrange point set for a geometry type
       *
       *  \note This method is not part of the DiscreteFunctionSpaceInterface. It
       *        is unique to the LagrangeDiscreteFunctionSpace.
       *
       *  \param[in]  type  type of geometry the Lagrange point set is requested for
       *
       *  \returns LagrangePointSetType
       */
      const LagrangePointSetType &lagrangePointSet ( const GeometryType type ) const
      {
        return lagrangePointSetContainer_.compiledLocalKey( type, polynomialOrder );
      }

    public:
      LagrangePointSetContainerType lagrangePointSetContainer_;
      MapperType *mapper_;
      BlockMapperType *blockMapper_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_SPACE_HH
