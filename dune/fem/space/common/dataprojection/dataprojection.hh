#ifndef DUNE_FEM_HPDG_SPACE_COMMON_DATAPROJECTION_HH
#define DUNE_FEM_HPDG_SPACE_COMMON_DATAPROJECTION_HH

#include <cstddef>

#include <vector>

namespace Dune
{

  namespace Fem
  {

    namespace hpDG
    {

      // DataProjection
      // --------------

      /** \brief Abstract definition of the local restriction and prolongation of discrete functions
       *
       *  \tparam  DiscreteFunctionSpace  an adaptive discrete function space
       *  \tparam  Implementation  type of an implementation of this interface
       *
       *  \ingroup DiscreteFunctionSpace_RestrictProlong
       */
      template< class DiscreteFunctionSpace, class Implementation >
      class DataProjection
      {
        using ThisType = DataProjection< DiscreteFunctionSpace, Implementation >;

      public:
        /** \brief discrete function space type */
        using DiscreteFunctionSpaceType = DiscreteFunctionSpace;
        /** \brief basis function set type */
        using BasisFunctionSetType = typename DiscreteFunctionSpaceType::BasisFunctionSetType;
        /** \brief entity type */
        using EntityType = typename BasisFunctionSetType::EntityType;

      protected:
        DataProjection () {}

      public:
        /** \name Move operations
         *  \{
         */

        DataProjection ( DataProjection && ) = default;

        DataProjection &operator= ( DataProjection && ) = default;

        /** \} */

        /** \name Deleted methods
         *  \{
         */

        DataProjection ( const ThisType & ) = delete;

        ThisType &operator= ( const ThisType & ) = delete;

        /** \} */

        /** \name Operation
         *  \{
         */

        /** \brief please doc me
         *
         *  \param[in]  entity  a grid part entity
         *  \param[in]  prior  basis function previously associated to entity
         *  \param[in]  present  basis function associated to entity
         *  \param[in]  origin  blocks previously associated to entity
         *  \param[in]  destination  blocks associated to entity
         */
        void operator() ( const EntityType &entity,
                          const BasisFunctionSetType &prior,
                          const BasisFunctionSetType &present,
                          const std::vector< std::size_t > &origin,
                          const std::vector< std::size_t > &destination )
        {
          asImp()( entity, prior, present, origin, destination );
        }

        /** \brief projection during space adapt
         *
         *  \param[in]  tmp  discrete function for temporary storage
         */
        template <class TemporaryStorage>
        void operator () ( TemporaryStorage& tmp )
        {
          asImp()( tmp );
        }

        /** \brief add discrete function to communicator
         *
         *  \param[in]  comm  communicator
         */
        template< class Communicator >
        void addToList ( Communicator &comm )
        {
          asImp().addToList( comm );
        }

        /** \} */

      protected:
        Implementation &asImp () { return static_cast< Implementation & >( *this ); }

        const Implementation &asImp () const { return static_cast< const Implementation & >( *this ); }
      };

    } // namespace hpDG

    // forward types to Fem namespace for convenience
    using hpDG::DataProjection ;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_HPDG_SPACE_COMMON_DATAPROJECTION_HH
