#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_BASISFUNCTIONSETS_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_BASISFUNCTIONSETS_HH

#include <utility>

#include <dune/common/documentation.hh>

#include <dune/fem/space/basisfunctionset/default.hh>

namespace Dune
{

  namespace Fem
  {

    // BasisFunctionSets
    // -----------------

    /** \class BasisFunctionSets
     *
     *  \brief interface class representing a family of basis function sets
     */
    class BasisFunctionSets
    {
    public:
      /** \brief basis function set */
      typedef ImplementationDefined BasisFunctionSetType;
      /** \brief entity type */
      typedef ImplementationDefined EntityType;

      /** \name Move construction/assignment
       *  \{
       */

      /** \brief move constructor */
      BasisFunctionSets ( BasisFunctionSets && );

      /** \} */

      /** \name Deleted methods
       *  \{
       */

      /** \brief copy constructor */
      BasisFunctionSets ( const BasisFunctionSets & ) = delete;

      /** \brief assignment constructor */
      BasisFunctionSets &operator= ( const BasisFunctionSets & ) = delete;

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \brief return maximum order */
      int order () const;

      /** \brief return order for given grid part entity */
      int order ( const EntityType &entity ) const;

      /** \brief return basis function set for given entity */
      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const;

      /** \} */
    };



    // DefaultBasisFunctionSets
    // ------------------------

    /** \class DefaultBasisFunctionSets
     *
     *  \brief generate a set of default basis function sets from given set of
     *         shape function sets
     *
     *  \tparam  GridPart  grid part type
     *  \tparam  ShapeFunctionSets  shape function sets type
     */
    template< class GridPart, class ShapeFunctionSets >
    class DefaultBasisFunctionSets
    {
      typedef DefaultBasisFunctionSets< GridPart, ShapeFunctionSets > ThisType;

    public:
      /** \brief grid part type */
      typedef GridPart GridPartType;

      /** \brief shape function sets type */
      typedef ShapeFunctionSets ShapeFunctionSetsType;
      /** \brief shape function set type */
      typedef typename ShapeFunctionSetsType::ShapeFunctionSetType ShapeFunctionSetType;

    private:
      static const int dimension = GridPartType::dimension;
      static const int mydimension = ShapeFunctionSetType::FunctionSpaceType::dimDomain;
      static const int codimension = dimension - mydimension;

    public:
      /** \copydoc Dune::Fem::BasisFunctionSets::EntityType */
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;

      /** \copydoc Dune::Fem::BasisFunctionSets::EntityType */
      typedef Dune::Fem::DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetType;

    public:
      /** \name Construction
       *  \{
       */

      explicit DefaultBasisFunctionSets ( ShapeFunctionSetsType &&shapeFunctionSets )
        : shapeFunctionSets_( std::move( shapeFunctionSets ) )
      {}

      template< class... Args, std::enable_if_t< std::is_constructible< ShapeFunctionSetsType, Args &&... >::value, int > = 0 >
      explicit DefaultBasisFunctionSets ( Args &&... args )
        : shapeFunctionSets_( std::forward< Args >( args )... )
      {}

      /** \} */

      /** \name Copying and assignment
       *  \{
       */

      DefaultBasisFunctionSets ( const ThisType & ) = delete;
      DefaultBasisFunctionSets ( ThisType &&other ) = default;

      DefaultBasisFunctionSets &operator= ( const ThisType & ) = delete;
      DefaultBasisFunctionSets &operator= ( ThisType && ) = delete;

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::BasisFunctionSets::order */
      int order () const { return shapeFunctionSets().order(); }

      /** \copydoc Dune::Fem::BasisFunctionSets::order */
      int order ( const EntityType &entity ) const { return shapeFunctionSets().order( entity.type() ); }

      /** \copydoc Dune::Fem::BasisFunctionSets::basisFunctionSet */
      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        return BasisFunctionSetType( entity, shapeFunctionSets().shapeFunctionSet( entity.type() ) );
      }

      /** \} */

      const ShapeFunctionSetsType &shapeFunctionSets () const { return shapeFunctionSets_; }
      ShapeFunctionSetsType &shapeFunctionSets () { return shapeFunctionSets_; }

    private:
      ShapeFunctionSetsType shapeFunctionSets_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_BASISFUNCTIONSETS_HH
