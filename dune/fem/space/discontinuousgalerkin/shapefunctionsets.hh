#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SHAPEFUNCTIONSETS_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SHAPEFUNCTIONSETS_HH

#include <algorithm>
#include <utility>
#include <vector>

#include <dune/common/documentation.hh>

#include <dune/geometry/type.hh>

#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/space/common/basesetlocalkeystorage.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>
#include <dune/fem/storage/singletonlist.hh>

namespace Dune
{

  namespace Fem
  {

    // ShapeFunctionSets
    // -----------------

    /** \class ShapeFunctionSets
     *
     *  \brief interface class representing a family of shape function sets
     */
    class ShapeFunctionSets
    {
    public:
      /** \brief shape function set type */
      typedef ImplementationDefined ShapeFunctionSetType;

      /** \name Move construction/assignemnt
       *  \{
       */

      /** \brief move constructor */
      ShapeFunctionSets ( ShapeFunctionSets && );

      /** \} */

      /** \name Deleted methods
       *  \{
       */

      /** \brief copy constructor */
      ShapeFunctionSets ( const ShapeFunctionSets & ) = delete;

      /** \brief assignment operator */
      ShapeFunctionSets &operator= ( const ShapeFunctionSets & ) = delete;

      /** \} */

      /** \name Public member functions
       *  \{
       */
      /** \brief return vector of geometry types */
      const std::vector< Dune::GeometryType > &types () const;

      /** \brief return maximum order */
      int order () const;

      /** \brief return order per geometry type */
      int order ( Dune::GeometryType type ) const;

      /** \brief return shape function set
       *
       *  \param[in]  type  geometry type
       *
       *  \returns shape function set
       */
      ShapeFunctionSetType shapeFunctionSet ( const Dune::GeometryType &type ) const;

      /** \} */
    };



    // CachedShapeFunctionSets
    // -----------------------

    template< class GridPart, class ShapeFunctionSet,
              class Factory = Dune::Fem::DefaultSingletonFactory< Dune::GeometryType, ShapeFunctionSet > >
    class CachedShapeFunctionSets
    {
      typedef CachedShapeFunctionSets< GridPart, ShapeFunctionSet, Factory > ThisType;

    public:
      /** \brief grid part type */
      typedef GridPart GridPartType;
      /** \brief shape function set type */
      typedef Dune::Fem::ShapeFunctionSetProxy< ShapeFunctionSet > ShapeFunctionSetType;

    private:
      static const int dimension = GridPartType::dimension;
      static const int mydimension = ShapeFunctionSet::FunctionSpaceType::dimDomain;
      static const int codimension = dimension - mydimension;

      typedef Dune::Fem::SingletonList< Dune::GeometryType, ShapeFunctionSet, Factory > SingletonProviderType;
      typedef Dune::Fem::BaseSetLocalKeyStorage< ShapeFunctionSet > ShapeFunctionSetStorageType;

    public:
      /** \name Construction
       *  \{
       */

      explicit CachedShapeFunctionSets ( const GridPartType &gridPart )
        : types_( types( gridPart ) )
      {
        typedef typename std::vector< Dune::GeometryType >::const_iterator const_iterator;
        const const_iterator end = types_.end();
        for( const_iterator it = types_.begin(); it != end; ++it )
        {
          const Dune::GeometryType type = *it;
          shapeFunctionSets_.template insert< SingletonProviderType >( type );
        }
      }

      /** \} */

      /** \name Copying and assignment
       *  \{
       */

      CachedShapeFunctionSets ( const ThisType & ) = delete;

      CachedShapeFunctionSets ( ThisType &&other )
        : types_( std::move( other.types_ ) ),
          shapeFunctionSets_( std::move( other.shapeFunctionSets_ ) )
      {}

      CachedShapeFunctionSets &operator= ( const ThisType & ) = delete;

      /** \} */

      /** \name Public member functions
       *  \{
       */

      /** \copydoc Dune::Fem::ShapeFunctionSets::types */
      const std::vector< Dune::GeometryType > &types () const { return types_; }

      /** \copydoc Dune::Fem::ShapeFunctionSets::order */
      int order () const
      {
        int order = 0;

        typedef typename std::vector< Dune::GeometryType >::const_iterator const_iterator;
        const const_iterator end = types_.end();
        for( const_iterator it = types_.begin(); it != end; ++it )
        {
          const Dune::GeometryType type = *it;
          order = std::max( this->order( type ), order );
        }

        return order;
      }

      /** \copydoc Dune::Fem::ShapeFunctionSets::order */
      int order ( Dune::GeometryType type ) const
      {
        return shapeFunctionSet( type ).order();
      }

      /** \copydoc Dune::Fem::ShapeFunctionSets::shapeFunctionSet */
      ShapeFunctionSetType shapeFunctionSet ( const Dune::GeometryType &type ) const
      {
        return ShapeFunctionSetType( &shapeFunctionSets_[ type ] );
      }

      /** \} */

    private:
      static std::vector< Dune::GeometryType > types ( const GridPartType &gridPart )
      {
        typedef typename GridPartType::GridType GridType;
        typedef typename GridPartType::IndexSetType IndexSetType;
        return Dune::Fem::AllGeomTypes< IndexSetType, GridType >( gridPart.indexSet() ).geomTypes( codimension );
      }

      std::vector< Dune::GeometryType > types_;
      ShapeFunctionSetStorageType shapeFunctionSets_;
    };



    // SelectCachingShapeFunctionSets
    // ------------------------------

    template< class GridPart, class ShapeFunctionSet, class Storage >
    class SelectCachingShapeFunctionSets
    {
      typedef SelectCachingShapeFunctionSets< GridPart, ShapeFunctionSet, Storage > ThisType;

      typedef SelectCachingShapeFunctionSet< ShapeFunctionSet, Storage > CachedShapeFunctionSetType;

      struct Factory
      {
        static CachedShapeFunctionSetType *createObject ( const Dune::GeometryType &type )
        {
          typedef typename CachedShapeFunctionSetType::ImplementationType ImplementationType;
          return new CachedShapeFunctionSetType( type, ImplementationType( type ) );
        }

        static void deleteObject ( CachedShapeFunctionSetType *object ) { delete object; }
      };

      typedef CachedShapeFunctionSets< GridPart, CachedShapeFunctionSetType, Factory > Implementation;

    public:
      /** \copydoc Dune::Fem::ShapeFunctionSets::ShapeFunctionSetType */
      typedef typename Implementation::ShapeFunctionSetType ShapeFunctionSetType;

      static constexpr bool codegenShapeFunctionSet = detail::IsCodegenShapeFunctionSet< CachedShapeFunctionSetType >::value;

      /** \name Construction
       *  \{
       */

      explicit SelectCachingShapeFunctionSets ( const GridPart &gridPart )
        : impl_( gridPart )
      {}

      /** \} */

      /** \name Copying and assignment
       *  \{
       */

      SelectCachingShapeFunctionSets ( const ThisType & ) = delete;

      SelectCachingShapeFunctionSets ( ThisType &&other )
        : impl_( std::move( other.impl_ ) )
      {}

      SelectCachingShapeFunctionSets &operator= ( const ThisType & ) = delete;

      /** \} */

      /** \copydoc Dune::Fem::ShapeFunctionSets::types */
      const std::vector< Dune::GeometryType > &types () const { return impl_.types(); }

      /** \copydoc Dune::Fem::ShapeFunctionSets::order */
      int order () const { return impl_.order(); }

      /** \copydoc Dune::Fem::ShapeFunctionSets::order */
      int order ( Dune::GeometryType type ) const { return impl_.order( type ); }

      /** \copydoc Dune::Fem::ShapeFunctionSets::shapeFunctionSet */
      ShapeFunctionSetType shapeFunctionSet ( const Dune::GeometryType &type ) const
      {
        return impl_.shapeFunctionSet( type );
      }

    private:
      Implementation impl_;
    };



    // VectorialShapeFunctionSets
    // --------------------------

    template< class Implementation, class Range >
    class VectorialShapeFunctionSets
    {
      typedef VectorialShapeFunctionSets< Implementation, Range > ThisType;

    public:
      static constexpr bool codegenShapeFunctionSet = detail::IsCodegenShapeFunctionSet< Implementation >::value;

      /** \brief shape function set type */
      typedef VectorialShapeFunctionSet< typename Implementation::ShapeFunctionSetType, Range > ShapeFunctionSetType;

      /** \name Construction
       *  \{
       */

      explicit VectorialShapeFunctionSets ( Implementation &&impl )
        : impl_( impl )
      {}

      template< class... Args >
      explicit VectorialShapeFunctionSets ( Args &&...args )
        : impl_( std::forward< Args >( args )... )
      {}

      /** \} */

      /** \name Copying and assignment
       *  \{
       */

      VectorialShapeFunctionSets ( const ThisType & ) = delete;

      VectorialShapeFunctionSets ( ThisType & ) = delete;

      VectorialShapeFunctionSets ( ThisType &&other )
        : impl_( std::move( other.impl_ ) )
      {}

      VectorialShapeFunctionSets &operator= ( const ThisType & ) = delete;

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::ShapeFunctionSets::types */
      const std::vector< Dune::GeometryType > &types () const { return impl().types(); }

      /** \copydoc Dune::Fem::ShapeFunctionSets::order */
      int order () const { return impl().order(); }

      /** \copydoc Dune::Fem::ShapeFunctionSets::order */
      int order ( Dune::GeometryType type ) const { return impl().order( type ); }

      /** \copydoc Dune::Fem::ShapeFunctionSets::shapeFunctionSet */
      ShapeFunctionSetType shapeFunctionSet ( const Dune::GeometryType &type ) const
      {
        return ShapeFunctionSetType( impl().shapeFunctionSet( type ) );
      }

      /** \} */

    private:
      const Implementation &impl () const { return impl_; }

      Implementation impl_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_SHAPEFUNCTIONSETS_HH
