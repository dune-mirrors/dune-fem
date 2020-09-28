#ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_SELECTCACHING_HH
#define DUNE_FEM_SPACE_SHAPEFUNCTIONSET_SELECTCACHING_HH

// dune-fem includes
#include <dune/fem/space/shapefunctionset/caching.hh>

namespace Dune
{

  namespace Fem
  {

    // External Forward Declarations
    // -----------------------------

    class CachingStorage {};

    class SimpleStorage {};

    class CodegenStorage {};

    // SelectCachingShapeFunctionSet
    // -----------------------------

    template< class ShapeFunctionSet, class Storage >
    class SelectCachingShapeFunctionSet;

    template< class ShapeFunctionSet >
    class SelectCachingShapeFunctionSet< ShapeFunctionSet, CachingStorage >
    : public CachingShapeFunctionSet< ShapeFunctionSet >
    {
      typedef CachingShapeFunctionSet< ShapeFunctionSet > BaseType;

    public:
      typedef ShapeFunctionSet ImplementationType;

      explicit SelectCachingShapeFunctionSet ( const GeometryType &type,
                                               const ShapeFunctionSet &shapeFunctionSet = ShapeFunctionSet() )
      : BaseType( type, shapeFunctionSet )
      {}
    };

    template< class ShapeFunctionSet >
    class SelectCachingShapeFunctionSet< ShapeFunctionSet, SimpleStorage >
    : public ShapeFunctionSet
    {
      typedef ShapeFunctionSet BaseType;

    public:
      typedef ShapeFunctionSet ImplementationType;

      explicit SelectCachingShapeFunctionSet ( const GeometryType &type,
                                               const ShapeFunctionSet &shapeFunctionSet = ShapeFunctionSet() )
      : BaseType( shapeFunctionSet )
      {}
    };

    template< class ShapeFunctionSet >
    class SelectCachingShapeFunctionSet< ShapeFunctionSet, CodegenStorage >
    : public CachingShapeFunctionSet< ShapeFunctionSet >
    {
      typedef CachingShapeFunctionSet< ShapeFunctionSet > BaseType;

    public:
      typedef ShapeFunctionSet ImplementationType;

      //! this indicates that generated codes for evaluate and axpy is available
      static constexpr bool codegenShapeFunctionSet = true ;

      explicit SelectCachingShapeFunctionSet ( const GeometryType &type,
                                               const ShapeFunctionSet &shapeFunctionSet = ShapeFunctionSet() )
      : BaseType( type, shapeFunctionSet )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_SELECTCACHING_HH
