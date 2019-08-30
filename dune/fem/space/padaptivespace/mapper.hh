#ifndef DUNE_FEM_SPACE_PADAPTIVE_MAPPER_HH
#define DUNE_FEM_SPACE_PADAPTIVE_MAPPER_HH

#include <dune/common/exceptions.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/grid/utility/persistentcontainer.hh>

#include <dune/fem/misc/capabilities.hh>
#include <dune/fem/misc/metaprogramming.hh>
#include <dune/fem/space/common/basesetlocalkeystorage.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/lagrange/lagrangepoints.hh>
#include <dune/fem/space/mapper/codimensionmapper.hh>
#include <dune/fem/space/mapper/dofmapper.hh>
#include <dune/fem/space/mapper/genericadaptivedofmapper.hh>


namespace Dune
{

  namespace Fem
  {

    // Internal forward declaration
    // ----------------------------
    template< class GridPart, int polOrder >
    class PAdaptiveLagrangeMapper;



    // PAdaptiveLagrangeMapperTraits
    // -----------------------------

    template< class GridPart, int polOrder >
    struct PAdaptiveLagrangeMapperTraits
    {
      typedef GridPart GridPartType;

      static const int polynomialOrder = polOrder;
      // if this is set to true the mapper behaves like a DG mapper
      static const bool discontinuousMapper = false;

      typedef typename GridPartType::template Codim< 0 >::EntityType ElementType;
      typedef PAdaptiveLagrangeMapper< GridPartType, polynomialOrder > DofMapperType;

      //! type of the compiled local key
      typedef LagrangePointSet< GridPartType, polynomialOrder >  CompiledLocalKeyType;
      typedef BaseSetLocalKeyStorage< CompiledLocalKeyType > BaseSetLocalKeyStorageType;

      typedef std::vector< BaseSetLocalKeyStorageType > CompiledLocalKeyVectorType ;

      typedef int SizeType ;
      typedef int GlobalKeyType ;
    };



    // Higher Order Lagrange Mapper
    // ----------------------------

    template< class GridPart, int polOrder >
    class PAdaptiveLagrangeMapper
      : public GenericAdaptiveDofMapper< PAdaptiveLagrangeMapperTraits< GridPart, polOrder > >
    {
    public:
      // my traits class
      typedef PAdaptiveLagrangeMapperTraits< GridPart, polOrder > Traits;

    private:
      typedef PAdaptiveLagrangeMapper< GridPart, polOrder > ThisType;
      typedef GenericAdaptiveDofMapper< Traits > BaseType;

    public:
      //! type of the grid part
      typedef typename Traits::GridPartType GridPartType;

      //! type of compiled local keys vector
      typedef typename Traits :: CompiledLocalKeyVectorType  CompiledLocalKeyVectorType;

      //! constructor
      PAdaptiveLagrangeMapper ( const GridPartType &gridPart,
                                const int order,
                                CompiledLocalKeyVectorType &compiledLocalKeys )
        : BaseType( gridPart, order, compiledLocalKeys )
      {}

      //! sort of copy constructor
      PAdaptiveLagrangeMapper ( const ThisType& other,
                                const int order,
                                CompiledLocalKeyVectorType &compiledLocalKeys )
        : BaseType( other, order, compiledLocalKeys )
      {}
    };

    template< class GridPart, int polOrder >
    class PAdaptiveDGMapper;

    template< class GridPart, int polOrder >
    struct PAdaptiveDGMapperTraits
      : public PAdaptiveLagrangeMapperTraits< GridPart, polOrder >
    {
      // this is a mapper for DG
      static const bool discontinuousMapper = true ;

      typedef typename GridPart::template Codim< 0 >::EntityType ElementType;
      typedef PAdaptiveDGMapper< GridPart, polOrder > DofMapperType;
      typedef int SizeType ;
      typedef int GlobalKeyType ;
    };


    // Higher Order Adaptive DG Mapper
    // -------------------------------

    template< class GridPart, int polOrder >
    class PAdaptiveDGMapper
      : public GenericAdaptiveDofMapper< PAdaptiveDGMapperTraits< GridPart, polOrder > >
    {
    public:
      // my traits class
      typedef PAdaptiveDGMapperTraits< GridPart, polOrder > Traits;

    private:
      typedef PAdaptiveDGMapper< GridPart, polOrder > ThisType;
      typedef GenericAdaptiveDofMapper< Traits > BaseType;

    public:
      //! type of the grid part
      typedef typename Traits::GridPartType GridPartType;

      //! type of compiled local keys vector
      typedef typename Traits :: CompiledLocalKeyVectorType  CompiledLocalKeyVectorType;

      //! constructor
      PAdaptiveDGMapper ( const GridPartType &gridPart,
                          const int order,
                          CompiledLocalKeyVectorType &compiledLocalKeys )
        : BaseType( gridPart, order, compiledLocalKeys )
      {}

      //! sort of copy constructor
      PAdaptiveDGMapper ( const ThisType& other,
                          const int order,
                          CompiledLocalKeyVectorType &compiledLocalKeys )
        : BaseType( other, order, compiledLocalKeys )
      {}
    };

    namespace Capabilities
    {
      // isConsecutiveIndexSet
      // ---------------------

      template< class GridPart, int polOrder >
      struct isConsecutiveIndexSet< PAdaptiveDGMapper< GridPart, polOrder > >
      {
        static const bool v = true;
      };

      template< class GridPart, int polOrder >
      struct isConsecutiveIndexSet< PAdaptiveLagrangeMapper< GridPart, polOrder > >
      {
        static const bool v = true;
      };

      template< class GridPart, int polOrder >
      struct isAdaptiveDofMapper< PAdaptiveLagrangeMapper< GridPart, polOrder > >
      {
        static const bool v = true;
      };

      template< class GridPart, int polOrder >
      struct isAdaptiveDofMapper< PAdaptiveDGMapper< GridPart, polOrder > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_PADAPTIVE_MAPPER_HH
