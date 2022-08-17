// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FEM_GRIDPART_COMMON_ENTITY_HH
#define DUNE_FEM_GRIDPART_COMMON_ENTITY_HH

#include <dune/grid/common/entity.hh>

namespace Dune
{

  /**
   @brief Wrapper class for wrapped entities added a cast operator to the host
   entity.

   \tparam cd Codimension of the entity
   \tparam dim Dimension of the grid
   \tparam GridImp Type that is a model of Dune::Grid
   \tparam EntityImp Class template that is a model of Dune::Entity

  */
  template<int cd, int dim, class GridImp, template<int,int,class> class EntityImp>
  class ExtendedEntity : public Dune::Entity< cd, dim, GridImp, EntityImp >
  {
    typedef Dune::Entity< cd, dim, GridImp, EntityImp > BaseType;
  public:
    using BaseType::BaseType;

  protected:
    template<typename T>
    struct ToVoid
    {
      typedef void type;
    };

    template <typename T, typename dummy = void>
    struct checkHostEntity : std::false_type {};

    template <typename T>
    struct checkHostEntity<T, typename ToVoid<typename T::HostEntityType>::type > : std::true_type{};


    template <class Impl, bool>
    struct HE {
      typedef Impl      HostEntity;
      typedef BaseType  GridEntity;
    };

    template <class Impl>
    struct HE< Impl, true >
    {
    private:
      // type of grid entity
      typedef typename Impl::HostGridPartType::GridType::template Codim<cd>::Entity __GEType;
    public:
      typedef typename Impl::HostEntityType HostEntity;
      typedef typename
        std::conditional< std::is_same<HostEntity, __GEType>::value,
                          BaseType, __GEType> :: type GridEntity;
    };

  public:
    typedef typename BaseType::Implementation  ImplementationType;
    static constexpr bool hasHostEntity = checkHostEntity< ImplementationType >::value;
    typedef typename HE< ImplementationType, hasHostEntity >::HostEntity   HostEntityType;
    typedef typename HE< ImplementationType, hasHostEntity >::GridEntity  GridEntityType;

    /** \brief cast to HostEntityType if such type exists otherwise return implementation */
    operator const HostEntityType& () const
    {
      if constexpr ( hasHostEntity )
        return this->impl().hostEntity();
      else
        return this->impl();
    }

    /** \brief cast to HostEntityType if such type exists otherwise return implementation */
    operator const GridEntityType& () const
    {
      if constexpr ( std::is_same< BaseType, GridEntityType > :: value )
        return *this;
      else
      {
        return gridEntity(*this);
      }
    }
  };

} // end namespace Dune
#endif
