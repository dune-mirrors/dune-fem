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
    struct HE { typedef Impl type; };

    template <class Impl>
    struct HE< Impl, true >  {typedef typename Impl::HostEntityType type; };

  public:
    typedef typename BaseType::Implementation  ImplementationType;
    static constexpr bool hasHostEntity = checkHostEntity< ImplementationType >::value;
    typedef typename HE< ImplementationType, hasHostEntity >::type  HostEntityType;

    /** \brief cast to HostEntityType if such type exists otherwise return implementation */
    operator const HostEntityType& () const
    {
      if constexpr ( hasHostEntity )
        return this->impl().hostEntity();
      else
        return this->impl();
    }
  };

} // end namespace Dune
#endif
