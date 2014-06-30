/**@file
 *
 * Define a default RestrictProlongTuple by chaining sufficiently many
 * Dune::ACFem::RestrictProlongPair templates together.
 */
#ifndef __DUNE_FEM_RESTRICT_PROLONG_TUPLE_HH__
#define __DUNE_FEM_RESTRICT_PROLONG_TUPLE_HH__

#include "restrictprolonginterface.hh"
#include "../../function/common/function.hh"
#include <dune/common/tuples.hh>

namespace Dune {

  namespace Fem {

    /**@addtogroup RestrictProlongInterface
     * @{
     */    

    /**Chain two RestrictProlong implementations together. Note that
     * RP1 and RP2 must be copy-constructible. This
     * RestrictProlongPair forms the base for forming RestrictProlong
     * tuples which chain the RestrictProlongDefault implementation of
     * arbitrary many DiscreteFunction-types together.
     */
    template<class RP1, class RP2>
    class RestrictProlongPair
      : public RestrictProlongInterface<
      RestrictProlongTraits<RestrictProlongPair<RP1, RP2>,
                                       typename RP1::DomainFieldType> >,
        public std::pair<RP1, RP2>
    {
      typedef std::pair<RP1, RP2> StorageType;
      typedef RestrictProlongInterface<
        RestrictProlongTraits<RestrictProlongPair<RP1, RP2>,
                                         typename RP1::DomainFieldType> >
      InterfaceType;

      using StorageType::first;
      using StorageType::second;

     public:
      typedef typename RP1::DomainFieldType DomainFieldType; // do we really need to check that ...

      RestrictProlongPair(const RP1& rp1, const RP2& rp2)
        : StorageType(rp1, rp2)
      {}

      void setFatherChildWeight(const DomainFieldType &weight) const {
        first.setFatherChildWeight(weight);
        second.setFatherChildWeight(weight);
      }

      template<class Entity>
      void restrictLocal(const Entity &father, const Entity &son, bool initialize) const {
        first.restrictLocal(father, son, initialize);
        second.restrictLocal(father, son, initialize);
      }
  
      template<class Entity, class LocalGeometry>
      void restrictLocal(const Entity &father, const Entity &son,
                         const LocalGeometry &geometryInFather,
                         bool initialize) const {
        first.restrictLocal(father, son, geometryInFather, initialize);
        second.restrictLocal(father, son, geometryInFather, initialize);
      }

      template<class Entity>
      void prolongLocal(const Entity &father, const Entity &son, bool initialize) const {
        first.prolongLocal(father, son, initialize);
        second.prolongLocal(father, son, initialize);
      }

      template<class Entity, class LocalGeometry>
      void prolongLocal(const Entity &father, const Entity &son,
                        const LocalGeometry &geometryInFather,
                        bool initialize) const {
        first.prolongLocal(father, son, geometryInFather, initialize);
        second.prolongLocal(father, son, geometryInFather, initialize);
      }

      template<class Communicator>
      void addToList(Communicator &comm) {
        first.addToList(comm);
        second.addToList(comm);
      }

      template<class LoadBalancer>
      void addToLoadBalancer(LoadBalancer &lb) {
        first.addToLoadBalancer(lb);
        second.addToLoadBalancer(lb);
      }
    };

    /**A helper class which defines the proper RestrictProlong
     * compound data-type for tuples of DiscreteFunction
     * implementations. This works by template recursion. This 
     */
    template<class... All>
    struct RestrictProlongDefaultTraits;

    //!@copydoc RestrictProlongDefaultTraits
    template<class DiscreteFunction, class... Rest>
    struct RestrictProlongDefaultTraits< DiscreteFunction, Rest... >
    {
      typedef DiscreteFunction FirstDiscreteFunctionType;
      typedef RestrictProlongDefault<FirstDiscreteFunctionType> FirstRestrictProlongType;
      typedef RestrictProlongPair<FirstRestrictProlongType, typename RestrictProlongDefaultTraits<Rest...>::Type> Type;
    };

    //!@copydoc RestrictProlongDefaultTraits
    template<class DiscreteFunction>
    struct RestrictProlongDefaultTraits<DiscreteFunction>
    {
      typedef DiscreteFunction FirstDiscreteFunctionType;
      typedef RestrictProlongDefault<FirstDiscreteFunctionType> FirstRestrictProlongType;
      typedef FirstRestrictProlongType Type;
    };

    /**Endpoint for the makeRestrictProlongDefault recursion. Define a
     * version for a single argument.
     */
    template<class DF>
    static inline
    RestrictProlongDefault<DF>
    makeRestrictProlongDefault(Function<typename DF::FunctionSpaceType, DF>& df_)
    {
      DF& df(static_cast<DF&>(df_));

      return RestrictProlongDefault<DF>(df);
    }

    /**Take arbitrary many discrete functions of potentially different
     * type and generate a suitable RestrictProlong implementation for
     * use with the AdaptationManager.
     */
    template<class DF1, class DF2, class... Rest>
    static inline
    typename RestrictProlongDefaultTraits<DF1, DF2, Rest...>::Type
    makeRestrictProlongDefault(DF1& df1, DF2& df2, Rest&... rest)
    {
      typedef typename RestrictProlongDefaultTraits<DF1, DF2, Rest...>::Type ResultType;

      return ResultType(makeRestrictProlongDefault(df1), makeRestrictProlongDefault(df2, rest...));
    }

    /**AFAIK, C++11 is not capable of expanding tuples into parameter
     * packs, so this recursive helper function is used to unpack a
     * tuple and perform the necessary constructions in order to
     * finally have a compound RestrictProlong type. This involved
     * O(N*N) copy constructions. maybe std::forward should be used
     * here ...
     */
    template<class Tuple, class Index = std::integral_constant<size_t, 0> >
    struct RestrictProlongDefaultTupleHelper
    {
      enum { index = Index::value };
      typedef std::integral_constant<size_t, index+1> NextIndexType;
      typedef typename remove_reference<typename tuple_element<index, Tuple>::type>::type DiscreteFunctionType;
      typedef RestrictProlongDefault<DiscreteFunctionType> RestrictProlongType;
      typedef RestrictProlongDefaultTupleHelper<Tuple, NextIndexType> NextHelperType;
      typedef RestrictProlongPair<RestrictProlongType, typename NextHelperType::Type> Type;

      static Type construct(const Tuple& arg)
      {
        return Type(RestrictProlongType(get<index>(arg)), NextHelperType::construct(arg));
      }
    };

    /**Recursion end-point to access the last argument of the tuple. */
    template<class Tuple>
    struct RestrictProlongDefaultTupleHelper<Tuple, std::integral_constant<size_t, std::tuple_size<Tuple>::value-1> >
    {
      enum { index = tuple_size<Tuple>::value - 1 };
      typedef typename remove_reference<typename tuple_element<index, Tuple>::type>::type DiscreteFunctionType;
      typedef RestrictProlongDefault<DiscreteFunctionType> RestrictProlongType;
      typedef RestrictProlongType Type;
  
      static Type construct(const Tuple& arg)
      {
        return Type(get<index>(arg));
      }
    };

    template<class... Args>
    struct RestrictProlongDefaultTraits<tuple<Args...> >
    {
      typedef
      typename RestrictProlongDefaultTupleHelper<tuple<Args...> >::Type
      Type;
    };

    /**Conveniently form a compound RestrictProlong-implementation
     * which conforms to RestrictProlongInterface. This
     * functions simply takes a std::tuple of references to
     * DiscreteFunction instances and glues them together. The result
     * is a RestrictProlong instance which can be added to a
     * AdaptationManager.
     *
     * @bug Rather a remark: C++14 will add index-sequences which can
     * be used to unpack tuples into template parameter packs.
     */
    template<class... Rest>
    static inline
    typename RestrictProlongDefaultTraits<Rest...>::Type
    makeRestrictProlongDefault(const tuple<Rest&...>& arg)
    {
      typedef tuple<Rest&...> TupleType;
      return RestrictProlongDefaultTupleHelper<TupleType>::construct(arg);
    }

    //!@} RestrictProlongInterface

  } // Fem

} // Dune

#endif // __DUNE_FEM_RESTRICT_PROLONG_TUPLE_HH__

