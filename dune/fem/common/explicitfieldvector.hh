#ifndef DUNE_FEM_COMMON_EXPLICITFIELDVECTOR_HH
#define DUNE_FEM_COMMON_EXPLICITFIELDVECTOR_HH

#include <type_traits>
#include <utility>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/typeutilities.hh>

namespace Dune
{

  namespace Fem
  {
    /**A variant of FieldVector which does not allow for implicit
     * type-conversion from an to everything.
     */
    template<class T, int N>
    class ExplicitFieldVector;
  }

  template<class T, int N>
  struct DenseMatVecTraits< Fem::ExplicitFieldVector<T, N> >
    : DenseMatVecTraits<FieldVector<T, N> >
  {
    typedef Fem::ExplicitFieldVector<T, N> derived_type;
  };

  template<class T, int N>
  struct FieldTraits<Fem::ExplicitFieldVector<T, N> >
    : FieldTraits<FieldVector<T, N> >
  {};

  namespace Fem {

    /**Accept implicit type conversion from any DenseVector to a field
     * vector only if both vectors contain field elements of some type
     * or both vectors do NOT contain field elements of some kind and
     * the elements of the DenseVector are convertible to the field
     * vector's elements.
     *
     * This inhibits the initialization of a FieldVector of complex
     * objects like matrices from a field vector of scalars.
     */
    template<class C, class T>
    struct AcceptElementImplicitConstruction
    {
      static constexpr bool value =
        (std::is_same<typename FieldTraits<typename DenseMatVecTraits<C>::value_type>::field_type,
         typename DenseMatVecTraits<C>::value_type
         >::value
         ==
         std::is_same<typename FieldTraits<T>::field_type, T>::value)
        &&
        std::is_convertible<typename DenseMatVecTraits<C>::value_type, T>::value;
    };

    template<class T, int N>
    class ExplicitFieldVector
      : public Dune::FieldVector<T, N>
    {
      typedef ExplicitFieldVector< T, N > ThisType;
      typedef Dune::FieldVector<T, N> BaseType;
     public:
      //! Inherit assignment
      using BaseType::operator=;

      /**Redirect any general construction to the base class during
       * explicit conversion
       */
      template< class... Args, disableCopyMove< ThisType, Args... > = 0, std::enable_if_t< std::is_constructible< BaseType, Args &&... >::value, int > = 0 >
      explicit ExplicitFieldVector ( Args &&... args )
        : BaseType( std::forward< Args >( args )... )
      {}

      ExplicitFieldVector ( const std::initializer_list< T > &values )
        : BaseType( values )
      {}

      /**Allow implicit conversion if bothe vector are either
       * composed of field-elements of some fields which can be
       * converted into each other or if both vectors are composed of
       * more complicated elements (which can be converted into each
       * other), but do not allow implicit conversion of a FieldVector
       * of scalars into a FieldVector composed of more complicated
       * stuff. In particalar, FunctionSpace::RangeType cannot be
       * implicitly converted to FunctionSpace::HessianRangeType.
       */
      template<class C>
      ExplicitFieldVector(const DenseVector<C>& x,
                          typename std::enable_if<(
                            IsFieldVectorSizeCorrect<C, N>::value
                            &&
                            AcceptElementImplicitConstruction<C, T>::value)
                          >::type* dummy=0 )
        : BaseType(x)
      {}

    };

  } // Fem

} // Dune


#endif // DUNE_FEM_COMMON_EXPLICITFIELDVECTOR_HH
