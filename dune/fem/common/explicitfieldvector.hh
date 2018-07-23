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

  template<typename T, int N, int M>
  struct IsFieldVectorSizeCorrect<Fem::ExplicitFieldVector<T, N>, M>
    : IsFieldVectorSizeCorrect<FieldVector<T, N>, M>
  {};

  namespace Fem {

    /**std::true_type for containers containing field elements.*/
    template<class T>
    struct IsFieldType
      : std::is_same<typename FieldTraits<T>::field_type, T>
    {};

    /**Accept implicit type conversion from any DenseVector to a field
     * vector only if both vectors contain field elements of some type
     * or both vectors do NOT contain field elements of some kind and
     * the elements of the DenseVector are convertible to the field
     * vector's elements.
     *
     * This inhibits the initialization of a FieldVector of complex
     * objects like matrices from a field vector of scalars.
     */
    template<class C, class T, class SFINAE = void>
    struct AcceptElementImplicitConstruction
      : IsFieldType<C>
    {};

    template<class C, class T>
    struct AcceptElementImplicitConstruction<
      C, T,
      std::enable_if_t<((IsFieldType<typename DenseMatVecTraits<C>::value_type>::value
                         ==
                         IsFieldType<T>::value)
        )> >
      : std::true_type
    {};

    template<class T, int N>
    class ExplicitFieldVector
      : public Dune::FieldVector<T, N>
    {
      typedef ExplicitFieldVector< T, N > ThisType;
      typedef Dune::FieldVector<T, N> BaseType;
     public:
      //! Constructor making default-initialized vector
      constexpr ExplicitFieldVector()
        : BaseType()
      {}

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

      /**Allow implicit conversion if both vectors are either
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

      //! Assignment operator for scalar
      template<typename C,
               std::enable_if_t<(
                 N == 1 &&
                 AcceptElementImplicitConstruction<C, T>::value &&
                 std::is_assignable<T, C>::value &&
                 ! std::is_base_of<DenseVector<typename FieldTraits<T>::field_type>, T
                                   >::value
                 ), int> = 0
               >
      ExplicitFieldVector& operator=(const C& c)
      {
        (*this)[0] = c;
        return *this;
      }

      //! Inherit assignment
      // using BaseType::operator=; <- give ambiguous overloads
      using DenseVector<FieldVector<T, N> >::operator=;

      //! copy assignment operator
      ExplicitFieldVector& operator=(const ExplicitFieldVector& other)
      {
        static_cast<BaseType&>(*this) = static_cast<const BaseType&>(other);
        return *this;
      }

      template <typename C, std::enable_if_t<std::is_assignable<T, C>::value, int> = 0>
      ExplicitFieldVector& operator=(const FieldVector<C, N>& other)
      {
        static_cast<BaseType&>(*this) = other;
        return *this;
      }

      template <typename C, std::enable_if_t<std::is_assignable<T, C>::value, int> = 0>
      ExplicitFieldVector& operator=(const ExplicitFieldVector<C, N>& other)
      {
        static_cast<BaseType&>(*this) = other;
        return *this;
      }

    };

    template<class FV>
    struct MakeExplicit
    {
      using Type = FV;
    };

    template<class Field, int Size>
    struct MakeExplicit<FieldVector<Field, Size> >
    {
      using Type = ExplicitFieldVector<Field, Size>;
    };

    template<class FV>
    using Explicit = typename MakeExplicit<FV>::Type;

  } // Fem

} // Dune


#endif // DUNE_FEM_COMMON_EXPLICITFIELDVECTOR_HH
