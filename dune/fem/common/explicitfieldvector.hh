#ifndef DUNE_FEM_COMMON_EXPLICITFIELDVECTOR_HH
#define DUNE_FEM_COMMON_EXPLICITFIELDVECTOR_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

namespace Dune
{

  namespace Fem
  {
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

    template<class T, int N>
    class ExplicitFieldVector
      : public Dune::FieldVector<T, N>
    {
      typedef Dune::FieldVector<T, N> BaseType;
     public:
      using BaseType::FieldVector;
      using BaseType::operator=;
      operator const BaseType&() const { return static_cast<const BaseType>(*this); }
      operator BaseType&() { return static_cast<BaseType>(*this); }
    };

  } // Fem

} // Dune


#endif // DUNE_FEM_COMMON_EXPLICITFIELDVECTOR_HH
