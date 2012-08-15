#ifndef DUNE_FEM_MISC_SELECTVALUE_HH
#define DUNE_FEM_MISC_SELECTVALUE_HH

namespace Dune
{

  namespace Fem
  {

    template< bool first, int V1, int V2 >
    struct SelectIntegerValue
    {
      typedef int Type;
      static const Type Value = V1;
    };

    template< int V1, int V2 >
    struct SelectIntegerValue< false, V1, V2 >
    {
      typedef int Type;
      static const Type Value = V2;
    };


    template< bool first, unsigned int V1, unsigned int V2 >
    struct SelectUnsignedValue
    {
      typedef unsigned int Type;
      static const Type Value = V1;
    };

    template< unsigned int V1, unsigned int V2 >
    struct SelectUnsignedValue< false, V1, V2 >
    {
      typedef unsigned int Type;
      static const Type Value = V2;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MISC_SELECTVALUE_HH
