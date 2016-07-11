#ifndef DUNE_FEM_COMMON_SIMPLETENSOR_HH
#define DUNE_FEM_COMMON_SIMPLETENSOR_HH

#include <dune/common/fmatrix.hh>

namespace Dune
{

  namespace Fem
  {

    template< class Field, class Arg >
    Arg scale ( const Field & a, const Arg & arg ) { return a * arg; }

    template< class Field, int i, int j >
    FieldMatrix< Field, i, j > scale ( const Field &a, const FieldMatrix< Field, i, j > & m )
    {
      FieldMatrix< Field, i, j > ret( m );
      return ret *= a;
    }


    // IdTensor
    // --------

    template< class Arg >
    struct IdTensor
    {
      Arg operator() ( const Arg &arg ) const { return arg; }
    };


    // ScaledIdTensor
    // --------------

    template< class Arg, class Tensor = IdTensor< Arg > >
    struct ScaledTensor
    {
      typedef typename FieldTraits< Arg >::field_type FieldType;
      ScaledTensor ( FieldType alpha, const Tensor &tensor = Tensor() ) : alpha_( alpha ), tensor_( tensor ) {}
      ScaledTensor ( const ScaledTensor & ) = default;
      ScaledTensor ( ScaledTensor && ) = default;

      Arg operator() ( const Arg &arg ) const { return scale( alpha_, tensor_( arg ) ); }

    protected:
      FieldType alpha_;
      Tensor tensor_;
    };



    // scaledTensor
    // ------------
    template< class Arg, class Tensor >
    ScaledTensor< Arg, Tensor > scaledTensor ( typename FieldTraits< Arg >::field_type alpha, Tensor tensor )
    {
      return ScaledTensor< Arg, Tensor >( alpha, tensor );
    }

    // DiagonalTensor
    // --------------
    template< class Arg >
    struct DiagonalTensor
    {
      DiagonalTensor ( const Arg &diag ) : diag_( diag ) {}
      DiagonalTensor ( const DiagonalTensor & ) = default;
      DiagonalTensor ( DiagonalTensor && ) = default;

      Arg operator() ( const Arg &arg ) const
      {
        Arg result;
        for (int i=0;i<arg.size();++i)
          result[i] = arg[i]*diag_[i];
        return result;
      }

    protected:
      Arg diag_;
    };

    // diagonalTensor
    // ------------
    template< class Arg >
    DiagonalTensor< Arg > diagonalTensor( const Arg &diag )
    {
      return DiagonalTensor< Arg >( diag );
    }

    // ZeroTensor
    // ----------

    template< class Dest >
    struct ZeroTensor
    {
      template< class ... Args >
      Dest operator() ( Args && ... ) const { return Dest(0); }
    };


    template< class Dest >
    ZeroTensor< Dest > zeroDest () {  return ZeroTensor< Dest >(); }


  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_COMMON_SIMPLETENSOR_HH
