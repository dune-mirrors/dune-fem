#ifndef DUNE_FEM_COMMON_SIMPLETENSOR_HH
#define DUNE_FEM_COMMON_SIMPLETENSOR_HH


namespace Dune
{

  namespace Fem
  {

    // SimpleTensor
    // ------------

    template< class Arg, class Dest >
    struct SimpleTensor
    {
      SimpleTensor () = default;
      SimpleTensor ( const SimpleTensor & ) = default;
      SimpleTensor ( SimpleTensor && ) = default;

      Dest operator() ( const Arg &arg ) const
      {
        return Dest();
      }
    };


    // IdTensor
    // --------

    template< class Arg >
    struct IdTensor
      : public SimpleTensor< Arg, Arg >
    {
      IdTensor () = default;
      IdTensor ( const IdTensor & ) = default;
      IdTensor ( IdTensor && ) = default;

      Arg operator() ( const Arg &arg ) const { return arg; }
    };


    // ScaledIdTensor
    // --------------

    template< class Arg >
    struct ScaledIdTensor
      : public SimpleTensor< Arg, Arg >
    {
      typedef typename FieldTraits< Arg >::field_type FieldType;
      ScaledIdTensor ( FieldType alpha ) : alpha_( alpha ) {}
      ScaledIdTensor ( const ScaledIdTensor & ) = default;
      ScaledIdTensor ( ScaledIdTensor && ) = default;

      Arg operator() ( const Arg &arg ) const
      {
        return alpha_ * arg;
      }

    protected:
      FieldType alpha_;
    };


    template< class Arg, class Dest >
    struct ZeroTensor
      : public SimpleTensor< Arg, Dest >
    {
      ZeroTensor () {}
      ZeroTensor ( const ZeroTensor & ) = default;
      ZeroTensor ( ZeroTensor && ) = default;

      Dest operator() ( const Arg &arg ) const
      {
        return Dest();
      }
    };


  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_COMMON_SIMPLETENSOR_HH
