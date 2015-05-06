#ifndef DUNE_FEM_ENVELOPE_HH
#define DUNE_FEM_ENVELOPE_HH

namespace Dune
{
  namespace Fem
  {

    template< class Object >
    class Envelope
    {
    public:
      typedef Object  ObjectType;

    protected:
      ObjectType object_;

    public:
      template< class ParamType >
      inline explicit Envelope ( ParamType param )
      : object_( param )
      {}

      inline Envelope ( const Envelope &other )
      : object_( other.object_ )
      {}

      inline const ObjectType& operator* () const
      {
        return object_;
      }

      inline ObjectType& operator* ()
      {
        return object_;
      }

      inline const ObjectType* operator-> () const
      {
        return &object_;
      }

      inline ObjectType* operator-> ()
      {
        return &object_;
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ENVELOPE_HH
