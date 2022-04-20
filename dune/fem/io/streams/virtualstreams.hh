#ifndef DUNE_FEM_VIRTUALSTREAMS_HH
#define DUNE_FEM_VIRTUALSTREAMS_HH

#include <dune/fem/io/streams/streams.hh>

namespace Dune
{

  namespace Fem
  {

    // Forward Declarations
    // --------------------

    class VirtualOutStream;

    class VirtualInStream;

    template< class Traits >
    VirtualOutStream virtualize ( OutStreamInterface< Traits > &stream );

    template< class Traits >
    VirtualInStream virtualize ( InStreamInterface< Traits > &stream );



    // VirtualOutStreamObject
    // ----------------------

    class VirtualOutStreamObject
    {
      typedef VirtualOutStreamObject ThisType;

      friend class VirtualOutStream;

    private:
      unsigned int refCount;

    protected:
      inline VirtualOutStreamObject ()
      : refCount( 0 )
      {}
      virtual ~VirtualOutStreamObject() {}

    public:
      virtual void flush () = 0;
      virtual void writeDouble ( double value ) = 0;
      virtual void writeFloat ( float value ) = 0;
      virtual void writeInt ( int value ) = 0;
      virtual void writeSignedInt64 ( int64_t value ) = 0;
      virtual void writeString ( const std::string &s ) = 0;
      virtual void writeUnsignedInt ( unsigned int value ) = 0;
      virtual void writeUnsignedInt64 ( uint64_t value ) = 0;
    };



    // VirtualOutStreamTraits
    // ----------------------

    struct VirtualOutStreamTraits
    {
      typedef VirtualOutStream OutStreamType;
    };



    // VirtualOutStream
    // ----------------

    class VirtualOutStream
    : public OutStreamInterface< VirtualOutStreamTraits >
    {
      typedef VirtualOutStream ThisType;
      typedef OutStreamInterface< VirtualOutStreamTraits > BaseType;

      template< class T >
      friend VirtualOutStream virtualize ( OutStreamInterface< T > & );

    private:
      VirtualOutStreamObject *const stream_;

    private:
      explicit VirtualOutStream ( VirtualOutStreamObject *stream )
      : stream_( stream )
      {
        ++stream_->refCount;
      }

    public:
      VirtualOutStream ( const ThisType &other )
      : stream_( other.stream_ )
      {
        ++stream_->refCount;
      }

      virtual ~VirtualOutStream ()
      {
        if( --stream_->refCount == 0 )
          delete stream_;
      }

    private:
      ThisType &operator= ( const ThisType & );

    public:
      void flush ()
      {
        stream_->flush();
      }

      void writeDouble ( double value )
      {
        stream_->writeDouble( value );
      }

      void writeFloat ( float value )
      {
        stream_->writeFloat( value );
      }

      void writeInt ( int value )
      {
        stream_->writeInt( value );
      }

      void writeSignedInt64 ( int64_t value )
      {
        stream_->writeSignedInt64( value );
      }

      void writeString ( const std :: string &s )
      {
        stream_->writeString( s );
      }

      void writeUnsignedInt ( unsigned int value )
      {
        stream_->writeUnsignedInt( value );
      }

      void writeUnsignedInt64 ( unsigned int value )
      {
        stream_->writeUnsignedInt64( value );
      }
    };



    // VirtualInStreamObject
    // ---------------------

    class VirtualInStreamObject
    {
      typedef VirtualInStream ThisType;

      friend class VirtualInStream;

    private:
      unsigned int refCount;

    protected:
      VirtualInStreamObject ()
      : refCount( 0 )
      {}

      virtual ~VirtualInStreamObject() {}

    public:
      virtual void readDouble ( double &value ) = 0;
      virtual void readFloat ( float &value ) = 0;
      virtual void readInt ( int &value ) = 0;
      virtual void readSignedInt64 ( int64_t &value ) = 0;
      virtual void readString ( std::string &s ) = 0;
      virtual void readUnsignedInt ( unsigned int &value ) = 0;
      virtual void readUnsignedInt64 ( uint64_t &value ) = 0;
    };



    // VirtualInStreamTraits
    // ---------------------

    struct VirtualInStreamTraits
    {
      typedef VirtualInStream InStreamType;
    };



    // VirtualInStream
    // ---------------

    class VirtualInStream
    : public InStreamInterface< VirtualInStreamTraits >
    {
      typedef VirtualInStream ThisType;
      typedef InStreamInterface< VirtualInStreamTraits > BaseType;

      template< class T >
      friend VirtualInStream virtualize ( InStreamInterface< T > & );

    private:
      VirtualInStreamObject *const stream_;

    private:
      explicit VirtualInStream ( VirtualInStreamObject *stream )
      : stream_( stream )
      {
        ++stream_->refCount;
      }

    public:
      VirtualInStream ( const ThisType &other )
      : stream_( other.stream_ )
      {
        ++stream_->refCount;
      }

      virtual ~VirtualInStream ()
      {
        if( --stream_->refCount == 0 )
          delete stream_;
      }

    private:
      ThisType &operator= ( const ThisType & );

    public:
      void readDouble ( double &value )
      {
        stream_->readDouble( value );
      }

      void readFloat ( float &value )
      {
        stream_->readFloat( value );
      }

      void readInt ( int &value )
      {
        stream_->readInt( value );
      }

      void readSignedInt64 ( int64_t &value )
      {
        stream_->readSignedInt64( value );
      }

      void readString ( std :: string &s )
      {
        stream_->readString( s );
      }

      void readUnsignedInt ( unsigned int &value )
      {
        stream_->readUnsignedInt( value );
      }

      void readUnsignedInt64 ( uint64_t &value )
      {
        stream_->readUnsignedInt64( value );
      }
    };



    // VirtualOutStreamWrapper
    // -----------------------

    template< class Traits >
    class VirtualOutStreamWrapper
    : public VirtualOutStreamObject
    {
      typedef VirtualOutStreamWrapper< Traits > ThisType;
      typedef VirtualOutStreamObject BaseType;

      template< class T >
      friend VirtualOutStream virtualize ( OutStreamInterface< T > & );

    public:
      typedef OutStreamInterface< Traits > StreamType;

    private:
      StreamType &stream_;

    private:
      explicit VirtualOutStreamWrapper ( StreamType &stream )
      : stream_( stream )
      {}

      VirtualOutStreamWrapper ( const ThisType & );
      ThisType &operator= ( const ThisType & );

    public:
      virtual void flush ()
      {
        stream_.flush();
      }

      virtual void writeDouble ( double value )
      {
        stream_.writeDouble( value );
      }

      virtual void writeFloat ( float value )
      {
        stream_.writeFloat( value );
      }

      virtual void writeInt ( int value )
      {
        stream_.writeInt( value );
      }

      virtual void writeSignedInt64 ( int64_t value )
      {
        stream_.writeSignedInt64( value );
      }

      virtual void writeString ( const std :: string &s )
      {
        stream_.writeString( s );
      }

      virtual void writeUnsignedInt( unsigned int value )
      {
        stream_.writeUnsignedInt( value );
      }

      virtual void writeUnsignedInt64( uint64_t value )
      {
        stream_.writeUnsignedInt64( value );
      }
    };



    // VirtualInStreamWrapper
    // ----------------------

    template< class Traits >
    class VirtualInStreamWrapper
    : public VirtualInStreamObject
    {
      typedef VirtualInStreamWrapper< Traits > ThisType;
      typedef VirtualInStreamObject BaseType;

      template< class T >
      friend VirtualInStream virtualize ( InStreamInterface< T > & );

    public:
      typedef InStreamInterface< Traits > StreamType;

    private:
      StreamType &stream_;

    private:
      explicit VirtualInStreamWrapper ( StreamType &stream )
      : stream_( stream )
      {}

      VirtualInStreamWrapper ( const ThisType & );
      ThisType &operator= ( const ThisType & );

    public:
      virtual void readDouble ( double &value )
      {
        stream_.readDouble( value );
      }

      virtual void readFloat ( float &value )
      {
        stream_.readFloat( value );
      }

      virtual void readInt ( int &value )
      {
        stream_.readInt( value );
      }

      virtual void readSignedInt64 ( int64_t &value )
      {
        stream_.readSignedInt64( value );
      }

      virtual void readString ( std :: string &s )
      {
        stream_.readString( s );
      }

      virtual void readUnsignedInt( unsigned int &value )
      {
        stream_.readUnsignedInt( value );
      }

      virtual void readUnsignedInt64( uint64_t &value )
      {
        stream_.readUnsignedInt64( value );
      }
    };



    // virtualize
    // ----------

    template< class Traits >
    VirtualOutStream virtualize ( OutStreamInterface< Traits > &stream )
    {
      return VirtualOutStream( new VirtualOutStreamWrapper< Traits >( stream ) );
    }

    template< class Traits >
    VirtualInStream virtualize ( InStreamInterface< Traits > &stream )
    {
      return VirtualInStream( new VirtualInStreamWrapper< Traits >( stream ) );
    }

  } // namespave Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_VIRTUALSTREAMS_HH
