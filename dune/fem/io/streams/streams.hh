#ifndef DUNE_FEM_STREAMS_HH
#define DUNE_FEM_STREAMS_HH

#include <string>
#include <fstream>

// we would use cstdint,
// if it would be available for all compilers, e.g. clang
//#include <stdint.h>
#include <cstdint>

#include <dune/common/exceptions.hh>

#include <dune/fem/misc/bartonnackmaninterface.hh>
#include <dune/fem/misc/mpimanager.hh>

namespace Dune
{

  namespace Fem
  {

    class StreamError : public Exception {};


    /** \class OutStreamInterface
     *  \ingroup InOutStreams
     *  \brief abstract interface for an output stream
     *
     *  An output stream provides methods to write the basic C++ types into the
     *  stream. Based on this information, more complicated types can be written
     *  to the stream by using these basic output operations.
     *
     *  Normally, the output methods of the stream are not used directly, but the
     *  operator << ist used to write information into the stream. This operator
     *  should also be overloaded for all types that should be writeable.
     *
     *  Unlike STL streams, dune-fem output streams throw a StreamError exception
     *  when a writing operation fails. Since the program is automatically
     *  aborted, if the exception is not caught, careless programming will not
     *  result in corrupted data files.
     *
     *  \interfaceclass
     */
    template< class TraitsImp >
    class OutStreamInterface
    : public BartonNackmanInterface< OutStreamInterface< TraitsImp >, typename TraitsImp::OutStreamType >
    {
      typedef OutStreamInterface< TraitsImp > ThisType;
      typedef BartonNackmanInterface< ThisType, typename TraitsImp::OutStreamType > BaseType;

    public:
      //! type of the traits
      typedef TraitsImp Traits;

      //! type of the implementation (Barton-Nackman)
      typedef typename Traits::OutStreamType OutStreamType;

      //! type of the interface
      typedef ThisType OutStreamInterfaceType;

    protected:
      using BaseType::asImp;

    public:
      /** \brief flush the stream
       *
       *  By calling the flush method, the user can ensure that the stream is
       *  actually transferred (e.g., written to disk)
       */
      void flush ()
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().flush() );
      }

      /** \brief write a double to the stream
       *
       * \param[in]  value  value to write to the stream
       */
      void writeDouble ( const double value )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().writeDouble( value ) );
      }

      /** \brief write a float to the stream
       *
       * \param[in]  value  value to write to the stream
       */
      void writeFloat ( const float value )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().writeFloat( value ) );
      }

      /** \brief write an int to the stream
       *
       * \param[in]  value  value to write to the stream
       */
      void writeInt ( const int value )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().writeInt( value ) );
      }

      /** \brief write an int64_t to the stream
       *
       * \param[in]  value  value to write to the stream
       */
      void writeSignedInt64 ( int64_t value )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().writeSignedInt64( value ) );
      }

      /** \brief write a char to the stream
       *
       * \param[in]  value  value to write to the stream
       */
      void writeChar ( const char value )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().writeChar( value ) );
      }

      /** \brief write a bool to the stream
       *
       * \param[in]  value  value to write to the stream
       */
      void writeBool ( const bool value )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().writeBool( value ) );
      }

      /** \brief write a string to the stream
       *
       * \param[in]  s  string to write to the stream
       */
      void writeString ( const std::string &s )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().writeString( s ) );
      }

      /** \brief write an unsigned int to the stream
       *
       * \param[in]  value  value to write to the stream
       */
      void writeUnsignedInt ( unsigned int value )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().writeUnsignedInt( value ) );
      }

      /** \brief write an uint64_t to the stream
       *
       * \param[in]  value  value to write to the stream
       */
      void writeUnsignedInt64 ( uint64_t value )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().writeUnsignedInt64( value ) );
      }

    protected:
      void writeError () const
      {
        DUNE_THROW( StreamError, "Unable to write to stream." );
      }
    };



    /** \class InStreamInterface
     *  \ingroup InOutStreams
     *  \brief abstract interface for an input stream
     *
     *  An input stream provides methods to read the basic C++ types from the
     *  stream. Based on this information, more complicated types can be read
     *  from the stream by using these basic input operations.
     *
     *  Normally, the input methods of the stream are not used directly, but the
     *  operator >> ist used to read information from the stream. This operator
     *  should also be overloaded for all types that should be readable.
     *
     *  Unlike STL streams, dune-fem input streams throw a StreamError exception
     *  when a reading operation fails. Since the program is automatically
     *  aborted, if the exception is not caught, careless programming will not
     *  result in uninitialized objects. If the user catches the exception, he
     *  may not assume the object, that should be read, to be in a defined or
     *  even useful state.
     *
     *  \interfaceclass
     */
    template< class TraitsImp >
    class InStreamInterface
    : public BartonNackmanInterface< InStreamInterface< TraitsImp >, typename TraitsImp::InStreamType >
    {
      typedef InStreamInterface< TraitsImp > ThisType;
      typedef BartonNackmanInterface< ThisType, typename TraitsImp::InStreamType > BaseType;

    public:
      //! type of the traits
      typedef TraitsImp Traits;

      //! type of the implementation (Barton-Nackman)
      typedef typename Traits::InStreamType InStreamType;

      //! type of the interface
      typedef ThisType InStreamInterfaceType;

    protected:
      using BaseType::asImp;

    public:
      /** \brief read a double from the stream
       *
       *  \param[out]  value  reference to the variable to read from the stream
       */
      void readDouble ( double &value )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readDouble( value ) );
      }

      /** \brief read a double from the stream
       *
       *  \returns a double read from the stream
       */
      double readDouble ()
      {
        double value;
        readDouble( value );
        return value;
      }

      /** \brief read a float from the stream
       *
       *  \param[out]  value  reference to the variable to read from the stream
       */
      void readFloat ( float &value )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readFloat( value ) );
      }

      /** \brief read a double from the stream
       *
       *  \returns a double read from the stream
       */
      float readFloat ()
      {
        float value;
        readFloat( value );
        return value;
      }

      /** \brief read an int from the stream
       *
       *  \param[out]  value  reference to the variable to read from the stream
       */
      void readInt ( int &value )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readInt( value ) );
      }

      /** \brief read an int from the stream
       *
       *  \returns an int read from the stream
       */
      int readInt ()
      {
        int value;
        readInt( value );
        return value;
      }

      /** \brief read an int64_t from the stream
       *
       *  \param[out]  value  reference to the variable to read from the stream
       */
      void readSignedInt64 ( int64_t &value )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readSignedInt64( value ) );
      }

      /** \brief read an int64_t from the stream
       *
       *  \param[out]  value  reference to the variable to read from the stream
       */
      int64_t readSignedInt64 ()
      {
        int64_t value;
        readSignedInt64( value );
        return value;
      }

      /** \brief read a char from the stream
       *
       *  \param[out]  value  reference to the variable to read from the stream
       */
      void readChar ( char &value )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readChar( value ) );
      }

      /** \brief read a char from the stream
       *
       *  \returns a char read from the stream
       */
      int readChar ()
      {
        char value;
        readChar( value );
        return value;
      }

      /** \brief read a bool from the stream
       *
       *  \param[out]  value  reference to the variable to read from the stream
       *
       */
      void readBool ( bool &value )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readBool( value ) );
      }

      /** \brief read a bool from the stream
       *
       *  \returns a bool read from the stream
       */
      bool readBool ()
      {
        bool value;
        readBool( value );
        return value;
      }

      /** \brief read a string from the stream
       *
       *  \param[out]  s  reference to the string to read from the stream
       */
      void readString ( std::string &s )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readString( s ) );
      }

      /** \brief read an unsigned int from the stream
       *
       *  \param[out]  value  reference to the variable to read from the stream
       */
      void readUnsignedInt ( unsigned int &value )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readUnsignedInt( value ) );
      }

      /** \brief read an unsigned int from the stream
       *
       *  \returns an unsigned int read from the stream
       */
      unsigned int readUnsignedInt ()
      {
        unsigned int value;
        readUnsignedInt( value );
        return value;
      }

      /** \brief read an uint64_t from the stream
       *
       *  \param[out]  value  reference to the variable to read from the stream
       */
      void readUnsignedInt64 ( uint64_t &value )
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readUnsignedInt64( value ) );
      }

      /** \brief read an uint64_t from the stream
       *
       *  \returns an uint64_t read from the stream
       */
      uint64_t readUnsignedInt64 ()
      {
        uint64_t value;
        readUnsignedInt64( value );
        return value;
      }

    protected:
      void readError () const
      {
        DUNE_THROW( StreamError, "Unable to read from stream." );
      }
    };

    /** \brief Factory class for Fem Streams to deal with different constructor
     *         parameters.
     */
    template <class StreamImpl>
    struct StreamFactory
    {
      //! type of MPI communicator
      typedef typename MPIHelper :: MPICommunicator MPICommunicatorType;

      /** \brief return pointer to stream object created by new.
       *
       *  \param[in] filename  name of file that the stream read/writes
       *  \param[in] rank      rank of process data is read/written (defaults to MPIManager::rank())
       *  \param[in] mpiComm   MPI communicator (defaults to MPIHelper :: getCommunicator())
       */
      static StreamImpl* create( const std::string& filename,
                                 const int rank = MPIManager::rank(),
                                 const MPICommunicatorType& mpiComm = MPIHelper :: getCommunicator() )
      {
        return new StreamImpl( filename );
      }
    };

  } // namespace Fem

} // end namespace Dune

#include "streams_inline.hh"

#endif // #ifndef DUNE_FEM_STREAMS_HH
