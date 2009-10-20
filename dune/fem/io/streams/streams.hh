#ifndef DUNE_FEM_STREAMS_HH
#define DUNE_FEM_STREAMS_HH

#include <string>
#include <fstream>

#include <dune/common/exceptions.hh>

#include <dune/fem/misc/bartonnackmaninterface.hh>

namespace Dune
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
  : public BartonNackmanInterface< OutStreamInterface< TraitsImp >,
                                   typename TraitsImp :: OutStreamType >
  {
  public:
    //! type of the traits
    typedef TraitsImp Traits;

    //! type of the implementation (Barton-Nackman)
    typedef typename Traits :: OutStreamType OutStreamType;

  private:
    typedef OutStreamInterface< Traits > ThisType;
    typedef BartonNackmanInterface< ThisType, OutStreamType > BaseType;

  public:
    //! type of the interface
    typedef ThisType OutStreamInterfaceType;

  protected:
    using BaseType :: asImp;

  public:
    /** \brief flush the stream
     * 
     *  By calling the flush method, the user can ensure that the stream is
     *  actually transferred (e.g., written to disk)
     */
    inline void flush ()
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().flush() );
    }

    /** \brief write a double to the stream
     *
     * \param[in]  value  value to write to the stream
     */
    inline void writeDouble ( double value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().writeDouble( value ) );
    }
    
    /** \brief write a float to the stream
     *
     * \param[in]  value  value to write to the stream
     */
    inline void writeFloat ( float value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().writeFloat( value ) );
    }
    
    /** \brief write an int to the stream
     *
     * \param[in]  value  value to write to the stream
     */
    inline void writeInt ( int value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().writeInt( value ) );
    }

    /** \brief write a string to the stream
     *
     * \param[in]  s  string to write to the stream
     */
    inline void writeString ( const std :: string &s )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().writeString( s ) );
    }

    /** \brief write an unsigned int to the stream
     *
     * \param[in]  value  value to write to the stream
     */
    inline void writeUnsignedInt ( unsigned int value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().writeUnsignedInt( value ) );
    }

  protected:
    inline void writeError () const
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
  : public BartonNackmanInterface< InStreamInterface< TraitsImp >,
                                   typename TraitsImp :: InStreamType >
  {
  public:
    //! type of the traits
    typedef TraitsImp Traits;

    //! type of the implementation (Barton-Nackman)
    typedef typename Traits :: InStreamType InStreamType;

  private:
    typedef InStreamInterface< Traits > ThisType;
    typedef BartonNackmanInterface< ThisType, InStreamType > BaseType;

  public:
    //! type of the interface
    typedef ThisType InStreamInterfaceType;

  protected:
    using BaseType :: asImp;
    
  public:
    /** \brief read a double from the stream
     *
     *  \param[out]  value  reference to the variable to read from the stream
     */
    inline void readDouble ( double &value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readDouble( value ) );
    }

    /** \brief read a double from the stream
     *
     *  \returns a double read from the stream
     */
    inline double readDouble ()
    {
      double value;
      readDouble( value );
      return value;
    }
    
    /** \brief read a float from the stream
     *
     *  \param[out]  value  reference to the variable to read from the stream
     */
    inline void readFloat ( float &value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readFloat( value ) );
    }
    
    /** \brief read a double from the stream
     *
     *  \returns a double read from the stream
     */
    inline float readFloat ()
    {
      float value;
      readFloat( value );
      return value;
    }
    
    /** \brief read an int from the stream
     *
     *  \param[out]  value  reference to the variable to read from the stream
     */
    inline void readInt ( int &value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readInt( value ) );
    }
    
    /** \brief read an int from the stream
     *
     *  \returns an int read from the stream
     */
    inline int readInt ()
    {
      int value;
      readInt( value );
      return value;
    }

    /** \brief read a string from the stream
     *
     *  \param[out]  s  reference to the string to read from the stream
     */
    inline void readString ( std :: string &s )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readString( s ) );
    }

    /** \brief read an unsigned int from the stream
     *
     *  \param[out]  value  reference to the variable to read from the stream
     */
    inline void readUnsignedInt ( unsigned int &value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readUnsignedInt( value ) );
    }

    /** \brief read an unsigned int from the stream
     *
     *  \returns an unsigned int read from the stream
     */
    inline unsigned int readUnsignedInt ()
    {
      unsigned int value;
      readUnsignedInt( value );
      return value;
    }

  protected:
    inline void readError () const
    {
      DUNE_THROW( StreamError, "Unable to read from stream." );
    }
  };

}

#include "streams_inline.hh"

#endif
