#ifndef DUNE_FEM_STREAMS_HH
#define DUNE_FEM_STREAMS_HH

#include <string>

#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune
{

  /** \class OutStreamInterface
   *  \brief abstract interface for an output stream
   */
  template< class TraitsImp >
  class OutStreamInterface
  {
  public:
    //! type of the traits
    typedef TraitsImp Traits;

  private:
    typedef OutStreamInterface< Traits > ThisType;

  public:
    //! type of the interface
    typedef ThisType OutStreamInterfaceType;

    //! type of the implementation (Barton-Nackman)
    typedef typename Traits :: OutStreamType OutStreamType;

  public:
    /** \brief returns \b true, if an error occured
     *
     *  This method is implemented to return the negated return value
     *  of valid.
     */
    inline bool operator! () const
    {
      return !valid();
    }

    /** \brief flush the stream
     * 
     *  By calling the flush method, the user can ensure that the stream is
     *  actually transferred (e.g., written to disk)
     */
    inline void flush ()
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().flush() );
    }
    
    /** \brief return \b false, if an error occured
     */
    inline bool valid () const
    {
      CHECK_INTERFACE_IMPLEMENTATON( asImp().valid() );
      return asImp().valid();
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
    // Barton-Nackman trick
    inline const OutStreamType &asImp () const
    {
      return static_cast< const OutStreamType & >( *this );
    }

    // Barton-Nackman trick
    inline OutStreamType &asImp ()
    {
      return static_cast< OutStreamType & >( *this );
    }
  };


  
  /** \class InStreamInterface
   *  \brief abstract interface for an input stream
   */
  template< class TraitsImp >
  class InStreamInterface
  {
  public:
    //! type of the traits
    typedef TraitsImp Traits;

  private:
    typedef InStreamInterface< Traits > ThisType;

  public:
    //! type of the interface
    typedef ThisType InStreamInterfaceType;

    //! type of the implementation (Barton-Nackman)
    typedef typename Traits :: InStreamType InStreamType;

  public:
    /** \brief returns \b true, if an error occured
     *
     *  This method is implemented to return the negated return value
     *  of valid.
     */
    inline bool operator! () const
    {
      return !valid();
    }

    /** \brief return \b false, if an error occured
     */
    inline bool valid () const
    {
      CHECK_INTERFACE_IMPLEMENTATON( asImp().valid() );
      return asImp().valid();
    }
   
    /** \brief read a double from the stream
     *
     *  \param[out]  value  reference to the variable to read from the stream
     */
    inline void readDouble ( double &value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readDouble( value ) );
    }
    
    /** \brief read a float from the stream
     *
     *  \param[out]  value  reference to the variable to read from the stream
     */
    inline void readFloat ( float &value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readFloat( value ) );
    }
    
    /** \brief read an int from the stream
     *
     *  \param[out]  value  reference to the variable to read from the stream
     */
    inline void readInt ( int &value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().readInt( value ) );
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

  protected:
    // Barton-Nackman trick
    inline const InStreamType &asImp () const
    {
      return static_cast< const InStreamType & >( *this );
    }

    // Barton-Nackman trick
    inline InStreamType &asImp ()
    {
      return static_cast< InStreamType & >( *this );
    }
  };
 
}

#include "streams_inline.hh"

#endif
