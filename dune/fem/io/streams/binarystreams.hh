#ifndef DUNE_FEM_BINARYSTREAMS_HH
#define DUNE_FEM_BINARYSTREAMS_HH

#include <dune/fem/io/streams/standardstreams.hh>

namespace Dune 
{
  namespace Fem 
  {
    /** \class BinaryFileOutStream
     *  \copydoc StandardOutStream 
     */
    class BinaryFileOutStream : public StandardOutStream 
    {
      typedef BinaryFileOutStream ThisType;
      typedef StandardOutStream   BaseType;

      // make copy constructor private because of the file pointers 
      BinaryFileOutStream( const BinaryFileOutStream& );
    public:
      /** \brief constructor
       *
       *  \param[in]  filename  name of a file to write to
       */
      explicit BinaryFileOutStream ( const std::string &filename )
        : BaseType( openFile( filename ) )
      {
      }

      /** \brief destructor deleteing file stream */
      ~BinaryFileOutStream() 
      {
        delete file_; file_ = 0;
      }
    protected:  
      std::ostream& openFile( const std::string& filename ) 
      {
        // init file 
        file_ = new std::ofstream( filename.c_str(), std::ios::binary );

        if( ! (*file_) ) 
          DUNE_THROW( Dune::IOError, "Unable to open file: " << filename );

        return *file_;  
      }

      //! standard file stream 
      std::ofstream* file_;  
    };

    /** \class BinaryFileInStream
     *  \copydoc StandardInStream 
     */
    class BinaryFileInStream: public StandardInStream 
    {
      typedef BinaryFileInStream ThisType;
      typedef StandardInStream   BaseType;

      // make copy constructor private because of the pointers 
      BinaryFileInStream ( const BinaryFileInStream& );

    public:
      /** \brief constructor
       *
       *  \param[in]  filename  name of a file to write to
       */
      explicit BinaryFileInStream ( const std::string &filename )
      : BaseType( openFile( filename ) )
      {
      }

      /** \brief destructor deleteing file stream */
      ~BinaryFileInStream() 
      {
        delete file_; file_ = 0;
      }

    protected:
      std::istream& openFile( const std::string& filename ) 
      {
        file_ = (new std::ifstream( filename.c_str(), std::ios::binary ));

        if( ! (*file_) ) 
          DUNE_THROW( Dune::IOError, "Unable to open file: " << filename );

        return *file_;  
      }

      //! standard file stream 
      std::ifstream* file_;
    };

  } // end namespace Fem   

  // #if DUNE_FEM_COMPATIBILITY  
  // put this in next version 1.4 
   
  using Fem :: BinaryFileOutStream ;
  using Fem :: BinaryFileInStream ;
  // #endif // DUNE_FEM_COMPATIBILITY

} // end namespace Dune

#endif // #ifndef DUNE_FEM_BINARYSTREAMS_HH
