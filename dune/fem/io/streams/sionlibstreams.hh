#ifndef DUNE_FEM_SIONLIBSTREAMS_HH
#define DUNE_FEM_SIONLIBSTREAMS_HH

#include <dune/fem/io/streams/standardstreams.hh>

#include <dune/fem/misc/mpimanager.hh>

#include <dune/common/collectivecommunication.hh>
#include <dune/common/mpicollectivecommunication.hh>

#if HAVE_SIONLIB

#if HAVE_MPI
#define SION_MPI
#endif

#include <sion.h>
#endif

namespace Dune
{

  namespace Fem 
  {
    /** \class SIONlibOutStream
     *  \ingroup InOutStreams
     *  \brief output stream writing into a single file with the SIONlib
     *  (http://www2.fz-juelich.de/jsc/sionlib/)
     *
     *  \note This stream directly stores the binary representation of the data.
     *        The binary representation of the stored data is always that of the
     *        current machine. On read the data is converted accordingly on machines 
     *        with different endianess.
     *
     *  \newimplementation
     */
    //template <class Communicator>
    class SIONlibOutStream : public StandardOutStream 
    {
      //typedef SIONlibOutStream< Communicator > ThisType;
      typedef SIONlibOutStream   ThisType;
      typedef StandardOutStream  BaseType;

      typedef MPIHelper  :: MPICommunicator          MPICommunicatorType;

      // don't allow copying because of interal pointers 
      SIONlibOutStream( const SIONlibOutStream& ) ;
    public:
      /** \brief constructor
       *
       *  \param[in]  filename  name of a global file to write to
       *  \param[in]  mpiComm   MPI communicator (defaults to MPIHelper :: getCommunicator() )
       *
       *  \note The filename must be the same on all ranks. 
       */
      SIONlibOutStream ( const std::string &filename,
                         MPICommunicatorType mpiComm = MPIHelper :: getCommunicator() )
        : BaseType( dataStream() ),
          mpiComm_( mpiComm ),
          filename_( filename )
      {
      }

      /** \brief destructor writing internal data buffer to the file via SIONlib */
      ~SIONlibOutStream () 
      { 
        writeFile();
        delete data_; data_ = 0;
      }

    protected:  
      std::ostream& dataStream() 
      {
        // init file 
        data_ = new std::stringstream();
        assert( data_ );
        return *data_;  
      }

      void writeFile() 
      {
#if HAVE_SIONLIB && HAVE_MPI
        if( data_ ) 
        {
          assert( data_ );
          // get data buffer 
          std::string data ( data_->str() );

          // get chunk size  
          size_t maxSize = data.size() + sizeof( int ); 
          // maxSize = MPIManager::comm().max( maxSize );

          sion_int64 chunkSize = maxSize;
          std::string fileMode( "wb" );

          int numFiles = 1;
          int blockSize = -1;
          int rank = MPIManager::rank();
          FILE* file = 0;

          // open sion file 
          int sid = 
            sion_paropen_mpi( (char *) filename_.c_str(),
                              (char *) fileMode.c_str(), 
                              &numFiles, // numFiles (0 means 1 file)
                              mpiComm_, // global comm 
                              &mpiComm_, // local comm 
                              &chunkSize, // maximal size of data to be written 
                              &blockSize, // default block size
                              &rank, // my rank 
                              &file, // file pointer that is set by sion lib
                              NULL 
                            ); 
          if( sid == -1 ) 
            DUNE_THROW( IOError, "opening sion_paropen_mpi for writing failed!" << filename_ );

          assert( file );

          const int dataSize = data.size();
          // write size 
          fprintf(file, "%d", dataSize );

          const char* buffer = data.c_str();
          // write data 
          for( int i=0; i<dataSize; ++i ) 
          {
            fprintf(file,"%c", buffer[ i ] ); 
          }

          // close file 
          sion_parclose_mpi( sid );
        }
#else 
        std::cerr << "WARNING: SIONlib or MPI not available" << std::endl;
#endif
      }

      MPICommunicatorType mpiComm_;
      const std::string filename_;

      //! standard file stream 
      std::stringstream* data_;  
    };

    /** \class SIONlibInStream
     *  \ingroup InOutStreams
     *  \brief input stream reading from a file in binary form
     *
     *  \note This stream directly stores the binary representation of the data.
     *        The binary representation might differ between different machines
     *        (e.g., little endian vs. big endian).
     *
     *  \newimplementation
     */
    //template <class Communicator>
    class SIONlibInStream : public StandardInStream 
    {
      typedef SIONlibInStream  ThisType;
      typedef StandardInStream  BaseType;

      typedef MPIHelper  :: MPICommunicator          MPICommunicatorType;

      // don't allow copying because of interal pointers 
      SIONlibInStream( const SIONlibInStream& ) ;
    public:
      /** \brief constructor
       *
       *  \param[in]  filename  name of a file to read from 
       *  \param[in]  mpiComm   MPI communicator (defaults to MPIHelper :: getCommunicator() )
       *
       *  \note The filename must be the same on all ranks. 
       */
      SIONlibInStream ( const std::string &filename,
                        MPICommunicatorType mpiComm = MPIHelper :: getCommunicator(),
                        const int rank = MPIManager :: rank() )
        : BaseType( readFile( filename , mpiComm, rank ) )
      {
      }

      /** \brief destructor deleting interal data buffer */
      ~SIONlibInStream () 
      { 
        delete data_; data_ = 0;
      }

    protected:  
      std::istream& readFile( const std::string& filename, 
                              MPICommunicatorType mpiComm,
                              int rank )
      {
        data_ = new std::stringstream();
        assert( data_ );
#if HAVE_SIONLIB 
        // chunkSize is set by sion_paropen_mpi 
        sion_int64 chunkSize = 0;
        std::string fileMode( "rb" );

        int numFiles = 1;
        int blockSize = -1;
        FILE* file = 0;

        // if MPI is avaialbe use sion_paropen_mpi
        const int mpiSize = MPIManager :: size () ;
        // sion file handle 
        int sid = 0;
#if HAVE_MPI
        if( mpiSize > 1 )
        {
          // open sion file 
          sid = sion_paropen_mpi( (char *) filename.c_str(),
                                  (char *) fileMode.c_str(), 
                                  &numFiles, // numFiles (0 means 1 file)
                                  mpiComm, // global comm 
                                  &mpiComm, // local comm 
                                  &chunkSize, // is set by library
                                  &blockSize, // default block size
                                  &rank, // my rank 
                                  &file, // file pointer that is set by sion lib
                                  NULL 
                                ); 

          if( sid == -1 ) 
            DUNE_THROW( IOError, "opening sion_paropen_mpi for reading failed!" << filename );
        }
        else 
#endif
        // serial open of file
        {
          // open sion file for reading only rank information 
          sid = sion_open_rank( (char *) filename.c_str(),
                                (char *) fileMode.c_str(), 
                                &chunkSize, // is set by library 
                                &blockSize, // default block size
                                &rank, // my rank 
                                &file  // file pointer that is set by sion lib
                              ); 

          if( sid == -1 ) 
            DUNE_THROW( IOError, "opening sion_ropen for reading failed!" << filename );
        }

        assert( file );

        // write size 
        int dataSize = chunkSize - sizeof( int );
        fscanf(file, "%d", &dataSize );
        // write data 
        std::string data ; 
        data.resize( dataSize );
        char* buffer = (char *) data.c_str();
        for( int i=0; i<dataSize; ++i ) 
        {
          fscanf(file, "%c", &buffer[ i ] ); 
        }

        data_ = new std::stringstream();
        (*data_) << data ;

#if HAVE_MPI
        // close file 
        if( mpiSize > 1 ) 
        {
          sion_parclose_mpi( sid );
        }
        else 
#endif
        {
          sion_close( sid );
        }
#endif // HAVE_SIONLIB

        return *data_;
      }
      //! standard file stream 
      std::stringstream* data_;  
    };

  } // end namespace Fem   

  // #if DUNE_FEM_COMPATIBILITY  
  // put this in next version 1.4 

  using Fem :: SIONlibOutStream ;
  using Fem :: SIONlibInStream ;
  // #endif // DUNE_FEM_COMPATIBILITY

} // end namespace Dune

#endif // #ifndef DUNE_FEM_BINARYSTREAMS_HH
