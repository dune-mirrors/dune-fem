#ifndef DUNE_FEM_SIONLIBSTREAMS_HH
#define DUNE_FEM_SIONLIBSTREAMS_HH

#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/streams/standardstreams.hh>

#include <dune/fem/misc/mpimanager.hh>

#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/mpicommunication.hh>

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
       *  \param[in]  rank      process rank (defaults to MPIManager::rank())
       *  \param[in]  mpiComm   MPI communicator (defaults to MPIHelper :: getCommunicator() )
       *
       *  \note The filename must be the same on all ranks.
       */
      SIONlibOutStream ( const std::string &filename, const int rank,
                         MPICommunicatorType mpiComm, ParameterReader parameter = Parameter::container() )
        : BaseType( dataStream() ),
          filename_( filename ),
          mpiComm_( mpiComm ),
          rank_( rank ),
          numFiles_( std::min( MPIManager::size(), parameter.getValue< int >( "fem.io.sionlib.numfiles", 1 ) ) ),
          blockSize_( parameter.getValue< int >( "fem.io.sionlib.blocksize", -1 ) )
      {}

      SIONlibOutStream ( const std::string &filename, ParameterReader parameter = Parameter::container() )
        : SIONlibOutStream( filename, MPIManager::rank(), MPIHelper :: getCommunicator(), parameter )
      {}

      SIONlibOutStream ( const std::string &filename, const int rank, ParameterReader parameter = Parameter::container() )
        : SIONlibOutStream( filename, rank, MPIHelper :: getCommunicator(), parameter )
      {}

      /** \brief destructor writing internal data buffer to the file via SIONlib */
      ~SIONlibOutStream ()
      {
        writeFile();
        delete data_; data_ = 0;
      }

    protected:
      std::ostream& dataStream()
      {
        // init data stream
        data_ = new std::stringstream();
        assert( data_ );
        return *data_;
      }

      void writeFile()
      {
#if HAVE_SIONLIB && HAVE_MPI
        if( data_ )
        {
          // get data buffer
          std::string data ( data_->str() );

          // get chunk size for this process
          // use sionlib int64
          sion_int64 chunkSize = data.size();

          // file mode is: write byte
          const char* fileMode = "wb";

          // number of physical files to be created
          int numFiles = numFiles_;

          // block size of filesystem, -1 use system default
          int blockSize = blockSize_;

          // my rank
          int rank = rank_;
          // file pointer
          FILE* file = 0;

          // open sion file
          int sid =
            sion_paropen_mpi( (char *) filename_.c_str(),
                              (char *) fileMode,
                              &numFiles, // number of physical files
                              mpiComm_, // global comm
                              &mpiComm_, // local comm
                              &chunkSize, // maximal size of data to be written
                              &blockSize, // filesystem block size
                              &rank, // my rank
                              &file, // file pointer that is set by sion lib
                              NULL
                            );
          if( sid == -1 )
            DUNE_THROW( IOError, "opening sion_paropen_mpi for writing failed!" << filename_ );

          assert( file );

          // get pointer to buffer
          const char* buffer = data.c_str();
          // write data
          assert( sizeof(char) == 1 );
          sion_fwrite( buffer, 1, chunkSize, sid);

          // close file
          sion_parclose_mpi( sid );
        }
#else
        std::cerr << "WARNING: SIONlib or MPI not available" << std::endl;
#endif
      }

      //! standard file stream
      std::stringstream* data_;

      const std::string filename_;
      MPICommunicatorType mpiComm_;

      const int rank_;
      const int numFiles_;
      const int blockSize_;
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
       *  \param[in]  rank      process rank data is read for
       *  \param[in]  mpiComm   MPI communicator (defaults to MPIHelper :: getCommunicator() )
       *
       *  \note The filename must be the same on all ranks.
       */
      SIONlibInStream ( const std::string &filename,
                        const int rank = MPIManager :: rank(),
                        MPICommunicatorType mpiComm = MPIHelper :: getCommunicator() )
        : BaseType( readFile( filename , rank, mpiComm ) )
      {
      }

      /** \brief destructor deleting interal data buffer */
      ~SIONlibInStream ()
      {
        delete data_; data_ = 0;
      }

    protected:
      std::istream& readFile( const std::string& filename,
                              int rank,
                              MPICommunicatorType mpiComm )
      {
        data_ = new std::stringstream();
        assert( data_ );
#if HAVE_SIONLIB
        // chunkSize is set by sion_paropen_mpi
        sion_int64 chunkSize = 0;
        // file mode is: read byte
        const char* fileMode = "rb";

        // blockSize, is recovered from stored files
        int blockSize = -1;

        // file handle
        FILE* file = 0;

        // sion file handle
        int sid = 0;

#if HAVE_MPI
        // number of files, is overwritten by sion_open
        int numFiles = 1;

        // if MPI is avaialbe use sion_paropen_mpi
        const int mpiSize = MPIManager :: size () ;

        if( mpiSize > 1 )
        {
          // open sion file
          sid = sion_paropen_mpi( (char *) filename.c_str(),
                                  (char *) fileMode,
                                  &numFiles, // numFiles
                                  mpiComm, // global comm
                                  &mpiComm, // local comm
                                  &chunkSize, // is set by library
                                  &blockSize, // block size
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
                                (char *) fileMode,
                                &chunkSize, // is set by library
                                &blockSize, // block size
                                &rank, // my rank
                                &file  // file pointer that is set by sion lib
                              );

          if( sid == -1 )
            DUNE_THROW( IOError, "opening sion_ropen for reading failed!" << filename );
        }

        assert( file );

        // get bytes available for reading (might differ from total chunkSize)
        chunkSize = sion_bytes_avail_in_block( sid );

        // create buffer
        std::string data;
        data.resize( chunkSize );

        // get pointer to buffer
        char* buffer = (char *) data.c_str();
        assert( sizeof(char) == 1 );
        // read data
        sion_fread( buffer, 1, chunkSize, sid );

        // write data to stream
        data_ = new std::stringstream();
        assert( data_ );
        data_->write( buffer, chunkSize );

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

    /** \brief Factory class for Fem Streams to deal with different constructor
     *         parameters.
     */
    template <>
    struct StreamFactory< SIONlibInStream >
    {
      //! type of MPI communicator
      typedef MPIHelper :: MPICommunicator MPICommunicatorType;

      /** \brief return pointer to stream object created by new.
       *
       *  \param[in] filename  name of file that the stream read/writes
       *  \param[in] rank      rank of process data is read/written (defaults to MPIManager::rank())
       *  \param[in] mpiComm   MPI communicator (defaults to MPIHelper :: getCommunicator())
       */
      static SIONlibInStream* create( const std::string& filename,
                                      const int rank = MPIManager::rank(),
                                      const MPICommunicatorType& mpiComm = MPIHelper :: getCommunicator() )
      {
        return new SIONlibInStream( filename, rank, mpiComm );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BINARYSTREAMS_HH
