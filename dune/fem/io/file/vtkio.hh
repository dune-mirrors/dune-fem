#ifndef DUNE_FEM_VTKIO_HH
#define DUNE_FEM_VTKIO_HH

#include <type_traits>

#include <dune/common/deprecated.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/fem/version.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem/misc/threads/domainthreaditerator.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class GridPart, bool subsampling = false >
    class VTKIO;



    // VTKFunctionWrapper
    // ------------------

    template< class DF >
    class VTKFunctionWrapper
    : public VTKFunction< typename DF::GridPartType::GridViewType >
    {
      typedef VTKFunctionWrapper< DF > ThisType;

    public:
      enum TypeOfField { real, complex_real, complex_imag };
      typedef DF DiscreteFunctionType;

      typedef ConstLocalFunction< DF > LocalFunctionType;
      typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;

      static const int dimRange = FunctionSpaceType::dimRange;
      static const int dimDomain = FunctionSpaceType::dimDomain;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;

      typedef typename DiscreteFunctionType::GridPartType GridPartType;
      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

      typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;

      //! constructor taking discrete function
      VTKFunctionWrapper ( const DiscreteFunctionType& df,
                           const std::string& dataName,
                           int component, bool vector, TypeOfField typeOfField )
      : localFunction_( df ),
        name_( ( dataName.size() > 0 ) ? dataName : df.name() ),
        vector_( vector ),
        typeOfField_( typeOfField ),
        component_( component )
      {}

      //! virtual destructor
      virtual ~VTKFunctionWrapper ()
      {}

      //! return number of components
      virtual int ncomps () const
      {
        return (!vector_) ? 1 : 3; // dimDomain;
      }

      //! evaluate single component comp in
      //! the entity
      virtual double evaluate ( int comp, const EntityType &e, const LocalCoordinateType &xi ) const
      {
        localFunction_.bind( e );
        typedef typename LocalFunctionType::RangeFieldType RangeFieldType;
        RangeType val;
        localFunction_.evaluate(xi,val);
        localFunction_.unbind();

        RangeFieldType outVal( 0 );
        if (vector_)
        {
          if( comp <= dimDomain )
            outVal = val[ comp + component_ ] ;
        }
        else
          outVal = val[ component_ ] ;
        if (typeOfField_ == TypeOfField::real || typeOfField_ == TypeOfField::complex_real )
          return std::real( outVal );
        else
          return std::imag( outVal );
      }

      //! get name
      virtual std::string name () const
      {
        std::stringstream ret;
        ret << name_;
        if (typeOfField_ == TypeOfField::complex_real)
          ret << "_real_";
        if (typeOfField_ == TypeOfField::complex_imag)
          ret << "_imag_";
        if (vector_)
          ret << "_vec" << component_;
        else
          ret << component_;
        return ret.str();
      }

    private:
      mutable LocalFunctionType localFunction_;
      const std::string name_ ;
      const bool vector_;
      const TypeOfField typeOfField_;
      const int component_;
    };



    // VTKIOBase
    // ---------

    //! \brief Output using VTK
    template< class GridPart, bool subsampling >
    class VTKIOBase
    {
      typedef VTKIOBase< GridPart, subsampling > ThisType;

    protected:
      typedef typename GridPart::GridViewType GridViewType;

      class VTKWriter;
      class SubsamplingVTKWriter;

      typedef typename std::conditional< subsampling, SubsamplingVTKWriter, VTKWriter >::type
        VTKWriterType;

    public:
      typedef GridPart GridPartType;

      typedef typename GridPartType::GridType GridType;
      typedef typename GridPartType::IndexSetType IndexSetType;

    protected:
      class PartitioningData
        : public VTKFunction< GridViewType >
      {
        typedef PartitioningData   ThisType;

      public:
        typedef typename GridViewType :: template Codim< 0 >::Entity EntityType;
        typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;

        typedef DomainDecomposedIteratorStorage< GridPartType >  ThreadIteratorType;

        //! constructor taking discrete function
        PartitioningData( const GridPartType& gridPart,
                          const std::string& name,
                          const int rank, const int nThreads )
          : iterators_( gridPart ), name_( name ), rank_( rank ), nThreads_( nThreads ) {}

        //! virtual destructor
        virtual ~PartitioningData () {}

        //! return number of components
        virtual int ncomps () const { return 1; }

        //! evaluate single component comp in
        //! the entity
        virtual double evaluate ( int comp, const EntityType &e, const LocalCoordinateType &xi ) const
        {
          const int thread = iterators_.thread( e );
          return (nThreads_ < 0) ? double( rank_ ) : double( rank_ * nThreads_ + thread );
        }

        //! get name
        virtual std::string name () const
        {
          return name_;
        }

      private:
        ThreadIteratorType iterators_;
        const std::string name_;
        const int rank_;
        const int nThreads_;

      };

      class VolumeData
        : public VTKFunction< GridViewType >
      {
        typedef PartitioningData   ThisType;

      public:
        typedef typename GridViewType :: template Codim< 0 >::Entity EntityType;
        typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;

        typedef DomainDecomposedIteratorStorage< GridPartType >  ThreadIteratorType;

        //! constructor taking discrete function
        VolumeData() {}

        //! virtual destructor
        virtual ~VolumeData () {}

        //! return number of components
        virtual int ncomps () const { return 1; }

        //! evaluate single component comp in
        //! the entity
        virtual double evaluate ( int comp, const EntityType &e, const LocalCoordinateType &xi ) const
        {
          return e.geometry().volume();
        }

        //! get name
        virtual std::string name () const
        {
          return std::string("volume");
        }
      };

      int getPartitionParameter ( const ParameterReader &parameter = Parameter::container() ) const
      {
        // 0 = none, 1 = MPI ranks only, 2 = ranks + threads, 3 = like 1 and also threads only
        const std::string names[] = { "none", "rank", "rank+thread", "rank/thread" };
        return parameter.getEnum( "fem.io.partitioning", names, 0 );
      }

    protected :
      VTKIOBase ( const GridPartType &gridPart, VTKWriterType *vtkWriter, const ParameterReader &parameter = Parameter::container() )
      : gridPart_( gridPart ),
        vtkWriter_( vtkWriter ),
        addPartition_( getPartitionParameter( parameter ) )
      {
        static const std::string typeTable[] = { "ascii", "base64", "appended-raw", "appended-base64" };
        static const VTK::OutputType typeValue[] = { VTK::ascii, VTK::base64, VTK::appendedraw, VTK::appendedbase64 };
        type_ = typeValue[ parameter.getEnum( "fem.io.vtk.type", typeTable, 2 ) ];
      }

      void addPartitionData( const int myRank = -1 )
      {
        if( addPartition_ > 0 )
        {
          std::shared_ptr<VolumeData> volumePtr( std::make_shared<VolumeData>() );
          vtkWriter_->addCellData( volumePtr );

          const int rank = ( myRank < 0 ) ? gridPart_.comm().rank() : myRank ;
          const int nThreads = ( addPartition_ > 1 ) ? MPIManager::maxThreads() : 1 ;
          if( addPartition_ <= 2 )
          {
            std::shared_ptr<PartitioningData> dataRankPtr( std::make_shared<PartitioningData>(gridPart_, "rank", rank, nThreads) );
            vtkWriter_->addCellData( dataRankPtr );
          }
          else
          {
            // rank only visualization
            std::shared_ptr<PartitioningData> dataRankPtr( std::make_shared<PartitioningData>(gridPart_, "rank", rank, -1) );
            vtkWriter_->addCellData( dataRankPtr );
            // thread only visualization
            std::shared_ptr<PartitioningData> dataThreadPtr( std::make_shared<PartitioningData>(gridPart_, "thread", 0, nThreads) );
            vtkWriter_->addCellData( dataThreadPtr );
          }
          addPartition_ = 0 ;
        }
      }

      template < class DF >
      static bool notComplex()
      {
        typedef typename DF::RangeFieldType RangeFieldType;
        typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;
        return ! std::is_same< typename std::remove_cv<RangeFieldType>::type, std::complex<RealType> >::value;
      }

    public:
      ~VTKIOBase ()
      {
        delete vtkWriter_;
      }

      //! return grid part
      const GridPartType &gridPart () const
      {
        return gridPart_;
      }

      template< class DF >
      void addCellData( DF &df , const std::string& dataName = "" )
      {
        static const int dimRange = DF::FunctionSpaceType::dimRange;
        for( int i = 0;i < dimRange; ++i )
        {
          if ( notComplex<DF>() )
          {
            std::shared_ptr<VTKFunctionWrapper< DF > > ptr( std::make_shared<VTKFunctionWrapper< DF > >( df, dataName, i,
                  false, VTKFunctionWrapper<DF>::TypeOfField::real) );
            vtkWriter_->addCellData( ptr );
          }
          else
          {
            std::shared_ptr<VTKFunctionWrapper< DF > > ptrR( std::make_shared<VTKFunctionWrapper< DF > >( df, dataName, i,
                  false, VTKFunctionWrapper<DF>::TypeOfField::complex_real) );
            vtkWriter_->addCellData( ptrR );
            std::shared_ptr<VTKFunctionWrapper< DF > > ptrI( std::make_shared<VTKFunctionWrapper< DF > >( df, dataName, i,
                  false, VTKFunctionWrapper<DF>::TypeOfField::complex_imag) );
            vtkWriter_->addCellData( ptrI );
          }
        }
      }

      template< class DF >
      void addVectorCellData( DF &df,
                              const std::string& dataName = "" ,
                              int startPoint = 0 )
      {
        if ( notComplex<DF>() )
        {
          std::shared_ptr<VTKFunctionWrapper< DF > > ptr( std::make_shared<VTKFunctionWrapper< DF > >( df, dataName, startPoint,
                  true, VTKFunctionWrapper<DF>::TypeOfField::real) );
          vtkWriter_->addCellData( ptr );
        }
        else
        {
          std::shared_ptr<VTKFunctionWrapper< DF > > ptrR( std::make_shared<VTKFunctionWrapper< DF > >( df, dataName, startPoint,
                  true, VTKFunctionWrapper<DF>::TypeOfField::complex_real) );
          vtkWriter_->addCellData( ptrR );
          std::shared_ptr<VTKFunctionWrapper< DF > > ptrI( std::make_shared<VTKFunctionWrapper< DF > >( df, dataName, startPoint,
                  true, VTKFunctionWrapper<DF>::TypeOfField::complex_imag) );
          vtkWriter_->addCellData( ptrI );
        }
      }

      template< class DF >
      void addVertexData( DF &df, const std::string& dataName = "" )
      {
        static const int dimRange = DF::FunctionSpaceType::dimRange;
        std::string name = ( dataName.size() > 0 ) ? dataName : df.name() ;
        for( int i = 0;i < dimRange; ++i )
        {
          if ( notComplex<DF>() )
          {
            std::shared_ptr<VTKFunctionWrapper< DF > > ptr( std::make_shared<VTKFunctionWrapper< DF > >( df, dataName, i,
                  false, VTKFunctionWrapper<DF>::TypeOfField::real) );
            vtkWriter_->addVertexData( ptr );
          }
          else
          {
            std::shared_ptr<VTKFunctionWrapper< DF > > ptrR( std::make_shared<VTKFunctionWrapper< DF > >( df, dataName, i,
                  false, VTKFunctionWrapper<DF>::TypeOfField::complex_real) );
            vtkWriter_->addVertexData( ptrR );
            std::shared_ptr<VTKFunctionWrapper< DF > > ptrI( std::make_shared<VTKFunctionWrapper< DF > >( df, dataName, i,
                  false, VTKFunctionWrapper<DF>::TypeOfField::complex_imag) );
            vtkWriter_->addVertexData( ptrI );
          }
        }
      }

      template< class DF >
      void addVectorVertexData( DF &df,
                                const std::string& dataName = "" ,
                                int startPoint = 0 )
      {
        if ( notComplex<DF>() )
        {
          std::shared_ptr<VTKFunctionWrapper< DF > > ptr( std::make_shared<VTKFunctionWrapper< DF > >( df, dataName, startPoint,
                  true, VTKFunctionWrapper<DF>::TypeOfField::real) );
          vtkWriter_->addVertexData( ptr );
        }
        else
        {
          std::shared_ptr<VTKFunctionWrapper< DF > > ptrR( std::make_shared<VTKFunctionWrapper< DF > >( df, dataName, startPoint,
                  true, VTKFunctionWrapper<DF>::TypeOfField::complex_real) );
          vtkWriter_->addVertexData( ptrR );
          std::shared_ptr<VTKFunctionWrapper< DF > > ptrI( std::make_shared<VTKFunctionWrapper< DF > >( df, dataName, startPoint,
                  true, VTKFunctionWrapper<DF>::TypeOfField::complex_imag) );
          vtkWriter_->addVertexData( ptrI );
        }
      }

      void clear ()
      {
        vtkWriter_->clear();
      }

      std::string write ( const std::string &name, VTK::OutputType type )
      {
        addPartitionData();
        size_t pos = name.find_last_of( '/' );
        if( pos != name.npos )
          return vtkWriter_->pwrite( name.substr( pos+1, name.npos ), name.substr( 0, pos ), "", type );
        else
          return vtkWriter_->write( name, type );
      }

      std::string write ( const std::string &name )
      {
        return write( name, type_ );
      }

      std::string pwrite ( const std::string &name,
                           const std::string &path,
                           const std::string &extendpath,
                           VTK::OutputType type )
      {
        addPartitionData();
        return vtkWriter_->pwrite( name, path, extendpath, type );
      }

      std::string pwrite ( const std::string &name,
                           const std::string &path,
                           const std::string &extendpath )
      {
        return pwrite( name, path, extendpath, type_ );
      }

      std::string write ( const std::string &name,
                          VTK::OutputType type,
                          const int rank,
                          const int size )
      {
        addPartitionData( rank );
        return vtkWriter_->write( name, type, rank, size );
      }

      std::string write ( const std::string &name,
                          const int rank,
                          const int size )
      {
        return write( name, type_, rank, size );
      }

    protected:
      const GridPartType &gridPart_;
      VTKWriterType *vtkWriter_;
      int addPartition_;
      VTK::OutputType type_;
    };



    // VTKIOBase::VTKWriter
    // --------------------

    template< class GridPart, bool subsampling >
    class VTKIOBase< GridPart, subsampling >::VTKWriter
    : public Dune::VTKWriter< GridViewType >
    {
      typedef Dune::VTKWriter< GridViewType > BaseType;

    public:
      // make all write methods public for data convert
      using BaseType::write;
      using BaseType::pwrite;

      //! constructor
      VTKWriter( const GridPartType &gridPart,
                 VTK::DataMode dm = VTK::conforming )
      : BaseType( gridPart, dm )
      {}
    };



    // VTKIOBase::SubSamplingVTKWriter
    // -------------------------------

    template< class GridPart, bool subsampling >
    class VTKIOBase< GridPart, subsampling >::SubsamplingVTKWriter
      : public Dune::SubsamplingVTKWriter< GridViewType >
    {
      typedef Dune::SubsamplingVTKWriter< GridViewType > BaseType;

    public:
      // make all write methods public for data convert
      using BaseType::write;
      using BaseType::pwrite;

      //! constructor
      SubsamplingVTKWriter( const GridPartType &gridPart,
                            Dune::RefinementIntervals intervals,
                            bool coerceToSimplex = false )
        : BaseType( gridPart, intervals, coerceToSimplex )
      {}
    };



    // VTKIO (without subsampling)
    // ---------------------------

    template< class GridPart >
    class VTKIO< GridPart, false >
      : public VTKIOBase< GridPart, false >
    {
      typedef VTKIO< GridPart > ThisType;
      typedef VTKIOBase< GridPart, false > BaseType;

      typedef typename BaseType::VTKWriterType VTKWriterType;

    public:
      typedef typename BaseType::GridPartType GridPartType;

      VTKIO ( const GridPartType &gridPart, VTK::DataMode dm, const ParameterReader &parameter = Parameter::container() )
        : BaseType( gridPart, new VTKWriterType( gridPart, dm ), parameter )
      {}

      explicit VTKIO ( const GridPartType &gridPart, const ParameterReader &parameter = Parameter::container() )
        : VTKIO( gridPart, VTK::conforming, parameter )
      {}
    };



    // VTKIO (with subsampling)
    // ------------------------

    template< class GridPart >
    class VTKIO< GridPart, true >
      : public VTKIOBase< GridPart , true >
    {
      typedef VTKIO< GridPart, true > ThisType;
      typedef VTKIOBase< GridPart, true > BaseType;

       typedef typename BaseType::VTKWriterType VTKWriterType;

    public:
      typedef typename BaseType::GridPartType GridPartType;

      explicit VTKIO ( const GridPartType &gridPart, Dune::RefinementIntervals intervals, bool coerceToSimplex, const ParameterReader &parameter = Parameter::container() )
        : BaseType( gridPart, new VTKWriterType( gridPart, intervals, coerceToSimplex ), parameter )
      {}

      explicit VTKIO ( const GridPartType &gridPart, unsigned int level, bool coerceToSimplex, const ParameterReader &parameter = Parameter::container() )
        : VTKIO( gridPart, Dune::refinementLevels(level), coerceToSimplex, parameter )
      {}

      VTKIO ( const GridPartType &gridPart, unsigned int level, const ParameterReader &parameter = Parameter::container() )
        : VTKIO( gridPart, level, false, parameter )
      {}

      [[deprecated( "pass level as unsigned int" )]]
      VTKIO ( const GridPartType &gridPart, int level, const ParameterReader &parameter = Parameter::container() )
        : VTKIO( gridPart, level, false, parameter )
      {}

      VTKIO ( const GridPartType &gridPart, const ParameterReader &parameter = Parameter::container() )
        : VTKIO( gridPart, 0, false, parameter )
      {}

      VTKIO ( const GridPartType &gridPart, bool coerceToSimplex, const ParameterReader &parameter = Parameter::container() )
        : VTKIO( gridPart, 0, coerceToSimplex, parameter )
      {}
    };



    // SubsamplingVTKIO
    // ----------------

    template< class GridPart >
    using SubsamplingVTKIO = VTKIO< GridPart, true >;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_VTKIO_HH
