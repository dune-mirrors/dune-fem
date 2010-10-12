#ifndef DUNE_DATAOUTPUT_HH
#define DUNE_DATAOUTPUT_HH

#ifndef USE_VTKWRITER
#define USE_VTKWRITER 1
#endif

//- Dune includes
#include <dune/fem/io/file/iointerface.hh>
#include <dune/fem/io/file/iotuple.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/space/common/loadbalancer.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

#if USE_VTKWRITER
#include <dune/fem/io/file/vtkio.hh>
#include <dune/fem/operator/projection/vtxprojection.hh>
#endif

#ifndef USE_GRAPE 
// define whether to use grape of not 
#define USE_GRAPE HAVE_GRAPE
#endif

#if USE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

namespace Dune
{

  /** \brief Parameter class for Dune::DataOutput

      Structure providing the main parameters used to setup the Dune::DataOutput.
      By default these parameters are set through the Dune::Parameter class, i.e.,
      can be set in a parameter file. To override this a user can derive from this
      structure and modify any subset of the parameters. An instance of this modified
      class can then be passed in the construction of the Dune::DataOutput.
   */
  struct DataOutputParameters
  : public LocalParameter< DataOutputParameters, DataOutputParameters >
  {
    //! path where the data is stored (path are always relative to fem.commonOutputPath)
    virtual std::string path() const
    {
      return Parameter::getValue< std::string >( "fem.io.path", "" );
    }

    //! base of file name for data file (fem.io.datafileprefix)
    virtual std::string prefix () const
    {
      return Parameter::getValue< std::string >( "fem.io.datafileprefix" );
    }

    //! format of output (fem.io.outputformat)
    virtual int outputformat () const
    {
      static const std::string formatTable[]
        = { "binary", "vtk-cell", "vtk-vertex", "gnuplot" , "sub-vtk-cell", "none" };
      int format = Parameter::getEnum( "fem.io.outputformat", formatTable, 1 );
      return format;
    }

    //! use online grape display (fem.io.grapedisplay)
    virtual bool grapedisplay () const
    {
#if USE_GRAPE
      return (Parameter::getValue( "fem.io.grapedisplay", 0 ) == 1);
#else
      return false;
#endif
    }

    //! save data every savestep interval (fem.io.savestep)
    virtual double savestep () const
    {
      return Parameter::getValue< double >( "fem.io.savestep", 0 );
    }

    //! save data every savecount calls to write method (fem.io.savecount)
    virtual int savecount () const
    {
      return Parameter::getValue< int >( "fem.io.savecount", 0 );
    }

    //! save data every savecount calls to write method (fem.io.savecount)
    virtual int subsamplingLevel() const
    {
      return Parameter::getValue< int >( "fem.io.subsamplinglevel", 1 );
    }

    //! number for first data file (no parameter available)
    virtual int startcounter () const
    {
      return 0;
    }

    //! method used for conditional data output - default
    //! value passed as argument.
    virtual bool willWrite ( bool write ) const
    {
      return write;
    }
  };



  /** @ingroup DiscFuncIO 
     \brief Implementation of the Dune::IOInterface. 
     This class manages data output.
     Available output formats are GRAPE, VTK and VTK Vertex projected
     using the VtxProjection operator. Details can be
     found in \ref DiscFuncIO.
  */    
  template< class GridImp, class DataImp >
  class DataOutput
  : public IOInterface
  {
    typedef DataOutput< GridImp, DataImp > ThisType;

  protected:  
  #if USE_VTKWRITER
    template< class Grid, class OutputTuple >
    struct GridPartGetter;

    template< class VTKIOType >
    struct VTKListEntry 
    {
      virtual ~VTKListEntry () {}
      virtual void add ( VTKIOType & ) const = 0;

    protected: 
      VTKListEntry () {}
    };

    template< class VTKIOType, class DFType >
    class VTKFunc;

    template< class VTKOut >
    struct VTKOutputerDG;

    template< class VTKOut >
    struct VTKOutputerLagrange;
  #endif

    template< class GridPartType >
    class GnuplotOutputer;

  protected:  
    enum OutputFormat { binary = 0, vtk = 1, vtkvtx = 2, gnuplot = 3, subvtk = 4 , none = 5 };

    //! \brief type of grid used 
    typedef GridImp GridType;
    //! \brief type of data tuple 
    typedef DataImp OutPutDataType; 

  public: 
   /** \brief Constructor creating data output class 
      \param grid corresponding grid 
      \param data Tuple containing discrete functions to write 
      \param parameter structure for tuning the behavior of the Dune::DataOutput 
                       defaults to Dune::DataOutputParameters
    */
    DataOutput ( const GridType &grid, OutPutDataType &data,
                 const DataOutputParameters &parameter = DataOutputParameters() );

    /** \brief Constructor creating data writer 
      \param grid corresponding grid 
      \param data Tuple containing discrete functions to write 
      \param tp   a time provider to set time (e.g. for restart)
      \param parameter structure for tuning the behavior of the Dune::DataOutput 
                       defaults to Dune::DataOutputParameters
    */ 
    DataOutput ( const GridType &grid, OutPutDataType &data,
                 const TimeProviderBase &tp,
                 const DataOutputParameters &parameter = DataOutputParameters() );

    void consistentSaveStep ( const TimeProviderBase &tp ) const;

    //! destructor 
    virtual ~DataOutput() { delete param_; }

  protected:  
    //! initialize data writer 
    void init ( const DataOutputParameters &parameter );

  public:
    /** \brief returns true if data will be written on next write call
    */
    virtual bool willWrite ( const TimeProviderBase &tp ) const
    {
      // make save step consistent
      consistentSaveStep( tp );

      const bool writeStep  = (saveStep_ > 0) && (tp.time() >= saveTime_);
      const bool writeCount = (saveCount_ > 0) && (writeCalls_ % saveCount_ == 0);
      return param_->willWrite( writeStep || writeCount );
    }

    /** \brief returns true if data will be written on next write call
    */
    virtual bool willWrite () const 
    {
      return param_->willWrite( (saveCount_ > 0) && (writeCalls_ % saveCount_ == 0) );
    }

    /** \brief write given data to disc, evaluates parameter savecount
    */
    void write( const std::string& outstring="" ) const
    {
      if( willWrite() )
        writeData( writeCalls_, outstring );
      ++writeCalls_;    
    }

    /** \brief write given data to disc, evaluates parameter savecount and savestep
    */
    void write(const TimeProviderBase& tp, const std::string& outstring="") const
    {
      if( willWrite(tp) )
        writeData( tp.time(), outstring );
      ++writeCalls_;    
    }

    void writeData ( double sequenceStamp, const std::string &outstring = "" ) const;

    //! print class name 
    virtual const char* myClassName () const { return "DataOutput"; }

    /** Return output path name */
    const std::string &path () const { return path_; }

  protected:  
  #if USE_VTKWRITER
    std::string writeVTKOutput () const;
  #endif

    std::string writeGnuPlotOutput () const;

    //! write binary data 
    virtual void writeBinaryData ( const double ) const
    {
      DUNE_THROW(NotImplemented, "Format 'binary' not supported by DataOutput, use DataWriter instead." );
    }

    //! display data with grape 
    virtual void display () const 
    {
      if( grapeDisplay_ )
        grapeDisplay( data_ );
    }

  protected:  
    void grapeDisplay ( Dune::Nil &data ) const {}

    //! display data with grape 
    template< class OutputTupleType >
    void grapeDisplay ( OutputTupleType &data ) const;

    //! \brief type of this class 
    const GridType& grid_;
    mutable OutPutDataType& data_; 
    const int myRank_;

    // name for data output 
    std::string path_;
    std::string datapref_;
    // use grape display
    bool grapeDisplay_;
    
    // counter for sequence of files
    mutable int writeStep_;
    // counter for call to write
    mutable int writeCalls_;
    // next point in time to save data
    mutable double saveTime_;
    // time interval to save files
    double saveStep_;
    // number of write call between writting files
    int saveCount_;
    // grape, vtk or ...
    OutputFormat outputFormat_;
    mutable std::ofstream sequence_;
    const DataOutputParameters* param_;
  }; // end class DataOutput



  // DataOutput::GridPartGetter
  // --------------------------

#if USE_VTKWRITER
  template< class GridImp, class DataImp >
  template< class Grid, class OutputTuple >
  struct DataOutput< GridImp, DataImp >::GridPartGetter
  {
    typedef typename TypeTraits< typename Dune::ElementType< 0, OutputTuple >::Type >::PointeeType DFType;
    typedef typename DFType :: DiscreteFunctionSpaceType :: GridPartType GridPartType;

    GridPartGetter ( const Grid &, const OutputTuple &data )
    : gridPart_( getGridPart( data ) )
    {}

    const GridPartType &gridPart () const { return gridPart_; }

  protected:
    static const GridPartType &getGridPart( const OutputTuple& data )
    {
      const DFType *df = Element< 0 >::get( data );
      assert( df );
      return df->space().gridPart();
    }

    const GridPartType &gridPart_;
  };


  template< class GridImp, class DataImp >
  template< class Grid >
  struct DataOutput< GridImp, DataImp >::GridPartGetter< Grid, Dune::Nil >
  {
    typedef HierarchicGridPart< Grid > GridPartType;

    GridPartGetter ( const Grid &grid, const Dune::Nil & )
    : gridPart_( const_cast< Grid & >( grid ) )
    {}

    const GridPartType &gridPart () const { return gridPart_; }

  protected:
    const GridPartType gridPart_;
  };
#endif // #if USE_VTKWRITER



  // DataOutput::VTKFunc
  // -------------------

#if USE_VTKWRITER
  template< class GridImp, class DataImp >
  template< class VTKIOType, class DFType >
  class DataOutput< GridImp, DataImp >::VTKFunc
  : public VTKListEntry< VTKIOType >
  {
    typedef typename VTKIOType :: GridPartType GridPartType ;
    typedef typename DFType :: DiscreteFunctionSpaceType :: FunctionSpaceType FunctionSpaceType;

    typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType,
            GridPartType, 1> LagrangeSpaceType;
    typedef AdaptiveDiscreteFunction< LagrangeSpaceType >  NewFunctionType;

    const DFType& df_;
    LagrangeSpaceType space_;
    mutable NewFunctionType* func_;
    
    VTKFunc( const VTKFunc&) ; 

  public:
    VTKFunc( const GridPartType& gridPart, const DFType& df ) 
      : df_ (df) 
      , space_( const_cast<GridPartType&> (gridPart) )
      , func_(0)
    {
      // space_.setDescription(df.space().getDescription());
    }

    ~VTKFunc () 
    {
      delete func_;
    }

    virtual void add(VTKIOType& vtkio) const 
    {
      if( df_.space().continuous() ) 
      {
        vtkio.addVertexData( df_ );
      }
      else 
      {
        func_ = new NewFunctionType ( df_.name()+"vtx-prj" , space_ );
        WeightDefault<GridPartType> weight;
        VtxProjectionImpl::project( df_, *func_, weight );
        vtkio.addVertexData( *func_ );
      }
      // vtkio.addVectorVertexData( *func_ );
    }
  };
#endif // #if USE_VTKWRITER



  // DataOutput::VTKOutputerDG
  // -------------------------

#if USE_VTKWRITER
  template< class GridImp, class DataImp >
  template< class VTKOut >
  struct DataOutput< GridImp, DataImp >::VTKOutputerDG
  {
    //! Constructor
    explicit VTKOutputerDG ( VTKOut &vtkOut )
    : vtkOut_( vtkOut )
    {}

    //! Applies the setting on every DiscreteFunction/LocalFunction pair.
    template< class DFType >
    void visit ( DFType *df )
    {
      if(df->space().order() > 0)
      {
        vtkOut_.addVertexData( *df );
        // vtkOut_.addVectorVertexData( *df );
      }
      else 
      {
        vtkOut_.addCellData( *df ); 
        // vtkOut_.addVectorCellData( *df );
      }
    }

    void forEach(Dune::Nil& ) {} 

    template <class OutputTuple>
    void forEach(OutputTuple& data) 
    {
      ForEachValue<OutputTuple> forEach(data); 
      forEach.apply( *this );
    }

  private:
    VTKOut &vtkOut_;
  };
#endif // #if USE_VTKWRITER



  // DataOutput::VTKOutputerLagrange
  // -------------------------------

#if USE_VTKWRITER
  template< class GridImp, class DataImp >
  template< class VTKOut >
  struct DataOutput< GridImp, DataImp >::VTKOutputerLagrange
  {
    //! Constructor
    explicit VTKOutputerLagrange ( VTKOut &vtkOut )
    : vtkOut_( vtkOut ), vec_()
    {}

    //! destructor, delete entries of list 
    ~VTKOutputerLagrange () 
    {
      for(size_t i=0; i<vec_.size(); ++i)
      {
        delete vec_[i];
        vec_[i] = 0;
      }
    }

    //! Applies the setting on every DiscreteFunction/LocalFunction pair.
    template< class DFType >
    void visit ( DFType *df )
    {
      typedef VTKFunc<VTKOut,DFType> EntryType; 
      EntryType* entry = new EntryType(vtkOut_.gridPart(), *df);
      entry->add( vtkOut_ );
      vec_.push_back( entry );
    }

    void forEach(Dune::Nil& ) {} 

    template <class OutputTuple>
    void forEach(OutputTuple& data) 
    {
      ForEachValue<OutputTuple> forEach(data); 
      forEach.apply( *this );
    }

  private:
    VTKOut &vtkOut_;
    typedef VTKListEntry<VTKOut> VTKListEntryType;
    std::vector< VTKListEntryType * > vec_;
  };
#endif // #if USE_VTKWRITER



  // DataOutput::GnuplotOutputer
  // ---------------------------

  template< class GridImp, class DataImp >
  template< class GridPartType >
  class DataOutput< GridImp, DataImp >::GnuplotOutputer
  {
    typedef typename GridPartType::template Codim< 0 >::IteratorType::Entity Entity;
    std::ostream& out_;
    CachingQuadrature<GridPartType,0> &quad_;
    int i_;
    const Entity& en_;

  public:
    //! Constructor
    GnuplotOutputer ( std::ostream& out,
                      CachingQuadrature<GridPartType,0> &quad,
                      int i,
                      const Entity &en )
    : out_(out), quad_(quad), i_(i), en_(en)
    {}

    //! Applies the setting on every DiscreteFunction/LocalFunction pair.
    template <class DFType>
    void visit(DFType* df) 
    {
      typedef typename DFType :: Traits Traits;
      typedef typename Traits :: LocalFunctionType LocalFunctionType;
      typedef typename Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
      typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

      typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
      typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
 
      enum{ dimDomain = DiscreteFunctionSpaceType :: dimDomain };
      enum{ dimRange = DiscreteFunctionSpaceType :: dimRange };

      LocalFunctionType lf = df->localFunction(en_);
      {
        RangeType u;
        lf.evaluate( quad_[ i_ ], u );
        for( int k = 0; k < dimRange; ++k )
          out_ << "  " << u[ k ];
      }
    }

    void forEach(Dune::Nil& ) {} 

    template <class OutputTuple>
    void forEach(OutputTuple& data) 
    {
      ForEachValue<OutputTuple> forEach(data); 
      forEach.apply( *this );
    }
  };



  // Implementation of DataOutput
  // ----------------------------

  template< class GridImp, class DataImp >
  inline DataOutput< GridImp, DataImp >
    ::DataOutput ( const GridType &grid, OutPutDataType &data,
                  const DataOutputParameters& parameter )
  : grid_( grid ),
    data_( data ), 
    myRank_(grid_.comm().rank()),
    writeStep_(0),
    writeCalls_(0),
    saveTime_(0.0),  // why 0.0?
    saveStep_(-1),
    saveCount_(-1),
    outputFormat_(vtkvtx),
    param_(parameter.clone())
  {
    // initialize class 
    init( parameter );
  }


  template< class GridImp, class DataImp >
  inline DataOutput< GridImp, DataImp >
    ::DataOutput ( const GridType &grid, OutPutDataType &data,
                   const TimeProviderBase &tp,
                   const DataOutputParameters &parameter )
    : grid_(grid),
      data_(data),
      myRank_(grid_.comm().rank()),
      writeStep_(0),
      writeCalls_(0),
      saveTime_(0),
      saveStep_(-1),
      saveCount_(-1),
      outputFormat_(vtkvtx),
      param_(parameter.clone())
  {
    // initialize class 
    init( parameter );

    // make save step consistent
    consistentSaveStep( tp );
  }


  template< class GridImp, class DataImp >
  inline void DataOutput< GridImp, DataImp >
    ::consistentSaveStep ( const TimeProviderBase &tp ) const
  {
    const double oldTime = tp.time() - saveStep_;
    // set old values according to new time 
    if( saveStep_ > 0 )
    {
      while( saveTime_ <= oldTime ) 
      {
        ++writeStep_;
        saveTime_ += saveStep_;
      }
    }
  }


  template< class GridImp, class DataImp >
  inline void DataOutput< GridImp, DataImp >::init ( const DataOutputParameters &parameter )
  {
    IOInterface :: createGlobalPath ( grid_.comm(), Parameter::commonOutputPath() );
    path_ = Parameter::commonOutputPath()+"/"+parameter.path();
    // create path if not already exists 
    IOInterface :: createGlobalPath ( grid_.comm(), path_ );

    
    // add prefix for data file
    datapref_ += parameter.prefix();

    int outputFormat = parameter.outputformat();
    switch( outputFormat ) 
    {
      case 0: outputFormat_ = binary; break;
      case 1: outputFormat_ = vtk; break;
      case 2: outputFormat_ = vtkvtx; break;
      case 3: outputFormat_ = gnuplot; break;
      case 4: outputFormat_ = subvtk; break;
      case 5: outputFormat_ = none; break;
      default:
        DUNE_THROW(NotImplemented,"DataOutput::init: wrong output format");
    }

    grapeDisplay_ = parameter.grapedisplay();

    // get parameters for data writing
    saveStep_ = parameter.savestep();
    saveCount_ = parameter.savecount();

    writeStep_ = parameter.startcounter();

    if (myRank_==0) {
      std::string name = path_ + "/" + datapref_;
      name += ".series";
      std::cout << "opening file: " << name << std::endl;
      sequence_.open(name.c_str());
      if (!sequence_)
        std::cout << "could not write sequence file" << std::endl;
    }

    // write parameter file 
    Parameter::write("parameter.log");
  }


  template< class GridImp, class DataImp >
  inline void DataOutput< GridImp, DataImp >
    ::writeData ( double sequenceStamp, const std::string &outstring ) const
  {
    std::string filename;
    // check online display 
    display(); 
    switch( outputFormat_ )
    {
    // if no output was chosen just return  
    case none: return ;
    case binary:   
      writeBinaryData( sequenceStamp );
      break;
#if USE_VTKWRITER
    case vtk : 
    case vtkvtx :
    case subvtk :
      // write data in vtk output format 
      filename = writeVTKOutput();
      break;
#endif // #if USE_VTKWRITER
    case gnuplot :
      filename = writeGnuPlotOutput();
      break;
    default:
      DUNE_THROW(NotImplemented,"DataOutput::write: wrong output format");
    }

    if (sequence_)
      sequence_ << writeStep_ << " "
                << filename << " "
                << sequenceStamp
                << outstring
                << std::endl;

    if(myRank_ <= 0)
    {
      // only write info for proc 0, otherwise on large number of procs
      // this is to much output 
      std::cout << myClassName() << "[" << myRank_ << "]::write data"
                << " writestep=" << writeStep_
                << " sequenceStamp=" << sequenceStamp
                << outstring
                << std::endl;
    }

    saveTime_ += saveStep_;
    ++writeStep_;
  }


#if USE_VTKWRITER
  template< class GridImp, class DataImp >
  inline std::string DataOutput< GridImp, DataImp >::writeVTKOutput () const
  {
    std::string filename;
    // check whether to use vertex data of discontinuous data 
    const bool vertexData = (outputFormat_ == vtkvtx);

    // check whether we have parallel run  
    const bool parallel = (grid_.comm().size() > 1);

    // generate filename, with path only for serial run  
    std::string name = genFilename( (parallel) ? "" : path_, datapref_, writeStep_ );

    if( vertexData ) 
    {
      // generate adaptive leaf grid part 
      // do not use leaf grid part since this will 
      // create the grids leaf index set, which might not be wanted. 
      typedef AdaptiveLeafGridPart< GridType > GridPartType; 
      GridPartType gridPart( const_cast<GridType&> (grid_) );

      // create vtk output handler 
      typedef VTKIO < GridPartType > VTKIOType; 
      VTKIOType vtkio ( gridPart, VTKOptions::conforming );

      // add all functions 
      VTKOutputerLagrange< VTKIOType > io( vtkio );
      io.forEach( data_ );

      if( parallel )
      {
        // write all data for parallel runs  
        filename = vtkio.pwrite( name, path_, "." , Dune::VTKOptions::binaryappended );
      }
      else
      {
        // write all data serial 
        filename = vtkio.write( name, Dune::VTKOptions::binaryappended );
      }
    }
    else if ( outputFormat_ == vtk )
    {
      typedef GridPartGetter< GridType, OutPutDataType> GridPartGetterType;
      GridPartGetterType gp( grid_, data_ );

      // create vtk output handler 
      typedef VTKIO < typename GridPartGetterType :: GridPartType > VTKIOType; 
      VTKIOType vtkio ( gp.gridPart() , VTKOptions::nonconforming );

      // add all functions 
      VTKOutputerDG< VTKIOType > io( vtkio );
      io.forEach( data_ );

      // write all data 
      if( parallel )
      {
        // write all data for parallel runs  
        filename = vtkio.pwrite( name, path_, "." , Dune::VTKOptions::binaryappended );
      }
      else
      {
        // write all data serial 
        filename = vtkio.write( name, Dune::VTKOptions::binaryappended );
      }
    }
    else if ( outputFormat_ == subvtk )
    {
      typedef GridPartGetter< GridType, OutPutDataType> GridPartGetterType;
      GridPartGetterType gp( grid_, data_ );

      // create vtk output handler 
      typedef SubsamplingVTKIO < typename GridPartGetterType :: GridPartType > VTKIOType; 
      VTKIOType vtkio ( gp.gridPart(), param_->subsamplingLevel() );

      // add all functions 
      VTKOutputerDG< VTKIOType > io( vtkio );
      io.forEach( data_ );

      // write all data 
      if( parallel )
      {
        // write all data for parallel runs  
        filename = vtkio.pwrite( name, path_, "." , Dune::VTKOptions::binaryappended );
      }
      else
      {
        // write all data serial 
        filename = vtkio.write( name, Dune::VTKOptions::binaryappended );
      }
    }
    return filename;
  }
#endif // #if USE_VTKWRITER


  template< class GridImp, class DataImp >
  inline std::string DataOutput< GridImp, DataImp >::writeGnuPlotOutput () const
  {
    // generate filename
    std::string name = genFilename( path_, datapref_, writeStep_ );
    name += ".gnu";
    std::ofstream gnuout(name.c_str());
    gnuout << std::scientific << std::setprecision( 16 );

    // create HierarchicGridPart
    // ! do not use LeafGridPart since this creates the grid's LeafIndexSet !
    typedef HierarchicGridPart< GridType > GridPartType;
    GridPartType gridPart( const_cast< GridType & >( grid_ ) );

    // start iteration
    typedef typename GridPartType::template Codim<0>::IteratorType IteratorType;
    IteratorType end = gridPart.template end< 0 >();
    for( IteratorType it = gridPart.template begin< 0 >(); it != end; ++it )
    {
      CachingQuadrature< GridPartType, 0 > quad( *it, 1 );
      for( size_t i = 0; i < quad.nop(); ++i )
      {
        typedef typename IteratorType::Entity Entity;
        const Entity &entity = *it;
        const typename Entity::Geometry::GlobalCoordinate x
          = entity.geometry().global( quad.point( i ) );
        for( int k = 0; k < x.size; ++k )
          gnuout << (k > 0 ? " " : "") << x[ k ];
        GnuplotOutputer< GridPartType > io( gnuout, quad, i, entity );
        io.forEach( data_ );
        gnuout << std::endl;
      }
    }
    return name;
  }


#if USE_GRAPE
  template< class GridImp, class DataImp >
  template< class OutputTupleType >
  inline void DataOutput< GridImp, DataImp >::grapeDisplay ( OutputTupleType &data ) const
  {
    GrapeDataDisplay<GridType> grape(grid_);
    IOTuple<OutPutDataType>::addToDisplay(grape,data);
    grape.display();
  }
#else // #if USE_GRAPE
  template< class GridImp, class DataImp >
  template< class OutputTupleType >
  inline void DataOutput< GridImp, DataImp >::grapeDisplay ( OutputTupleType &data ) const
  {
    std::cerr << "WARNING: HAVE_GRAPE == 0, but grapeDisplay == true, recompile with grape support for online display!" << std::endl;
  }
#endif // #else // #if USE_GRAPE

} // end namespace Dune 

#endif // #ifndef DUNE_DATAOUTPUT_HH
