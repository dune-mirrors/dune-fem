#ifndef DUNE_FEM_VTKIO_HH
#define DUNE_FEM_VTKIO_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/fem/version.hh>
#include <dune/fem/misc/field.hh>
#include <dune/fem/gridpart/common/gridpartview.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem/misc/threads/threadmanager.hh>
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
    : public VTKFunction< typename DF::DiscreteFunctionSpaceType::GridPartType::GridViewType >
    {
      typedef VTKFunctionWrapper< DF > ThisType;

    public:
      typedef DF DiscreteFunctionType;

      typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
        DiscreteFunctionSpaceType;

      static const int dimRange = DiscreteFunctionSpaceType::dimRange;
      static const int dimDomain = DiscreteFunctionSpaceType::dimDomain;

      typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
      typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
      typedef typename DiscreteFunctionSpaceType::GridPartType::template Codim< 0 >::EntityType EntityType;

      typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;


      //! constructor taking discrete function 
      VTKFunctionWrapper ( const DiscreteFunctionType& df,
                           const std::string& dataName,
                           int component, bool vector )
      : discFunc_( df ),
        name_( ( dataName.size() > 0 ) ? dataName : df.name() ),
        vector_( vector ),
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
        const LocalFunctionType lf = discFunc_.localFunction(e);
        RangeType val;
        lf.evaluate(xi,val);
        if (vector_)
          return (comp > dimDomain ? 0.0 : Dune::field_cast< double >( val[ comp + component_ ] ));
        else 
          return Dune::field_cast< double >( val[component_] );
      }

      //! get name
      virtual std::string name () const
      {
        if (vector_) {
          std::stringstream ret_vec;
          ret_vec << name_ << "_vec" << component_;
          return ret_vec.str(); 
        }
        std::stringstream ret;
        ret << name_ << component_;
        return ret.str();
      }

    private:
      const DiscreteFunctionType &discFunc_;
      const std::string name_ ;
      const bool vector_;
      const int component_;
    };



    // VTKIOBase
    // ---------

    //! /brief Output using VTK
    template< class GridPart, bool subsampling >
    class VTKIOBase
    {
      typedef VTKIOBase< GridPart, subsampling > ThisType;

    protected:
      typedef typename GridPart::GridViewType GridViewType;

      class VTKWriter;
      class SubsamplingVTKWriter;

      typedef typename SelectType< subsampling, SubsamplingVTKWriter, VTKWriter >::Type
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

      int getPartitionParameter() const 
      {
        // 0 = none, 1 = MPI ranks only, 2 = ranks + threads, 3 = like 1 and also threads only
        const std::string names[] = { "none", "rank", "rank+thread", "rank/thread" };
        return Parameter :: getEnum ("fem.io.partitioning", names, 0 );
      }

    protected :
      VTKIOBase ( const GridPartType &gridPart, VTKWriterType *vtkWriter )
      : gridPart_( gridPart ),
        vtkWriter_( vtkWriter ),
        addPartition_( getPartitionParameter() )
      {
      }

      void addPartitionData( const int myRank = -1 ) 
      {
        if( addPartition_ > 0 ) 
        {
          vtkWriter_->addCellData( new VolumeData() );

          const int rank = ( myRank < 0 ) ? gridPart_.grid().comm().rank() : myRank ;
          const int nThreads = ( addPartition_ > 1 ) ? ThreadManager::maxThreads() : 1 ;
          if( addPartition_ <= 2 ) 
            vtkWriter_->addCellData( new PartitioningData( gridPart_, "rank", rank, nThreads ) );
          else 
          {
            // rank only visualization 
            vtkWriter_->addCellData( new PartitioningData( gridPart_, "rank", rank, -1 ) );
            // thread only visualization 
            vtkWriter_->addCellData( new PartitioningData( gridPart_, "thread", 0, nThreads ) );
          }
          addPartition_ = 0 ;
        }
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
          vtkWriter_->addCellData( new VTKFunctionWrapper< DF >( df, dataName, i, false ) );
      }

      template< class DF >
      void addVectorCellData( DF &df, 
                              const std::string& dataName = "" , 
                              int startPoint = 0 )
      {
        vtkWriter_->addCellData( new VTKFunctionWrapper< DF >( df, dataName, startPoint, true ) );
      }

      template< class DF >
      void addVertexData( DF &df, const std::string& dataName = "" )
      {
        static const int dimRange = DF::FunctionSpaceType::dimRange;
        std::string name = ( dataName.size() > 0 ) ? dataName : df.name() ;
        for( int i = 0;i < dimRange; ++i )
          vtkWriter_->addVertexData( new VTKFunctionWrapper< DF >( df, dataName, i, false ) );
      }

      template< class DF >
      void addVectorVertexData( DF &df, 
                                const std::string& dataName = "" , 
                                int startPoint = 0 )
      {
        vtkWriter_->addVertexData( new VTKFunctionWrapper< DF >( df, dataName, startPoint, true ) );
      }

      void clear ()
      {
        vtkWriter_->clear();
      }

      std::string write ( const std::string &name, VTK::OutputType type = VTK::ascii )
      {
        addPartitionData();
        size_t pos = name.find_last_of( '/' );
        if( pos != name.npos )
          return vtkWriter_->pwrite( name.substr( pos+1, name.npos ), name.substr( 0, pos ), "", type );
        else
          return vtkWriter_->write( name, type );
      }

      std::string pwrite ( const std::string &name,
                           const std::string &path,
                           const std::string &extendpath,
                           VTK::OutputType type = VTK::ascii )
      {
        addPartitionData();
        return vtkWriter_->pwrite( name, path, extendpath, type );
      }

      std::string write ( const std::string &name,
                          VTK::OutputType type,
                          const int rank, 
                          const int size )
      {
        addPartitionData( rank );
        return vtkWriter_->write( name, type, rank, size );
      }

    protected:
      const GridPartType &gridPart_;
      VTKWriterType *vtkWriter_;
      int addPartition_;
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
      : BaseType( gridPart.gridView(), dm )
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
                            const int level, 
                            bool coerceToSimplex = false )
      : BaseType( gridPart.gridView(), level, coerceToSimplex )
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

      explicit VTKIO ( const GridPartType &gridPart,
                       VTK::DataMode dm = VTK::conforming )
      : BaseType( gridPart, new VTKWriterType( gridPart, dm ) )
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

      explicit VTKIO ( const GridPartType &gridPart,
                       unsigned int level = 0,
                       bool coerceToSimplex = false )
      : BaseType( gridPart, new VTKWriterType( gridPart, level, coerceToSimplex ) )
      {}
    };



    // SubsamplingVTKIO
    // ----------------

    template< class GridPart >
    class SubsamplingVTKIO
    : public VTKIO< GridPart, true >
    {
      typedef SubsamplingVTKIO< GridPart > ThisType;
      typedef VTKIO< GridPart, true > BaseType;

    public:
      typedef typename BaseType::GridPartType GridPartType;

      explicit SubsamplingVTKIO ( const GridPartType &gridPart,
                                  unsigned int level = 0,
                                  bool coerceToSimplex = false )
      : BaseType( gridPart, level, coerceToSimplex )
      {}
    };

  } // namespace Fem

#if DUNE_FEM_COMPATIBILITY  
  // put this in next version 1.4 

  using Fem :: VTKIO  ;
  using Fem :: SubsamplingVTKIO  ;
#endif // DUNE_FEM_COMPATIBILITY

} // namespace Dune

#endif // #ifndef DUNE_FEM_VTKIO_HH
