/**************************************************************************

  The dune-fem module is a module of DUNE (see www.dune-project.org).
  It is based on the dune-grid interface library
  extending the grid interface by a number of discretization algorithms
  for solving non-linear systems of partial differential equations.

  Copyright (C) 2003 - 2015 Robert Kloefkorn
  Copyright (C) 2003 - 2010 Mario Ohlberger
  Copyright (C) 2004 - 2015 Andreas Dedner
  Copyright (C) 2005        Adrian Burri
  Copyright (C) 2005 - 2015 Mirko Kraenkel
  Copyright (C) 2006 - 2015 Christoph Gersbacher
  Copyright (C) 2006 - 2015 Martin Nolte
  Copyright (C) 2011 - 2015 Tobias Malkmus
  Copyright (C) 2012 - 2015 Stefan Girke
  Copyright (C) 2013 - 2015 Claus-Justus Heine
  Copyright (C) 2013 - 2014 Janick Gerstenberger
  Copyright (C) 2013        Sven Kaulman
  Copyright (C) 2013        Tom Ranner
  Copyright (C) 2015        Marco Agnese
  Copyright (C) 2015        Martin Alkaemper


  The dune-fem module is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The dune-fem module is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

**************************************************************************/
#ifndef DUNE_DIRICHLETCONSTRAINTS_HH
#define DUNE_DIRICHLETCONSTRAINTS_HH

#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/operator/common/temporarylocalmatrix.hh>
#include <dune/fem/common/bindguard.hh>
#include <dune/fem/common/coordinate.hh>
#include <dune/fem/common/localcontribution.hh>
#include <dune/fem/space/common/localinterpolation.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>

namespace Dune {

namespace Fem {
class HasLocalFunction;
}

template < class Model, class DiscreteFunctionSpace >
class DirichletConstraints
{
public:
  enum Operation { set = 0, sub = 1 };
  typedef Model ModelType;
  typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename DiscreteFunctionSpaceType::HessianRangeType HessianRangeType;

  //! type of grid partition
  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  //! type of grid
  typedef typename DiscreteFunctionSpaceType :: GridType GridType;

  // types for boundary treatment
  // ----------------------------
  static const int localBlockSize = DiscreteFunctionSpaceType :: localBlockSize ;
  static const int dimRange = DiscreteFunctionSpaceType::FunctionSpaceType::dimRange;
  static_assert( localBlockSize <= dimRange,
       "local block size of the space must be less than or equal to the dimension of the range of the function space.");
  typedef std::array<int,localBlockSize> DirichletBlock;
  typedef std::vector< DirichletBlock > DirichletBlockVector;

  class BoundaryWrapper
  {
    public:
    typedef typename DiscreteFunctionSpaceType::EntityType EntityType;
    typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpace::RangeType RangeType;
    typedef typename DiscreteFunctionSpace::JacobianRangeType JacobianRangeType;
    typedef typename DiscreteFunctionSpace::HessianRangeType HessianRangeType;

    private:
    const ModelType& impl_;
    const EntityType& entity_;
    const int order_;
    int bndId_;

    public:
    static const int dimRange = RangeType::dimension;
    BoundaryWrapper( const ModelType& impl, const EntityType& entity, const int order, int bndId )
    : impl_( impl ), entity_(entity), order_(order), bndId_(bndId) {}
    const EntityType& entity() const { return entity_; }
    const int order () const { return order_; }
    template <class Point>
    void evaluate( const Point& x, RangeType& ret ) const
    { impl_.dirichlet(bndId_,Dune::Fem::coordinate(x),ret); }
    template <class Point>
    void jacobian( const Point& x, JacobianRangeType& ret ) const
    { DUNE_THROW(Dune::NotImplemented,"rhs jacobian not implemented"); }
  };

  DirichletConstraints( ModelType &model, const DiscreteFunctionSpaceType& space )
    : model_(model),
      space_( space ),
      dirichletBlocks_(),
      // mark DoFs on the Dirichlet boundary
      hasDirichletDofs_( false ),
      sequence_( -1 )
  {}

  /*! treatment of Dirichlet-const DoFs for given discrete function
   *
   *   \note A LagrangeDiscreteFunctionSpace is implicitly assumed.
   *
   *   \param[in]  u   discrete function providing the constraints
   *   \param[out] w   discrete function the constraints are applied to
   */
  template < class DiscreteFunctionType >
  void operator ()( const DiscreteFunctionType& u, DiscreteFunctionType& w ) const
  {
    updateDirichletDofs();

    // if Dirichlet Dofs have been found, treat them
    if( hasDirichletDofs_ )
    {
      typedef typename DiscreteFunctionType :: DofIteratorType DofIteratorType ;
      typedef typename DiscreteFunctionType :: ConstDofIteratorType ConstDofIteratorType ;

      ConstDofIteratorType uIt = u.dbegin();
      DofIteratorType wIt = w.dbegin();

      // loop over all blocks
      const unsigned int blocks = space_.blockMapper().size();
      for( unsigned int blockDof = 0; blockDof < blocks ; ++ blockDof )
      {
        for( int l = 0; l < localBlockSize ; ++ l, ++ wIt, ++ uIt )
        {
          if( dirichletBlocks_[ blockDof ][l] )
          {
            // copy dofs of the block
            assert( uIt != u.dend() );
            assert( wIt != w.dend() );
            (*wIt) = (*uIt);
          }
        }
      }
    }
  }
  /*! treatment of Dirichlet-DoFs for given discrete function
   *
   *   \note A LagrangeDiscreteFunctionSpace is implicitly assumed.
   *
   *   \param[in]  value   a range vector
   *   \param[out] w       discrete function the constraints are applied to
   */
  template < class DiscreteFunctionType >
  void operator ()( const typename DiscreteFunctionType::RangeType& value, DiscreteFunctionType& w ) const
  {
    updateDirichletDofs();

    // if Dirichlet Dofs have been found, treat them
    if( hasDirichletDofs_ )
    {
      typedef typename DiscreteFunctionType :: DofIteratorType DofIteratorType ;

      DofIteratorType wIt = w.dbegin();

      // loop over all blocks
      const unsigned int blocks = space_.blockMapper().size();
      for( unsigned int blockDof = 0; blockDof < blocks ; ++ blockDof )
      {
        for( int l = 0; l < localBlockSize ; ++ l, ++ wIt )
        {
          if( dirichletBlocks_[ blockDof ][l] )
          {
            // copy dofs of the block
            assert( wIt != w.dend() );
            (*wIt) = value[l];
          }
        }
      }
    }
  }

  /*! treatment of Dirichlet-DoFs for given discrete function
   *
   *   \note A LagrangeDiscreteFunctionSpace is implicitly assumed.
   *
   *   \param[in]  u   discrete function providing the constraints
   *   \param[out] w   discrete function the constraints are applied to
   */
  template < class DiscreteFunctionType >
  void operator ()( DiscreteFunctionType& w ) const
  {
    updateDirichletDofs();

    if( hasDirichletDofs_ )
    {
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
      typedef typename IteratorType :: Entity EntityType;

      Dune::Fem::SetSelectedLocalContribution< DiscreteFunctionType > wLocal( w );
      Dune::Fem::LocalInterpolation< DiscreteFunctionSpaceType > interpolation( space_ );
      for( const EntityType &entity : space_ )
      {
        model_.init(entity);

        // bind local contribution to entity
        auto wGuard = Dune::Fem::bindGuard( wLocal, entity );
        // bind interpolation to entity
        auto iGuard = bindGuard( interpolation, entity );

        // interpolate dirichlet dofs
        dirichletDofTreatment( interpolation, wLocal );
      }
    }
  }
  template < class GridFunctionType, class DiscreteFunctionType,
           typename = std::enable_if_t< std::is_base_of<Dune::Fem::HasLocalFunction, GridFunctionType>::value > >
  void operator ()( const GridFunctionType &u,
                    DiscreteFunctionType& w, Operation op=Operation::setDF ) const
  {
    updateDirichletDofs();

    if( hasDirichletDofs_ )
    {
      Dune::Fem::ConstLocalFunction< GridFunctionType > uLocal( u );
      Dune::Fem::SetSelectedLocalContribution< DiscreteFunctionType > wLocal( w );
      Dune::Fem::LocalInterpolation< DiscreteFunctionSpaceType > interpolation( space_ );

      for( const auto &entity : space_ )
      {
        auto guard = Dune::Fem::bindGuard( std::tie(uLocal,wLocal), entity );
        // bind interpolation to entity
        auto iGuard = bindGuard( interpolation, entity );

        // interpolate dirichlet dofs
        if (op == Operation::sub)
          model_.init(entity);
        dirichletDofTreatment( interpolation, uLocal, wLocal, op );
      }
    }
  }

  /*! treatment of Dirichlet-DoFs for solution and right-hand-side
   *
   *   delete rows for dirichlet-DoFs, setting diagonal element to 1.
   *
   *   \note A LagrangeDiscreteFunctionSpace is implicitly assumed.
   *
   *   \param[out] linearOperator  linear operator to be adjusted
   */
  template <class LinearOperator>
  void applyToOperator( LinearOperator& linearOperator ) const
  {
    updateDirichletDofs();

    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
    typedef typename IteratorType :: Entity EntityType;

    // if Dirichlet Dofs have been found, treat them
    if( hasDirichletDofs_ )
    {
      typedef typename DiscreteFunctionSpaceType :: BlockMapperType BlockMapperType;
      Dune::Fem::NonBlockMapper< BlockMapperType, localBlockSize > mapper( space_.blockMapper() );

      std::vector<std::size_t> globalBlockDofs;
      std::vector<std::size_t> globalDofs;

      // storage for unit rows and auxiliary rows, should be unique and sorted
      std::set<std::size_t> unitRows;
      std::set<std::size_t> auxRows;

      const auto& space = space_; // linearOperator.rangeSpace();
      // get auxiliary dof structure (for parallel runs)   /*@LST0S@*/
      const auto &auxiliaryDofs = space.auxiliaryDofs();

      const IteratorType end = space_.end();
      for( IteratorType it = space_.begin(); it != end; ++it )
      {
        const EntityType &entity = *it;

        // get number of basis functions
        const int localBlocks = space.blockMapper().numDofs( entity );

        // map local to global dofs
        globalBlockDofs.resize(localBlocks);
        // obtain all DofBlocks for this element
        space.blockMapper().map( entity, globalBlockDofs );

        // obtain all non-blocked dofs
        globalDofs.resize(localBlocks * localBlockSize);
        mapper.map( entity, globalDofs );

        /*
        unitRows.reserve( globalDofs.size() );
        unitRows.clear();
        auxRows.reserve( globalDofs.size() );
        auxRows.clear();
        */

        // counter for all local dofs (i.e. localBlockDof * localBlockSize + ... )
        int localDof = 0;
        // iterate over face dofs and set unit row
        for( int localBlockDof = 0 ; localBlockDof < localBlocks; ++ localBlockDof )
        {
          int global = globalBlockDofs[localBlockDof];
          std::set<std::size_t>& rows = auxiliaryDofs.contains( global ) ? auxRows : unitRows;
          for( int l = 0; l < localBlockSize; ++ l, ++ localDof )
          {
            if( dirichletBlocks_[global][l] )
            {
              // push non-blocked dof
              rows.insert( globalDofs[ localDof ] );
            }
          }
        }
      } // end for elements

      // set unit and auxiliary rows at once
      linearOperator.setUnitRows( unitRows, auxRows );
    }
  }

  const DirichletBlockVector &dirichletBlocks() const
  {
    updateDirichletDofs();
    return dirichletBlocks_;
  }

protected:
  //! set the Dirichlet points to exact values
  template< class LocalInterpolationType, class LocalFunctionType >
  void dirichletDofTreatment( const LocalInterpolationType& interpolation, LocalFunctionType &wLocal ) const
  {
    // get entity
    const typename LocalFunctionType::EntityType &entity = wLocal.entity();

    // get number of Lagrange Points
    const int localBlocks = space_.blockMapper().numDofs( entity );

    // map local to global BlockDofs
    std::vector<std::size_t> globalBlockDofs(localBlocks);
    space_.blockMapper().map( entity, globalBlockDofs );
    std::vector<typename LocalFunctionType::RangeFieldType> values( localBlocks*localBlockSize );

    int localDof = 0;

    // iterate over face dofs and set unit row
    for( int localBlock = 0 ; localBlock < localBlocks; ++ localBlock )
    {
      // store result to dof vector
      int global = globalBlockDofs[localBlock];
      for( int l = 0; l < localBlockSize ; ++ l, ++localDof )
      {
        if( dirichletBlocks_[ global ][l] )
        {
          interpolation( BoundaryWrapper(model_, entity, wLocal.order(), dirichletBlocks_[global][l]), values );
          // store result
          assert( (unsigned int)localDof < wLocal.size() );
          wLocal[ localDof ] = values[ localDof ];
        }
      }
    }
  }

  template< class LocalInterpolationType, class GridLocalFunctionType, class LocalFunctionType >
  void dirichletDofTreatment( const LocalInterpolationType& interpolation,
                              const GridLocalFunctionType &uLocal,
                              LocalFunctionType &wLocal, Operation op ) const
  {
    // get entity
    const typename LocalFunctionType::EntityType &entity = wLocal.entity();

    // get number of Lagrange Points
    const int localBlocks = space_.blockMapper().numDofs( entity );

    typedef typename DiscreteFunctionSpaceType::BlockMapperType::GlobalKeyType  GlobalKeyType;

    // map local to global BlockDofs
    std::vector< GlobalKeyType > globalBlockDofs(localBlocks);
    space_.blockMapper().map(entity,globalBlockDofs);
    std::vector<double> values( localBlocks*localBlockSize );
    std::vector<double> valuesModel( localBlocks*localBlockSize );

    interpolation( uLocal, values );

    int localDof = 0;

    // iterate over face dofs and set unit row
    for( int localBlock = 0 ; localBlock < localBlocks; ++ localBlock )
    {
      // store result to dof vector
      int global = globalBlockDofs[localBlock];
      for( int l = 0; l < localBlockSize ; ++ l, ++localDof )
      {
        if( dirichletBlocks_[ global ][l] )
        {
          if (op == Operation::sub)
          {
            interpolation(BoundaryWrapper(model_, entity, wLocal.order(), dirichletBlocks_[global][l]), valuesModel);
            values[ localDof ] -= valuesModel[ localDof ];
          }
          // store result
          assert( (unsigned int)localDof < wLocal.size() );
          wLocal[ localDof ] = values[ localDof ];
        }
      }
    }
  }

protected:
  // detect all DoFs on the Dirichlet boundary
  void updateDirichletDofs() const
  {
    if( sequence_ != space_.sequence() )
    {
      // only start search if Dirichlet boundary is present
      if( ! model_.hasDirichletBoundary() )
      {
        hasDirichletDofs_ = false ;
        return ;
      }

      // resize flag vector with number of blocks and reset flags
      const int blocks = space_.blockMapper().size() ;
      dirichletBlocks_.resize( blocks );
      for( int i=0; i<blocks; ++i )
        dirichletBlocks_[ i ].fill(0);

      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
      typedef typename IteratorType :: Entity EntityType;

      bool hasDirichletBoundary = false;
      const IteratorType end = space_.end();
      for( IteratorType it = space_.begin(); it != end; ++it )
      {
        const EntityType &entity = *it;
        // if entity has boundary intersections
        if( entity.hasBoundaryIntersections() )
        {
          hasDirichletBoundary |= searchEntityDirichletDofs( entity, model_ );
        }
      }

      // update sequence number
      sequence_ = space_.sequence();
      if( space_.gridPart().comm().size() > 1 )
      {
        try
        {
          DirichletBuilder handle( *this, space_ , space_.blockMapper() );
          space_.gridPart().communicate
            ( handle, GridPartType::indexSetInterfaceType, ForwardCommunication );
        }
        // catch possible exceptions here to have a clue where it happend
        catch( const Exception &e )
        {
          std::cerr << e << std::endl;
          std::cerr << "Exception thrown in: " << __FILE__ << " line:" << __LINE__ << std::endl;
          abort();
        }
        hasDirichletDofs_ = space_.gridPart().grid().comm().max( hasDirichletBoundary );
      }
      else
      {
        hasDirichletDofs_ = hasDirichletBoundary;
      }
    }
  }

  // detect all DoFs on the Dirichlet boundary of the given entity
  template< class EntityType >
  bool searchEntityDirichletDofs( const EntityType &entity, ModelType& model ) const
  {
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

    typedef typename GridPartType :: IntersectionIteratorType
      IntersectionIteratorType;

    const GridPartType &gridPart = space_.gridPart();

    // default is false
    bool hasDirichletBoundary = false;

    //map local to global BlockDofs
    std::vector<size_t> globalBlockDofs(space_.blockMapper().numDofs(entity));
    space_.blockMapper().map(entity,globalBlockDofs);

    std::vector<bool> globalBlockDofsFilter(space_.blockMapper().numDofs(entity));

    IntersectionIteratorType it = gridPart.ibegin( entity );
    const IntersectionIteratorType endit = gridPart.iend( entity );
    for( ; it != endit; ++it )
    {
      typedef typename IntersectionIteratorType :: Intersection IntersectionType;
      const IntersectionType& intersection = *it;

      // if intersection is with boundary, adjust data
      if( intersection.boundary() )
      {
        // get dirichlet information from model
        DirichletBlock block;
        block.fill(0);
        std::array<int,dimRange> dblock;

        model.init(intersection);
        const bool isDirichletIntersection = model.isDirichletIntersection( intersection, dblock );
        if (isDirichletIntersection)
        {
          for ( unsigned int i=0;i<localBlockSize;++i)
            block[i] = dblock[i];
          // get face number of boundary intersection
          const int face = intersection.indexInInside();
          space_.blockMapper().onSubEntity(entity,face,1,globalBlockDofsFilter);
          for( unsigned int i=0;i<globalBlockDofs.size();++i)
          {
            assert( i < globalBlockDofsFilter.size());
            if ( !globalBlockDofsFilter[i] ) continue;
            // mark global DoF number
            for(int k = 0; k < localBlockSize; ++k)
              dirichletBlocks_[ globalBlockDofs[ i ] ][k] = block [k];
            // we have Dirichlet values
            hasDirichletBoundary = true ;
          }
        }
      }
    }

    return hasDirichletBoundary;
  }

  ModelType &model_;
  const DiscreteFunctionSpaceType& space_;
  mutable DirichletBlockVector dirichletBlocks_;
  mutable bool hasDirichletDofs_ ;
  mutable int sequence_ ;

  class DirichletBuilder;
};

template < class Model, class Space >
class DirichletConstraints< Model,Space > :: DirichletBuilder
    : public CommDataHandleIF< DirichletBuilder, int >
{
public:
  typedef Space SpaceType;
  typedef typename SpaceType::BlockMapperType MapperType;

  enum { nCodim = SpaceType :: GridType :: dimension + 1 };

public:
  typedef int DataType;

  const int myRank_;
  const int mySize_;

  typedef DirichletConstraints< Model,Space > DirichletType;
  const DirichletType &dirichlet_;

  const SpaceType &space_;
  const MapperType &mapper_;

  static const int blockSize = SpaceType::localBlockSize;

public:
  DirichletBuilder( const DirichletType &dirichlet,
                    const SpaceType &space,
                    const MapperType &mapper )
  : myRank_( space.gridPart().comm().rank() ),
    mySize_( space.gridPart().comm().size() ),
    dirichlet_( dirichlet ),
    space_( space ),
    mapper_( mapper )
  {
  }
  bool contains ( int dim, int codim ) const
  {
    return mapper_.contains( codim );
  }

  bool fixedSize ( int dim, int codim ) const
  {
    return false;
  }

  //! read buffer and apply operation
  template< class MessageBuffer, class Entity >
  inline void gather ( MessageBuffer &buffer,
                       const Entity &entity ) const
  {
    unsigned int localBlocks = mapper_.numEntityDofs( entity );
    std::vector<std::size_t> globalBlockDofs(localBlocks);
    mapper_.mapEntityDofs( entity, globalBlockDofs );
    assert( size(entity) == globalBlockDofs.size()*blockSize );
    for( unsigned int localBlock = 0 ; localBlock < globalBlockDofs.size(); ++ localBlock )
    {
      int global = globalBlockDofs[localBlock];
      for (int r=0;r<blockSize;++r)
        if (dirichlet_.dirichletBlocks_[ global ][r] )
          buffer.write( 1 );
        else
          buffer.write( 0 );
    }
  }

  //! read buffer and apply operation
  //! scatter is called for one every entity
  //! several times depending on how much data
  //! was gathered
  template< class MessageBuffer, class EntityType >
  inline void scatter ( MessageBuffer &buffer,
                        const EntityType &entity,
                        size_t n )
  {
    unsigned int localBlocks = mapper_.numEntityDofs( entity );
    std::vector<std::size_t> globalBlockDofs(localBlocks);
    mapper_.mapEntityDofs( entity, globalBlockDofs );
    assert( n == globalBlockDofs.size()*blockSize );
    assert( n == size(entity) );
    for( unsigned int localBlock = 0 ; localBlock < globalBlockDofs.size(); ++ localBlock )
    {
      int global = globalBlockDofs[localBlock];
      for (int r=0;r<blockSize;++r)
      {
        int val;
        buffer.read(val);
        if ( !dirichlet_.dirichletBlocks_[ global ][r] && val == 1)
          dirichlet_.dirichletBlocks_[ global ][r] = true;
      }
    }
  }
  //! return local dof size to be communicated
  template< class Entity >
  size_t size ( const Entity &entity ) const
  {
    return blockSize * mapper_.numEntityDofs( entity );
  }
};

} // end namespace Dune
#endif
