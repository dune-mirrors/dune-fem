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
#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_MINIELEMENT_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_MINIELEMENT_HH

#include <dune/grid/common/gridenums.hh>

// dune-fem includes
#include <dune/fem/common/hybrid.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>

#include <dune/fem/space/basisfunctionset/default.hh>

#include <dune/fem/space/shapefunctionset/selectcaching.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>

#include <dune/fem/space/mapper/indexsetdofmapper.hh>

#include <dune/fem/operator/matrix/istlmatrixadapter.hh>

#include "localinterpolation.hh"
#include "localkeymap.hh"
#include "shapefunctions.hh"

namespace Dune
{

  namespace Fem
  {

    // Forward declaration
    // -------------------

    template< class FunctionSpace, class GridPart, template< class > class Storage >
    class BubbleElementSpace;

    // BubbleElementSpaceTraits
    // ----------------------

    template< class FunctionSpace, class GridPart, template< class > class Storage >
    struct BubbleElementSpaceTraits
    {
      typedef BubbleElementSpace< FunctionSpace, GridPart, Storage > DiscreteFunctionSpaceType;

      typedef FunctionSpace FunctionSpaceType;
      typedef GridPart GridPartType;

      static const int codimension = 0;

    private:
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;
    public:
      typedef typename FunctionSpaceType :: ScalarFunctionSpaceType ScalarFunctionSpaceType;

      // defined in shapefunctionset.hh
      // First check if we have a grid with only one element type
      // (hybrid // grids are not supported with this space yet)
      // If it only has one get the topologyId of that element type
      // (note simplex=0 and cube=2^dim-1)
      static_assert( Dune::Capabilities::hasSingleGeometryType<typename GridPartType::GridType>::v,
                     "bubble elements only implemented for grids with a single geometry type" );
      static const unsigned int topologyId =
        Dune::Capabilities::hasSingleGeometryType<typename GridPartType::GridType>::topologyId;
      // now choose either the simplex or the cube implementation of the
      // bubble space
      typedef SimplexBubbleElementShapeFunctionSet< ScalarFunctionSpaceType > ScalarSimplexShapeFunctionSetType;
      typedef Cube2DBubbleElementShapeFunctionSet< ScalarFunctionSpaceType > ScalarCubeShapeFunctionSetType;
      typedef typename std::conditional< topologyId == 0,
              ScalarSimplexShapeFunctionSetType, ScalarCubeShapeFunctionSetType >::type
                ScalarShapeFunctionSetType;
      // now extend it to a vector valued funcion space
      typedef VectorialShapeFunctionSet< ScalarShapeFunctionSetType, typename FunctionSpaceType::RangeType > ShapeFunctionSetType;
      // and add some default functionality
      typedef DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetType;
      // finished with the shape function set

      static const int localBlockSize = FunctionSpaceType::dimRange;
      static const int polynomialOrder = ScalarShapeFunctionSetType::polynomialOrder;

      typedef IndexSetDofMapper< GridPartType > BlockMapperType;
      typedef Hybrid::IndexRange< int, FunctionSpaceType::dimRange > LocalBlockIndices;

      template <class DiscreteFunction, class Operation = DFCommunicationOperation::Add >
      struct CommDataHandle
      {
        typedef Operation OperationType;
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
      };
    };

    // BubbleElementSpace
    // ----------------

    //! [Class definition for new space]
    template< class FunctionSpace, class GridPart, template< class > class Storage = CachingStorage >
    class BubbleElementSpace
    : public DiscreteFunctionSpaceDefault< BubbleElementSpaceTraits< FunctionSpace, GridPart, Storage > >
    //! [Class definition for new space]
    {
      typedef BubbleElementSpace< FunctionSpace, GridPart, Storage > ThisType;
      typedef DiscreteFunctionSpaceDefault< BubbleElementSpaceTraits< FunctionSpace, GridPart, Storage > > BaseType;

    public:
      typedef BubbleElementSpaceTraits< FunctionSpace, GridPart, Storage > Traits;
      typedef typename Traits::ShapeFunctionSetType ShapeFunctionSetType;
      static const int polynomialOrder = Traits::polynomialOrder;

      typedef typename BaseType::GridPartType GridPartType;

      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;
      typedef typename BaseType::IntersectionType IntersectionType;
      typedef typename BaseType::EntityType EntityType;

      typedef typename BaseType::BlockMapperType BlockMapperType;

      // type of local interpolation
      typedef LocalBubbleElementInterpolation< FunctionSpace > InterpolationType;

      // static const InterfaceType defaultInterface = InteriorBorder_InteriorBorder_Interface;
      static const InterfaceType defaultInterface = GridPartType::indexSetInterfaceType;
      static const CommunicationDirection defaultDirection = ForwardCommunication;

      BubbleElementSpace ( GridPartType &gridPart,
                         const InterfaceType commInterface = defaultInterface,
                         const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        blockMapper_(  gridPart, BubbleElementDofMapperCodeFactory() )
      {
      }

      ~BubbleElementSpace ()
      {
      }

      bool contains ( const int codim ) const
      {
        // forward to mapper since this information is held there
        return blockMapper().contains( codim );
      }

      bool continuous () const { return true; }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous ( const IntersectionType & intersection ) const
      {
        // forward to the subsapces
        return true;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order () const
      {
        return polynomialOrder;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      template<class Entity>
      int order ( const Entity &entity ) const
      {
        return polynomialOrder;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::basisFunctionSet(const EntityType &entity) const */
      template< class EntityType >
      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        return BasisFunctionSetType( entity, ShapeFunctionSetType( entity.geometry().type() ) );
      }

      /** \brief obtain the DoF block mapper of this space
          \return BlockMapperType
      **/
      BlockMapperType &blockMapper () const { return blockMapper_; }

      // Non-interface methods

      /** \brief return local interpolation for given entity
       *
       *  \param[in]  entity  grid part entity
       *  \note this method is needed to call the global function inteprolate( ... )
       */
      InterpolationType interpolation ( const EntityType &entity ) const
      {
        return InterpolationType();
      }

    protected:
      mutable BlockMapperType blockMapper_;
    };

    template< class FunctionSpace,
              class GridPart,
              template<class> class Storage,
              class NewFunctionSpace >
    struct DifferentDiscreteFunctionSpace< BubbleElementSpace < FunctionSpace, GridPart, Storage >, NewFunctionSpace >
    {
      typedef BubbleElementSpace< NewFunctionSpace, GridPart, Storage > Type;
    };

    template< class FunctionSpace, class GridPart, template< class > class Storage >
    class DefaultLocalRestrictProlong < BubbleElementSpace< FunctionSpace, GridPart, Storage > >
    : public EmptyLocalRestrictProlong< BubbleElementSpace< FunctionSpace, GridPart, Storage > >
    {
      typedef EmptyLocalRestrictProlong< BubbleElementSpace< FunctionSpace, GridPart, Storage > > BaseType;
      public:
      DefaultLocalRestrictProlong( const BubbleElementSpace< FunctionSpace, GridPart, Storage > &space )
        : BaseType()
      {}
    };

    template< class MatrixImp,
              class FunctionSpace, class GridPart, template< class > class Storage>
    struct ISTLParallelMatrixAdapter< MatrixImp, BubbleElementSpace< FunctionSpace,GridPart,Storage> >
    {
      typedef LagrangeParallelMatrixAdapter<MatrixImp> Type;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_MINIELEMENT_HH
