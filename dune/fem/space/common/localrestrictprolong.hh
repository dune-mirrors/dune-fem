#ifndef DUNE_FEM_LOCALRESTRICTPROLONG_HH
#define DUNE_FEM_LOCALRESTRICTPROLONG_HH

#include <dune/fem/function/localfunction/localfunction.hh>

namespace Dune
{

  namespace Fem
  {

    // DefaultLocalRestrictProlong
    // ---------------------------

    template< class DiscreteFunctionSpace >
    class DefaultLocalRestrictProlong;



    // ConstantLocalRestrictProlong
    // ----------------------------

    template< class DiscreteFunctionSpace >
    class ConstantLocalRestrictProlong
    {
      typedef ConstantLocalRestrictProlong< DiscreteFunctionSpace > ThisType;

    public:
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

      typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;

      ConstantLocalRestrictProlong ()
      : weight_( -1 )
      {}

      /** \brief explicit set volume ratio of son and father
       *
       *  \param[in]  weight  volume of son / volume of father
       *
       *  \note If this ratio is set, it is assume to be constant.
       */
      void setFatherChildWeight ( const DomainFieldType &weight )
      {
        weight_ = weight;
      }

      //! restrict data to father
      template< class LFFather, class LFSon, class LocalGeometry >
      void restrictLocal ( LFFather &lfFather, const LFSon &lfSon,
                           const LocalGeometry &geometryInFather, bool initialize ) const
      {
        const DomainFieldType weight = (weight_ < DomainFieldType( 0 ) ? calcWeight( lfFather.entity(), lfSon.entity() ) : weight_);

        assert( weight > 0.0 );
        //assert( std::abs( geometryInFather.volume() - weight ) < 1e-8 );

        const int size = lfFather.size();
        assert( lfFather.size() == lfSon.size() );
        if( initialize )
        {
          for( int i = 0; i < size; ++i )
            lfFather[ i ] = weight * lfSon[ i ];
        }
        else
        {
          for( int i = 0; i < size; ++i )
            lfFather[ i ] += weight * lfSon[ i ];
        }
      }
      template< class LFFather >
      void restrictFinalize ( LFFather &lfFather ) const
      {}

      //! prolong data to children
      template< class LFFather, class LFSon, class LocalGeometry >
      void prolongLocal ( const LFFather &lfFather, LFSon &lfSon,
                          const LocalGeometry &geometryInFather, bool initialize ) const
      {
        const int size = lfFather.size();
        assert( lfFather.size() == lfSon.size() );
        for( int i = 0; i < size; ++i )
          lfSon[ i ] = lfFather[ i ];
      }

      //! do discrete functions need a communication after restriction / prolongation?
      bool needCommunication () const { return true; }


      template< class Entity >
      static DomainFieldType calcWeight ( const Entity &father, const Entity &son )
      {
        return son.geometry().volume() / father.geometry().volume();
      }

    protected:
      DomainFieldType weight_;
    };



    // EmptyLocalRestrictProlong
    // -------------------------

    template< class DiscreteFunctionSpace >
    class EmptyLocalRestrictProlong
    {
      typedef EmptyLocalRestrictProlong< DiscreteFunctionSpace > ThisType;

    public:
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

      typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;

      /** \brief explicit set volume ratio of son and father
       *
       *  \param[in]  weight  volume of son / volume of father
       *
       *  \note If this ratio is set, it is assume to be constant.
       */
      void setFatherChildWeight ( const DomainFieldType &weight ) {}

      //! restrict data to father
      template< class LFFather, class LFSon, class LocalGeometry >
      void restrictLocal ( LFFather &lfFather, const LFSon &lfSon,
                           const LocalGeometry &geometryInFather, bool initialize ) const
      {}
      template< class LFFather >
      void restrictFinalize ( LFFather &lfFather ) const
      {}

      //! prolong data to children
      template< class LFFather, class LFSon, class LocalGeometry >
      void prolongLocal ( const LFFather &lfFather, LFSon &lfSon,
                          const LocalGeometry &geometryInFather, bool initialize ) const
      {}

      //! do discrete functions need a communication after restriction / prolongation?
      bool needCommunication () const { return false; }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_LOCALRESTRICTPROLONG_HH
