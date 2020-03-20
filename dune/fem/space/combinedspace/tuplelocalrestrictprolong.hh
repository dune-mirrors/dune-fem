#ifndef DUNE_FEM_SPACE_COMBINEDSPACE_TUPLELOCALRESTIRCTPROLONG_HH
#define DUNE_FEM_SPACE_COMBINEDSPACE_TUPLELOCALRESTIRCTPROLONG_HH

#include <algorithm>
#include <array>
#include <tuple>
#include <vector>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/fem/common/forloop.hh>
#include <dune/fem/common/utility.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/function/localfunction/localfunction.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>
#include <dune/fem/storage/subvector.hh>

namespace Dune
{

  namespace Fem
  {

    // TupleLocalRestricProlong
    // ------------------------

    template< class ... DiscreteFunctionSpaces >
    class TupleLocalRestrictProlong
    {
      typedef TupleLocalRestrictProlong< DiscreteFunctionSpaces ... > ThisType;

      typedef std::tuple< DefaultLocalRestrictProlong< DiscreteFunctionSpaces > ... > LocalRestrictProlongTupleType;

      static const int setSize = sizeof...( DiscreteFunctionSpaces )-1;

      // helper structs
      template< int > struct RestrictLocal;
      template< int > struct RestrictFinalize;
      template< int > struct ProlongLocal;

      template< std::size_t  ... i >
      static LocalRestrictProlongTupleType localRestrictProlongTuple ( std::tuple< const DiscreteFunctionSpaces &... > tuple, std::index_sequence< i ... > )
      {
        return std::make_tuple( typename std::tuple_element< i, LocalRestrictProlongTupleType >::type( std::get< i >( tuple ) ) ...);
      }

    public:
      static_assert( Std::are_all_same< typename DiscreteFunctionSpaces::DomainFieldType ... >::value,
          "TupleLocalRestrictProlong needs common DomainFieldType in the Spaces!" );

      typedef std::tuple_element_t< 0, std::tuple< typename DiscreteFunctionSpaces::DomainFieldType... > > DomainFieldType;

      TupleLocalRestrictProlong ( std::tuple< const DiscreteFunctionSpaces & ... > tuple )
        : localRestrictProlongTuple_( localRestrictProlongTuple( tuple, std::index_sequence_for< DiscreteFunctionSpaces ... >() ) )
      {}

      void setFatherChildWeight ( const DomainFieldType &weight )
      {
        Hybrid::forEach( std::make_index_sequence< sizeof ... ( DiscreteFunctionSpaces ) >{},
          [ & ]( auto i ){ std::get< i >( localRestrictProlongTuple_ ).setFatherChildWeight( weight ); } );
      }

      //! restrict data to father
      template< class LFFather, class LFSon, class LocalGeometry >
      void restrictLocal ( LFFather &lfFather, const LFSon &lfSon,
                           const LocalGeometry &geometryInFather, bool initialize ) const
      {
        Fem::ForLoop< RestrictLocal, 0, setSize >::apply( lfFather, lfSon, geometryInFather, initialize, localRestrictProlongTuple_ );
      }
      template< class LFFather >
      void restrictFinalize ( LFFather &lfFather ) const
      {
        Fem::ForLoop< RestrictFinalize, 0, setSize >::apply( lfFather, localRestrictProlongTuple_ );
      }

      template< class LFFather, class LFSon, class LocalGeometry >
      void prolongLocal ( const LFFather &lfFather, LFSon &lfSon,
                          const LocalGeometry &geometryInFather, bool initialize ) const
      {
        Fem::ForLoop< ProlongLocal, 0, setSize >::apply( lfFather, lfSon, geometryInFather, initialize, localRestrictProlongTuple_ );
      }

      bool needCommunication () const
      {
        return needCommunication( std::index_sequence_for< DiscreteFunctionSpaces ... >() );
      }

    protected:
      template< std::size_t ... i >
      bool needCommunication ( std::index_sequence< i...> ) const
      {
        return Std::Or( std::get< i >( localRestrictProlongTuple_ ).needCommunication() ... );
      }

    private:
      LocalRestrictProlongTupleType localRestrictProlongTuple_;
    };



    // ProlongLocal
    // ------------

    template< class ... DiscreteFunctionSpaces >
    template< int i >
    struct TupleLocalRestrictProlong< DiscreteFunctionSpaces ... >::
    ProlongLocal
    {
      template< class LFFather, class LFSon, class LocalGeometry, class Tuple >
      static void apply ( const LFFather &lfFather, LFSon &lfSon, const LocalGeometry &geometryInFather, bool initialize,
          const Tuple &tuple )
      {
        typedef SubVector< const typename LFFather::LocalDofVectorType, OffsetSubMapper > SubDofVectorTypeFather;
        typedef SubVector< typename LFSon::LocalDofVectorType, OffsetSubMapper > SubDofVectorTypeSon;

        typedef typename LFFather::BasisFunctionSetType::template SubBasisFunctionSet< i >::type SubFatherBasisFunctionSetType;
        typedef typename LFSon::BasisFunctionSetType::template SubBasisFunctionSet< i >::type SubSonBasisFunctionSetType;

        SubFatherBasisFunctionSetType subFatherBasisFunctionSet = lfFather.basisFunctionSet().template subBasisFunctionSet< i >();
        SubSonBasisFunctionSetType subSonBasisFunctionSet = lfSon.basisFunctionSet().template subBasisFunctionSet< i >();

        std::size_t fatherBasisSetOffset = lfFather.basisFunctionSet().offset( i );
        std::size_t sonBasisSetOffset = lfSon.basisFunctionSet().offset(i);

        SubDofVectorTypeSon sonSubDofVector( lfSon.localDofVector(), OffsetSubMapper( subSonBasisFunctionSet.size(), sonBasisSetOffset ) );
        SubDofVectorTypeFather fatherSubDofVector( lfFather.localDofVector(), OffsetSubMapper( subFatherBasisFunctionSet.size(),
                                                                                               fatherBasisSetOffset ) );

        BasicConstLocalFunction< SubFatherBasisFunctionSetType, SubDofVectorTypeFather > subLFFather( subFatherBasisFunctionSet,
                                                                                                      fatherSubDofVector );
        LocalFunction< SubSonBasisFunctionSetType, SubDofVectorTypeSon > subLFSon( subSonBasisFunctionSet, sonSubDofVector );

        std::get< i >( tuple ).prolongLocal( subLFFather, subLFSon, geometryInFather, initialize );
      }
    };


    // RestrictFinalize
    // ----------------

    template< class ... DiscreteFunctionSpaces >
    template< int i >
    struct TupleLocalRestrictProlong< DiscreteFunctionSpaces ... >::
    RestrictFinalize
    {
      template< class LFFather, class Tuple >
      static void apply ( LFFather &lfFather,
          const Tuple &tuple )
      {
        typedef SubVector< typename LFFather::LocalDofVectorType, OffsetSubMapper > SubDofVectorTypeFather;

        typedef typename LFFather::BasisFunctionSetType::template SubBasisFunctionSet< i >::type SubFatherBasisFunctionSetType;

        SubFatherBasisFunctionSetType subFatherBasisFunctionSet = lfFather.basisFunctionSet().template subBasisFunctionSet< i >();

        std::size_t fatherBasisSetOffset = lfFather.basisFunctionSet().offset(i);

        SubDofVectorTypeFather fatherSubDofVector( lfFather.localDofVector(), OffsetSubMapper( subFatherBasisFunctionSet.size(),
                                                                                               fatherBasisSetOffset ) );
        LocalFunction< SubFatherBasisFunctionSetType, SubDofVectorTypeFather > subLFFather( subFatherBasisFunctionSet, fatherSubDofVector );

        std::get< i >( tuple ).restrictFinalize( subLFFather );
      }
    };


    // RestrictLocal
    // -------------

    template< class ... DiscreteFunctionSpaces >
    template< int i >
    struct TupleLocalRestrictProlong< DiscreteFunctionSpaces ... >::
    RestrictLocal
    {
      template< class LFFather, class LFSon, class LocalGeometry, class Tuple >
      static void apply ( LFFather &lfFather, const LFSon &lfSon, const LocalGeometry &geometryInFather, bool initialize,
          const Tuple &tuple )
      {
        typedef SubVector< typename LFFather::LocalDofVectorType, OffsetSubMapper > SubDofVectorTypeFather;
        typedef SubVector< const typename LFSon::LocalDofVectorType, OffsetSubMapper > SubDofVectorTypeSon;

        typedef typename LFFather::BasisFunctionSetType::template SubBasisFunctionSet< i >::type SubFatherBasisFunctionSetType;
        typedef typename LFSon::BasisFunctionSetType::template SubBasisFunctionSet< i >::type SubSonBasisFunctionSetType;

        SubFatherBasisFunctionSetType subFatherBasisFunctionSet = lfFather.basisFunctionSet().template subBasisFunctionSet< i >();
        SubSonBasisFunctionSetType subSonBasisFunctionSet = lfSon.basisFunctionSet().template subBasisFunctionSet< i >();

        std::size_t fatherBasisSetOffset = lfFather.basisFunctionSet().offset(i);
        std::size_t sonBasisSetOffset = lfSon.basisFunctionSet().offset(i);

        SubDofVectorTypeSon sonSubDofVector( lfSon.localDofVector(), OffsetSubMapper( subSonBasisFunctionSet.size(), sonBasisSetOffset ) );
        SubDofVectorTypeFather fatherSubDofVector( lfFather.localDofVector(), OffsetSubMapper( subFatherBasisFunctionSet.size(),
                                                                                               fatherBasisSetOffset ) );

        LocalFunction< SubFatherBasisFunctionSetType, SubDofVectorTypeFather > subLFFather( subFatherBasisFunctionSet, fatherSubDofVector );
        BasicConstLocalFunction< SubSonBasisFunctionSetType, SubDofVectorTypeSon > subLFSon( subSonBasisFunctionSet, sonSubDofVector );

        std::get< i >( tuple ).restrictLocal( subLFFather, subLFSon, geometryInFather, initialize );
      }
    };


  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMBINEDSPACE_TUPLELOCALRESTIRCTPROLONG_HH
