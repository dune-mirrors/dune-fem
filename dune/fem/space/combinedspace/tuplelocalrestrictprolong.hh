#ifndef DUNE_FEM_SPACE_COMBINEDSPACE_TUPLELOCALRESTIRCTPROLONG_HH
#define DUNE_FEM_SPACE_COMBINEDSPACE_TUPLELOCALRESTIRCTPROLONG_HH

#include <algorithm>
#include <array>
#include <tuple>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/forloop.hh>
#include <dune/common/std/utility.hh>
#include <dune/common/tuples.hh>

#include <dune/fem/common/utility.hh>
#include <dune/fem/common/subvector.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/function/localfunction/localfunction.hh>

#include <dune/fem/space/common/localrestrictprolong.hh>


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

      typedef Dune::tuple< DefaultLocalRestrictProlong< DiscreteFunctionSpaces > ... > LocalRestrictProlongTupleType;

      static const int setSize = sizeof...( DiscreteFunctionSpaces )-1;

      // helper structs
      template< int > struct SetFatherChildWeight;
      template< int > struct RestrictLocal;
      template< int > struct ProlongLocal;

    public:
      static_assert( Std::are_all_same< typename DiscreteFunctionSpaces::DomainFieldType ... >::value, 
          "TupleLocalRestrictProlong needs common DomainFieldType in the Spaces!" );

      typedef typename Dune::tuple_element< 0, LocalRestrictProlongTupleType >::type::DomainFieldType DomainFieldType;

      TupleLocalRestrictProlong ( const DiscreteFunctionSpaces & ...spaces )
        : localRestrictProlongTuple_( spaces ... )
      {}

      void setFatherChildWeight ( const DomainFieldType &weight )
      {
        ForLoop< SetFatherChildWeight, 0, setSize >::apply( weight, localRestrictProlongTuple_ );
      }

      //! restrict data to father
      template< class LFFather, class LFSon, class LocalGeometry >
      void restrictLocal ( LFFather &lfFather, const LFSon &lfSon,
                           const LocalGeometry &geometryInFather, bool initialize ) const
      {
        ForLoop< RestrictLocal, 0, setSize >::apply( lfFather, lfSon, geometryInFather, initialize, localRestrictProlongTuple_ );
      }


      template< class LFFather, class LFSon, class LocalGeometry >
      void prolongLocal ( const LFFather &lfFather, LFSon &lfSon,
                          const LocalGeometry &geometryInFather, bool initialize ) const
      {
        ForLoop< ProlongLocal, 0, setSize >::apply( lfFather, lfSon, geometryInFather, initialize, localRestrictProlongTuple_ );
      }

      bool needCommunication () const 
      {
        return needCommunication( Std::index_sequence_for< DiscreteFunctionSpaces ... >() ); 
      }

    protected:
      template< std::size_t ... i >
      bool needCommunication ( Std::index_sequence< i...> ) const
      {
        return Std::Or( std::get< i >( localRestrictProlongTuple_ ).needCommunication() ... );
      }

    private:
      LocalRestrictProlongTupleType localRestrictProlongTuple_;
    };


    // SetFatherChildWeight
    // --------------------

    template< class ... DiscreteFunctionSpaces >
    template< int i >
    struct TupleLocalRestrictProlong< DiscreteFunctionSpaces ... >::
    SetFatherChildWeight
    {
      template< class Tuple >
      static void apply( const DomainFieldType &weight, Tuple &tuple )
      {
        std::get< i >( tuple ).setFatherChildWeight( weight );
      }
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
        typedef DenseSubVector< const typename LFFather::LocalDofVectorType > SubDofVectorTypeFather;
        typedef DenseSubVector< typename LFSon::LocalDofVectorType > SubDofVectorTypeSon;
        
        typedef typename LFFather::BasisFunctionSetType::template SubBasisFunctionSet< i >::type SubFatherBasisFunctionSetType;
        typedef typename LFSon::BasisFunctionSetType::template SubBasisFunctionSet< i >::type SubSonBasisFunctionSetType;

        SubFatherBasisFunctionSetType subFatherBasisFunctionSet = lfFather.basisFunctionSet().template subBasisFunctionSet< i >();
        SubSonBasisFunctionSetType subSonBasisFunctionSet = lfSon.basisFunctionSet().template subBasisFunctionSet< i >();

        std::size_t fatherBasisSetOffset = lfFather.basisFunctionSet().offset( i );
        std::size_t sonBasisSetOffset = lfSon.basisFunctionSet().offset(i);

        SubDofVectorTypeFather fatherSubDofVector( lfFather.localDofVector(), subFatherBasisFunctionSet.size(), fatherBasisSetOffset );
        SubDofVectorTypeSon sonSubDofVector( lfSon.localDofVector(), subSonBasisFunctionSet.size(), sonBasisSetOffset );

        BasicConstLocalFunction< SubFatherBasisFunctionSetType, SubDofVectorTypeFather > 
          subLFFather( subFatherBasisFunctionSet, fatherSubDofVector );
        LocalFunction< SubSonBasisFunctionSetType, SubDofVectorTypeSon > 
          subLFSon( subSonBasisFunctionSet, sonSubDofVector );

        std::get< i >( tuple ).prolongLocal( subLFFather, subLFSon, geometryInFather, initialize );
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
        typedef DenseSubVector< typename LFFather::LocalDofVectorType > SubDofVectorTypeFather;
        typedef DenseSubVector< const typename LFSon::LocalDofVectorType > SubDofVectorTypeSon;
        
        typedef typename LFFather::BasisFunctionSetType::template SubBasisFunctionSet< i >::type SubFatherBasisFunctionSetType;
        typedef typename LFSon::BasisFunctionSetType::template SubBasisFunctionSet< i >::type SubSonBasisFunctionSetType;

        SubFatherBasisFunctionSetType subFatherBasisFunctionSet = lfFather.basisFunctionSet().template subBasisFunctionSet< i >();
        SubSonBasisFunctionSetType subSonBasisFunctionSet = lfSon.basisFunctionSet().template subBasisFunctionSet< i >();

        std::size_t fatherBasisSetOffset = lfFather.basisFunctionSet().offset(i);
        std::size_t sonBasisSetOffset = lfSon.basisFunctionSet().offset(i);

        SubDofVectorTypeFather fatherSubDofVector( lfFather.localDofVector(), subFatherBasisFunctionSet.size(), fatherBasisSetOffset );
        SubDofVectorTypeSon sonSubDofVector( lfSon.localDofVector(), subSonBasisFunctionSet.size(), sonBasisSetOffset );

        LocalFunction< SubFatherBasisFunctionSetType, SubDofVectorTypeFather > subLFFather( subFatherBasisFunctionSet, fatherSubDofVector );
        BasicConstLocalFunction< SubSonBasisFunctionSetType, SubDofVectorTypeSon > subLFSon( subSonBasisFunctionSet, sonSubDofVector );

        std::get< i >( tuple ).restrictLocal( subLFFather, subLFSon, geometryInFather, initialize );
      }
    };


  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMBINEDSPACE_TUPLELOCALRESTIRCTPROLONG_HH
