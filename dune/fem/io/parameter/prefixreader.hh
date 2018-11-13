#ifndef DUNE_FEM_IO_PARAMETER_PREFIXREADER_HH
#define DUNE_FEM_IO_PARAMETER_PREFIXREADER_HH

#include <functional>
#include <string>
#include <utility>

#include <dune/fem/io/parameter/reader.hh>

namespace Dune
{

  namespace Fem
  {

    // PrefixParameter
    // ------------

    template< class Parameter >
    struct PrefixParameter
    {
      PrefixParameter ( std::string prefix, Parameter parameter )
        : prefix_( std::move( prefix ) ), parameter_( std::move( parameter ) )
      {}

      const std::string *operator() ( const std::string &key, const std::string *defaultValue ) const
      {
        return parameter()( prefix() + key, defaultValue );
      }

      const std::string &prefix () const { return prefix_; }
      const Parameter &parameter () const { return parameter_; }

    private:
      std::string prefix_;
      Parameter parameter_;
    };



    // PrefixParameterReader
    // ------------------

    template< class Parameter >
    class PrefixParameterReader
      : public BasicParameterReader< PrefixParameter< Parameter > >
    {
      typedef BasicParameterReader< PrefixParameter< Parameter > > BaseType;

    public:
      using BaseType::parameter;

      PrefixParameterReader ( std::string prefix, Parameter parameter )
        : BaseType( PrefixParameter< Parameter >( std::move( prefix ), std::move( parameter ) ) )
      {}

      PrefixParameterReader ( std::string prefix, const PrefixParameterReader &reader )
        : BaseType( PrefixParameter< Parameter >( std::move( prefix ) + reader.parameter().prefix(), reader.parameter().parameter() ) )
      {}

      operator ParameterReader () const { return ParameterReader( parameter() ); }
    };



    // PrefixParameterReader
    // ------------------

    inline static PrefixParameterReader< std::reference_wrapper< const ParameterContainerData > >
    prefixParameterReader ( std::string prefix, const ParameterContainer &container = Parameter::container() )
    {
      return PrefixParameterReader< std::reference_wrapper< const ParameterContainerData > >( std::move( prefix ), std::ref( container.parameter() ) );
    }

    inline static PrefixParameterReader< std::function< const std::string *( const std::string &, const std::string * ) > >
    prefixParameterReader ( std::string prefix, const ParameterReader &reader )
    {
      return PrefixParameterReader< std::function< const std::string *( const std::string &, const std::string * ) > >( std::move( prefix ), reader.parameter() );
    }

    template< class Parameter >
    inline static PrefixParameterReader< Parameter >
    prefixParameterReader ( std::string prefix, const PrefixParameterReader< Parameter > &reader )
    {
      return PrefixParameterReader< Parameter >( std::move( prefix ), reader );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_IO_PARAMETER_PREFIXREADER_HH
