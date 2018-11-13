#ifndef DUNE_FEM_IO_PARAMETER_SUBREADER_HH
#define DUNE_FEM_IO_PARAMETER_SUBREADER_HH

#include <functional>
#include <string>
#include <utility>

#include <dune/fem/io/parameter/reader.hh>

namespace Dune
{

  namespace Fem
  {

    // SubParameter
    // ------------

    template< class Parameter >
    struct SubParameter
    {
      SubParameter ( std::string prefix, Parameter parameter )
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



    // SubParameterReader
    // ------------------

    template< class Parameter >
    class SubParameterReader
      : public BasicParameterReader< SubParameter< Parameter > >
    {
      typedef BasicParameterReader< SubParameter< Parameter > > BaseType;

    public:
      using BaseType::parameter;

      SubParameterReader ( std::string prefix, Parameter parameter )
        : BaseType( SubParameter< Parameter >( std::move( prefix ), std::move( parameter ) ) )
      {}

      SubParameterReader ( std::string prefix, const SubParameterReader &reader )
        : BaseType( SubParameter< Parameter >( std::move( prefix ) + reader.parameter().prefix(), reader.parameter().parameter() ) )
      {}

      operator ParameterReader () const { return ParameterReader( parameter() ); }
    };



    // subParameterReader
    // ------------------

    inline static SubParameterReader< std::reference_wrapper< const ParameterContainerData > >
    subParameterReader ( std::string prefix, const ParameterContainer &container = Parameter::container() )
    {
      return SubParameterReader< std::reference_wrapper< const ParameterContainerData > >( std::move( prefix ), std::ref( container.parameter() ) );
    }

    inline static SubParameterReader< std::function< const std::string *( const std::string &, const std::string * ) > >
    subParameterReader ( std::string prefix, const ParameterReader &reader )
    {
      return SubParameterReader< std::function< const std::string *( const std::string &, const std::string * ) > >( std::move( prefix ), reader.parameter() );
    }

    template< class Parameter >
    inline static SubParameterReader< Parameter >
    subParameterReader ( std::string prefix, const SubParameterReader< Parameter > &reader )
    {
      return SubParameterReader< Parameter >( std::move( prefix ), reader );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_IO_PARAMETER_SUBREADER_HH
