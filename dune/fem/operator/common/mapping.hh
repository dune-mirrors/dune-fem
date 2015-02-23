#ifndef DUNE_FEM_MAPPING_HH
#define DUNE_FEM_MAPPING_HH

#include <iostream>
#include <vector>

namespace Dune
{

  namespace Fem
  {


    /** @addtogroup Mapping

        Mappings in Dune always map from one vector space into another vector space.
        Mapping are the base class for Operators and Functions.
        Operators work on vector spaces containing Functions (i.e. Domain and Range are Functions). In
        contrary Functions work on vector spaces containing real numbers (i.e.
        \f$R^n\f$). For both Mapping the base interface class. Furthermore,
        Mapping provided a machinery to combine different mapping linearly.

        \remarks
        The interface for Mappings is defined by the class Mapping.

      @{
     */

    /** \brief A mapping from one vector space into another
        This class describes a general mapping from the domain vector space into
        the range vector space.
        It can also be used to construct linear combinations of mappings.

        This two-sided character has the following consequence: when you address
        an object of type mapping or any of its descendants through a reference
        or pointer of type Mapping, the linear combination defined for that mapping
        is evaluated. On the other hand, if you address through a reference of the
        type of any of its descendants (notably Operator and Function), you'll
        get the functionality specific for that type.

        \note The Domain type and Range type must have operator *=, operator +=,
        and operator -= (e.g. FieldVector and DiscreteFunction fit that
        interface).
    */
    template<typename DFieldType,typename RFieldType, class DType, class RType>
    class Mapping
    {
    protected:
      //! type of mapping
      typedef Mapping<DFieldType,RFieldType,DType,RType> MappingType;

    public:
      /** \brief domain vector space (for usage in derived classes)
          This can either be for example a discrete function (in case of
          operators) or a FieldVector (in case of discrete functions)
      */
      typedef DType DomainType;

      /** \brief range vector space (for usage in derived classes)
          This can either be for example a discrete function (in case of
          operators) or a FieldVector (in case of discrete functions)
      */
      typedef RType  RangeType;

      /** \brief type of field the domain vector space, i.e. double */
      typedef DFieldType DomainFieldType;

      /** \brief type of field the range vector space, i.e. double */
      typedef RFieldType RangeFieldType;

      //! create Mappiung with empty linear combination
      Mapping()
        : lincomb_(),
          loopDetected_( false )
      {
        // store pointer to myself
        lincomb_.push_back( term( *this, 1. ) );
      }

      //! delete linear combination if necessary
      virtual ~Mapping ()
      {
      }

      /** \brief assignment of mapping mapping
          \param[in] mapping which is copied
          \returns reference to mapping
      */
      MappingType& operator = (const MappingType &mapping)
      {
        const MappingType &m = dynamic_cast<const MappingType& >( mapping );

        lincomb_.erase( lincomb_.begin(), lincomb_.end() );
        typedef typename std::vector<term>::const_iterator iterator;

        // reserve memory
        lincomb_.reserve( m.lincomb_.size() );

        iterator end = m.lincomb_.end();
        for (iterator it = m.lincomb_.begin(); it != end; it++ )
        {
          lincomb_.push_back( *it );
        }
        return *this;
      }

      /** \brief Application operator that applies all operators in the
          linear combination stack.
          \param[in] arg argument
          \param[out] dest destination
      */
      void operator() (const DomainType &arg, RangeType &dest ) const
      {
        int count = 0;
        typedef typename std::vector<term>::const_iterator const_iterator;
        const_iterator end = lincomb_.end();
        for ( const_iterator it = lincomb_.begin(); it != end; ++it )
        {
          const term& current = (*it);
          if ( count == 0 )
          {
            // call apply of term, see below
            current.apply( arg, dest );
          }
          else
          {
            // call applyAdd of term, see below
            current.applyAdd( arg, dest );
          }
          ++count;
        }

        loopDetected_ = false ;
      }

    private:
      //! apply operators
      virtual void apply (const DomainType &arg, RangeType &dest) const
      {
        // do noting if apply method leads to myself
        if( loopDetected_ )
        {
#ifndef NDEBUG
          std::cerr << "WARNING: Mapping::apply: Loop detected! Check overloading of operator() and apply! " << std::endl;
#endif
          return ;
        }

        // set applyLoop to true again
        loopDetected_ = true ;

        operator()(arg, dest);
      }

      //! linear comnination object
      struct term {
        term() : v_(NULL), scalar_(1.0), scaleIt_(false) { }
        term(const term& other)
          : v_(other.v_), scalar_(other.scalar_), scaleIt_(other.scaleIt_) { }

        term(const MappingType &mapping, RangeFieldType scalar ) : v_(&mapping), scalar_(scalar), scaleIt_( true ) {
          if ( scalar_ == 1. ) {
            scaleIt_ = false;
          }
        }

        void apply(const DomainType &arg, RangeType &dest) const
        {
          v_->apply( arg, dest );
          if ( scaleIt_ )
          {
            dest *= scalar_;
          }
        }

        void applyAdd(const DomainType &arg, RangeType &dest) const
        {
          // note, copying here might be costly
          RangeType tmp( dest );

          v_->apply( arg, tmp );
          if ( scalar_ == 1. )
          {
            dest += tmp;
          }
          else if ( scalar_ == -1. )
          {
            dest -= tmp;
          }
          else
          {
            tmp *= scalar_;
            dest += tmp;
          }
        }

      protected:
        const MappingType *v_;
        RangeFieldType scalar_;
        bool scaleIt_;

        // friendship for operations
        friend struct MappingOperators;
      };


      //! vector holding linear combination factors
      std::vector<term> lincomb_;

      //! true if we are stuck in an apply - operator () loop
      mutable bool loopDetected_ ;

      // friendship for operations
      friend struct MappingOperators;
    };


    /** \brief Implementation of Mapping +, -, *, / operations. */
    struct MappingOperators
    {
      //! \brief copy mapping
      template<typename DFieldType,typename RFieldType, class DType, class RType>
      static inline void
      copyMapping(const Mapping<DFieldType,RFieldType,DType,RType>& org,
                  Mapping<DFieldType,RFieldType,DType,RType>& copy)
      {
        typedef Mapping<DFieldType,RFieldType,DType,RType> MappingType;
        typedef typename std::vector< typename MappingType :: term > :: const_iterator iterator;

        // clear mapping entries
        copy.lincomb_.clear();

        {
          iterator end = org.lincomb_.end();
          for ( iterator it = org.lincomb_.begin(); it != end; ++it )
          {
            copy.lincomb_.push_back( *it );
          }
        }
      }

      //! \brief add mappings
      template<typename DFieldType,typename RFieldType, class DType, class RType>
      static inline Mapping<DFieldType,RFieldType,DType,RType>
      addMappings(const Mapping<DFieldType,RFieldType,DType,RType>& a,
                  const Mapping<DFieldType,RFieldType,DType,RType>& b)
      {
        typedef Mapping<DFieldType,RFieldType,DType,RType> MappingType;
        // new mapping
        MappingType newMapping;

        // copy mapping
        copyMapping(a, newMapping);

        typedef typename std::vector< typename MappingType :: term > :: const_iterator iterator;

        iterator end = b.lincomb_.end();
        for ( iterator it = b.lincomb_.begin(); it != end; ++it )
        {
          newMapping.lincomb_.push_back( *it );
        }
        return newMapping;
      }

      //! \brief substract mappings
      template<typename DFieldType,typename RFieldType, class DType, class RType>
      static inline Mapping<DFieldType,RFieldType,DType,RType>
      substractMappings(const Mapping<DFieldType,RFieldType,DType,RType>& a,
                        const Mapping<DFieldType,RFieldType,DType,RType>& b)
      {
        typedef Mapping<DFieldType,RFieldType,DType,RType> MappingType;
        typedef typename  MappingType :: term term;
        typedef typename std::vector< term > :: const_iterator iterator;
        // new mapping
        MappingType newMapping;

        // copy mapping
        copyMapping(a, newMapping);

        iterator end = b.lincomb_.end();
        for ( iterator it = b.lincomb_.begin(); it != end; ++it )
        {
          newMapping.lincomb_.push_back( term( *it->v_, -it->scalar_ ) );
        }
        return newMapping;
      }

      //! \brief multiply mapping
      template<typename DFieldType,typename RFieldType, class DType, class RType>
      static inline Mapping<DFieldType,RFieldType,DType,RType>
      multiplyMapping(const Mapping<DFieldType,RFieldType,DType,RType>& a,
                      const RFieldType& scalar)
      {
        typedef Mapping<DFieldType,RFieldType,DType,RType> MappingType;
        typedef typename  MappingType :: term term;
        typedef typename std::vector< term > :: iterator iterator;
        // new mapping
        MappingType newMapping;

        // copy mapping
        copyMapping(a, newMapping);

        iterator end = newMapping.lincomb_.end();
        for ( iterator it = newMapping.lincomb_.begin(); it != end; ++it )
        {
          it->scalar_ *= scalar;
        }
        return newMapping;
      }

      //! \brief divide mapping
      template<typename DFieldType,typename RFieldType, class DType, class RType>
      static inline Mapping<DFieldType,RFieldType,DType,RType>
      divideMapping(const Mapping<DFieldType,RFieldType,DType,RType>& a,
                    const RFieldType& scalar)
      {
        RFieldType factor = RFieldType(1)/scalar;
        return multiplyMapping(a,factor);
      }
    };

    /** @} end documentation group */

    /** \relates Mapping
        \brief add two mappings
        \param[in] a mapping 1
        \param[in] b mapping 2
        \returns new object mapping
    */
    template<class DFieldType, class RFieldType, class DType, class RType>
    static inline Mapping<DFieldType,RFieldType,DType,RType>
    operator +(const Mapping<DFieldType,RFieldType,DType,RType>& a,
               const Mapping<DFieldType,RFieldType,DType,RType>& b)
    {
      return MappingOperators::addMappings(a,b);
    }

    /** \relates Mapping
        \brief substract two mappings
        \param[in] a mapping 1
        \param[in] b mapping 2
        \returns new object mapping
    */
    template<class DFieldType, class RFieldType, class DType, class RType>
    static inline Mapping<DFieldType,RFieldType,DType,RType>
    operator -(const Mapping<DFieldType,RFieldType,DType,RType>& a,
               const Mapping<DFieldType,RFieldType,DType,RType>& b)
    {
      return MappingOperators::substractMappings(a,b);
    }

    /** \relates Mapping
        \brief scale mapping with factor
        \param[in] mapping Mapping which is scaled
        \param[in] factor  factor with which mapping is scaled
        \returns new object mapping
    */
    template<class DFieldType, class RFieldType, class DType, class RType>
    static inline Mapping<DFieldType,RFieldType,DType,RType>
    operator *(const Mapping<DFieldType,RFieldType,DType,RType>& mapping,
               const RFieldType& factor)
    {
      return MappingOperators::multiplyMapping(mapping,factor);
    }

    /** \relates Mapping
        \brief scale mapping with factor
        \param[in] factor  factor with which mapping is scaled
        \param[in] mapping Mapping which is scaled
        \returns new object mapping
    */
    template<class DFieldType, class RFieldType, class DType, class RType>
    static inline Mapping<DFieldType,RFieldType,DType,RType>
    operator *(const RFieldType& factor,
               const Mapping<DFieldType,RFieldType,DType,RType>& mapping)
    {
      return MappingOperators::multiplyMapping(mapping,factor);
    }

    /** \relates Mapping
        \brief operator / for mappings
        \param[in] mapping mapping which is divided
        \param[in] factor f factor by which result of mapping is divided
        \returns new object mapping
    */
    template<class DFieldType, class RFieldType, class DType, class RType>
    static inline Mapping<DFieldType,RFieldType,DType,RType>
    operator /(const Mapping<DFieldType,RFieldType,DType,RType>& mapping,
               const RFieldType& factor)
    {
      return MappingOperators::divideMapping(mapping,factor);
    }

    /** \relates Mapping
        \brief operator / for mappings
        \param[in] factor by which result of mapping is divided
        \param[in] mapping which is divided
        \returns new object mapping
    */
    template<class DFieldType, class RFieldType, class DType, class RType>
    static inline Mapping<DFieldType,RFieldType,DType,RType>
    operator /(const RFieldType& factor,
               const Mapping<DFieldType,RFieldType,DType,RType>& mapping)
    {
      return MappingOperators::divideMapping(mapping,factor);
    }

  } // namespace Fem

} // namespace Dune
#endif // #ifndef DUNE_FEM_MAPPING_HH
