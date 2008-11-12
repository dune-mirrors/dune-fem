namespace Dune
{
  
  template <class DiscreteFunctionSpaceImp>
  AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>
    :: AdaptiveFunctionImplementation ( const std :: string &name,
                                        const DiscreteFunctionSpaceType &spc )
  : spc_(spc),
    memObject_( 0 ),
    dofVec_( allocateDofStorage( name ) )
  {
  }

  // create discrete function with vector 
  template <class DiscreteFunctionSpaceImp>
  template <class VectorPointerType> 
  AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>
    :: AdaptiveFunctionImplementation ( const std :: string &name,
                                        const DiscreteFunctionSpaceType &spc,
                                        VectorPointerType *vector )
  : spc_(spc),
    memObject_( 0 ),
    dofVec_( allocateDofStorageWrapper( spc_.mapper(), name, vector ) )
  {}

  template <class DiscreteFunctionSpaceImp>
  AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>
    :: AdaptiveFunctionImplementation( const std :: string &name,
                                       const DiscreteFunctionSpaceType &spc,
                                       DofStorageType &dofVec)
  : spc_(spc),
    memObject_( 0 ),
    dofVec_(dofVec)
  {}

  template <class DiscreteFunctionSpaceImp>
  AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  AdaptiveFunctionImplementation( const std :: string &name,
                                  const ThisType &other )
  : spc_( other.spc_ ),
    memObject_( 0 ),
    dofVec_( allocateDofStorage( name ) )
  {
    // copy values
    dofVec_ = other.dofVec_;
  }

  template <class DiscreteFunctionSpaceImp>
  AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  ~AdaptiveFunctionImplementation() 
  {
    if( memObject_ ) 
    {
      delete memObject_;
    }
  }
 
  template <class DiscreteFunctionSpaceImp>
  inline void 
  AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  clear () 
  {
    dofVec_.clear();
  }

  template <class DiscreteFunctionSpaceImp>
  inline void 
  AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  addScaled(const ThisType& org, const RangeFieldType& c) 
  {
    dofVec_.axpy(org.dofVec_ , c);
  }

  // operator=
  template <class DiscreteFunctionSpaceImp>
  inline void 
  AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  assignFunction(const ThisType& org)
  {
    assert(this->size() == org.size());
    dofVec_ = org.dofVec_;
  }

  // operator +=
  template <class DiscreteFunctionSpaceImp>
  inline void 
  AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  addFunction(const ThisType& org)
  {
    assert(this->size() == org.size());
    dofVec_ += org.dofVec_;
  }

  // operator +=
  template <class DiscreteFunctionSpaceImp>
  inline void 
  AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  substractFunction(const ThisType& org)
  {
    assert(this->size() == org.size());
    dofVec_ -= org.dofVec_;
  }

  template <class DiscreteFunctionSpaceImp>
  int AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  size() const 
  {
    return dofVec_.size();
  }

  template <class DiscreteFunctionSpaceImp>
  typename AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  DofIteratorType
  AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  dbegin() 
  {
    return dofStorage().begin();
  }

  template <class DiscreteFunctionSpaceImp>
  typename AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  DofIteratorType
  AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  dend() 
  {
    return dofStorage().end();
  }

  template <class DiscreteFunctionSpaceImp>
  typename AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  ConstDofIteratorType
  AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  dbegin() const 
  {
    return dofStorage().begin();
  }

  template <class DiscreteFunctionSpaceImp>
  typename AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  ConstDofIteratorType
  AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  dend() const
  {
    return dofStorage().end();
  }

#if DUNE_FEM_COMPATIBILITY
  //- Read/write methods
  template<class DiscreteFunctionSpaceImp>
  bool AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  write_xdr(const std::string fn) const
  {
    // create write stream 
    XDRWriteStream xdr(fn);
    // make sure data is only written in compressed state. 
    if( dofVec_.size() != spc_.size() )
    {
      DUNE_THROW(InvalidStateException,"DofVector not in compressed state while writing!"); 
    }

    // write data 
    return dofVec_.processXdr(xdr);
  }

  template <class DiscreteFunctionSpaceImp>
  bool AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  read_xdr(const std::string fn)
  {
    XDRReadStream xdr(fn);
    // make sure data is only written in compressed state. 
    if( dofVec_.size() != spc_.size() )
    {
      DUNE_THROW(InvalidStateException,"DofVector size does not match while reading!"); 
    }

    // read data 
    return dofVec_.processXdr(xdr);
  }

  template <class DiscreteFunctionSpaceImp>
  bool AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  write_ascii(const std::string fn) const
  {
    std::fstream outfile( fn.c_str() , std::ios::out );
    if (!outfile)
    { 
      printf( "\aERROR in AdaptiveDiscreteFunction::write_ascii(..): could not open <%s>!\n", fn.c_str());
      fflush(stderr);
      return false;
    }

    {
      int length = spc_.size();
      outfile << length << "\n";
      ConstDofIteratorType enddof = dend ( );
      for(ConstDofIteratorType itdof = dbegin ( );itdof != enddof; ++itdof) 
      {
        outfile << (*itdof) << " ";
      }
      outfile << "\n";
    }

    outfile.close();
    return true;
  }


  template<class DiscreteFunctionSpaceImp>
  bool AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  read_ascii(const std::string fn)
  {
    std::ifstream infile( fn.c_str() ); 
    if ( infile )
    {
      int length;
      infile >> length; 

      if ( length != spc_.size( ))
      {
        DUNE_THROW(InvalidStateException,"Length to read has wrong value!");
      } 

      DofIteratorType enddof = dend ( );
      for(DofIteratorType itdof = dbegin ( );itdof != enddof; ++itdof) 
      {
        infile >> (*itdof); 
      }
      return true;
    }
    return false;
  }
#endif

  template<class DiscreteFunctionSpaceImp>
  void AdaptiveFunctionImplementation<DiscreteFunctionSpaceImp>::
  enableDofCompression() 
  {
    if( memObject_ )
      memObject_->enableDofCompression();
  }
  
} // end namespace Dune
