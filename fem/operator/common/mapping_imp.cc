
template<typename DFieldType,typename RFieldType, class DType, class RType> 
Mapping<DFieldType,RFieldType,DType,RType> Mapping<DFieldType,RFieldType,DType,RType>::operator+(const MappingType &mapping) const {
  const Mapping<DFieldType,RFieldType,DType,RType> &m = dynamic_cast<const Mapping<DFieldType,RFieldType,DType,RType>& >( mapping );

  Mapping<DFieldType,RFieldType,DType,RType> newMapping = *this;

  for ( typename std::vector<term>::const_iterator it = m.lincomb_.begin(); it != m.lincomb_.end(); it++ ) {
    newMapping.lincomb_.push_back( *it );
  }

  return newMapping;
}

template<typename DFieldType,typename RFieldType, class DType, class RType> 
Mapping<DFieldType,RFieldType,DType,RType> Mapping<DFieldType,RFieldType,DType,RType>::operator-(const MappingType &mapping) const {
  const Mapping<DFieldType,RFieldType,DType,RType> &m = dynamic_cast<const Mapping<DFieldType,RFieldType,DType,RType>& >( mapping );

  Mapping<DFieldType,RFieldType,DType,RType> newMapping = *this;

  for ( typename std::vector<term>::const_iterator it = m.lincomb_.begin(); it != m.lincomb_.end(); it++ ) {
    newMapping.lincomb_.push_back( term( *it->v_, -it->scalar_ ) );
  }

  return newMapping;
}

template<typename DFieldType,typename RFieldType, class DType, class RType> 
Mapping<DFieldType,RFieldType,DType,RType>& Mapping<DFieldType,RFieldType,DType,RType>::operator=(const MappingType &mapping)  {
  const Mapping<DFieldType,RFieldType,DType,RType> &m = dynamic_cast<const Mapping<DFieldType,RFieldType,DType,RType>& >( mapping );

  lincomb_.erase( lincomb_.begin(), lincomb_.end() );

  for ( typename std::vector<term>::const_iterator it = m.lincomb_.begin(); it != m.lincomb_.end(); it++ ) {
    lincomb_.push_back( term( *it->v_, -it->scalar_ ) );
  }

  return *this;
}

template<typename DFieldType,typename RFieldType, class DType, class RType> 
Mapping<DFieldType,RFieldType,DType,RType> Mapping<DFieldType,RFieldType,DType,RType>::operator*(const RangeFieldType &factor) const 
{
  Mapping<DFieldType,RFieldType,DType,RType> newMapping = *this;

  for ( typename std::vector<term>::iterator it = newMapping.lincomb_.begin(); it != newMapping.lincomb_.end(); it++ ) {
    it->scalar_ *= factor;
  }

  return newMapping;
}

template<typename DFieldType,typename RFieldType, class DType, class RType> 
Mapping<DFieldType,RFieldType,DType,RType> Mapping<DFieldType,RFieldType,DType,RType>::operator/(const RangeFieldType &divisor) const {
  Mapping<DFieldType,RFieldType,DType,RType> newMapping = *this;
  for ( typename std::vector<term>::iterator it = newMapping.lincomb_.begin(); it != newMapping.lincomb_.end(); it++ ) {
    it->scalar_ /= divisor;
  }
  return newMapping;
}
