#ifndef DYNAMICAL_OBJECT_HPP
#define DYNAMICAL_OBJECT_HPP


#include <cassert>
#include <iostream>


class DynamicalObject
{
public:
  void set_eta(double eta_lo, double eta_hi);
  void set_output(std::ostream &os);

protected:
  DynamicalObject(const char name[], int id, int components = 1);
  virtual ~DynamicalObject();

  int id() const;

  int size(int component = 0) const;
  int new_size(int requested_new_size, int component = 0);

  virtual void resize(int new_size, int component) = 0;

private:
  const int components;  
  int *_size;
  char *name;
  const int _id;
  double eta_lo, eta_hi;

  std::ostream *os;
};




// ======== inline implementation
inline
DynamicalObject::DynamicalObject(const char name[], int id, int components) : 
  components(components), _size(new int[components]), _id(id), os(NULL) 
{
  assert(_size);
  for(int i=0; i<components; i++) _size[i] = 0;

  this->name = new char[strlen(name) + 1];
  strcpy(this->name, name);

  // set this to some useful values
  eta_lo = 0.2;
  eta_hi = 0.1;
}



inline 
DynamicalObject::~DynamicalObject() 
{
  delete[] _size;
  delete[] name;
}


inline
int DynamicalObject::id() const
{
  return _id;
}


inline
int DynamicalObject::size(int component) const
{
  return _size[component];
}


// eta_lo = 0.0, eta_hi = 0.0: exact size
// eta_lo = 1.0, eta_hi = 0.0: enlarge
inline
void DynamicalObject::set_eta(double eta_lo, double eta_hi)
{
  assert(eta_lo >= 0.0 && eta_lo <= 1.0 && eta_hi >= 0.0);
  this->eta_lo = eta_lo;
  this->eta_hi = eta_hi;
}


inline
void DynamicalObject::set_output(std::ostream &os) 
{ 
  this->os = &os;
}


// todo: eta_hi not used! why? -> check
// new_size >= requested_new_size
inline 
int DynamicalObject::new_size(int requested_new_size, int component)
{
  const int add_size = static_cast<int>(eta_lo * requested_new_size);

  if (requested_new_size > _size[component]
      || requested_new_size < (1.0-eta_lo)*_size[component]){

    resize(requested_new_size + add_size, component);

    if (os){
      if (_id >= 0){
	*os << name << " " << _id << ":   " 
	    << "component: " << component << "   "
	    << _size[component] << " -> " << requested_new_size+add_size
	    << std::endl;
      }
      else{
	*os << name << " " << "unknown id:   " 
	    << "component: " << component << "   "
	    << _size[component] << " -> " << requested_new_size+add_size
	    << std::endl;
      }
    }

    _size[component] = requested_new_size + add_size;
  }
  return _size[component];
}


#endif
