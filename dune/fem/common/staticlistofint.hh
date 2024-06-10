#ifndef DUNE_FEM_STATICLISTOFINTS_HH
#define DUNE_FEM_STATICLISTOFINTS_HH

#include <iostream>
#include <sstream>
#include <map>
#include <cassert>

// macro that defines a struct containing a list of static const ints
// and the names of those ints as strings. This is mainly intended for
// parameter identifiers, e.g. in dune/fem/solver/parameter.hh
#define LIST_OF_INT(ListName, ...) \
struct ListName {\
static constexpr int __VA_ARGS__; \
typedef std::pair< std::map<int, int>, std::vector<std::string> > EntriesType; \
static inline EntriesType& entries() { \
static EntriesType idMap; \
static bool initialized = false;\
if( !initialized ){\
std::string str = #__VA_ARGS__; \
int len = str.length(); \
std::ostringstream temp; \
int id = 0;\
int n = 0;\
for(int i = 0; i < len; ++i) { \
  if(str[i] == '=') { \
    size_t c = str.find(',', i+1);\
    std::string number = str.substr(i+1, c-(i+1));\
    id = stoi(number);\
    while( i < len && str[i] != ',') \
     ++i;\
  }\
  if(isspace(str[i])) continue; \
  else if(str[i] == ',') { \
  idMap.first[ id ] = n++; \
  idMap.second.push_back(temp.str()); \
  temp.str(std::string());\
  } \
  else if (str[i] == '_') {\
    temp << "-";\
  }\
  else temp<< str[i]; \
} \
idMap.first[id] = n; \
idMap.second.push_back( temp.str() );\
idMap.second.back().pop_back();\
initialized = true;\
}\
return idMap;} \
static inline int to_id(int value){\
  for( const auto& item : entries().first ){\
    if( item.second == value )\
      return item.first;\
  }\
  assert(false); \
  return -1;\
}\
static inline const std::vector<std::string>& names() { \
  return entries().second;}\
static int length() { return entries().second.size(); }\
static std::string to_string(int value){\
  auto it = entries().first.find( value );\
  assert( it != entries().first.end() );\
  assert( it->second < int(entries().second.size()) );\
  return entries().second[ it->second ];\
}\
};

// same as LIST_OF_INT. In addition the integers are forwarded to the
// current namespace. This is needed for SolverParameters.
#define LIST_OF_INT_FORWARDED(ListName, ...) \
  LIST_OF_INT(ListName, __VA_ARGS__);\
  static const int __VA_ARGS__;

#endif
