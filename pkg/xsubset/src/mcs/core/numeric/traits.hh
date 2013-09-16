#ifndef _MCS_CORE_NUMERIC_TRAITS_HH_
#define _MCS_CORE_NUMERIC_TRAITS_HH_


#include <cstddef>


namespace mcs     {
namespace core    {
namespace numeric {


// forward declarations
template<typename Value>
class vector;

template<typename Value>
class const_vector;

template<typename Value>
class vector_reference;

template<typename Value>
class const_vector_reference;

template<typename Value>
class matrix;

template<typename Value>
class const_matrix;

template<typename Value>
class matrix_reference;

template<typename Value>
class const_matrix_reference;

template<typename Size>
class subscript;


// traits
template<typename Value>
struct traits
{
  typedef std::size_t size_type;
  typedef Value value_type;
  typedef Value& reference;
  typedef const Value& const_reference;
  typedef Value* pointer;
  typedef const Value* const_pointer;

  typedef numeric::vector<Value> vector;
  typedef numeric::const_vector<Value> const_vector;
  typedef numeric::vector_reference<Value> vector_reference;
  typedef numeric::const_vector_reference<Value> const_vector_reference;

  typedef numeric::matrix<Value> matrix;
  typedef numeric::const_matrix<Value> const_matrix;
  typedef numeric::matrix_reference<Value> matrix_reference;
  typedef numeric::const_matrix_reference<Value> const_matrix_reference;

  typedef numeric::subscript<size_type> subscript;
};


template<typename Value>
struct traits<vector<Value> > : traits<Value>  { };

template<typename Value>
struct traits<const_vector<Value> > : traits<Value>  { };

template<typename Value>
struct traits<vector_reference<Value> > : traits<Value>  { };

template<typename Value>
struct traits<const_vector_reference<Value> > : traits<Value>  { };

template<typename Value>
struct traits<matrix<Value> > : traits<Value>  { };

template<typename Value>
struct traits<const_matrix<Value> > : traits<Value>  { };

template<typename Value>
struct traits<matrix_reference<Value> > : traits<Value>  { };

template<typename Value>
struct traits<const_matrix_reference<Value> > : traits<Value>  { };


}  // namespace numeric
}  // namespace core
}  // namespace mcs


#endif
