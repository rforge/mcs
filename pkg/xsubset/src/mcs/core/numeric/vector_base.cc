/**
 * @file vector_base.cc
 */
#ifndef MCS_CORE_NUMERIC_VECTOR_BASE_CC
#define MCS_CORE_NUMERIC_VECTOR_BASE_CC


#include <utility>

#include "../../mcs.hh"

#include "slice.hh"
#include "vector_base.hh"
#include "slice_vector.hh"


#define VECTOR_BASE vector_base<Value,          \
                                Size,           \
                                Derived>


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      VECTOR_BASE::vector_base() :
        base_(0),
        len_(0)
      {
        //std::cout << "vector_base::vector_base()" << std::endl;
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      VECTOR_BASE::vector_base(Value* const base,
                               const Size len) :
        base_(base),
        len_(len)
      {
        //std::cout << "vector_base::vector_base(Value*,Size)" << std::endl;
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      VECTOR_BASE::vector_base(const vector_base<Value, Size, Derived>& x) :
        base_(x.base_),
        len_(x.len_)
      {
        //std::cout << "vector_base::vector_base(const vector_base&)" << std::endl;
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      VECTOR_BASE::~vector_base()
      {
        //std::cout << "vector_base::~vector_base()" << std::endl;

        // TODO
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      Value*
      VECTOR_BASE::base()
      {
        return base_;
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      const Value*
      VECTOR_BASE::base() const
      {
        return base_;
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      Size
      VECTOR_BASE::inc() const
      {
        return static_cast<const Derived<Value, Size>*>(this)->inc();
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      Size
      VECTOR_BASE::len() const
      {
        return len_;
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      Value&
      VECTOR_BASE::at(const Size pos)
      {
        return static_cast<Derived<Value, Size>*>(this)->at(pos);
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      const Value&
      VECTOR_BASE::at(const Size pos) const
      {
        return static_cast<const Derived<Value, Size>*>(this)->at(pos);
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      slice_vector<Value, Size>
      VECTOR_BASE::at(const slice<Size>& s)
      {
        return static_cast<Derived<Value, Size>*>(this)->at(s);
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      const slice_vector<Value, Size>
      VECTOR_BASE::at(const slice<Size>& s) const
      {
        return static_cast<const Derived<Value, Size>*>(this)->at(s);
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      template<template<typename V,
                        typename S>
               class D>
      void
      VECTOR_BASE::copy(const vector_base<Value, Size, D>& x)
      {
        //std::cout << "vector_base::copy(const vector_base&)" << std::endl;

        static_cast<Derived<Value, Size>*>(this)->copy(x);
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      template<template<typename V,
                        typename S>
               class D>
      void
      VECTOR_BASE::swap(vector_base<Value, Size, D>&& x)
      {
        //std::cout << "vector_base::swap(vector_base&&)" << std::endl;

        static_cast<Derived<Value, Size>*>(this)->swap(x);
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      void
      VECTOR_BASE::fill(const Value x)
      {
        //std::cout << "vector_base::fill(Value)" << std::endl;

        static_cast<Derived<Value, Size>*>(this)->fill(x);
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      vector_base<Value, Size, Derived>&
      VECTOR_BASE::operator =(const vector_base<Value, Size, Derived>& x)
      {
        //std::cout << "vector_base::operator=(const vector_base&) 1" << std::endl;

        return static_cast<Derived<Value, Size>*>(this)->operator =(*static_cast<const Derived<Value, Size>*>(&x));
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      template<template<typename V,
                        typename S>
               class D>
      vector_base<Value, Size, Derived>&
      VECTOR_BASE::operator =(const vector_base<Value, Size, D>& x)
      {
        //std::cout << "vector_base::operator=(const vector_base&) 2" << std::endl;

        return static_cast<Derived<Value, Size>*>(this)->operator =(*static_cast<const D<Value, Size>*>(&x));
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      template<template<typename V,
                        typename S>
               class D>
      vector_base<Value, Size, Derived>&
      VECTOR_BASE::operator =(vector_base<Value, Size, D>&& x)
      {
        //std::cout << "vector_base::operator=(vector_base&&)" << std::endl;

        return static_cast<Derived<Value, Size>*>(this)->operator =(std::move(*static_cast<const D<Value, Size>*>(&x)));
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      vector_base<Value, Size, Derived>&
      VECTOR_BASE::operator =(const Value x)
      {
        //std::cout << "vector_base::operator=(Value)" << std::endl;

        fill(x);
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      Value&
      VECTOR_BASE::operator ()(const Size pos)
      {
        return at(pos);
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      const Value&
      VECTOR_BASE::operator ()(const Size pos) const
      {
        return at(pos);
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      slice_vector<Value, Size>
      VECTOR_BASE::operator ()(const numeric::slice<Size>& pos)
      {
        return at(pos);
      }


      template<typename Value,
               typename Size,
               template<typename V,
                        typename S>
               class Derived>
      const slice_vector<Value, Size>
      VECTOR_BASE::operator ()(const numeric::slice<Size>& pos) const
      {
        return at(pos);
      }


    }

  }

}


#undef VECTOR_BASE


namespace std
{


  template<typename Value,
           typename Size,
           template<typename V,
                    typename S>
           class Derived>
  void
  swap(mcs::core::numeric::vector_base<Value, Size, Derived>&& x,
       mcs::core::numeric::vector_base<Value, Size, Derived>&& y)
  {
    x.swap(y);
  }


}


#endif
