/**
 * @file subset.cc
 *
 * @author Marc Hofmann
 */

#ifndef MCS_SUBSET_SUBSET_CC
#define MCS_SUBSET_SUBSET_CC


#include <limits>
#include <iterator>
#include <numeric>

#include "subset.hh"


#define SUBSET subset<Variable, Value>


namespace mcs
{

  namespace subset
  {


    template<typename Variable,
	     typename Value>
    SUBSET::subset() :
      variables_(),
      value_()
    {
    }


    template<typename Variable,
	     typename Value>
    SUBSET::subset(const size_type size) :
      variables_(),
      value_(std::numeric_limits<value_type>::max())
    {
      variables_.reserve(size);
    }


    template<typename Variable,
	     typename Value>
    SUBSET::subset(const size_type size,
		   const value_type val) :
      variables_(size),
      value_(val)
    {
      std::iota(variables_.begin(), variables_.end(), 0);
    }


    template<typename Variable,
	     typename Value>
    typename SUBSET::size_type
    SUBSET::get_max_size() const
    {
      return variables_.capacity();
    }


    template<typename Variable,
	     typename Value>
    typename SUBSET::size_type
    SUBSET::get_size() const
    {
      return variables_.size();
    }


    template<typename Variable,
	     typename Value>
    template<typename Iterator>
    void
    SUBSET::assign(const Iterator first,
		   const Iterator last)
    {
      variables_.assign(first, last);
    }


    template<typename Variable,
	     typename Value>
    template<typename Iterator>
    void
    SUBSET::append(const Iterator first,
		   const Iterator last)
    {
      std::copy(first, last, std::back_inserter(variables_));
    }


    template<typename Variable, typename Value>
    typename SUBSET::iterator
    SUBSET::begin()
    {
      return variables_.begin();
    }


    template<typename Variable,
	     typename Value>
    typename SUBSET::const_iterator
    SUBSET::begin() const
    {
      return variables_.begin();
    }


    template<typename Variable,
	     typename Value>
    typename SUBSET::iterator
    SUBSET::end()
    {
      return variables_.end();
    }


    template<typename Variable,
	     typename Value>
    typename SUBSET::const_iterator
    SUBSET::end() const
    {
      return variables_.end();
    }


    template<typename Variable,
	     typename Value>
    void
    SUBSET::set_value(const value_type value)
    {
      value_ = value;
    }


    template<typename Variable,
	     typename Value>
    typename SUBSET::value_type
    SUBSET::get_value() const
    {
      return value_;
    }


    template<typename Variable,
	     typename Value>
    void
    SUBSET::swap(subset& other)
    {
      std::swap(variables_, other.variables_);
      std::swap(value_, other.value_);
    }


  }

}


#undef SUBSET


#define SUBSET mcs::subset::subset<Variable, Value>


namespace std
{


  template<typename Variable,
	   typename Value>
  void
  swap(SUBSET& x, SUBSET& y)
  {
    x.swap(y);
  }


}


#undef SUBSET


#endif
