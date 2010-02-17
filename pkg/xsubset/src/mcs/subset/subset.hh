/**
 * @file subset.hh
 *
 * @author Marc Hofmann
 */

#ifndef MCS_SUBSET_SUBSET_HH
#define MCS_SUBSET_SUBSET_HH


#include <vector>


namespace mcs
{

  namespace subset
  {


    template<typename Variable,
	     typename Value>
    class subset
    {

    private:

      typedef std::vector<Variable> Container;


    public:

      typedef Variable variable_type;

      typedef Value value_type;

      typedef typename Container::size_type size_type;

      typedef typename Container::iterator iterator;

      typedef typename Container::const_iterator const_iterator;


    private:

      Container variables_;

      value_type value_;


    public:

      subset();

      subset(size_type max_size);

      subset(size_type size,
	     value_type val);

      size_type
      get_max_size() const;

      size_type
      get_size() const;

      template<typename Iterator>
      void
      assign(Iterator first,
	     Iterator last);

      template<typename Iterator>
      void
      append(Iterator first,
	     Iterator last);

      iterator
      begin();

      const_iterator
      begin() const;

      iterator
      end();

      const_iterator
      end() const;

      void
      set_value(value_type value);

      value_type
      get_value() const;

      void
      swap(subset& other);

    };


  }

}


#define SUBSET mcs::subset::subset<Variable, Value>


namespace std
{


  template<typename Variable,
	   typename Value>
  void
  swap(SUBSET& x, SUBSET& y);


}


#undef SUBSET


#include "subset.cc"
#endif
