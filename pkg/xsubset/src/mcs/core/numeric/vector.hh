/**
 * @file vector.hh
 *
 * @author Marc Hofmann
 */

#ifndef MCS_CORE_NUMERIC_VECTOR_HH
#define MCS_CORE_NUMERIC_VECTOR_HH


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Size>
      class range;


      template<typename Value,
	       typename Alloc>
      class subvector;


      template<typename Value,
               typename Alloc = std::allocator<Value>>
      class vector
      {

        friend class subvector<Value, Alloc>;


      public:

        typedef Value value_type;

        typedef Alloc allocator_type;

        typedef typename allocator_type::size_type size_type;

        typedef typename allocator_type::difference_type difference_type;

        typedef typename allocator_type::pointer pointer;

        typedef typename allocator_type::const_pointer const_pointer;

        typedef typename allocator_type::reference reference;

        typedef typename allocator_type::const_reference const_reference;

        typedef range<size_type> range_type;

        typedef subvector<value_type, allocator_type> subvector_type;


      private:

        allocator_type alloc_;

        pointer stor_;

        pointer base_;

        size_type inc_;

        size_type len_;


      public:

        vector();

        explicit
        vector(size_type len);

        template<typename V>
        vector(size_type len,
               V val);

        template<typename V>
        vector(size_type len,
               const V* data);

        vector(const vector& vec);

        vector(vector&& vec);

        vector(subvector_type&& vec);


      private:

        vector(pointer stor,
               size_type len);

        vector(pointer base,
               size_type inc,
               size_type len);


      public:

        ~vector();

        vector&
        operator =(const vector& vec);

        vector&
        operator =(vector&& vec);

        vector&
        operator =(subvector_type&& vec);

        reference
        operator ()(size_type pos);

        const_reference
        operator ()(size_type pos) const;

        subvector_type
        operator ()(range_type rng);

        const subvector_type
        operator ()(range_type rng) const;

        size_type
        inc() const;

        size_type
        len() const;

        template<typename V>
        void
        fill(V val);

        template<typename V>
        void
        fill(const V* data);

        void
        copy(const vector& vec);

        void
        swap(vector& vec);


      private:

        reference
        at(size_type pos) const;

      };


      template<typename Value,
	       typename Alloc>
      class matrix;


      template<typename Value,
               typename Alloc>
      class subvector : public vector<Value, Alloc>
      {

        friend class vector<Value, Alloc>;
	friend class matrix<Value, Alloc>;


      private:

        typedef vector<Value, Alloc> vector_type;

        typedef typename vector_type::size_type size_type;

        typedef typename vector_type::pointer pointer;


      private:

        using vector_type::len_;


      private:

        subvector(pointer base,
                  size_type inc,
                  size_type len);


      public:

        ~subvector();

        subvector&
        operator =(const vector_type& vec);

        subvector&
        operator =(vector_type&& vec);

        using vector_type::copy;

      };


    }

  }

}


#define VECTOR mcs::core::numeric::vector<Value, Alloc>


namespace std
{


  template<typename Value,
	   typename Alloc>
  void
  swap(VECTOR& x, VECTOR& y);


}


#undef VECTOR


#include "vector.cc"
#endif
