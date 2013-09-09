#ifndef _MCS_CORE_NUMERIC_SUBSCRIPT_HH_
#define _MCS_CORE_NUMERIC_SUBSCRIPT_HH_


namespace mcs     {
namespace core    {
namespace numeric {


template<typename Size>
struct subscript
{
  typedef Size size_type;

  size_type pos;
  size_type len;
  size_type inc;

  explicit
  subscript(size_type pos) : subscript(pos, 1)
  {
  }
  subscript(size_type pos, size_type len) : subscript(pos, len, 1)
  {
  }
  subscript(size_type pos, size_type len, size_type inc) :
    pos(pos),
    len(len),
    inc(inc)
  {
    assert(0 <= pos);
    assert(0 <= len);
    assert(0 <= inc);
  }
};


}  // namespace numeric
}  // namespace core
}  // namespace mcs


#endif
