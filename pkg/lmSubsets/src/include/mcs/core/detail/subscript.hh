#ifndef MCS_CORE_DETAIL_SUBSCRIPT_HH
#define MCS_CORE_DETAIL_SUBSCRIPT_HH



namespace mcs    {
namespace core   {
namespace detail {



struct subscript
{

    int off;

    int len;



    constexpr
    subscript(
        const int off,
        const int len
    ) noexcept :
        off(off),
        len(len)
    {
    }

};



}  // end namespace detail
}  // end namespace core
}  // end namespace mcs



#endif
