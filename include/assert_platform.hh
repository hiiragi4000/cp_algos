#ifndef ASSERT_PLATFORM_HH
#define ASSERT_PLATFORM_HH

#include<climits>
#include<cstdlib>

#if CHAR_BIT != 8
#error "We're not ready for platforms whose char size is other than 8 QQ"
#endif

#if (-1&3) != 3
#error "We're not ready for platforms whose negative integers aren't of 2's complement QQ"
#endif

#if -1>>1 != -1
#error "We're not ready for platforms whose >> isn't SAR for signed integers QQ"
#endif

#define I16 short
#define I64 long long

static_assert(
   sizeof(I16)==2 && sizeof(int)==4 && (sizeof(long)==4 || sizeof(long)==8) && sizeof(I64)==8,
   "We're not ready for platforms with exotic sizes of integers QQ"
);

#undef I64
#undef I16

#define U8 unsigned char
#define U32 unsigned

namespace impl{
struct AssertLittleEndian{
   AssertLittleEndian(){
      U32 buf = 0x01020304;
      U8 *p = reinterpret_cast<U8*>(&buf);
      if(p[0]!=4 || p[1]!=3 || p[2]!=2 || p[3]!=1){
         std::abort();
      }
   }
};
inline AssertLittleEndian assert_little_endian;
}

#undef U32
#undef U8

#endif // ASSERT_PLATFORM_HH
