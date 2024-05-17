#ifndef ASSERT_PLATFORM_HH
#define ASSERT_PLATFORM_HH

#include<climits>

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

#endif // ASSERT_PLATFORM_HH
