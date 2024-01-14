#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#define SYM_ASSERT(c) (!(c) && fprintf(stderr, "assertion failed: " #c "\nat " __FILE__ ":%d\n", __LINE__) && (abort(), 0))

#ifdef __cplusplus
}
#endif
