#include "types.h"

#include <math.h>

void rot3_from_tangent(f64* rot_tangent, f64* rot_storage, f64 epsilon) {
    f64 tmp0 = sqrt(epsilon * epsilon + rot_tangent[0] * rot_tangent[0] + rot_tangent[1] * rot_tangent[1] + rot_tangent[2] * rot_tangent[2]);
    f64 tmp1 = 0.5 * tmp0;
    f64 tmp2 = sin(tmp1) / tmp0;

    rot_storage[0] = tmp2 * rot_tangent[0];
    rot_storage[1] = tmp2 * rot_tangent[1];
    rot_storage[2] = tmp2 * rot_tangent[2];
    rot_storage[3] = cos(tmp1);
}
