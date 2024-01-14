#include "types.h"

#include <math.h>

void sym_pose3_retract_in_place(f64* pose_storage, f64* pose_delta, f64 epsilon) {
    f64 a0 = pose_storage[0];
    f64 a1 = pose_storage[1];
    f64 a2 = pose_storage[2];
    f64 a3 = pose_storage[3];
    f64 a4 = pose_storage[4];
    f64 a5 = pose_storage[5];
    f64 a6 = pose_storage[6];

    f64 tmp0 = sqrt(epsilon * epsilon + pose_delta[0] * pose_delta[0] + pose_delta[1] * pose_delta[1] + pose_delta[2] * pose_delta[2]);
    f64 tmp1 = 0.5 * tmp0;
    f64 tmp2 = sin(tmp1) / tmp0;
    f64 tmp3 = tmp2 * pose_delta[1];
    f64 tmp4 = tmp2 * pose_delta[2];
    f64 tmp5 = tmp2 * pose_delta[0];
    f64 tmp6 = cos(tmp1);
    f64 tmp7 = a0 * tmp2; 

    pose_storage[0] = a0 * tmp6 + a1 * tmp4 - a2 * tmp3 + a3 * tmp5;
    pose_storage[1] = a1 * tmp6 + a2 * tmp5 + a3 * tmp3 - tmp7 * pose_delta[2];
    pose_storage[2] = -a1 * tmp5 + a2 * tmp6 + a3 * tmp4 + tmp7 * pose_delta[1];
    pose_storage[3] = -a1 * tmp3 - a2 * tmp4 + a3 * tmp6 - tmp7 * pose_delta[0];
    pose_storage[4] = a4 + pose_delta[3];
    pose_storage[5] = a5 + pose_delta[4];
    pose_storage[6] = a6 + pose_delta[5];
}
