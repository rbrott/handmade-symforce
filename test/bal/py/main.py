import numpy as np
import bal
import sym

import open3d as o3d
import sys


# So slow! Takes 5-10 seconds to load even a relatively small problem.
def read_problem(path):
    with open(path, "r") as file:
        num_cameras, num_points, num_observations = map(int, file.readline().split())

        camera_indices = []
        point_indices = []
        pixels = []

        for _ in range(num_observations):
            camera, point, px, py = map(float, file.readline().split())
            camera_indices.append(camera)
            point_indices.append(point)
            pixels.append(np.array([px, py]))

        camera_poses = []
        intrinsics = []

        for _ in range(num_cameras):
            rx = float(file.readline())
            ry = float(file.readline())
            rz = float(file.readline())
            tx = float(file.readline())
            ty = float(file.readline())
            tz = float(file.readline())
            f = float(file.readline())
            k1 = float(file.readline())
            k2 = float(file.readline())
            camera_poses.append((np.array([rx, ry, rz]), np.array([tx, ty, tz])))
            intrinsics.append(np.array([f, k1, k2]))

        points = []

        for _ in range(num_points):
            x = float(file.readline())
            y = float(file.readline())
            z = float(file.readline())
            points.append(np.array([x, y, z]))

    return {
        "num_cameras": num_cameras,
        "num_points": num_points,
        "num_observations": num_observations,
        "camera_indices": camera_indices,
        "point_indices": point_indices,
        "pixels": pixels,
        "camera_poses": camera_poses,
        "intrinsics": intrinsics,
        "points": points,
    }


def draw_scene():
    cam_T_world, intrinsics, point, camera_indices, point_indices, pixels = bal.read_problem(sys.argv[1])

    # values = np.fromfile("/tmp/values.bin", dtype=np.float64) 
    # camera_values = values[:10 * cam_T_world.shape[0]].reshape(-1, 10)
    # point_values = values[10 * cam_T_world.shape[0]:].reshape(-1, 3)

    # cam_T_world = camera_values[:, :7]
    # intrinsics = camera_values[:, 7:]
    # point = point_values

    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(point)

    pcd2 = o3d.geometry.PointCloud()
    pcd2.points = o3d.utility.Vector3dVector(cam_T_world[:, 4:])

    # selected_camera = 1
    # camera_lines = o3d.geometry.LineSet()
    # camera_lines.points = o3d.utility.Vector3dVector(
    #     np.vstack([
    #         cam_T_world[selected_camera, 4:],
    #         point[point_indices[camera_indices == selected_camera]]
    #     ])
    # )
    # camera_lines.lines = o3d.utility.Vector2iVector([
    #     [0, i] for i in range(1, len(camera_lines.points))
    # ])
    
    pcd.paint_uniform_color([1.0, 0, 0])
    pcd2.paint_uniform_color([0, 1.0, 0])

    o3d.visualization.draw_geometries([pcd, pcd2])


def debug_bad_camera():
    cam_T_world, intrinsics, point, camera_indices, point_indices, pixels = bal.read_problem(sys.argv[1])
    
    norms = np.linalg.norm(cam_T_world[:, 4:], axis=1)
    sort_perm = np.argsort(norms)

    # print(norms)

    # print(sort_perm[::-1])

    # print(intrinsics[238])

    return

    # print(cam_T_world[0])
    # print(intrinsics[0])
    # return

    res_norm_by_camera = np.zeros(cam_T_world.shape[0])
    for i in range(len(camera_indices)):
        # if camera_indices[i] != 72:
        #     continue

        res, H_lower, rhs = bal.snavely_reprojection_factor(cam_T_world[camera_indices[i]], intrinsics[camera_indices[i]], point[point_indices[i]], pixels[i])
        # H = H_lower + H_lower.T - np.diag(np.diag(H_lower))
        # H += np.eye(12) * np.max(np.diag(H)) * 1e-12
        # print(H)
        
        # print(np.linalg.norm(res))
        res_norm_by_camera[camera_indices[i]] += res[0] ** 2 + res[1] ** 2
        
        # try:
        #     np.linalg.solve(H, rhs)
        # except np.linalg.LinAlgError:
        #     print("singular matrix")
        #     print(H)
        
        # print(res)

        # H = H_lower + H_lower.T - np.diag(np.diag(H_lower)) + np.eye(12) * 1e-10
        # if np.linalg.cond(H) > 1e10:
        #     print(i)

    print(res_norm_by_camera)
    print(res_norm_by_camera[72])

    # res, H_lower, rhs = bal.snavely_reprojection_factor(cam_T_world[camera_indices[0]], intrinsics[camera_indices[0]], point[point_indices[0]], pixels[0])
    # H = H_lower + H_lower.T - np.diag(np.diag(H_lower)) + np.eye(12) * 1e-10
    # print(H)
    # print(rhs)
    # print(np.linalg.solve(H, rhs))

    # // 100+ have problems with camera at index 72


if __name__ == "__main__":
    # debug_bad_camera()
    # draw_scene()

    cam_T_world, intrinsics, point, camera_indices, point_indices, pixels = bal.read_problem(sys.argv[1])
    cam_T_world = cam_T_world[camera_indices[0]]
    intrinsics = intrinsics[camera_indices[0]]
    point = point[point_indices[0]]
    pixels = pixels[0]

    # print(cam_T_world)
    # print(intrinsics)
    # print(point)
    # print(pixels)

    # res, H_lower, rhs = bal.snavely_reprojection_factor(cam_T_world[camera_indices[0]], intrinsics[camera_indices[0]], point[point_indices[0]], pixels[0])
    # res2, H_lower2, rhs2 = bal.snavely_reprojection_factor2(cam_T_world[camera_indices[0]], intrinsics[camera_indices[0]], point[point_indices[0]], pixels[0])

    # cam_T_world = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    # intrinsics = np.array([1.0, 0.0, 0.0])
    # point = np.array([1.0, 2.0, 10.0])
    # pixels = np.array([1.0, 0.0]) 

    focal_length, k1, k2 = intrinsics
    point_cam = sym.Pose3.from_storage(cam_T_world) * point

    p = -point_cam[:2] / point_cam[2]
    n = p.T @ p
    r = 1 + k1 * n + k2 * n ** 2
    pixel_projected = focal_length * r * p
    res = pixel_projected - pixels

    res3, H_lower, rhs = bal.snavely_reprojection_factor(cam_T_world, intrinsics, point, pixels)

    res2, H_lower2, rhs2 = bal.snavely_reprojection_factor2(cam_T_world, intrinsics, point, pixels)

    # print(res)
    # print(res2, rhs)
    # print(res3, rhs2)

    print(H_lower)
    print(H_lower2)

