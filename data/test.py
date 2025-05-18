'''
import numpy as np

def generate_uniform_sphere_mesh(center, radius, num_points, output_file):

    # 生成球面上的均匀分布点（Fibonacci Sphere）
    vertices = []
    for i in range(num_points):
        z = 1 - (2 * i + 1) / num_points  # 均匀分布的 z 坐标
        r = np.sqrt(1 - z**2)            # 对应 z 的半径
        theta = 2 * np.pi * i / ((1 + np.sqrt(5)) / 2)  # 黄金角度分布
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        vertices.append((x * radius + center[0], 
                         y * radius + center[1], 
                         z * radius + center[2]))

    # 简单生成三角形索引（基于 Delaunay 三角化）
    from scipy.spatial import ConvexHull
    hull = ConvexHull(vertices)  # 计算凸包
    faces = hull.simplices       # 凸包的三角形面

    # 保存为 .obj 文件
    with open(output_file, 'w') as f:
        f.write("# Sphere Mesh with Uniformly Distributed Points\n")
        for v in vertices:
            f.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        for face in faces:
            f.write(f"f {face[0] + 1} {face[1] + 1} {face[2] + 1}\n")  # OBJ 顶点索引从 1 开始

    print(f"Mesh saved to {output_file}")


# 参数设置
center = (0.0, 0.0, 0.0)   # 球心
radius = 0.5               # 半径
num_points = 196            # 球表面均匀分布的点数量#156
output_file = "196.obj"  # 输出文件名

# 生成球体网格
generate_uniform_sphere_mesh(center, radius, num_points, output_file)
'''

import numpy as np

def generate_cuboid_mesh_equal_triangles(dimensions, subdivisions, output_file):
    """
    生成一个长方体网格，表面均匀分布的三角形面片，确保每个三角形的两条直角边尽量等长。

    Args:
        dimensions (tuple): 长方体的长、宽、高 (L, W, H)。
        subdivisions (tuple): 每个面在长、宽方向的划分数 (nx, ny)。
        output_file (str): 输出的 .obj 文件路径。
    """
    # 获取长方体的长、宽、高
    L, W, H = dimensions

    # 定义长方体的 8 个顶点坐标
    vertices = [
        (0, 0, 0),       # v1
        (L, 0, 0),       # v2
        (L, W, 0),       # v3
        (0, W, 0),       # v4
        (0, 0, H),       # v5
        (L, 0, H),       # v6
        (L, W, H),       # v7
        (0, W, H),       # v8
    ]

    # 定义长方体的 6 个矩形面（每个面由 4 个顶点组成）
    faces = {
        "bottom": [(0, 1, 2, 3)],  # 底面
        "top": [(4, 5, 6, 7)],     # 顶面
        "front": [(0, 1, 5, 4)],   # 前面
        "back": [(3, 2, 6, 7)],    # 后面
        "left": [(0, 3, 7, 4)],    # 左面
        "right": [(1, 2, 6, 5)],   # 右面
    }

    # 初始化全局顶点和面列表
    global_vertices = []
    global_faces = []

    # 对每个面进行细分
    for face_name in faces:
        for quad in faces[face_name]:
            # 获取矩形面四个顶点
            v0, v1, v2, v3 = [vertices[idx] for idx in quad]
            # 生成均匀网格
            subdivided_vertices, subdivided_faces = subdivide_rectangle(v0, v1, v2, v3, subdivisions)
            # 更新全局顶点和面
            offset = len(global_vertices)
            global_vertices.extend(subdivided_vertices)
            global_faces.extend([(f[0] + offset, f[1] + offset, f[2] + offset) for f in subdivided_faces])

    # 将所有顶点和面写入 .obj 文件
    with open(output_file, 'w') as f:
        f.write("# Cuboid Mesh with Uniformly Distributed Triangles\n")
        for v in global_vertices:
            f.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        for face in global_faces:
            f.write(f"f {face[0] + 1} {face[1] + 1} {face[2] + 1}\n")  # OBJ 文件索引从 1 开始

    print(f"Mesh saved to {output_file}")


def subdivide_rectangle(v0, v1, v2, v3, subdivisions):
    """
    将矩形面细分为均匀网格，并生成三角形面片。

    Args:
        v0, v1, v2, v3 (tuple): 矩形面四个顶点 (x, y, z)。
        subdivisions (tuple): 在两个方向的划分数 (nx, ny)。

    Returns:
        subdivided_vertices (list): 细分后的顶点列表。
        subdivided_faces (list): 细分后的三角形面片列表。
    """
    nx, ny = subdivisions  # 获取两个方向的划分数
    subdivided_vertices = []
    subdivided_faces = []

    # 生成均匀网格顶点
    for i in range(nx + 1):
        for j in range(ny + 1):
            alpha = i / nx
            beta = j / ny
            # 插值计算顶点位置
            new_vertex = (
                (1 - alpha) * (1 - beta) * np.array(v0) +
                alpha * (1 - beta) * np.array(v1) +
                alpha * beta * np.array(v2) +
                (1 - alpha) * beta * np.array(v3)
            )
            subdivided_vertices.append(tuple(new_vertex))

    # 生成三角形面片
    for i in range(nx):
        for j in range(ny):
            # 当前网格单元的四个顶点索引
            v00 = i * (ny + 1) + j
            v10 = (i + 1) * (ny + 1) + j
            v01 = i * (ny + 1) + (j + 1)
            v11 = (i + 1) * (ny + 1) + (j + 1)

            # 将网格单元划分为两个三角形
            subdivided_faces.append((v00, v10, v11))  # 第一三角形
            subdivided_faces.append((v00, v11, v01))  # 第二三角形

    return subdivided_vertices, subdivided_faces


# 参数设置
dimensions = (1.0, 0.1, 0.1)  # 长方体的长、宽、高
subdivisions = (10, 2)        # 每个面在长、宽方向的划分数
output_file = "cuboid.obj"  # 输出文件路径

# 生成长方体网格
generate_cuboid_mesh_equal_triangles(dimensions, subdivisions, output_file)
