import torch
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
from argparse import Namespace
import pyPBAD


def get_test_config():
    # DATA 配置
    data_cfg = Namespace()
    data_cfg.DT = 0.01
    data_cfg.FREQ = -1
    data_cfg.GRAVITY = np.array([0, 9.81, 0])
    data_cfg.GROUND_Y = -10.0

    # 主配置
    cfg = Namespace()
    cfg.N_STEP = 5
    cfg.DATA = data_cfg
    return cfg


def generate_cube_points():
    """生成立方体点"""
    vertices = np.array([
        [-1, -1, -1], [1, -1, -1], [-1, 1, -1], [1, 1, -1],
        [-1, -1, 1], [1, -1, 1], [-1, 1, 1], [1, 1, 1]
    ], dtype=np.float64)

    # 重塑为[m,1]格式
    vertices_reshaped = vertices.reshape(-1, 1)
    return vertices_reshaped


def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")

    # 获取配置
    cfg = get_test_config()

    # 初始化PBD仿真器
    sim = pyPBAD.ConvHullPBDSimulator(cfg.DATA.DT)
    #render = pyPBAD.SimulatorVisualizer(sim)

    # 创建场景
    # 设置地面
    floor = pyPBAD.BBoxExact(
        np.array([-15, cfg.DATA.GROUND_Y, -15], dtype=np.float64),
        np.array([15, cfg.DATA.GROUND_Y + 0.1, 15], dtype=np.float64)
    )
    
    # 设置场景
    body=pyPBAD.ArticulatedLoader.createInitJoints(8,1000)
    sim.setArticulatedBody(body)
    sim.addShape(floor)
    sim.gravity=cfg.DATA.GRAVITY
    
    # 获取点云数据
    points = generate_cube_points()
    print(f"Points shape: {points.shape}")
    print(points)
    sim.updateConvexPoints(points)
    points = torch.from_numpy(sim.getConvexPoints()).to(device)
    points.requires_grad_(False)
    print(f"Generated points shape: {points.shape}")

    # 仿真主循环
    n_frames = 1000
    traj = []
    sim.reset()

    # 初始状态
    pos = sim.pos()
    pos = torch.from_numpy(pos).to(device)
    last_pos = pos.clone()
    vel = torch.zeros_like(pos)

    for frame_id in range(n_frames):
        # 仿真步进
        for _ in range(cfg.N_STEP):
            new_pos = torch.from_numpy(sim.step(pos.cpu().numpy(), last_pos.cpu().numpy())).to(device)
            last_pos = pos
            pos = new_pos
            vel = (pos - last_pos) / cfg.DATA.DT

        if frame_id % 10 == 0:
            print(f"Frame {frame_id}: ")
            print(f"Max velocity: {vel.abs().max().item():.4f}")
            print(f"Mean position: {pos.mean().item():.4f}")

            # 可视化当前帧
            plot_points(pos.detach().cpu().numpy(), frame_id, cfg)
            traj.append(pos.detach().cpu().numpy())

    # 生成动画
    render.visualize(traj, 0)


def plot_points(points, frame_id, cfg):
    """绘制点云"""
    plt.clf()
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111, projection='3d')

    # 绘制顶点
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], c='b', marker='o', s=20)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Frame {frame_id}')

    # 坐标轴范围
    ax.set_xlim([-5, 5])
    ax.set_ylim([5, -15])
    ax.set_zlim([-5, 5])

    # 绘制地面
    xx, zz = np.meshgrid([-5, 5], [-5, 5])
    yy = np.ones_like(xx) * cfg.DATA.GROUND_Y
    ax.plot_surface(xx, yy, zz, alpha=0.2, color='gray')

    # 视角设置
    ax.view_init(elev=25, azim=45)
    ax.dist = 12
    ax.grid(True)

    # 保存
    os.makedirs('results', exist_ok=True)
    plt.savefig(f'results/frame_{frame_id:04d}.png', dpi=150)
    plt.close()


if __name__ == '__main__':
    main()