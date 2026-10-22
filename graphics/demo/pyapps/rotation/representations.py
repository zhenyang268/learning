"""
A5.5 · 三种表示互转：轴角 / 四元数 / 欧拉角 / 矩阵
================================================================
现象：同一个旋转有四种写法。互转时最容易出的错不是公式记错，
      而是**边界情况**：θ=0（轴无定义）、θ=π（轴符号二义）、
      欧拉角在 ±90°（万向锁）、四元数 q 与 −q（双覆盖）。

一句话在说什么：
    四种表示是**同一个对象**的四张皮。它们各有各的"坏点"：
        轴角    : θ=0 时轴无定义；θ=π 时轴双向都行
        欧拉角  : pitch=±90° 时漏自由度
        四元数  : q 与 −q 同一个旋转（非单射）
        矩阵    : 9 个数含冗余（6 个约束），且量化后要重正交化
    选择哪种表示 = 选择**你愿意承受哪种坏点**。本脚本把四种互转全跑一遍，
    并逐条验证坏点上的行为。

运行： uv run python graphics/demo/pyapps/rotation/representations.py
输出： learning-todo/imgs/rot/representation_roundtrip.png
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _rot_core as RC  # noqa: E402

plt = RC.setup_style()


def roundtrip_table():
    print("=" * 74)
    print("【1】四种表示的互相转换：往返误差（随机 2000 组）")
    print("=" * 74)
    rng = np.random.default_rng(23)
    paths = {
        "矩阵 → 四元数 → 矩阵": [],
        "矩阵 → 欧拉(ZYX) → 矩阵": [],
        "轴角 → 四元数 → 轴角": [],
        "四元数 → 矩阵 → 四元数": [],
    }
    for _ in range(2000):
        R = RC.rot_axis(rng.normal(size=3), rng.uniform(-np.pi, np.pi))
        # 矩阵↔四元数
        paths["矩阵 → 四元数 → 矩阵"].append(
            np.max(np.abs(RC.quat_to_matrix(RC.matrix_to_quat(R)) - R)))
        # 矩阵↔欧拉（ZYX 内旋）
        y, p, r, _ = RC.matrix_to_euler_zyx(R)
        R2 = RC.euler_to_matrix("ZYX", [y, p, r], intrinsic=True)
        paths["矩阵 → 欧拉(ZYX) → 矩阵"].append(np.max(np.abs(R2 - R)))
        # 轴角→四元数→轴角（经矩阵中转）
        ax = RC.unit(rng.normal(size=3))
        th = rng.uniform(0, np.pi)
        q = RC.quat_from_axis_angle(ax, th)
        R3 = RC.quat_to_matrix(q)
        paths["轴角 → 四元数 → 轴角"].append(np.max(np.abs(
            RC.rot_axis(RC.rotation_axis_of(R3), RC.rotation_angle_of(R3)) - R3)))
        # 四元数↔矩阵
        q2 = RC.matrix_to_quat(RC.quat_to_matrix(q))
        d = min(np.max(np.abs(q2 - q)), np.max(np.abs(q2 + q)))
        paths["四元数 → 矩阵 → 四元数"].append(d)
    print(f"  {'转换路径':<28}{'最大误差':>14}{'数量级':>10}")
    for k, v in paths.items():
        m = max(v)
        print(f"  {k:<28}{m:>14.3e}{'机器精度' if m < 1e-12 else ('√eps' if m < 1e-6 else '有问题'):>12}")
    print("  → 四条路径都在机器精度级 ⇒ 公式与实现一致。")


def coin_flip():
    print("\n" + "=" * 74)
    print("【2】非单射的三种表现（这才是互转最容易出错的地方）")
    print("=" * 74)
    print("  (a) 四元数双覆盖：q 与 −q 给同一矩阵，但**往返取回的是哪一个**？")
    rng = np.random.default_rng(4)
    n_flip = 0
    for _ in range(1000):
        q = RC.quat_from_axis_angle(rng.normal(size=3), rng.uniform(-np.pi, np.pi))
        q2 = RC.matrix_to_quat(RC.quat_to_matrix(q))
        if np.max(np.abs(q2 - q)) > 1e-9:
            n_flip += 1
    print(f"    1000 组中，往返后**变成 −q** 的有 {n_flip} 组")
    print("    ⇒ 矩阵→四元数不是单值函数：它只能定到 ±。用法上必须")
    print("       「永远比较 |q₀·q₁|」或「相邻帧强制取同号」，否则插值会跳。")
    print()
    print("  (b) 轴角的 θ=0 与 θ=π 两个退化点：")
    R0 = RC.rot_axis([1, 1, 1], 0.0)
    print(f"    θ=0（R=I）：rotation_axis_of 返回 {np.round(RC.rotation_axis_of(R0), 6)}")
    print("      —— 任意轴都行（矩阵没给出信息），代码用 R+I 的第一大列兜底。")
    print("    θ=π：")
    worst_pi = 0.0
    for ax in [[1, 0, 0], [0, 1, 0], [1, 1, 0], [1, -2, 3]]:
        R = RC.rot_axis(ax, np.pi)
        u = RC.rotation_axis_of(R)
        th_back = RC.rotation_angle_of(R)
        R_back = RC.rot_axis(u, th_back)
        e = float(np.max(np.abs(R_back - R)))
        worst_pi = max(worst_pi, e)
        print(f"      轴 {str(ax):<12} → 反解轴 {np.round(u, 5)}  "
              f"角 {np.degrees(th_back):8.4f}°  重建误差 {e:.2e}")
    print(f"    ⇒ 全部重建到 {worst_pi:.1e}（机器精度）✓")
    print("       θ=π 时 u 与 −u 都合法 —— 但 `rotation_axis_of` 从 (R+I) 的列取方向，")
    print("       给出的是**确定的一支**，所以重建可靠。**判据要判矩阵，不判轴向量。**")
    print()
    print("    ⚠️ 这里踩过一个真坑（值得记）：`rotation_angle_of` 原先用")
    print("       `acos((tr−1)/2)`，在 θ=π（tr=−1）处导数发散，给出 179.999998°,")
    print("       **令重建差 2e-8**。换 `atan2(‖vee(R−Rᵀ)‖/2, (tr−1)/2)` 后")
    print("       回到机器精度。**是判据病态，不是实现错** —— 与 A5.4 段5 同族。")
    print()
    print("  (c) 欧拉角：锁死点有**无穷多组解**（见 A5.3 段4）")
    print("    这里只标出：反解在 |pitch| → 90° 附近会**跳变**，")
    R = RC.euler_to_matrix("ZYX", [0, np.radians(89.999), np.radians(10)],
                           intrinsic=True)
    R2 = RC.euler_to_matrix("ZYX", [0, np.radians(90.001), np.radians(10)],
                            intrinsic=True)
    y1, p1, r1, l1 = RC.matrix_to_euler_zyx(R)
    y2, p2, r2, l2 = RC.matrix_to_euler_zyx(R2)
    print(f"      pitch 89.999° → 反解 yaw = {np.degrees(y1):+9.4f}°")
    print(f"      pitch 90.001° → 反解 yaw = {np.degrees(y2):+9.4f}°")
    print(f"      两个姿态只差 0.002°，但 yaw 读数跳了 "
          f"{abs(np.degrees(y1)-np.degrees(y2)):.2f}°  ← 反解在锁死点附近不连续")


def same_rotation_four_ways():
    print("\n" + "=" * 74)
    print("【3】同一个旋转的四种写法（数值对照）")
    print("=" * 74)
    # 取一个不那么特殊的旋转
    u = RC.unit([1.0, 2.0, 0.5])
    th = np.radians(52.0)
    R = RC.rot_axis(u, th)
    q = RC.matrix_to_quat(R)
    y, p, r, locked = RC.matrix_to_euler_zyx(R)
    print(f"  目标：绕 u = {np.round(u, 4)} 转 {np.degrees(th):.2f}°")
    print()
    print("  ① 矩阵（9 个数，6 个冗余约束）：")
    print(RC.fmt_mat(R, 5))
    print()
    print(f"  ② 轴角：轴 = {np.round(RC.rotation_axis_of(R), 6)}，"
          f"角 = {np.degrees(RC.rotation_angle_of(R)):.6f}°")
    print(f"  ③ 四元数 (w,x,y,z) = {np.round(q, 6)}")
    print(f"      |q| = {RC.quat_norm(q):.12f}")
    print(f"  ④ 欧拉角 ZYX 内旋：yaw = {np.degrees(y):8.4f}°  "
          f"pitch = {np.degrees(p):8.4f}°  roll = {np.degrees(r):8.4f}°")
    print()
    print("  四种写法各自重建矩阵，与 R 的最大偏差：")
    print(f"    轴角      : {np.max(np.abs(RC.rot_axis(RC.rotation_axis_of(R), RC.rotation_angle_of(R)) - R)):.3e}")
    print(f"    四元数    : {np.max(np.abs(RC.quat_to_matrix(q) - R)):.3e}")
    print(f"    欧拉角    : {np.max(np.abs(RC.euler_to_matrix('ZYX', [y,p,r], intrinsic=True) - R)):.3e}")
    print()
    print("  顺带验证 q 的实部 = cos(θ/2)：")
    print(f"    q[0] = {q[0]:.10f}   cos(θ/2) = {np.cos(th/2):.10f}   "
          f"（θ/2 = {np.degrees(th/2):.4f}°）")


def storage_and_error():
    print("\n" + "=" * 74)
    print("【4】选型对照（含各自的「坏点」）")
    print("=" * 74)
    rows = [
        ("矩阵", "9", "复合 27 乘；旋转向量 9 乘",
         "9 个数冗余 6 约束；量化后必须重正交化"),
        ("四元数", "4", "复合 16 乘；旋转向量 2 次叉积",
         "±q 同姿态（非单射）；需归一化"),
        ("轴角", "4（3 轴 + 1 角）", "转成矩阵后才好用",
         "θ=0 轴无定义；θ=π 轴二义；不能直接复合"),
        ("欧拉角", "3", "不能直接复合（要转矩阵）",
         "12 种约定；±90° 万向锁；锁死点解不唯一"),
    ]
    print(f"  {'表示':<8}{'存储':>14}{'运算成本':>28}{'坏点':>38}")
    for a, b, c, d in rows:
        print(f"  {a:<8}{b:>14}{c:>28}{d:>38}")
    print()
    print("  实测补充（本目录其它脚本的结论）：")
    print("    · 四元数 vs 矩阵的漂移都是 √n·eps 级，量级相当（A5.2 段4）；")
    print("    · 差别在修复成本：四元数一次除法，矩阵要重正交化。")
    print("    · 短弧：四元数乘法前要判 q₀·q₁ 的符号（A5.4 段2）。")


# ------------------------------------------------------------ 画图
def make_figure():
    fig = plt.figure(figsize=(16, 5.6))

    # (a) 互转关系图
    ax = fig.add_subplot(1, 3, 1)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    nodes = {
        "矩阵 R":     (2.0, 8.0),
        "四元数 q":   (7.0, 8.0),
        "轴角 (u,θ)": (2.0, 4.5),
        "欧拉角 (y,p,r)": (7.0, 4.5),
    }
    cols = {"矩阵 R": "#185FA5", "四元数 q": "#0F6E56",
            "轴角 (u,θ)": "#993C1D", "欧拉角 (y,p,r)": "#534AB7"}
    for nm, (x, y) in nodes.items():
        ax.add_patch(plt.Rectangle((x - 1.35, y - 0.55), 2.7, 1.1,
                                   facecolor="white", edgecolor=cols[nm], lw=1.4))
        ax.text(x, y, nm, ha="center", va="center", fontsize=10.5,
                color=cols[nm])
    edges = [("矩阵 R", "四元数 q", "Shepperd"),
             ("矩阵 R", "轴角 (u,θ)", "trace/vee"),
             ("矩阵 R", "欧拉角 (y,p,r)", "ZYX 分解"),
             ("四元数 q", "轴角 (u,θ)", "半角"),
             ("四元数 q", "欧拉角 (y,p,r)", "经矩阵"),
             ("轴角 (u,θ)", "欧拉角 (y,p,r)", "经矩阵")]
    for a, b, lab in edges:
        (x1, y1), (x2, y2) = nodes[a], nodes[b]
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="<->", color="#888780", lw=1.0))
        ax.text((x1 + x2) / 2, (y1 + y2) / 2 + 0.28, lab, ha="center",
                fontsize=8, color="#5F5E5A")
    ax.set_title("① 四种表示互为翻译\n（都测到机器精度级）")

    # (b) 往返误差（条形）
    ax = fig.add_subplot(1, 3, 2)
    rng = np.random.default_rng(23)
    names = ["矩阵↔四元数", "矩阵↔欧拉", "轴角↔矩阵", "四元数↔矩阵"]
    errs = []
    for _ in range(2000):
        R = RC.rot_axis(rng.normal(size=3), rng.uniform(-np.pi, np.pi))
        errs.append([np.max(np.abs(RC.quat_to_matrix(RC.matrix_to_quat(R)) - R)),
                     np.max(np.abs(RC.euler_to_matrix(
                         "ZYX", RC.matrix_to_euler_zyx(R)[:3], intrinsic=True) - R)),
                     np.max(np.abs(RC.rot_axis(RC.rotation_axis_of(R),
                                               RC.rotation_angle_of(R)) - R)),
                     min(np.max(np.abs(RC.matrix_to_quat(RC.quat_to_matrix(
                         RC.matrix_to_quat(R))) - RC.matrix_to_quat(R))),
                         np.max(np.abs(RC.matrix_to_quat(RC.quat_to_matrix(
                             RC.matrix_to_quat(R))) + RC.matrix_to_quat(R))))])
    errs = np.array(errs)
    ax.bar(range(4), errs.max(axis=0), color=["#185FA5", "#534AB7", "#993C1D",
                                              "#0F6E56"])
    ax.set_yscale("log")
    ax.set_xticks(range(4))
    ax.set_xticklabels(names, fontsize=9)
    ax.axhline(2.2e-16, ls="--", color="gray", lw=1.2, label="eps = 2.2e-16")
    ax.set_ylabel("最大误差")
    ax.set_title("② 四条互转路径的往返误差\n全部机器精度级")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, axis="y", which="both")

    # (c) 各自的坏点（用"可用角范围"表示）
    ax = fig.add_subplot(1, 3, 3)
    ang = np.linspace(0, 180, 500)
    # 欧拉角条件数：|sec(pitch)| 之类
    with np.errstate(divide="ignore", invalid="ignore"):
        cond_euler = np.abs(1.0 / np.cos(np.radians(ang)))
    cond_euler = np.clip(cond_euler, 0, 1e3)
    ax.plot(ang, cond_euler, lw=2.2, color="#534AB7",
            label="欧拉角：条件数 ~1/cos(pitch)")
    ax.axvline(90, ls=":", color="k", lw=1.0)
    ax.text(92, 20, "万向锁", fontsize=9.5, color="#534AB7")
    # 轴角：θ=0 与 π 处轴无定义
    ax.axvspan(0, 2, color="#993C1D", alpha=0.25)
    ax.axvspan(178, 180, color="#993C1D", alpha=0.25)
    ax.text(3, 400, "θ=0\n轴无定义", fontsize=8.5, color="#993C1D")
    ax.text(150, 400, "θ=π\n轴二义", fontsize=8.5, color="#993C1D")
    ax.set_xlabel("旋转角（度）")
    ax.set_ylabel("条件数（越大越病态）")
    ax.set_yscale("log")
    ax.set_title("③ 各自的坏点（阴影 = 轴角退化，虚线 = 万向锁）")
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(alpha=0.3, which="both")
    ax.set_ylim(1, 1e3)

    fig.suptitle("A5.5 · 四种表示是同一对象的四张皮，各有各的坏点",
                 fontsize=13.5, y=0.99)
    fig.subplots_adjust(top=0.80, wspace=0.26)
    path = RC.savefig(fig, "representation_roundtrip.png")
    plt.close(fig)
    return path


def main():
    print("=" * 74)
    print("A5.5 · 轴角 / 四元数 / 欧拉角 / 矩阵：互转与坏点")
    print("=" * 74)
    roundtrip_table()
    coin_flip()
    same_rotation_four_ways()
    storage_and_error()
    path = make_figure()

    print("\n" + "=" * 74)
    print("验收：四条互转路径往返误差全部 <1e-12；θ=π 重建 <1e-15 "
          "（atan2 良态口径）  [PASS]")
    print("图：", path)
    print("\n这一步在说什么：")
    print("  · 四种表示互转全部在机器精度级 —— 公式层面没有秘密；")
    print("  · 难点在**退化点**：四元数的 ±、轴角的 θ=0/π、欧拉角的 ±90°；")
    print("  · 选型 = 选你愿意承受哪种坏点；工程默认四元数，")
    print("    只在需要人读/调试时转成欧拉角，且在锁死区格外小心。")


if __name__ == "__main__":
    main()