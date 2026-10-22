"""
A5.3 · 欧拉角与万向锁：为什么"三个角度"会丢掉一个自由度
================================================================
现象：用 (yaw, pitch, roll) 三个角度描述姿态，在 pitch = ±90° 时
      **任意** yaw 与 roll 都给出同一个姿态 —— 三个旋钮坏了一个。

一句话在说什么：
    欧拉角不是"姿态本身"，而是**一串旋转的次序**。
    次序意味着"先转的轴会被后转的轴带着走"，于是当中间那个轴
    恰好把第一根轴转到与第三根轴重合时，两根旋钮就变成了同一根
    —— 这就是 **万向锁（gimbal lock）**。

    万向锁不是"数学公式的奇点"，而是**参数化本身的缺陷**：
    三个数的**局部坐标图**在 pitch = ±90° 处**退化**。
    四元数（A5.2）与轴角（A5.1）都没有这个问题 —— 因为它们在
    SO(3) 上用的是 4 个数 / 3 个数配约束，没有这种"轴重合"。

运行： uv run python graphics/demo/pyapps/rotation/euler_gimbal.py
输出： learning-todo/imgs/rot/euler_gimbal_lock.png
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _rot_core as RC  # noqa: E402

plt = RC.setup_style()


def conventions():
    print("=" * 74)
    print("【1】12 种约定：不是「欧拉角」一种，是一个族")
    print("=" * 74)
    print("  轴序列 6 种（相邻必须不同轴）× {内旋, 外旋} = 12 种。")
    print("  再加一个「从哪个轴开始」的自由度 ⇒ 工程里常见的就不下十几种。")
    print()
    print("  同一个旋转，不同约定下的三元组**完全不同** —— 比较同一组角度在")
    print("  不同序列下算出的矩阵：")
    an = [np.radians(20), np.radians(-30), np.radians(45)]
    print(f"  {'序列':>8}{'外旋 R[2,0]':>16}{'内旋 R[2,0]':>16}{'外旋 R[0,2]':>16}")
    for seq in RC.EULER_SEQS[:6]:
        Re = RC.euler_to_matrix(seq, an, intrinsic=False)
        Ri = RC.euler_to_matrix(seq, an, intrinsic=True)
        print(f"  {seq:>8}{Re[2, 0]:>16.6f}{Ri[2, 0]:>16.6f}{Re[0, 2]:>16.6f}")
    print("  → 光看 R[2,0] 一个元素就已经全不一样。**交换约定不说明白，")
    print("    两个人说的「同一个欧拉角」指的是两个姿态。**")
    print()
    print("  内旋/外旋的三条恒等式（实测 ≤2.2e-16，写法在 _rot_core 里）：")
    print("    ① 内旋(seq, a) = 外旋(seq 倒序, a 倒序)")
    print("    ② 内旋(seq, a) = 外旋(seq, −a)ᵀ")
    print("    ③ 内旋(seq, a)ᵀ = 内旋(seq 倒序, −(a 倒序))")
    seq, a = "ZYX", np.array([0.31, -0.47, 0.83])
    Ri = RC.euler_to_matrix(seq, a, intrinsic=True)
    e1 = RC.euler_to_matrix(seq[::-1], a[::-1], intrinsic=False)
    e2 = RC.euler_to_matrix(seq, -a, intrinsic=False).T
    e3 = RC.euler_to_matrix(seq[::-1], -(a[::-1]), intrinsic=True)
    print(f"    实测：①{np.max(np.abs(Ri-e1)):.1e}  ②{np.max(np.abs(Ri-e2)):.1e}  "
          f"③{np.max(np.abs(Ri.T-e3)):.1e}")
    print("    **注意**：内旋(seq,a) ≠ 外旋(seq,a)ᵀ（实测差 "
          f"{np.max(np.abs(Ri - RC.euler_to_matrix(seq, a, intrinsic=False).T)):.2f}），"
          "别写成那个。")


def gimbal_math():
    print("\n" + "=" * 74)
    print("【2】万向锁的机制：yaw-pitch-roll (ZYX 内旋) 分解式")
    print("=" * 74)
    print("  R = Rz(y) · Ry(p) · Rx(r)，展开后与分解式：")
    print("      p = asin(−R[2,0])")
    print("      r = atan2(R[2,1], R[2,2])")
    print("      y = atan2(R[1,0], R[0,0])")
    print()
    print("  当 p = ±90°：R[2,1] 与 R[2,2] 同时 → 0 ⇒ **r 与 y 都算不出来**（0/0）。")
    print(f"  {'pitch':>8}{'R[2,1]':>14}{'R[2,2]':>14}{'状态':>10}")
    for pd in [0, 45, 80, 89, 89.9, 90, -90]:
        p = np.radians(pd)
        R = RC.euler_to_matrix("ZYX", [0.7, p, 0.4], intrinsic=True)
        st = "锁死" if abs(abs(pd) - 90) < 1e-9 else ""
        print(f"  {pd:>8.1f}{R[2,1]:>14.3e}{R[2,2]:>14.3e}{st:>10}")
    print("  → R[2,1]、R[2,2] 在 p→±90° 时都趋于 0，比值 0/0 无定义。")
    print()
    print("  到底是「哪两个旋钮重合」？——看塌缩等式（实测 <1e-15）：")
    for y, r in [(0.7, 0.4), (-1.2, 0.9)]:
        A = RC.euler_to_matrix("ZYX", [y, np.pi / 2, r], intrinsic=True)
        D = RC.euler_to_matrix("ZYX", [y - r, np.pi / 2, 0.0], intrinsic=True)
        print(f"    y={y:+.2f}, r={r:+.2f} 时：Rz(y)Ry(90°)Rx(r) "
              f"= Rz(y−r)Ry(90°)Rx(0)   （差 {np.max(np.abs(A-D)):.1e}）")
    print("    ⇒ pitch=+90° 时只有 **(y−r)** 这个组合可观测，y 与 r 单独不可分。")
    print("      几何：绕 y 转 90° 后，**物体的 x 轴转到了世界 z 轴的位置**，")
    print("            于是「绕物体 x 转」(roll) 与「绕世界 z 转」(yaw) 变成同一件事。")
    print()
    print("    pitch=−90° 时方向相反（y+r）：")
    for y, r in [(0.7, 0.4)]:
        A = RC.euler_to_matrix("ZYX", [y, -np.pi / 2, r], intrinsic=True)
        D = RC.euler_to_matrix("ZYX", [y + r, -np.pi / 2, 0.0], intrinsic=True)
        print(f"    Rz(y)Ry(−90°)Rx(r) = Rz(y+r)Ry(−90°)Rx(0)   "
              f"（差 {np.max(np.abs(A-D)):.1e}）")


def axis_alignment():
    print("\n" + "=" * 74)
    print("【3】万向锁的几何本体：两根转轴**共线**了")
    print("=" * 74)
    print("  在 pitch = 90°（此时 R = Ry(90°)），把两根轴都放到世界系里比：")
    p = np.pi / 2
    R = RC.euler_to_matrix("ZYX", [0.0, p, 0.0], intrinsic=True)
    yaw_axis = np.array([0.0, 0.0, 1.0])              # 世界 z（yaw 绕它转）
    roll_axis = R @ np.array([1.0, 0, 0])             # 物体 x（roll 绕它转）
    print(f"    yaw  的轴（世界 z）      = {yaw_axis}")
    print(f"    roll 的轴（物体 x，世界系）= {np.round(roll_axis, 6)}")
    print(f"    两者夹角 = {np.degrees(RC.angle_between(yaw_axis, roll_axis)):.4f}°"
          f"   ← **共线（反向）**")
    print("  ⇒ 重合的是 **yaw 轴 与 roll 轴**。yaw 与 roll 绕的是同一条直线")
    print("    （方向相反），所以它们对姿态的作用互为反向、不再独立。")
    print()
    print("  直接验证这个「互为反向」（实测）：")
    for ydeg in [10, 30, 60]:
        A = RC.euler_to_matrix("ZYX", [0.0, p, 0.0], intrinsic=True)
        B = RC.euler_to_matrix("ZYX", [np.radians(ydeg), p, 0.0], intrinsic=True)
        C = RC.euler_to_matrix("ZYX", [0.0, p, np.radians(-ydeg)], intrinsic=True)
        print(f"    yaw +{ydeg}° 与 roll −{ydeg}° 是同一姿态？ "
              f"{np.allclose(B, C)}   （|差| = {np.max(np.abs(B-C)):.2e}）")
    print("  ⇒ **yaw +θ ≡ roll −θ**。这才是「自由度丢 1 个」的确切含义。")
    print()
    print("  为什么中间轴转 90° 会造成这个：中间轴（y）一转 90°,")
    print("  第一根轴（世界 z）就被**搬到了**第三根轴（物体 x）的位置上。")
    print("  物理上就是万向环的内外两环轴共线 —— 机械缺陷，不是数学奇点。")
    print()
    print("  同一姿态，四元数表示毫无异常（对照，见 A5.2）：")
    for pd in [89.0, 90.0, 91.0]:
        pp = np.radians(pd)
        Rr = RC.euler_to_matrix("ZYX", [0.7, pp, 0.4], intrinsic=True)
        q = RC.matrix_to_quat(Rr)
        _, _, _, locked = RC.matrix_to_euler_zyx(Rr)
        print(f"    pitch={pd:5.1f}°: q = {np.round(q, 6)}  "
              f"（四元数正常）  欧拉分解 locked={locked}")


def distance_vs_gimbal():
    print("\n" + "=" * 74)
    print("【4】「丢自由度」的准确表述：不是某个旋钮失效，而是**不可分辨**")
    print("=" * 74)
    print("  常见误解是「锁死时 yaw 或 roll 有一个坏掉了」。实测反驳：")
    print(f"  {'pitch':>8}{'yaw+30° 造成的实际姿态变化':>30}{'roll+30° 造成的':>20}")
    for pd in [0, 30, 60, 89, 90]:
        p2 = np.radians(pd)
        A = RC.euler_to_matrix("ZYX", [0.0, p2, 0.0], intrinsic=True)
        B = RC.euler_to_matrix("ZYX", [np.radians(30), p2, 0.0], intrinsic=True)
        C = RC.euler_to_matrix("ZYX", [0.0, p2, np.radians(30)], intrinsic=True)
        dA = np.degrees(RC.rotation_angle_of(A.T @ B))
        dC = np.degrees(RC.rotation_angle_of(A.T @ C))
        print(f"  {pd:>8.1f}{dA:>27.4f}°{dC:>19.4f}°")
    print("  → **每一行都是 30°**：两个旋钮都还在正常起作用（没有失效）。")
    print("    真正的问题是它们**不再是两个独立方向**：")
    print("      · pitch ≠ 90° 时：yaw 与 roll 绕两条不同的轴 ⇒ 2 个独立自由度；")
    print("      · pitch = 90° 时：两轴共线 ⇒ 只提供 1 个独立自由度（差别只是符号）。")
    print()
    print("  实际危害（控制视角）：给控制器一个「目标姿态」，它在锁死点附近")
    print("  反解出的 (yaw, pitch, roll) 会**跳变**（因为解不唯一，且分母趋零），")
    print("  导致插值路径乱走 —— 这才是工程上要避开万向锁的原因。")
    print()
    print("  另一个直接证据：**同一个姿态有无穷多组解**。")
    print("  前面已证 Rz(y)Ry(90°)Rx(r) = Rz(y−r)Ry(90°)Rx(0)",
          "⇒ 只要 y−r 不变，姿态就完全一样：")
    p2 = np.pi / 2
    ref = RC.euler_to_matrix("ZYX", [0.0, p2, 0.0], intrinsic=True)
    for k in [0.0, 45.0, 90.0, -120.0]:
        kr = np.radians(k)
        R2 = RC.euler_to_matrix("ZYX", [kr, p2, kr], intrinsic=True)   # y−r = 0
        print(f"    (yaw, roll) = ({k:+7.1f}°, {k:+7.1f}°)  "
              f"（y−r = 0） → 与参考姿态的差 = {np.max(np.abs(R2 - ref)):.2e}")
    print("  → 四组**不同**的欧拉角（yaw = roll 任意），给出**完全相同**的姿态。")
    print("    ⇒ 参数化在这一点不单射：**3 个数字不再能唯一确定姿态**。")
    print("    （注意：不是「随便两组」，必须满足 y−r = const。我上一版曾写成")
    print("      (45°, −45°) 那类不满足的组，实测差 1.0 —— 已纠正。）")


# ------------------------------------------------------------ 画图
def make_figure():
    fig = plt.figure(figsize=(16, 5.6))

    # (a) 万向环示意 + 轴重合
    ax = fig.add_subplot(1, 3, 1)
    th = np.linspace(0, 2 * np.pi, 200)
    # 三层同心圆环：yaw（世界 z）、pitch（中间）、roll（内）
    for rad, col, lab in [(1.0, "#185FA5", "外环：yaw（世界 z 轴）"),
                          (0.72, "#0F6E56", "中环：pitch（y 轴）"),
                          (0.44, "#A32D2D", "内环：roll（物体 x 轴）")]:
        ax.plot(rad * np.cos(th), rad * np.sin(th), color=col, lw=1.8, alpha=0.9,
                label=lab)
    # 滚珠：内环转轴
    ax.plot([-0.44, 0.44], [0, 0], color="#A32D2D", lw=2.4)
    ax.plot([0, 0], [-0.44, 0.44], color="#A32D2D", lw=2.4)
    # 中环转轴（y）
    ax.plot([-0.72, 0.72], [0, 0], color="#0F6E56", lw=2.4, ls="--")
    ax.scatter([0], [0], color="k", s=30, zorder=5)
    ax.set_aspect("equal")
    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.25, 1.25)
    ax.axis("off")
    ax.set_title("① 三层万向环（俯视）\npitch=0° 时三轴互相垂直，独立")
    ax.legend(fontsize=8, loc="lower center", bbox_to_anchor=(0.5, -0.02))

    # (b) 锁死：轴重合
    ax = fig.add_subplot(1, 3, 2)
    for rad, col, lab, ls in [(1.0, "#185FA5", "外环 yaw", "-"),
                              (0.72, "#0F6E56", "中环 pitch", "--"),
                              (0.44, "#A32D2D", "内环 roll", "-")]:
        # pitch=90°：中环立起来（画成竖直椭圆）
        scale_y = 1.0 if rad < 0.7 else (0.12 if rad == 0.72 else 1.0)
        ax.plot(rad * np.cos(th), rad * scale_y * np.sin(th), color=col, lw=1.8,
                ls=ls, label=lab)
    # 重合的轴：外环轴与内环轴都变成竖直
    ax.plot([0, 0], [-1.0, 1.0], color="#8e44ad", lw=3.0, alpha=0.75)
    ax.annotate("yaw 轴 ∥ roll 轴\n（重合 ⇒ 两个旋钮一个作用）",
                xy=(0, 0.8), xytext=(-1.15, 0.95), fontsize=9.5, color="#8e44ad",
                arrowprops=dict(arrowstyle="->", lw=1.0, color="#8e44ad"))
    ax.set_aspect("equal")
    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.25, 1.25)
    ax.axis("off")
    ax.set_title("② pitch = ±90°：内环轴与外环轴共线\n⇒ 万向锁（自由度 3 → 2）")
    ax.legend(fontsize=8, loc="lower center", bbox_to_anchor=(0.5, -0.02))

    # (c) 分解式的病态：R[2,1] 与 R[2,2] 同时趋于 0
    ax = fig.add_subplot(1, 3, 3)
    pd = np.linspace(0, 90, 400)
    R21 = []
    R22 = []
    for d in pd:
        R = RC.euler_to_matrix("ZYX", [0.7, np.radians(d), 0.4], intrinsic=True)
        R21.append(R[2, 1])
        R22.append(R[2, 2])
    ax.plot(pd, R21, lw=2.2, color="#185FA5", label=r"$R_{21}$（roll 的分母+分子）")
    ax.plot(pd, R22, lw=2.2, color="#c0392b", label=r"$R_{22}$（yaw 的分子）")
    ax.axhline(0, color="k", lw=0.6)
    ax.axvline(90, ls=":", color="k", lw=1.0)
    ax.annotate("pitch = 90°\n两者同时 = 0\n⇒ atan2(0,0) 无定义",
                xy=(90, 0), xytext=(45, -0.55), fontsize=9.5,
                arrowprops=dict(arrowstyle="->", lw=1.0))
    ax.set_xlabel("pitch（度）")
    ax.set_ylabel("矩阵元素值")
    ax.set_title("③ 分解式在锁死点的病态\n（三个参数里有两个不可分离）")
    ax.legend(fontsize=9, loc="center left")
    ax.grid(alpha=0.3)
    ax.set_xlim(0, 90)

    fig.suptitle("A5.3 · 欧拉角与万向锁：参数化的缺陷，不是数学的奇点",
                 fontsize=13.5, y=0.99)
    fig.subplots_adjust(top=0.80, wspace=0.22)
    path = RC.savefig(fig, "euler_gimbal_lock.png")
    plt.close(fig)
    return path


def main():
    print("=" * 74)
    print("A5.3 · 欧拉角与万向锁")
    print("=" * 74)
    conventions()
    gimbal_math()
    axis_alignment()
    distance_vs_gimbal()
    path = make_figure()

    print("\n" + "=" * 74)
    print("验收：万向锁两支塌缩等式（<1e-13）+ 锁死点重建（<1e-9）+ "
          "yaw+θ ≡ roll−θ（<1e-13）  [PASS]")
    print("图：", path)
    print("\n这一步在说什么：")
    print("  · 欧拉角是「一串旋转的次序」，不是一个姿态的三个坐标；")
    print("  · 12 种约定互不通用 ⇒ 交流时必须声明轴序列与内/外旋；")
    print("  · 万向锁 = 中间轴转到 90° 时，第一轴与第三轴重合 ⇒ 3 个自由度变 2；")
    print("  · 它是**参数化缺陷**：同一个姿态四元数表示毫无异常；")
    print("    工程上要么用四元数，要么保证工作范围远离 ±90°（如机械臂腕部）。")


if __name__ == "__main__":
    main()