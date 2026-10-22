"""
A5.2 · 四元数：转 θ 为什么要用 θ/2
================================================================
现象：单位四元数 q = cos(θ/2) + sin(θ/2) u 表示"绕 u 转 θ"。
      为什么是**半角**？为什么 q 和 −q 表示同一个旋转（双覆盖）？

一句话在说什么：
    四元数是把"旋转"解释成一个东西的**平方**。
    绕 u 转 θ 可以看成"绕 u 转 θ/2 两次" —— 而四元数乘法恰好实现这个复合，
    于是 q = (半角那个四元数)，平方一次才得到矩阵 R。
    这就是半角的来历，也是双覆盖（q ↔ −q）的来历：
    取负号的那个 q 是"转 θ/2 + 180°"的版本，平方后 180° 那部分被吃掉了。

四元数乘法的非交换性：
    (w₁ + v₁)(w₂ + v₂) = (w₁w₂ − v₁·v₂) + (w₁v₂ + w₂v₁ + v₁×v₂)
    —— 最后那一项 v₁×v₂ 就是非交换性的唯一来源。

运行： uv run python graphics/demo/pyapps/rotation/quaternion.py
输出： learning-todo/imgs/rot/quaternion_geometry.png
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _rot_core as RC  # noqa: E402

plt = RC.setup_style()


def half_angle_story():
    print("=" * 72)
    print("【1】半角：为什么是 θ/2 而不是 θ")
    print("=" * 72)
    print("  四元数的旋转公式（三明治积）：v' = q v q*，其中 q = cos(θ/2) + sin(θ/2)u")
    print("  如果按「转 θ 就是 θ」来猜，会写成 q = cosθ + sinθ·u，那就错了。")
    print()
    print("  实测：把两种候选都算出来，看谁给出正确的旋转：")
    for th in [np.pi / 3, np.pi / 2, 2.0]:
        u = RC.unit([0, 0, 1])
        v = np.array([1.0, 0, 0])
        # 正确的：半角
        q_ok = RC.quat_from_axis_angle(u, th)
        v_ok = RC.quat_rotate(q_ok, v)
        # 错误的：全角
        q_bad = np.array([np.cos(th), *(np.sin(th) * u)])
        v_bad = RC.quat_rotate(q_bad, v)
        # 参照：矩阵
        v_ref = RC.rot_axis(u, th) @ v
        print(f"    θ={th:6.3f}: 半角版 → {np.round(v_ok, 6)}  "
              f"（参照 {np.round(v_ref, 6)}，差 {np.max(np.abs(v_ok-v_ref)):.1e}）")
        print(f"             全角版 → {np.round(v_bad, 6)}  "
              f"（差 {np.max(np.abs(v_bad-v_ref)):.1e}）  ← 错")
    print()
    print("  为什么是半角？—— 四元数乘法是「复合两次」：")
    q_half = RC.quat_from_axis_angle([0, 0, 1], np.pi / 2)
    q_double = RC.quat_mul(q_half, q_half)
    q_full = RC.quat_from_axis_angle([0, 0, 1], np.pi)
    print(f"    绕 z 转 90° 的 q        = {np.round(q_half, 6)}")
    print(f"    自己乘自己 (q·q)        = {np.round(q_double, 6)}")
    print(f"    绕 z 转 180° 的 q       = {np.round(q_full, 6)}")
    print(f"    → |q·q − q(180°)| = {np.max(np.abs(q_double - q_full)):.2e}")
    print("    ⇒ **四元数相乘 = 旋转复合**（角相加）。所以半角两次才凑成一个 θ。")


def double_cover():
    print("\n" + "=" * 72)
    print("【2】双覆盖：q 与 −q 是同一个旋转（SO(3) ≅ SU(2)/{±1}）")
    print("=" * 72)
    rng = np.random.default_rng(5)
    worst = 0.0
    for _ in range(1000):
        q = RC.quat_from_axis_angle(rng.normal(size=3), rng.uniform(-np.pi, np.pi))
        worst = max(worst, np.max(np.abs(RC.quat_to_matrix(q) - RC.quat_to_matrix(-q))))
    print(f"  1000 组随机：max|M(q) − M(−q)| = {worst:.2e}  "
          f"[{'PASS' if worst < 1e-13 else 'FAIL'}]")
    print()
    print("  数值解释：M(q) 里每一项都是 q 的**二次型**（如 2(xy − wz)），")
    print("  整体取负后 (−q) 的二次型不变 ⇒ 矩阵不变。")
    print()
    print("  几何解释：q 转 θ ⟺ −q 也转 θ，但 −q 走了「另一条路」——把 q 看成")
    print("  cos(θ/2)+sin(θ/2)u，当 θ/2 走过 180° 时 q 的符号翻转。")
    print("  所以**绕任意轴转 360°**（不是 720°）时 q 会变成 −q：")
    for th in [np.pi, 2 * np.pi, 4 * np.pi]:
        q = RC.quat_from_axis_angle([0, 0, 1], th)
        print(f"    θ={np.degrees(th):5.0f}°: q = {np.round(q, 6)}"
              f"   （q 是 +1 还是 −1 的倍？{'原样' if q[0] > -0.99 else '变号'}）")
    print("  → θ=360° 时 q = −1·q₀（变成了负代表元），θ=720° 才回到 q₀ 本身。")
    print("    这就是「旋量」（spinor）需要转两圈才复原的数学原型。")


def noncommutativity():
    print("\n" + "=" * 72)
    print("【3】非交换：转动的次序不能换（这也是四元数不是复数的原因）")
    print("=" * 72)
    Rx90 = RC.rot_xyz("X", np.pi / 2)
    Rz90 = RC.rot_xyz("Z", np.pi / 2)
    qx = RC.quat_from_axis_angle([1, 0, 0], np.pi / 2)
    qz = RC.quat_from_axis_angle([0, 0, 1], np.pi / 2)
    print("  矩阵语言：Rx·Rz 与 Rz·Rx 的差 = "
          f"{np.max(np.abs(Rx90 @ Rz90 - Rz90 @ Rx90)):.4f}")
    print("  四元数语言：qx·qz 与 qz·qx 的差 = "
          f"{np.max(np.abs(RC.quat_mul(qx, qz) - RC.quat_mul(qz, qx))):.4f}")
    print("  一致性：M(qx qz) 是否 = Rx Rz ？  "
          f"{np.allclose(RC.quat_to_matrix(RC.quat_mul(qx, qz)), Rx90 @ Rz90)}")
    print()
    print("  直观：绕 x 转 90° 再绕 z 转 90°，与反过来，结果不同 ——")
    v = np.array([0.0, 0.0, 1.0])
    a = Rx90 @ (Rz90 @ v)
    b = Rz90 @ (Rx90 @ v)
    print(f"    (先 z 后 x) z 轴 → {np.round(a, 6)}")
    print(f"    (先 x 后 z) z 轴 → {np.round(b, 6)}")
    print("  复数乘法可交换（ab = ba），所以**复数表示不了三维旋转**。")
    print("  四元数非交换，正好与「旋转不交换」匹配 —— 这是它好用的根本原因。")


def norm_and_drift():
    print("\n" + "=" * 72)
    print("【4】四元数 vs 矩阵：漂移与修复成本（以实测为准，不背结论）")
    print("=" * 72)
    print(f"  {'复合次数':>10}{'四元数 |q|−1':>18}{'矩阵 |RᵀR−I|':>18}"
          f"{'√n·eps 参考':>16}")
    prev_q = prev_R = None
    for n in [1_000, 10_000, 100_000, 1_000_000]:
        rng = np.random.default_rng(9)
        q = np.array([1.0, 0, 0, 0])
        R = np.eye(3)
        for _ in range(n):
            axv = RC.unit(rng.normal(size=3))
            d = rng.uniform(-0.01, 0.01)
            q = RC.quat_mul(q, RC.quat_from_axis_angle(axv, d))
            R = R @ RC.rot_axis(axv, d)
        dq = abs(RC.quat_norm(q) - 1.0)
        dR = np.max(np.abs(R.T @ R - np.eye(3)))
        print(f"  {n:>10d}{dq:>18.3e}{dR:>18.3e}{np.sqrt(n)*2.22e-16:>16.2e}")
    print()
    print("  **实测结论（与常见说法不同，值得记住）**：")
    print("    · 两者的漂移都是 **√n·eps 级**（随机游走），量级**相当**，")
    print("      不存在「矩阵漂得快、四元数不漂」这回事。")
    print("    · 真正的差别在**修复成本与副作用**：")
    print("        四元数：q ← q/|q|。3 次乘法、O(1)，而且这个**新四元数仍是")
    print("                单位四元数**（数域内闭合），反复修不会累积几何误差。")
    print("        矩阵  ：必须 Gram–Schmidt 或 SVD 重正交化，代价 O(n³)→(这里是 O(27))")
    print("                且**改变了原矩阵**（最近的正交阵不是它本身），是投影而非修正。")
    print("    · 所以四元数的优势是「4 个数 + 廉价且无害的归一化」，")
    print("      不是「不漂移」。这点很多教程讲错了。")
    print()
    print("  运算成本（每次复合）：")
    print("    四元数：4 个数；16 乘 + 12 加")
    print("    矩阵  ：9 个数；27 乘 + 18 加")
    print("    旋转一个向量 v：四元数 2 次叉积（无需展开成矩阵）。")


# ------------------------------------------------------------ 画图
def make_figure():
    fig = plt.figure(figsize=(15.5, 5.6))

    # (a) 半角：q 的"角度"是旋转角的一半
    ax = fig.add_subplot(1, 3, 1)
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.degrees(th), np.cos(th / 2), lw=2.2, color="#185FA5",
            label=r"$w=\cos(\theta/2)$（q 的实部）")
    ax.plot(np.degrees(th), np.sin(th / 2), lw=2.2, color="#c0392b",
            label=r"$\|v\|=\sin(\theta/2)$（虚部模长）")
    ax.axhline(0, color="k", lw=0.6)
    for k, lab in [(1, "360°"), (2, "720°")]:
        ax.axvline(360 * k, ls=":", color="k", lw=0.7)
    ax.set_xlabel("实际旋转角 θ（度）")
    ax.set_title("① 半角：转 360° 时 q 才走到 −1\n转 720° 回到 +1")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(alpha=0.3)
    ax.set_xlim(0, 720)
    ax.annotate("q = −1\n（负代表元）", xy=(360, -1), xytext=(430, -0.55),
                fontsize=9, arrowprops=dict(arrowstyle="->", lw=0.9))
    ax.annotate("q = +1\n（回到自身）", xy=(720, 1), xytext=(600, 0.45),
                fontsize=9, arrowprops=dict(arrowstyle="->", lw=0.9))

    # (b) 双覆盖：q 空间是"球面 S³"，两条半圈对应同一个旋转
    ax = fig.add_subplot(1, 3, 2, projection="3d")
    ph = np.linspace(0, 2 * np.pi, 60)
    arc = np.linspace(0, np.pi, 60)
    for p0, col, alp, lab in [(0.0, "#185FA5", 0.9, "q 的上半圆"),
                              (np.pi, "#c0392b", 0.9, "−q 的下半圆")]:
        x = np.sin(arc) * np.cos(p0)
        y = np.sin(arc) * np.sin(p0)
        # 用 3D 画一个"大圆"的两个半条，示意 q 与 −q 两端
        ax.plot(np.cos(arc), np.sin(arc) * np.cos(p0), np.sin(arc) * np.sin(p0),
                color=col, lw=2.2, alpha=alp, label=lab if p0 == 0 else None)
    ax.scatter([1, -1], [0, 0], [0, 0], color="k", s=45)
    ax.text(1.1, 0, 0, "q", fontsize=11)
    ax.text(-1.3, 0, 0, "−q", fontsize=11)
    ax.set_box_aspect((1, 1, 1))
    ax.set_title("② 单位四元数住在 S³ 上\nq 与 −q 是对径点，却是同一个旋转")
    ax.set_axis_off()

    # (c) 归一化漂移对比
    ax = fig.add_subplot(1, 3, 3)
    rng = np.random.default_rng(9)
    ns = [10, 30, 100, 300, 1000, 3000, 10000]
    qs_drift, R_drift = [], []
    for n in ns:
        rng2 = np.random.default_rng(9)
        q = np.array([1.0, 0, 0, 0])
        R = np.eye(3)
        for _ in range(n):
            axv = RC.unit(rng2.normal(size=3))
            d = rng2.uniform(-0.01, 0.01)
            q = RC.quat_mul(q, RC.quat_from_axis_angle(axv, d))
            R = R @ RC.rot_axis(axv, d)
        qs_drift.append(abs(RC.quat_norm(q) - 1.0))
        R_drift.append(np.max(np.abs(R.T @ R - np.eye(3))))
    ax.loglog(ns, qs_drift, "o-", ms=4, lw=1.6, color="#185FA5",
              label="四元数 |q|−1（不归一）")
    ax.loglog(ns, R_drift, "s-", ms=4, lw=1.6, color="#c0392b",
              label=r"矩阵 $|R^TR-I|$（不重正交化）")
    ax.loglog(ns, np.sqrt(ns) * 1e-17, "--", color="gray", lw=1.2,
              label=r"参考 $\sqrt{n}$（随机游走）")
    ax.set_xlabel("复合次数 n")
    ax.set_ylabel("误差")
    ax.set_title("③ 漂移都是 √n 级（量级相当）\n差别在修复成本：除法 vs 重正交化")
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.3, which="both")

    fig.suptitle("A5.2 · 四元数：半角来自「平方」、双覆盖来自「对径点」",
                 fontsize=13.5, y=0.99)
    fig.subplots_adjust(top=0.80, wspace=0.28)
    path = RC.savefig(fig, "quaternion_geometry.png")
    plt.close(fig)
    return path


def main():
    print("=" * 72)
    print("A5.2 · 四元数：转 θ 为什么用 θ/2")
    print("=" * 72)
    half_angle_story()
    double_cover()
    noncommutativity()
    norm_and_drift()
    path = make_figure()

    print("\n图：", path)
    print("\n这一步在说什么：")
    print("  · 四元数乘法 = 旋转复合 ⇒ 半角，取两次才凑成一个 θ；")
    print("  · q 与 −q 给出同一个旋转 ⇒ SO(3) ≅ SU(2)/{±1}（双覆盖）；")
    print("  · 非交换（v₁×v₂ 项）不是缺陷，正好匹配「旋转不交换」这一事实；")
    print("  · 4 个数 + 一次除法就能抗漂移 ⇒ 工程首选。")


if __name__ == "__main__":
    main()