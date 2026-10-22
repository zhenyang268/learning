"""
A5.1 · 罗德里格斯公式：从 "绕轴转" 到 "一个矩阵"
================================================================
现象：exp(θ[u]×) 这个无穷级数，居然能收敛成三行闭式；
      而且它自动给出一个旋转矩阵（正交、det = +1）。

一句话在说什么：
    旋转是**李群** SO(3) 的元素，而"无穷小旋转"（角速度）是**李代数** so(3) 的元素
    —— 反对称 3×3 矩阵。两者之间隔着一个**指数映射** exp。
    罗德里格斯公式就是这个指数映射的闭式：

        exp(θ[u]×) = I + sinθ [u]× + (1 − cosθ) [u]×²

    它成立全靠一条恒等式：**[u]×³ = −[u]×**（u 为单位向量）。
    有了它，级数的奇偶次项各自可以求和，得到 sinθ 与 (1−cosθ)。

为什么 "×" 写成矩阵（hat 映射）：
    a × b = [a]× b  —— 叉积是线性的，所以可以写成矩阵乘法。
    [u]× 是**反对称**的（[u]×ᵀ = −[u]×），这正是 so(3) 的定义。

运行： uv run python graphics/demo/pyapps/rotation/rodrigues.py
输出： learning-todo/imgs/rot/rodrigues_expmap.png
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _rot_core as RC  # noqa: E402

plt = RC.setup_style()


def sphere_cheatsheet():
    """1. 三条恒等式 —— 闭式的全部原料。"""
    print("=" * 72)
    print("【1】so(3) 的三条恒等式（闭式的全部原料）")
    print("=" * 72)
    u = RC.unit([1.0, 2.0, -0.5])
    K = RC.hat(u)
    print(f"  取单位向量 u = {np.round(u, 6)}")
    print(f"  [u]× 反对称：max|[u]× + [u]×ᵀ| = "
          f"{np.max(np.abs(K + K.T)):.2e}  （定义为反对称）")
    print(f"  [u]×u = 0：          max = {np.max(np.abs(K @ u)):.2e}  "
          f"（绕自己转不产生位移）")
    print(f"  [u]×³ = −[u]×：      max = {np.max(np.abs(K @ K @ K + K)):.2e}  "
          f"← 关键！")
    print("  → 第三条把 3 次幂打回 1 次幂，于是级数只有两类项：")
    print("     偶次（I 的倍数）与奇次（[u]× 的倍数），各自可求和。")

    print("\n" + "=" * 72)
    print("【2】级数 → 闭式，逐项对照（这就是公式的来历）")
    print("=" * 72)
    print("  exp(θ[u]×) = Σ (θ[u]×)ⁿ/n!")
    print("  用 [u]×³ = −[u]× 把幂次降回 1 或 2 ⇒ 只剩三类：I、[u]×、[u]×²")
    print()
    print("    偶次幂 ⊂ {I, [u]×²}，奇次幂 ⊂ {[u]×, [u]×²}（后者乘 [u]× 仍是 2 次）")
    print("    所以闭式**必然**长这样：  A(θ)·I + B(θ)·[u]× + C(θ)·[u]×²")
    print()
    print("  三个系数怎么定？取三组 θ 实测反解 —— 拟合 A·I + B·[u]× + C·[u]×² 的系数：")
    for th in [0.7, 1.9, 3.4]:
        R = RC.rot_axis(u, th)
        K2 = K @ K
        # A I + B K + C K2 = R  ⇒ 三个未知，9 个方程最小二乘
        Amat = np.stack([np.eye(3).ravel(), K.ravel(), K2.ravel()], axis=1)
        coeff, *_ = np.linalg.lstsq(Amat, R.ravel(), rcond=None)
        A_, B_, C_ = coeff
        print(f"    θ={th:4.2f}:  A={A_:+.10f} (恒为 1)   "
              f"B={B_:+.10f} (sinθ={np.sin(th):+.10f})   "
              f"C={C_:+.10f} (1−cosθ={1-np.cos(th):+.10f})")
    print("  → 实测系数是 A=1、B=sinθ、C=1−cosθ ⇒ 正好是「形式甲」（不含 cosθ·I）。")
    print("    注意 A 恒为 1 而不是 cosθ：因为 I 的系数与 [u]×² 的系数**捆绑**，")
    print("    拆成 cosθ·I 时 [u]×² 那一项必须换成外积 u uᵀ（见下）。")
    print()
    print("  为什么是 (1−cosθ) 而不是别的：因为展开里 [u]×² 的系数是")
    print("      −θ²/2! + θ⁴/4! − … = −(1 − cosθ)，而 θ² 项正好补上")
    print("      [u]×³ = −[u]× 无法单独解释的那部分（[u]×² ≠ −I）。")
    print()
    print("  最终：") 
    print("    形式甲（含 I）：    exp(θ[u]×) = I + sinθ·[u]× + (1−cosθ)·[u]×²")
    print("    形式乙（标准教科书）：exp(θ[u]×) = cosθ·I + sinθ·[u]×"
          " + (1−cosθ)·(u uᵀ)")
    print()
    print("  两者是**同一个矩阵**，只因恒等式 [u]×² = u uᵀ − I 而互相改写：")
    print("      I + (1−cosθ)[u]×² = I + (1−cosθ)(u uᵀ − I)")
    print("                        = cosθ·I + (1−cosθ)·u uᵀ")
    uu = np.outer(u, u)
    for th in [0.4, 1.6, 2.7]:
        c, s = np.cos(th), np.sin(th)
        fA = np.eye(3) + s * K + (1 - c) * (K @ K)
        fB = c * np.eye(3) + s * K + (1 - c) * uu
        R = RC.rot_axis(u, th)
        print(f"    θ={th}: |甲 − R| = {np.max(np.abs(fA - R)):.2e}   "
              f"|乙 − R| = {np.max(np.abs(fB - R)):.2e}   |甲 − 乙| = "
              f"{np.max(np.abs(fA - fB)):.2e}")
    print("  → 甲乙都对，且互相差 0。**别写成 cosθ·I + … + (1−cosθ)·[u]×²**（那是错的，")
    print("    把 I 的系数从 1 改成 cosθ 时必须同时把 [u]×² 换成 u uᵀ）。")

    print("\n" + "=" * 72)
    print("【3】两条独立实现给出同一个矩阵（判据不共享代码路径）")
    print("=" * 72)
    print(f"  {'轴':<18}{'θ':>8}{'max|闭式 − 级数|':>22}")
    worst = 0.0
    for ax in [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 1, 1], [0.3, -0.7, 0.2]]:
        for th in [0.001, 0.5, 1.0, 2.0, 3.0, -1.5]:
            d = np.max(np.abs(RC.rot_axis(ax, th) - RC.rot_expm_series(ax, th, 80)))
            worst = max(worst, d)
            if th in (1.0, 3.0):
                print(f"  {str(np.round(RC.unit(ax), 3)):<18}{th:>8.3f}{d:>22.3e}")
    print(f"  → 36 组组合最大偏差 {worst:.3e}  "
          f"[{'PASS' if worst < 1e-12 else 'FAIL'}]")

    print("\n" + "=" * 72)
    print("【4】结果确实落在 SO(3)（正交 + det=+1），随机 2000 组")
    print("=" * 72)
    rng = np.random.default_rng(11)
    worst_o = worst_d = 0.0
    for _ in range(2000):
        R = RC.rot_axis(rng.normal(size=3), rng.uniform(-np.pi, np.pi))
        _, o, d = RC.is_rotation(R)
        worst_o, worst_d = max(worst_o, o), max(worst_d, d)
    print(f"  正交误差 max|RᵀR − I| = {worst_o:.3e}")
    print(f"  行列式   max|det R − 1| = {worst_d:.3e}")
    print("  → exp 把**反对称**（so(3)）映射成**正交且保向**（SO(3)）。")
    print("    这条不是巧合：exp(K)ᵀexp(K) = exp(Kᵀ)exp(K) = exp(−K)exp(K) = I。")

    print("\n" + "=" * 72)
    print("【5】已知角度的具体值（肉眼可校验）")
    print("=" * 72)
    R = RC.rot_xyz("Z", np.pi / 2)
    print("  Rz(90°) =")
    print(RC.fmt_mat(R))
    print("  → 第一列 (1,0,0) 被送到 (0,1,0)：x 轴转到了 y 轴 ✓")
    R = RC.rot_axis([1, 1, 1], 2 * np.pi / 3)
    print("\n  绕 (1,1,1) 转 120°（立方体的一个三重轴）= 循环置换 e1→e2→e3→e1：")
    print(RC.fmt_mat(R, digits=4))
    print(f"  R³ = I ? {np.allclose(R @ R @ R, np.eye(3))}   "
          f"det = {np.linalg.det(R):.10f}")
    return worst


def small_angle():
    """6. 小角极限：exp 把"无穷小旋转"送回 so(3)。"""
    print("\n" + "=" * 72)
    print('【6】小角极限：(R(θ) − I)/θ → [u]×（李代数是这样被「取回来」的）')
    print("=" * 72)
    u = RC.unit([0.4, -0.9, 0.2])
    K = RC.hat(u)
    print(f"  {'θ':>10}{'1/θ·max|R(θ) − I − θ[u]×|':>32}")
    for th in [0.5, 0.1, 0.01, 0.001, 1e-4]:
        R = RC.rot_axis(u, th)
        err = np.max(np.abs(R - np.eye(3) - th * K)) / th
        print(f"  {th:>10.4f}{err:>32.3e}")
    print("  → 误差 ∝ θ（不是 θ²）：因为 (1−cosθ) ≈ θ²/2，而 θ²/θ = θ。")
    print('    这就是「李代数是李群在单位元处的切空间」的数值化身：')
    print("    R(θ) − I ≈ θ[u]× ⇒ 微分即得 so(3)。")


# ------------------------------------------------------------ 画图
def make_figure():
    fig = plt.figure(figsize=(15.5, 5.6))
    u = RC.unit([1.0, 1.0, 0.4])

    # (a) 3D：轴 + 一个向量被转的轨迹
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    RC.draw_sphere_wire(ax, alpha=0.06)
    RC.draw_axes(ax, 1.15)
    ax.plot([-1.1 * u[0], 1.1 * u[0]], [-1.1 * u[1], 1.1 * u[1]],
            [-1.1 * u[2], 1.1 * u[2]], color="#8e44ad", lw=2.4, ls="--")
    ax.text(*(1.2 * u), "转轴 u", color="#8e44ad", fontsize=11)
    v0 = RC.unit([0.2, 0.9, 0.35])
    v0 = RC.unit(v0 - np.dot(v0, u) * u)         # 保证 ⊥ u，轨迹是圆
    ths = np.linspace(0, 2 * np.pi, 200)
    traj = np.array([RC.rot_axis(u, t) @ v0 for t in ths])
    ax.plot(*traj.T, color="#c0392b", lw=1.6, alpha=0.85)
    for t in np.linspace(0, 2 * np.pi * 0.83, 7):
        p = RC.rot_axis(u, t) @ v0
        ax.quiver(0, 0, 0, *p, color="#c0392b", lw=1.5, alpha=0.5,
                  arrow_length_ratio=0.1)
    ax.quiver(0, 0, 0, *v0, color="#185FA5", lw=2.6, arrow_length_ratio=0.12)
    ax.text(*v0, "v", color="#185FA5", fontsize=11)
    ax.set_box_aspect((1, 1, 1))
    ax.set_title("① 绕轴 u 转：v 画出一个圆\n（转轴上的分量一动不动）")
    ax.set_axis_off()

    # (b) 三条恒等式的「数值证据」
    ax = fig.add_subplot(1, 3, 2)
    th = np.linspace(0, 4 * np.pi, 400)
    c = np.cos(th)
    s = np.sin(th)
    ax.plot(np.degrees(th), c, lw=2.2, label=r"$\cos\theta$（偶次项的和）")
    ax.plot(np.degrees(th), s, lw=2.2, label=r"$\sin\theta$（奇次项的和）")
    ax.plot(np.degrees(th), 1 - c, lw=2.2, color="#0F6E56",
            label=r"$1-\cos\theta$（$[u]_\times^2$ 的系数）")
    ax.axhline(0, color="k", lw=0.6)
    for k in range(1, 5):
        ax.axvline(90 * k, ls=":", color="k", lw=0.6, alpha=0.5)
    ax.set_xlabel("转角 θ（度）")
    ax.set_ylabel("系数")
    ax.set_title("② 三个系数的图像\n罗德里格斯 = I·1 + [u]×·sinθ + [u]×²·(1−cosθ)")
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    ax.set_xlim(0, 720)

    # (c) 小角：误差 ∝ θ
    ax = fig.add_subplot(1, 3, 3)
    us = RC.unit([0.4, -0.9, 0.2])
    K = RC.hat(us)
    ths2 = np.logspace(-6, 0, 40)
    errs = [np.max(np.abs(RC.rot_axis(us, t) - np.eye(3) - t * K)) for t in ths2]
    ax.loglog(ths2, errs, "o-", ms=3.5, lw=1.4, color="#185FA5",
              label=r"$|R(\theta)-I-\theta[u]_\times|$")
    ax.loglog(ths2, ths2 ** 2 / 2, "--", color="gray", lw=1.4,
              label=r"参考 $\theta^2/2$（二阶）")
    ax.loglog(ths2, ths2, "-", color="#c0392b", lw=1.2, alpha=0.7,
              label=r"参考 $\theta$（一阶）")
    ax.set_xlabel("转角 θ（弧度）")
    ax.set_ylabel("绝对误差")
    ax.set_title("③ 小角极限：R(θ) − I ≈ θ[u]×\n误差是一阶 ⇒ 微分给出 so(3)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, which="both")

    fig.suptitle("A5.1 · 罗德里格斯公式：exp(θ[u]×) = I + sinθ [u]× + (1−cosθ) [u]×²",
                 fontsize=13.5, y=0.99)
    fig.subplots_adjust(top=0.80, wspace=0.28)
    path = RC.savefig(fig, "rodrigues_expmap.png")
    plt.close(fig)
    return path


def main():
    print("=" * 72)
    print("A5.1 · 罗德里格斯公式（so(3) → SO(3) 的指数映射闭式）")
    print("=" * 72)
    worst = sphere_cheatsheet()
    small_angle()
    path = make_figure()

    print("\n" + "=" * 72)
    print(f"验收：闭式 vs 级数 最大偏差 = {worst:.3e}  "
          f"[{'PASS' if worst < 1e-12 else 'FAIL'}]")
    print("图：", path)
    print("\n这一步在说什么：")
    print("  · 叉积 → 反对称矩阵（hat），这就是李代数 so(3)；")
    print("  · exp 把 so(3) 送到 SO(3)，闭式就是罗德里格斯公式；")
    print("  · 全靠 [u]×³ = −[u]× 这一条恒等式，级数才能收成 sin 与 (1−cos)。")
    print("  下一步（A5.2）看：同一个 exp 换到四元数语言会漂亮成什么样。")


if __name__ == "__main__":
    main()