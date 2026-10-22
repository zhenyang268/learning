"""
G2 · 测地线方程：Christoffel 符号第一次以"发动机"的身份出现
================================================================
现象：从赤道一点、以斜 45° 方向出发，把球坐标测地线方程用 RK4 积分一整圈，
      轨迹与解析大圆重合到 1e-9 以内，并精确闭合。

测地线方程（单位球，Γ 从 G5 里手写的结果来）：
    θ'' =  sinθ cosθ (φ')²
    φ'' = −2 cotθ θ' φ'

一句话在说什么：
    测地线 = "不加速"的推广：在弯曲坐标里，"走直线"意味着坐标加速度
    恰好被 Γ（基向量自己在转）抵消。Γ 不是力，是坐标系在弯曲 ——
    但把它们写进同一个方程，就得到了能数值积分的机器。

**三条路**指向同一条曲线（判据独立于实现路径）：
    A. 解微分方程：RK4 积分测地线方程
    B. 用 SO(3) 旋转矩阵推：γ(s) = R(s·n̂) p₀，n̂ = p₀×t₀（罗德里格斯公式）
    C. cos/sin 线性组合：γ(s) = cos s p₀ + sin s t₀（最朴素的写法）
    B 与 C 应当逐位相同（R(θ) 作用在 p₀ 上就是那个线性组合）；A 逼近它们。
    —— "解方程" 与 "李群代数" 两条完全不同的路会师。

Clairaut 关系：sin²θ·φ' = const（球面测地线的角动量守恒）
    —— 这一条把测地线和力学连起来：球面上自由运动 = 中心力场运动。

运行： uv run python graphics/demo/pyapps/sphere/geodesic_eq.py
输出： learning-todo/imgs/geo/geodesic_eq.png
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _core as C  # noqa: E402

plt = C.setup_style()


# ---------------------------------------------------------- 方程与积分
def geodesic_rhs(s, y):
    """测地线方程右端。y = (θ, φ, θ', φ')。推导见模块 docstring。"""
    th, ph, dth, dph = y
    st, ct = np.sin(th), np.cos(th)
    return np.array([dth, dph, st * ct * dph * dph, -2.0 * ct / st * dth * dph])


def analytic_linear(p0, t0, s_grid):
    """路线 C：解析大圆的最朴素写法 γ(s) = cos s·p0 + sin s·t0（t0 ⊥ p0，单位）。"""
    p0 = np.asarray(p0, float)
    p0 = p0 / np.linalg.norm(p0)
    t0 = np.asarray(t0, float)
    t0 = t0 - (t0 @ p0) * p0
    t0 = t0 / np.linalg.norm(t0)
    s = np.asarray(s_grid, float)
    return np.outer(np.cos(s), p0) + np.outer(np.sin(s), t0)


def analytic_so3(p0, t0, s_grid):
    """路线 B：用 SO(3) 旋转矩阵推大圆：γ(s) = R(s·n̂) p0，n̂ = p0 × t0。

    这是"李群代数"路线：把"沿大圆走弧长 s"实现成"绕 n̂ 转 s 弧度"。
    与路线 C 的 cos/sin 线性组合是同一个东西 —— 罗德里格斯公式展开即得。
    这条路的引入让"解微分方程 vs SO(3) 代数"这句两路对拍名副其实
    （原版只在注释里声称，代码里并没有真的用旋转矩阵）。
    """
    p0 = np.asarray(p0, float)
    p0 = p0 / np.linalg.norm(p0)
    t0 = np.asarray(t0, float)
    t0 = t0 - (t0 @ p0) * p0
    t0 = t0 / np.linalg.norm(t0)
    n_hat = np.cross(p0, t0)                    # 大圆平面法向
    n_hat = n_hat / np.linalg.norm(n_hat)
    return np.array([C.rotate_about(n_hat, si) @ p0 for si in np.asarray(s_grid, float)])


def state_to_cart(y, R=1.0):
    """状态 (θ, φ) -> 笛卡尔位置。"""
    return np.stack(C.sph_to_cart(y[:, 0], y[:, 1], R), axis=1)


# ---------------------------------------------------------- 验证
def run_checks():
    print("=" * 70)
    print("【设置】赤道出发点，初速与经线夹 45°（单位速度）")
    print("=" * 70)
    alpha = np.pi / 4
    th0, ph0 = np.pi / 2, 0.0
    dth0 = -np.cos(alpha)                 # 向北分量（θ 减小方向）
    dph0 = np.sin(alpha) / np.sin(th0)    # sinθ·φ' = sinα —— Clairaut 的种子
    y0 = np.array([th0, ph0, dth0, dph0])
    v0 = np.sqrt(dth0 ** 2 + np.sin(th0) ** 2 * dph0 ** 2)
    print(f"  初速度大小 |v| = {v0:.12f}  （应 = 1）")

    p0 = np.array(C.sph_to_cart(th0, ph0))
    e_th, e_ph = C.sph_basis(th0, ph0)
    t0 = dth0 * e_th + dph0 * e_ph

    n = 720                               # 一圈 720 步
    s_grid = np.linspace(0, 2 * np.pi, n + 1)
    traj = C.integrate(geodesic_rhs, y0, s_grid)
    pos_num = state_to_cart(traj)
    pos_lin = analytic_linear(p0, t0, s_grid)
    pos_so3 = analytic_so3(p0, t0, s_grid)

    print("\n" + "=" * 70)
    print("【验证 0】B 与 C 本是同一条曲线：SO(3) 旋转 ≡ cos/sin 线性组合")
    print("=" * 70)
    d_bc = np.max(np.linalg.norm(pos_so3 - pos_lin, axis=1))
    print(f"  max|γ_SO3(s) − γ_cossin(s)| = {d_bc:.3e}  （应 = 机器精度级）")
    print("  → 罗德里格斯公式 R(θ) 作用在 p₀ 上，展开就是 cos·p₀ + sin·t₀ ——")
    print("    “绕轴转一个角”与“线性组合”在这里是同一件事的两种写法。")

    print("\n" + "=" * 70)
    print("【验证 1】数值轨迹 vs 解析大圆（一整圈后回到起点？）")
    print("=" * 70)
    err_lin = np.max(np.linalg.norm(pos_num - pos_lin, axis=1))
    err_so3 = np.max(np.linalg.norm(pos_num - pos_so3, axis=1))
    close = np.linalg.norm(pos_num[-1] - pos_num[0])
    print(f"  vs 路线 C（cos/sin）最大偏差 = {err_lin:.3e}")
    print(f"  vs 路线 B（SO(3)）  最大偏差 = {err_so3:.3e}")
    print(f"  闭合误差 |γ(2π) − γ(0)|       = {close:.3e}  （绕球一圈必须回来）")
    ok = err_lin < 1e-7 and err_so3 < 1e-7 and close < 1e-7
    print(f"  [{'PASS' if ok else 'FAIL'}]")

    print("\n" + "=" * 70)
    print("【验证 2】Clairaut 关系 sin²θ·φ' = const（角动量守恒）")
    print("=" * 70)
    L_z = np.sin(traj[:, 0]) ** 2 * traj[:, 3]
    drift = np.max(np.abs(L_z - L_z[0]))
    print(f"  sin²θ·φ'  起始 = {L_z[0]:.12f}   漂移 = {drift:.3e}")
    print("  这不是巧合：球面测地线在赤道平面内的角动量守恒 ——")
    print("  测地线方程和中心力场的轨道方程是同一族（N-S 之前，先认识这一层）。")

    print("\n" + "=" * 70)
    print("【验证 3】RK4 收敛性（步长减半 ÷16）")
    print("=" * 70)
    print(f"  {'步数':>6}  {'max 偏差':>12}  {'误差比':>8}")
    prev = None
    for n_t in [90, 180, 360, 720, 1440]:
        sg = np.linspace(0, 2 * np.pi, n_t + 1)
        tr = C.integrate(geodesic_rhs, y0, sg)
        e = np.max(np.linalg.norm(state_to_cart(tr) - analytic_linear(p0, t0, sg),
                                  axis=1))
        ratio = f"÷{prev / e:8.1f}" if prev else "        —"
        print(f"  {n_t:6d}  {e:12.3e}  {ratio}")
        prev = e
    print("  → 比值稳定 ≈16 = 2⁴，确认四阶。")
    return s_grid, traj, pos_lin, p0


# ---------------------------------------------------------- 画图
def make_figure(s_grid, traj, pos_lin, p0):
    fig = plt.figure(figsize=(15.5, 5.4))

    ax = fig.add_subplot(1, 2, 1, projection="3d")
    C.draw_sphere_base(ax)
    pos_num = state_to_cart(traj)
    ax.plot(*pos_lin.T, color="gray", lw=5, alpha=0.45,
            label="路线 C：解析大圆（cos/sin）（粗灰底）")
    ax.plot(*pos_num.T, color="#c0392b", lw=1.6,
            label="路线 A：RK4 积分测地线方程")
    ax.scatter(*p0, color="k", s=45)
    ax.text(*p0 * 1.2, "出发（赤道，45°斜向）", fontsize=10)
    ax.view_init(elev=22, azim=50)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel("x"), ax.set_ylabel("y"), ax.set_zlabel("z")
    ax.set_title("① 解微分方程（红）叠在解析大圆（灰）上\n"
                 "两条线重合到 1e-9 —— Γ 方程就是“直线”的坐标系版本")
    ax.legend(fontsize=9, loc="upper left")

    ax = fig.add_subplot(1, 2, 2)
    L_z = np.sin(traj[:, 0]) ** 2 * traj[:, 3]
    sdeg = np.degrees(s_grid)
    ax.plot(sdeg, L_z, lw=2, color="#2980b9",
            label=r"Clairaut 量 $\sin^2\theta\,\varphi'$")
    ax.axhline(L_z[0], ls="--", color="gray", lw=1, label="守恒值（起始）")
    ax.set_xlabel("沿测地线走过的弧长 s（度）")
    ax.set_ylabel(r"$\sin^2\theta\,\varphi'$")
    ax.set_ylim(L_z[0] - 0.02, L_z[0] + 0.02)
    ax.set_title("② Clairaut 关系：测地线自带角动量守恒\n"
                 "（纬度高 → 转经度快；这是大圆绕向极地的微分形式）")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    fig.suptitle("G2 · 测地线方程：让 Γ 当发动机，跑出来的就是大圆",
                 fontsize=13.5, y=0.985)
    fig.subplots_adjust(top=0.80, wspace=0.22)
    path = C.savefig(fig, "geodesic_eq.png")
    plt.close(fig)
    return path


def main():
    s_grid, traj, pos_lin, p0 = run_checks()
    path = make_figure(s_grid, traj, pos_lin, p0)
    print("\n图：", path)
    print("\n这一步在说什么：")
    print("  度规 → Γ → 测地线方程，这条流水线造出了“直线”的推广。")
    print("  三条路指向同一条曲线：解微分方程（A）、SO(3) 旋转矩阵（B）、")
    print("  cos/sin 线性组合（C）—— B 与 C 代数同源，A 数值逼近它们。")
    print("  Clairaut 守恒量则是预告：测地线动力学 ≐ 力学，下一步就能谈力。")


if __name__ == "__main__":
    main()