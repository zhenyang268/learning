"""
G3 · 平行移动：沿球面曲线搬运向量，一圈之后回不到原方向
================================================================
现象：把向量沿纬线圈平行移动一整圈，它转过了 −2π·cos θ₀ ——
      赤道（θ₀=90°）转 0°，θ₀=60° 转 −180°，越靠近极点越接近 ∓360°。

一句话在说什么：
    "向量的方向"在曲面上**没有绝对意义** —— 不同点的切平面是不同的空间，
    比较两个向量必须先规定"怎么搬"。平行移动就是最自然的规定（方向尽量不变），
    但搬一圈回来发现**必须差一个角度**，这个角度由圈住的曲率决定 ——
    这就是曲率最可触摸的样子（holonomy）。

**注意符号与量法**（原版这里口径混乱，已统一）：
    · 方程给出的是"在**正交标架**下连续展开"的转角 = −2π cos θ₀。
    · 直接比较首末向量（欧氏内积）得到的是它的 mod 360° 代表元：
      θ₀=60° → −180° 与 +180° 等价；θ₀=10° → −354.5° 与 +5.47° 等价。
    · 两者**不是矛盾**，是同一个角的两种读法。本脚本两个都给，并校核 mod 360°。

两套实现互为对拍（这正是 G5 图里那个教学点的续集）：
    A. 球坐标 + Christoffel：dV^k/ds + Γ^k_{ij} V^i γ̇^j = 0
       → 纬线上化简为 dV^θ/dφ = sinθ₀cosθ₀ V^φ,  dV^φ/dφ = −cotθ₀ V^θ
    B. 笛卡尔嵌入：dV/ds = −(V·γ̇) γ
       （球面上 ∇ = "ℝ³ 导数的切向投影"；V·γ=0 保持不变 ⇒ 只剩法向修正项）
    B 在极点照样能跑 —— 坐标奇点 ≠ 几何奇点。

自检：|V| = √((V^θ)² + sin²θ₀ (V^φ)²) 漂移 < 1e-8；RK4 步长减半误差 ÷16。

运行： uv run python graphics/demo/pyapps/sphere/parallel_transport.py
输出： learning-todo/imgs/geo/transport_latitude.png
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _core as C  # noqa: E402

plt = C.setup_style()


# ---------------------------------------------------------- 纬线圈解析解
def latitude_transport_exact(theta0):
    """纬线 θ=θ₀ 平行移动的**解析解**（用矩阵指数的闭式）。

    球坐标方程（φ 为参数）：
        d/dφ [V^θ, V^φ]ᵀ = M [V^θ, V^φ]ᵀ,   M = [[0, sinθ₀cosθ₀], [−cotθ₀, 0]]
    M² = −cos²θ₀ · I  ⇒  exp(φM) = cos(cφ) I + sin(cφ) M / c,  c = cosθ₀
    （这正是"二维旋转"，holonomy = 转角 = −2π cosθ₀ —— 在正交标架下量。）

    返回 (V_end 球坐标分量, 转角_rad)。
    """
    c = np.cos(theta0)
    M = np.array([[0.0, np.sin(theta0) * c],
                  [-c / np.sin(theta0), 0.0]])
    phi = 2 * np.pi
    if abs(c) < 1e-14:                    # 赤道：M = 0，向量不动
        E = np.eye(2)
    else:
        E = np.cos(c * phi) * np.eye(2) + np.sin(c * phi) * M / c
    V_end = E @ np.array([1.0, 0.0])
    angle = -2 * np.pi * c                # 正交标架下的连续转角
    return V_end, angle


# ---------------------------------------------------------- 路线 A：球坐标
def transport_latitude_spherical(theta0, n_steps=2000):
    """球坐标 RK4：沿 θ=θ₀ 纬线平行移动 V₀=∂_θ 一圈。

    返回 (phi 网格, V 轨迹 (N,2))。方程推导见模块 docstring。
    """
    def f(phi, y):
        Vth, Vph = y
        return np.array([np.sin(theta0) * np.cos(theta0) * Vph,
                         -np.cos(theta0) / np.sin(theta0) * Vth])

    phi_grid = np.linspace(0, 2 * np.pi, n_steps + 1)
    traj = C.integrate(f, np.array([1.0, 0.0]), phi_grid)
    return phi_grid, traj


# ---------------------------------------------------------- 路线 B：笛卡尔
def latitude_gamma(theta0, phi0=0.0):
    """纬线的解析描述：返回 gamma_fn(s) -> (位置 γ, 单位切向 γ̇)。

    γ(s) = (sinθ₀ cos(φ₀+s/sinθ₀), sinθ₀ sin(φ₀+s/sinθ₀), cosθ₀)，单位速度。
    """
    st = np.sin(theta0)

    def gamma_fn(s):
        ph = phi0 + s / st
        gamma = np.array([st * np.cos(ph), st * np.sin(ph), np.cos(theta0)])
        gdot = np.array([-np.sin(ph), np.cos(ph), 0.0])
        return gamma, gdot

    return gamma_fn


# ---------------------------------------------------------- 自检与收敛
def checks():
    theta0 = np.pi / 3

    print("=" * 70)
    print("【自检 1】|V| 守恒（球坐标路线，θ₀=60°）")
    print("=" * 70)
    phi_grid, traj = transport_latitude_spherical(theta0)
    norms = np.sqrt(traj[:, 0] ** 2 + np.sin(theta0) ** 2 * traj[:, 1] ** 2)
    drift = np.max(np.abs(norms - 1.0))
    print(f"  起始 |V| = {norms[0]:.15f}   末尾 |V| = {norms[-1]:.15f}")
    print(f"  全程最大漂移 = {drift:.3e}   [{'PASS' if drift < 1e-8 else 'FAIL'}]")
    print("  → 平行移动是等距（保长度）：它只旋转方向，不改变大小。")

    print("\n" + "=" * 70)
    print("【自检 2】RK4 收敛性：步长减半，误差应 ÷16（四阶）")
    print("=" * 70)
    V_exact, _ = latitude_transport_exact(theta0)
    print(f"  {'步长 h':>10}  {'max|ΔV|':>12}  {'误差比':>8}")
    prev = None
    for n in [125, 250, 500, 1000, 2000]:
        _, tr = transport_latitude_spherical(theta0, n_steps=n)
        err = np.max(np.abs(tr[-1] - V_exact))
        ratio = f"{prev / err:8.1f}" if prev else "     —"
        print(f"  {2 * np.pi / n:10.5f}  {err:12.3e}  {ratio}")
        prev = err
    print("  → 比值稳定 ≈16 = 2⁴，确认四阶。")

    print("\n" + "=" * 70)
    print("【自检 3】两套实现 + 解析解三方对拍（holonomy，θ₀=60°）")
    print("=" * 70)
    V_exact, ang_exact = latitude_transport_exact(theta0)
    _, traj = transport_latitude_spherical(theta0)
    err_sph = np.max(np.abs(traj[-1] - V_exact))

    # 笛卡尔路线：同一条纬线（解析 γ，不用点列折线）
    e_th, _ = C.sph_basis(theta0, 0.0)
    s_grid = np.linspace(0, 2 * np.pi * np.sin(theta0), 2001)
    V_cart = C.transport_along(latitude_gamma(theta0), s_grid, e_th)
    e_th_e, e_ph_e = C.sph_basis(theta0, 0.0)
    V_exact_cart = V_exact[0] * e_th_e + V_exact[1] * e_ph_e
    err_cart = np.linalg.norm(V_cart[-1] - V_exact_cart)

    print(f"  解析转角（正交标架，连续展开） = {np.degrees(ang_exact):+.4f}°")
    print(f"  球坐标 RK4 误差 vs 解析 = {err_sph:.3e}")
    print(f"  笛卡尔 RK4 误差 vs 解析 = {err_cart:.3e}")
    ok = err_sph < 1e-8 and err_cart < 1e-8
    print(f"  → {'PASS' if ok else 'FAIL'}：两套完全不同的公式给出同一个向量。")
    return drift


# ---------------------------------------------------------- holonomy 扫描
def holonomy_scan():
    """对一串纬度扫描 holonomy，验证 = −2π cos θ₀（连续展开口径）。"""
    print("\n" + "=" * 70)
    print("【holonomy 扫描】沿纬线一圈的转角 vs 理论 −2π·cos θ₀")
    print("=" * 70)
    print(f"  {'θ₀(度)':>8}  {'数值(度)':>12}  {'理论(度)':>12}  {'差':>10}  {'mod360 校验':>12}")
    worst = 0.0
    rows = []
    for deg in [10, 30, 50, 70, 90, 110, 130, 150, 170]:
        t0 = np.radians(deg)
        _, traj = transport_latitude_spherical(t0, n_steps=4000)
        # 正交标架分量 W = (V^θ, V^φ sinθ₀)，角度用连续展开
        W = np.stack([traj[:, 0], np.sin(t0) * traj[:, 1]], axis=1)
        ang = np.unwrap(np.arctan2(W[:, 1], W[:, 0]))
        num = ang[-1] - ang[0]
        th_exact = -2 * np.pi * np.cos(t0)
        d = abs(num - th_exact)
        # mod 360 校验：连续展开角与"直接比首末向量"的那个角同余
        # 正交标架分量 = (V^θ, V^φ·sinθ₀)（这是**乘**不是除）
        v0 = np.array([1.0, 0.0])
        v_end = np.array([traj[-1][0], traj[-1][1] * np.sin(t0)])
        ang_dir = np.degrees(np.arctan2(v_end[1], v_end[0]))
        mod_ok = abs(((np.degrees(num) - ang_dir + 180) % 360) - 180) < 1e-6
        worst = max(worst, d)
        rows.append((deg, num, th_exact))
        print(f"  {deg:8.0f}  {np.degrees(num):12.4f}  "
              f"{np.degrees(th_exact):12.4f}  {np.degrees(d):10.2e}  "
              f"{'✓' if mod_ok else '✗':>12}")
    print(f"  → 最大偏差 {np.degrees(worst):.2e} 度")
    print("  注：θ₀→0（贴着极点的小圈）转角 → −360°≡0°：极点附近几何是平的。")
    print("      赤道（90°）转角 = 0：赤道是大圆=测地线，切向自己就是平行移动的。")
    return rows


# ---------------------------------------------------------- 画图
def make_figure(theta0=np.pi / 3, rows=None):
    fig = plt.figure(figsize=(15.5, 5.8))

    # ---- (a) 3D：纬线上的平行移动 ----
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    C.draw_sphere_base(ax)
    phi = np.linspace(0, 2 * np.pi, 400)
    pts = np.stack(C.sph_to_cart(np.full_like(phi, theta0), phi), axis=1)
    ax.plot(*pts.T, "k-", lw=1.5, label=f"纬线圈 θ₀={np.degrees(theta0):.0f}°")

    s_grid = np.linspace(0, 2 * np.pi * np.sin(theta0), 400)
    V_cart = C.transport_along(latitude_gamma(theta0), s_grid,
                               C.sph_basis(theta0, 0.0)[0])
    idx = np.linspace(0, len(pts) - 2, 9).astype(int)
    L = 0.28
    for i in idx:
        p = pts[i]
        ax.quiver(*p, *(L * V_cart[i]), color="#c0392b", lw=2,
                  arrow_length_ratio=0.25)
    ax.quiver(*pts[0], *(L * V_cart[0]), color="blue", lw=2.5,
              arrow_length_ratio=0.25, label="初始向量")
    ax.quiver(*pts[0], *(L * V_cart[-1]), color="green", lw=2.5,
              arrow_length_ratio=0.25, label="绕一圈后回到起点")
    ax.set_title(f"平行移动一圈：蓝→绿转了 "
                 f"{np.degrees(C.angle_between(V_cart[0], V_cart[-1])):.1f}°\n"
                 f"（红箭头 = 沿途的向量，方向“尽量不变”，却被迫整体转了）")
    ax.set_xlabel("x"), ax.set_ylabel("y"), ax.set_zlabel("z")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_box_aspect((1, 1, 1))

    # ---- (b) holonomy vs 纬度 ----
    ax = fig.add_subplot(1, 2, 2)
    if rows:
        degs = [r[0] for r in rows]
        nums = [np.degrees(r[1]) for r in rows]
        ax.plot(degs, nums, "o", ms=7, color="#c0392b",
                label="数值（球坐标 RK4 + 连续展开）")
    d = np.linspace(1, 179, 300)
    ax.plot(d, -360 * np.cos(np.radians(d)), "-", color="gray", lw=1.5,
            label=r"理论  $-2\pi\cos\theta_0$")
    ax.axhline(0, color="k", lw=0.6)
    ax.axvline(90, ls=":", color="k", lw=0.8)
    ax.annotate("赤道：转角 0\n（大圆 = 测地线）", xy=(90, 0), xytext=(100, 60),
                fontsize=9, arrowprops=dict(arrowstyle="->", lw=0.8))
    ax.annotate("贴北极的小圈：≈−360°≡0\n（局部近似平面）", xy=(10, -354),
                xytext=(28, -250), fontsize=9,
                arrowprops=dict(arrowstyle="->", lw=0.8))
    ax.set_xlabel("纬线圈所在极角 θ₀（度）")
    ax.set_ylabel("绕一圈的转角（度）")
    ax.set_title("holonomy = −2π cos θ₀\n曲率通过“搬一圈转多少”显形")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    fig.suptitle("G3 · 平行移动：方向没有绝对意义，搬一圈的差角就是曲率",
                 fontsize=13.5, y=0.99)
    fig.subplots_adjust(top=0.82, wspace=0.25)
    path = C.savefig(fig, "transport_latitude.png")
    plt.close(fig)
    return path


def main():
    drift = checks()
    rows = holonomy_scan()
    path = make_figure(rows=rows)

    print("\n" + "=" * 70)
    print(f"验收：|V| 漂移 = {drift:.3e}  "
          f"[{'PASS' if drift < 1e-8 else 'FAIL'}]（阈值 1e-8）")
    print("图：", path)
    print("\n这一步在说什么：")
    print("  平行移动 = 「向量只做被迫的转动」。球面上它有两条等价写法：")
    print("  坐标语言（Γ 方程）和嵌入语言（dV/ds = −(V·γ̇)γ）。")
    print("  搬一圈的转角不为零 —— 缺口就是曲率留下的指纹（holonomy）。")
    print("  θ₀=60° 转 −180°、θ₀=90° 转 0°、贴极点 ≈ −360°≡0°：")
    print("  纬度越靠极点，经线越密集，被迫的扭转越接近一整圈。")


if __name__ == "__main__":
    main()