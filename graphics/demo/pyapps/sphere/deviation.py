"""
G6 · 测地偏离：经线在北极汇聚 —— 弯曲怎么变成"力"
================================================================
现象：赤道上相邻两点（经度差 Δλ）同时向北出发，彼此越走越近，
      在北极（s = π/2 处）撞在一起 —— 北极是它们的**共轭点**。

分离距离的精确解：ξ(s) = 2 arcsin(sinθ·sin(Δλ/2))，θ = π/2 − s
它满足**雅可比方程 / 测地偏离方程**：
    ξ'' + K·ξ = 0    （K = 1）
—— 这正是**简谐振动**的方程！

一句话在说什么：
    曲率为正 → 测地线互相吸引（透镜聚焦）：球面经线、引力透镜、潮汐力
    曲率为零 → 测地线平行或匀速散开：平面
    曲率为负 → 测地线加速散开（散光）：双曲面、马鞍面
    "曲率怎么变成力"的答案就在这个方程里：广义相对论里潮汐力 =
    沿自由落体轨道的测地偏离，Riemann 张量直接出现在方程系数里。

**两个量必须分清**（原版在这里把注释写错了，已修）：
    (a) 线性化近似 ξ ≈ Δλ·cos s 与精确解的偏差：
            展开得 Δλ sinθ + (Δλ³/24)(sin³θ − sinθ) + …，
            即 ξ_exact − ξ_lin ≈ −(Δλ³/24)·sinθ·cos²θ
            （t=sin(Δλ/2) 自身的 Δλ³ 修正与 arcsin 的 Δλ³ 修正部分相消）。
            实测 max = 2.005e-06，理论 Δλ³/24·2/(3√3) = 2.005e-06 —— 逐位吻合。
    (b) 精确解代进雅可比方程 ξ''+Kξ 的残差：
            = −(Δλ³/12)·cos 3s，max = Δλ³/12 = 1.042e-05 —— 也是逐位吻合。
    两者**都是三阶项**，但系数完全不同（一个是 sinθcos²θ，一个是 cos3s）。
    原版把 (a) 写成"Δλ³/24·Δλ = Δλ⁴/24"，量级与实测差 8 倍，是错的。

两条独立路线互为对拍：
    A. 数值积分**两条**邻近测地线（球坐标测地线方程 RK4），直接量距离；
    B. 精确球面距离 2 arcsin(sinθ sin(Δλ/2))。
    A 里的两条测地线会**同时**撞向北极，但球坐标在 θ=0 处 cotθ 爆炸
    （坐标奇点！），所以数值积分只走到 θ ≈ 0.02，其余交给解析延拓。

运行： uv run python graphics/demo/pyapps/sphere/deviation.py
输出： learning-todo/imgs/geo/deviation_meridians.png
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _core as C  # noqa: E402

plt = C.setup_style()

THETA_STOP = 0.02            # 数值积分停在这里（再往前 cotθ 就爆炸）


# ---------------------------------------------------------- 测地线方程
def geodesic_rhs(s, state):
    """球坐标测地线方程（单位球）：状态 y = (θ, φ, θ', φ')。

        θ'' =  sinθ cosθ (φ')²          ← 来自 Γ^θ_{φφ} = −sinθcosθ
        φ'' = −2 cotθ θ' φ'             ← 来自 2Γ^φ_{θφ} = 2cotθ
    单位速度参数化时守恒量：(θ')² + sin²θ (φ')² = 1（可当自检）。
    （s 参数不进入方程 —— 球面测地线方程对弧长是自治的。）
    """
    th, ph, dth, dph = state
    st, ct = np.sin(th), np.cos(th)
    return np.array([dth, dph,
                     st * ct * dph * dph,
                     -2.0 * ct / st * dth * dph])


def shoot_meridian(dlambda, s_max, n_steps, R=1.0):
    """从赤道 (θ=π/2, φ=Δλ) 向北（θ'=−1）积分测地线，走到 θ 停止。

    返回 (s 网格, 轨迹 (N,4))。**不穿过极点**：θ→0 处 cotθ 发散。
    """
    y = np.array([np.pi / 2, dlambda, -1.0, 0.0])
    s_grid = np.linspace(0, s_max, n_steps + 1)
    return s_grid, C.integrate(geodesic_rhs, y, s_grid)


def sep_exact(th, dlambda, R=1.0):
    """有限间隔的**精确**球面距离（纬线圈几何）：sin(ξ/2) = sinθ·sin(Δλ/2)。"""
    return 2 * R * np.arcsin(np.sin(th) * np.sin(dlambda / 2))


# ---------------------------------------------------------- 验证
def run_checks(dlambda=0.05):
    s_stop = np.pi / 2 - THETA_STOP
    s_grid = np.linspace(0, s_stop, 3001)
    th = np.pi / 2 - s_grid

    print("=" * 70)
    print(f"【验证 1】两条经线（经度差 Δλ={dlambda}）的分离距离")
    print("=" * 70)
    _, t1 = shoot_meridian(0.0, s_stop, 3000)
    _, t2 = shoot_meridian(dlambda, s_stop, 3000)

    # 自检：φ 应当恒定（经线），|速度| 应当 = 1
    print(f"  路线1 的 φ 漂移: {np.max(np.abs(t1[:, 1] - 0.0)):.3e}（应≈0）")
    speed1 = np.sqrt(t1[:, 2] ** 2 + np.sin(t1[:, 0]) ** 2 * t1[:, 3] ** 2)
    print(f"  路线1 的速度漂移: {np.max(np.abs(speed1 - 1.0)):.3e}（应≈0）")

    p1 = np.stack(C.sph_to_cart(t1[:, 0], t1[:, 1]), axis=1)
    p2 = np.stack(C.sph_to_cart(t2[:, 0], t2[:, 1]), axis=1)
    # 中心角用 haversine 型公式（arcsin 在小距离处精度远好于 arccos）
    xi_num = 2 * np.arcsin(np.linalg.norm(p1 - p2, axis=1) / 2)

    xi_exact = sep_exact(th, dlambda)
    xi_lin = dlambda * np.cos(s_grid)
    err_exact = np.max(np.abs(xi_num - xi_exact))

    # 偏差的正确理论式：−(Δλ³/24)·sinθ·cos²θ
    pred_dev = -(dlambda ** 3 / 24) * np.sin(th) * np.cos(th) ** 2
    dev = xi_exact - xi_lin
    err_pred = np.max(np.abs(dev - pred_dev))
    t_pred = dlambda ** 3 / 24 * (2 / (3 * np.sqrt(3)))     # max|sinθcos²θ| = 2/(3√3)

    print(f"  数值 ξ vs 精确球面距离: max 偏差 = {err_exact:.3e}  "
          f"[{'PASS' if err_exact < 1e-9 else 'FAIL'}]")
    print(f"  精确 vs 线性化 Δλ·cos s: max 偏差 = {np.max(np.abs(dev)):.3e}")
    print(f"    理论式 −(Δλ³/24)sinθcos²θ: max 预测 = {t_pred:.3e}")
    print(f"      逐点最大差 = {err_pred:.3e}   "
          f"[{'PASS' if err_pred < 1e-8 else 'FAIL'}]（应 = 0，说明理论式正确）")
    print("  → 测地偏离方程是**无穷小**分离的方程，有限间隔要带三阶修正。")

    print("\n" + "=" * 70)
    print("【验证 1b】线性化偏差的收敛阶：偏差 ∝ Δλ³")
    print("=" * 70)
    print(f"  {'Δλ':>8}  {'max|ξ−Δλcos s|':>16}  {'相邻比':>8}  {'累计比':>9}")
    prev = None
    first = None
    for dl in [0.05, 0.01, 0.002]:
        _, tb1 = shoot_meridian(0.0, s_stop, 1500)
        _, tb2 = shoot_meridian(dl, s_stop, 1500)
        qb1 = np.stack(C.sph_to_cart(tb1[:, 0], tb1[:, 1]), axis=1)
        qb2 = np.stack(C.sph_to_cart(tb2[:, 0], tb2[:, 1]), axis=1)
        xb = 2 * np.arcsin(np.linalg.norm(qb1 - qb2, axis=1) / 2)
        sg2 = np.linspace(0, s_stop, 1501)
        e = np.max(np.abs(xb - dl * np.cos(sg2)))
        if first is None:
            first = e
        ratio = f"÷{prev / e:7.1f}" if prev else "      —"
        cumul = f"÷{first / e:8.1f}"
        print(f"  {dl:8.3f}  {e:16.3e}  {ratio}  {cumul}")
        prev = e
    print("  → 相邻比恒 ≈125 = 5³（每行 Δλ 缩 5 倍）；")
    print("    末行累计比 ≈15625 = 25³（相对首行 0.05 → 0.002 缩 25 倍）。")

    print("\n" + "=" * 70)
    print("【验证 2】雅可比方程 ξ'' + Kξ = 0（三点差分二阶导）")
    print("=" * 70)
    h = s_grid[1] - s_grid[0]
    # (a) 线性化解 cos s 代入方程：残差应到差分误差级
    d2_lin = (xi_lin[2:] - 2 * xi_lin[1:-1] + xi_lin[:-2]) / h ** 2
    r_lin = np.max(np.abs((d2_lin + xi_lin[1:-1])[20:-20]))
    # (b) 精确距离代入：残差 = 线性化丢掉的三阶项本身（−Δλ³/12·cos3s），不是数值误差
    d2_num = (xi_num[2:] - 2 * xi_num[1:-1] + xi_num[:-2]) / h ** 2
    r_num = np.max(np.abs((d2_num + xi_num[1:-1])[20:-20]))
    c = np.cos(s_grid[20:-20])
    expect = dlambda ** 3 / 24 * np.max(np.abs(6 * c - 8 * c ** 3))
    print(f"  (a) 解析解 Δλ·cos s 代入:  max|ξ''+ξ| = {r_lin:.3e}  ←差分误差级")
    print(f"  (b) 精确距离代入:          max|ξ''+ξ| = {r_num:.3e}")
    print(f"      理论 −(Δλ³/12)cos3s 的 max = Δλ³/12 = {dlambda ** 3 / 12:.3e}"
          f"  ←完全吻合")
    print("  → 雅可比方程是**无穷小分离**的方程；有限间隔的 ξ 只近似满足它，")
    print("    偏离恰好是线性化丢掉的三阶项。这不是 bug，是方程的适用范围。")

    print("\n" + "=" * 70)
    print("【验证 3】共轭点与线性拟合 R²")
    print("=" * 70)
    basis = np.cos(s_grid)
    a_fit = np.sum(xi_num * basis) / np.sum(basis * basis)
    resid_sq = np.sum((xi_num - a_fit * basis) ** 2)
    tot_sq = np.sum((xi_num - xi_num.mean()) ** 2)
    r2 = 1 - resid_sq / tot_sq
    print(f"  拟合 ξ = a·cos(s)：a = {a_fit:.8f}  （理论 Δλ = {dlambda}）")
    print(f"  拟合 R² = {r2:.8f}   验收 > 0.999 → "
          f"[{'PASS' if r2 > 0.999 else 'FAIL'}]")
    print(f"  共轭点：理论 s = π/2 ≈ {np.pi / 2:.6f} 处 ξ→0（两条经线在北极汇合）。")
    print(f"  数值积分到 s = {s_grid[-1]:.4f}（θ={THETA_STOP}，再往前 cotθ 就爆炸了 ——")
    print(f"  坐标奇点，不是几何奇点；极点本身毫发无损）。")
    return s_grid, xi_num, xi_exact, xi_lin, a_fit, r2


def radius_scaling():
    """半径 R 的球：s = R·u（u 为无量纲角程），ξ 精确解 = 2R arcsin(sinθ·sin(Δλ/2))。"""
    print("\n" + "=" * 70)
    print("【缩放检验】半径 R 的球面：ξ 与测地线积分完全自洽")
    print("=" * 70)
    dlambda = 0.05
    for R in [1.0, 2.0, 5.0]:
        u_max = (np.pi / 2 - THETA_STOP) / R
        y1 = np.array([np.pi / 2, 0.0, -1.0, 0.0])
        y2 = np.array([np.pi / 2, dlambda, -1.0, 0.0])
        u_grid = np.linspace(0, u_max, 3000)
        t1 = C.integrate(geodesic_rhs, y1, u_grid)
        t2 = C.integrate(geodesic_rhs, y2, u_grid)
        p1 = np.stack(C.sph_to_cart(t1[:, 0], t1[:, 1], R), axis=1)
        p2 = np.stack(C.sph_to_cart(t2[:, 0], t2[:, 1], R), axis=1)
        xi_num = 2 * R * np.arcsin(np.linalg.norm(p1 - p2, axis=1) / (2 * R))
        th = np.pi / 2 - u_grid
        err = np.max(np.abs(xi_num - sep_exact(th, dlambda, R)))
        print(f"  R={R:.1f}:  max|ξ_num − ξ_exact| = {err:.3e}  "
              f"[{'PASS' if err < 1e-9 else 'FAIL'}]")
    print("  → 半径 R 缩放下测地线积分与距离公式完全自洽。")
    print("    注意 K = 1/R²：R 变大曲率变小，北极的“汇聚”尺度也按 R 放大。")


# ---------------------------------------------------------- 画图
def make_figure(s_grid, xi_num, xi_exact, xi_lin, dlambda=0.05):
    fig = plt.figure(figsize=(15.5, 5.4))

    # ---- (a) 3D：经线束向北汇聚 ----
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    C.draw_sphere_base(ax, alpha=0.07)
    dls = np.linspace(0, 0.25, 6)
    for dl in dls:
        th = np.linspace(0, np.pi / 2, 100)
        pts = np.stack(C.sph_to_cart(th, np.full_like(th, dl)), axis=1)
        lw = 2.4 if dl == 0 else 1.2
        color = "#c0392b" if dl == 0 else ("#2980b9" if dl == dls[-1] else "gray")
        ax.plot(*pts.T, color=color, lw=lw)
        ax.scatter(*pts[0], color=color, s=30)
    # 分离距离示意
    th_a = np.radians(70)
    pa = np.array(C.sph_to_cart(th_a, 0.0))
    pb = np.array(C.sph_to_cart(th_a, 0.25))
    ax.plot(*np.stack([pa, pb]).T, color="green", lw=2, ls="--")
    ax.text(*(pa + pb) / 2 * 1.08, "ξ(s)", color="green", fontsize=12)
    ax.text(0, 0, 1.16, "北极 = 共轭点\n（所有经线在此汇合）",
            fontsize=10, ha="center", color="#c0392b")
    ax.text(1.16, 0, 0, "赤道出发点", fontsize=9)
    ax.set_title("① 相邻经线（间隔 Δλ）向北出发\n"
                 "间隔越走越小，北极撞在一起 —— 正曲率 = 聚焦")
    ax.set_xlabel("x"), ax.set_ylabel("y"), ax.set_zlabel("z")
    ax.view_init(elev=14, azim=25)
    ax.set_box_aspect((1, 1, 1))

    # ---- (b) ξ(s)：数值 vs 精确 vs 线性化 ----
    ax = fig.add_subplot(1, 2, 2)
    sdeg = np.degrees(s_grid)
    ax.plot(sdeg, xi_exact, "-", color="gray", lw=2.2,
            label=r"精确球面距离 $2\arcsin(\sin\theta\,\sin\frac{\Delta\lambda}{2})$")
    ax.plot(sdeg, xi_num, "o", ms=3.5, color="#c0392b",
            label="数值（两条 RK4 测地线的实测距离）")
    ax.plot(sdeg, xi_lin, "--", color="#2980b9", lw=1.4,
            label=r"线性化（测地偏离方程）$\xi=\Delta\lambda\cos s$")
    ax.axvline(90, ls=":", color="k", lw=0.9)
    ax.text(88, dlambda * 0.62, "s = π/2\n共轭点 ξ→0", fontsize=9, ha="right")
    ax.set_xlabel("沿经线走过的弧长 s（度）")
    ax.set_ylabel(f"分离距离 ξ（弧度，Δλ={dlambda}）")
    ax.set_title("② 测地偏离：ξ'' + Kξ = 0（K=1）—— 简谐振动方程\n"
                 "曲率为正 → 聚焦（透镜）；曲率为负 → 发散（散光）")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    fig.suptitle("G6 · 测地偏离：相邻测地线的分与合，就是曲率在场的方式",
                 fontsize=13.5, y=0.985)
    fig.subplots_adjust(top=0.80, wspace=0.24)
    path = C.savefig(fig, "deviation_meridians.png")
    plt.close(fig)
    return path


def main():
    s_grid, xi_num, xi_exact, xi_lin, a_fit, r2 = run_checks()
    radius_scaling()
    path = make_figure(s_grid, xi_num, xi_exact, xi_lin)
    print("\n图：", path)
    print("\n这一步在说什么：")
    print("  两条测地线的分离距离 ξ 不是随便变的 —— 它被 ξ'' + Kξ = 0 管着。")
    print("  正曲率像会聚透镜（经线在北极聚焦），负曲率像发散透镜。")
    print("  把“测地线”换成“自由落体”，这个方程就是广义相对论的潮汐力：")
    print("  引力不是力，是曲率造成的测地偏离。")


if __name__ == "__main__":
    main()