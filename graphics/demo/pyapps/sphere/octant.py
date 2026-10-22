"""
G4 · 球面八分仪三重验证：角超 = 面积 = holonomy
================================================================
取单位球面上三条大圆弧围成的三角形：
    A = (1,0,0),  B = (0,1,0),  C = (0,0,1)
它恰好是坐标参数区域 {θ,φ} ∈ [0,π/2]²，三条边分别是赤道弧和两条经线弧。

三个完全不同来源的数**同时等于 π/2**：
    1. 内角和 − π = 3·(π/2) − π = π/2        （量角器）
    2. 三角形面积 = 4π/8 = π/2               （几何：八分之一个球）
    3. 平行移动绕一圈的转角 = π/2            （RK4 积分 dV/ds = −(V·γ̇)γ）

一句话在说什么：
    这就是球面角超公式 (α+β+γ) − π = 面积/R² 的实例化，也是 Gauss–Bonnet 的
    测地多边形版本：**角超 = 面积 = holonomy**。平面三角形这三个数全是 0。
    "搬一圈转 90°"这个数字不是巧合，它就是被圈住的那块曲率。

绕行方向与符号（原版这里含糊，已补实测）：
    顶点顺序 A→B→C→A 是**逆时针**（左手法则朝外，实测 cross(AB,AC)·A = +1）。
    · 直接比较首末向量（两者都在 A 的切平面内）得 **+90°**。
    · 在沿途切向标架下连续展开得 **−180°**（向量相对标架一路落后）。
    两者 mod 360° 并不同余 —— 因为前者量的是"首/末向量夹角"，
    后者量的是"向量相对局部标架的落后量"，参照系不是一个。
    本脚本把两个都算出来并如实说明，不再凭记忆写一个数。

运行： uv run python graphics/demo/pyapps/sphere/octant.py
输出： learning-todo/imgs/geo/octant_triangle.png
       learning-todo/imgs/geo/octant_transport.gif
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _core as C  # noqa: E402

plt = C.setup_style()

A = np.array([1.0, 0.0, 0.0])
B = np.array([0.0, 1.0, 0.0])
Cc = np.array([0.0, 0.0, 1.0])
N_OUT = np.array([1.0, 1.0, 1.0]) / np.sqrt(3)   # 八分仪外法向（顶点和方向）
QUARTER = np.pi / 2          # 每条边都是四分之一大圆，弧长 π/2


# ---------------------------------------------------------- 三条边的解析描述
def edge_gamma(p_from, p_to):
    """大圆边的解析描述：γ(s) = cos(s)·from + sin(s)·(n×from)，s ∈ [0, 弧长]。

    n = normalize(from × to)：起终点唯一确定的大圆平面法向。
    单位速度、弧长参数化。
    """
    n = np.cross(p_from, p_to)
    n = n / np.linalg.norm(n)
    t0 = np.cross(n, p_from)                  # 起点单位切向（指向 to 一侧）
    L = C.arc_length(p_from, p_to)

    def gamma_fn(s):
        return (np.cos(s) * p_from + np.sin(s) * t0,
                -np.sin(s) * p_from + np.cos(s) * t0)

    return gamma_fn, L


EDGES = [edge_gamma(A, B), edge_gamma(B, Cc), edge_gamma(Cc, A)]


def silver_traj_pos(i, n):
    """第 i 个轨迹样本对应的球面位置（沿三段均分还原）。作图用。"""
    u = i / (n - 1) * 3 * QUARTER
    k = min(int(u // QUARTER), 2)
    s = u - k * QUARTER
    gamma_fn, _ = EDGES[k]
    return gamma_fn(s)[0]


# ---------------------------------------------------------- 三重验证
def check_angles():
    """内角和：每顶点处两条边（从该顶点出发的方向）的夹角。"""
    print("=" * 70)
    print("【验证 1】内角和")
    print("=" * 70)
    verts = [("A", A, B, Cc), ("B", B, A, Cc), ("C", Cc, A, B)]
    total = 0.0
    for name, P, Q, R in verts:
        # 从 P 出发沿大圆走向 Q / R 的初始切向 = 邻点去掉法向分量
        u = Q - (Q @ P) * P
        v = R - (R @ P) * P
        u, v = u / np.linalg.norm(u), v / np.linalg.norm(v)
        ang = np.arccos(np.clip(u @ v, -1.0, 1.0))
        total += ang
        print(f"  角 {name} = {np.degrees(ang):10.6f}°")
    excess = total - np.pi
    print(f"  内角和  = {np.degrees(total):10.6f}°   （平面三角形是 180°）")
    print(f"  角超    = {np.degrees(excess):10.6f}°   验收阈值 ±0.1° → "
          f"[{'PASS' if abs(np.degrees(excess) - 90) < 0.1 else 'FAIL'}]")
    return total, excess


def check_area():
    """面积：解析值（八分球）与数值双重积分互证。"""
    print("\n" + "=" * 70)
    print("【验证 2】面积")
    print("=" * 70)
    exact = 4 * np.pi / 8
    # 八分仪 = 参数区域 [0,π/2]²，面积元 dA = sinθ dθ dφ（R=1）
    n = 2000
    th = np.linspace(0, np.pi / 2, n)
    thm = (th[1:] + th[:-1]) / 2
    numeric = np.sum(np.sin(thm)) * (th[1] - th[0]) * (np.pi / 2)   # 中点规则
    print(f"  解析（八分球 4π/8）   = {exact:.10f}")
    print(f"  数值 ∫∫ sinθ dθdφ     = {numeric:.10f}   (中点规则, {n} 格)")
    print(f"  两者之差              = {abs(exact - numeric):.3e}   （O(h²) ≈ 4e-8）")
    print("  注：面积是纯几何量（弧度²=sr），本身没有“度”这个单位；")
    print("     下面汇总里把它换算成角度只是为了跟另两个**角**同轴对比——")
    print("     用的是球面恒等式“面积 = 角超”，那正是本脚本要证的结论，")
    print("     所以这个 90° 不是独立测量，而是同一恒等式的另一种读法。")
    return exact


def check_holonomy():
    """holonomy：沿 A→B→C→A 平行移动一圈，V 的转角（两种口径都算）。"""
    print("\n" + "=" * 70)
    print("【验证 3】平行移动绕一圈的转角（holonomy）")
    print("=" * 70)
    n_per = 1500
    V0 = EDGES[0][0](0.0)[1]                       # A 处沿 AB 的切向
    V0 = V0 / np.linalg.norm(V0)
    segs = []
    v = V0.copy()
    for gamma_fn, L in EDGES:
        s_grid = np.linspace(0, L, n_per + 1)
        Vt = C.transport_along(gamma_fn, s_grid, v)
        v = Vt[-1]
        segs.append(Vt)

    # 口径 1：直接比较首末向量（都在 A 的切平面内）
    ang = C.angle_between(V0, v, signed_about=N_OUT)

    # 口径 2：在沿途切向标架 [γ̇, γ×γ̇] 下连续跟踪（每段内标架连续）
    prof = []
    for gi, (gamma_fn, L) in enumerate(EDGES):
        for j, s in enumerate(np.linspace(0, L, n_per + 1)):
            pos, tan = gamma_fn(s)
            e1, e2 = tan, np.cross(pos, tan)
            prof.append(np.arctan2(segs[gi][j] @ e2, segs[gi][j] @ e1))
    prof = np.unwrap(np.array(prof))
    tot_ang = np.degrees(prof[-1] - prof[0])

    print(f"  初向量 V0   = {np.array2string(V0, precision=4)}")
    print(f"  末向量 Vend = {np.array2string(v, precision=4)}")
    print(f"  口径① 首/末向量夹角（A 切平面内，逆时针为正） = "
          f"{np.degrees(ang):+.6f}°   → [{'PASS' if abs(np.degrees(ang) - 90) < 0.1 else 'FAIL'}]")
    print(f"  口径② 沿途标架下连续展开的落后量             = "
          f"{tot_ang:+.6f}°")
    print("  → 两个数**不 mod 360 同余**：口径①量的是两向量夹角（几何量）；")
    print("    口径②量的是向量相对局部切向标架的落后量（含标架自身的转）。")
    print("    三重验证用的是口径①（它与角超、面积同源）。")
    return V0, v, ang, segs


# ---------------------------------------------------------- 静态图
def make_figure(segs, V0, V_end, ang, excess, area):
    fig = plt.figure(figsize=(15.5, 5.4))

    # ---- (a) 3D：三角形 + 沿途平行移动的向量 ----
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    C.draw_sphere_base(ax, alpha=0.07)
    for gamma_fn, L in EDGES:
        s = np.linspace(0, L, 100)
        pts = np.array([gamma_fn(si)[0] for si in s])
        ax.plot(*pts.T, "k-", lw=2)
    for P, name in [(A, "A"), (B, "B"), (Cc, "C")]:
        ax.scatter(*P, color="#c0392b", s=45)
        ax.text(*P * 1.14, name, fontsize=13, weight="bold")

    # 把三段轨迹拼成一条，沿途画向量
    traj = np.vstack(segs)
    n = len(traj)
    idx = np.linspace(0, n - 1, 12).astype(int)
    Lv = 0.22
    for i in idx:
        g = silver_traj_pos(i, n)
        ax.quiver(*g, *(Lv * traj[i]), color="#2980b9", lw=1.6,
                  arrow_length_ratio=0.3)
    ax.quiver(*A, *(Lv * V0), color="green", lw=2.5, arrow_length_ratio=0.3,
              label="初始向量（A 点）")
    ax.quiver(*A, *(Lv * V_end), color="red", lw=2.5, arrow_length_ratio=0.3,
              label="绕一圈后（回到 A 点）")
    ax.set_title(f"① 八分仪三角形：三边都是 π/2 的大圆弧\n"
                 f"蓝向量沿边平行移动，回到 A 时（红）与初始（绿）差 "
                 f"{abs(np.degrees(ang)):.1f}°")
    ax.set_xlabel("x"), ax.set_ylabel("y"), ax.set_zlabel("z")
    ax.legend(fontsize=8, loc="upper left")
    ax.view_init(elev=18, azim=35)
    ax.set_box_aspect((1, 1, 1))

    # ---- (b) 放大：A 点切平面上的错位 ----
    ax = fig.add_subplot(1, 2, 2)
    e1, e2 = np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0])  # A 处切平面基
    p0 = (V0 @ e1, V0 @ e2)
    p1 = (V_end @ e1, V_end @ e2)
    # 从 p0 转到 p1 的圆弧（带符号角）
    a0 = np.arctan2(p0[1], p0[0])
    th_arc = np.linspace(0, ang, 60) + a0
    ax.plot(np.cos(th_arc), np.sin(th_arc), "-", color="gray", lw=1.2)
    ax.annotate("", xy=p1, xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", color="red", lw=2.2))
    ax.annotate("", xy=p0, xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", color="green", lw=2.2))
    ax.text(p0[0] * 1.03, p0[1] * 1.10 - 0.06, "初始方向（沿 AB 出发）",
            fontsize=10, color="green")
    ax.text(p1[0] * 1.03 + 0.05, p1[1] - 0.02, "绕一圈回来（沿 CA 到达）",
            fontsize=10, color="red")
    ax.text(0.10, -0.32, f"转角 = {np.degrees(ang):.2f}°", fontsize=13,
            color="#c0392b", weight="bold")
    ax.axhline(0, color="k", lw=0.5), ax.axvline(0, color="k", lw=0.5)
    ax.set_xlim(-0.35, 1.35), ax.set_ylim(-0.45, 1.25)
    ax.set_aspect("equal")
    ax.set_xlabel("切平面基 e1 = (0,1,0)")
    ax.set_ylabel("切平面基 e2 = (0,0,1)")
    ax.set_title("② A 点切平面放大：起点与终点的方向错位\n"
                 "（终点回到同一点，切平面相同，比较才有意义）")
    ax.grid(alpha=0.3)

    fig.suptitle(
        "G4 · 三重验证：内角和−π = 面积 = 平行移动转角 = 90°",
        fontsize=12.5, y=1.01)
    path = C.savefig(fig, "octant_triangle.png")
    plt.close(fig)
    return path


# ---------------------------------------------------------- 动画
def make_gif(V0, n_frames=90):
    from matplotlib.animation import FuncAnimation, PillowWriter

    fig = plt.figure(figsize=(7.2, 7.2))
    ax = fig.add_subplot(projection="3d")
    L = 0.30

    # 重新跑一遍拿到每段的轨迹（与静态图独立，避免耦合）
    n_per = 200
    traj_parts = []
    v = V0.copy()
    for gamma_fn, Le in EDGES:
        sg = np.linspace(0, Le, n_per + 1)
        Vt = C.transport_along(gamma_fn, sg, v)
        v = Vt[-1]
        traj_parts.append(Vt)

    def update(k):
        ax.cla()
        X, Y, Z = C.sphere_surface(n_theta=36, n_phi=72)
        ax.plot_surface(X, Y, Z, alpha=0.10, color="gray", linewidth=0)
        for gamma_fn, Le in EDGES:
            s = np.linspace(0, Le, 60)
            pts = np.array([gamma_fn(si)[0] for si in s])
            ax.plot(*pts.T, "k-", lw=1.6)
        # 第 k 帧落在哪一段、段内第几个点
        seg_i = min(int(k / n_frames * 3), 2)
        loc = k / n_frames * 3 - seg_i
        j = min(int(loc * n_per), n_per)
        gamma_fn, Le = EDGES[seg_i]
        pos_now = gamma_fn(loc * Le)[0]
        # 已走过的路径
        for gi in range(seg_i):
            sg = np.linspace(0, EDGES[gi][1], n_per + 1)
            pts = np.array([EDGES[gi][0](si)[0] for si in sg])
            ax.plot(*pts.T, color="#c0392b", lw=2.4)
        sg = np.linspace(0, loc * Le, max(j, 1) + 1)
        pts = np.array([gamma_fn(si)[0] for si in sg])
        ax.plot(*pts.T, color="#c0392b", lw=2.4)
        ax.scatter(*pos_now, color="#c0392b", s=35)
        # 长矛（当前平行移动向量）与初始 ghost
        ax.quiver(*pos_now, *(L * traj_parts[seg_i][j]), color="#c0392b", lw=2.6,
                  arrow_length_ratio=0.28)
        ax.quiver(*A, *(L * V0), color="green", lw=2.0, ls="dashed",
                  arrow_length_ratio=0.28, alpha=0.85)
        ax.set_xlim(-1.1, 1.1), ax.set_ylim(-1.1, 1.1), ax.set_zlim(-1.1, 1.1)
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=16, azim=38)
        ax.set_axis_off()
        ax.set_title(f"平行移动沿八分仪走一圈   进度 {k / (n_frames - 1) * 100:4.1f}%",
                     fontsize=11)
        if k == n_frames - 1:
            ax.text2D(0.5, 0.06, "向量转过了 90.0°\n= 内角和 − π = 面积（holonomy）",
                      transform=ax.transAxes, ha="center", fontsize=13,
                      color="#c0392b", weight="bold")
        return []

    anim = FuncAnimation(fig, update, frames=n_frames, interval=80)
    out = os.path.join(C.FIG_DIR, "octant_transport.gif")
    anim.save(out, writer=PillowWriter(fps=12))
    plt.close(fig)
    return out


def main():
    total, excess = check_angles()
    area = check_area()
    V0, V_end, ang, segs = check_holonomy()

    print("\n" + "=" * 70)
    print("【三重验证汇总】三个数应当相等（球面角超 = Gauss–Bonnet 实例）")
    print("=" * 70)
    # 面积换算成“角当量”只是为了同轴展示；三者相等本身是恒等式，非独立测量
    print(f"  内角和 − π    = {np.degrees(excess):10.6f}°")
    print(f"  面积（角当量）= {np.degrees(area):10.6f}°   [面积 {area:.6f} sr]")
    print(f"  holonomy      = {abs(np.degrees(ang)):10.6f}°")
    spread = max(abs(np.degrees(excess) - 90),
                 abs(np.degrees(area) - 90),
                 abs(np.degrees(ang) - 90))
    print(f"  三者最大偏离 90° 的幅度 = {spread:.2e}°  "
          f"[{'PASS' if spread < 0.1 else 'FAIL'}]")

    path_png = make_figure(segs, V0, V_end, ang, excess, area)
    path_gif = make_gif(V0)
    print("\n图：", path_png)
    print("动画：", path_gif)
    print("\n这一步在说什么：")
    print("  平面三角形：内角和 = π，面积圈住的曲率 = 0，搬一圈向量回原方向。")
    print("  球面八分仪：三个数同时变成 90° —— 角超、面积、转角是**同一个东西**")
    print("  的三种量法。曲率不再是一个数，而是「绕一圈必须付出的角度债」。")


if __name__ == "__main__":
    main()