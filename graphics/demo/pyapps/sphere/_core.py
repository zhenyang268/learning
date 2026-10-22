"""
球面几何公共数学层 —— 被同目录下各演示脚本共享。

坐标约定（全目录统一，改动会连锁影响所有脚本）
--------------------------------------------------
theta = 极角 / colatitude，从北极量起，取值 [0, pi]
phi   = 经度 / longitude，取值 [0, 2pi)
R     = 球半径（默认 1 = 单位球）

笛卡尔嵌入：  x = R sin(theta) cos(phi)
              y = R sin(theta) sin(phi)
              z = R cos(theta)

度规（第一基本形式）：  ds^2 = R^2 dtheta^2 + R^2 sin^2(theta) dphi^2
即  g_theta theta = R^2,  g_phi phi = R^2 sin^2 theta,  g_theta phi = 0
**g 含 R^2，但 Christoffel 不含 R** —— 常系数缩放不改变联络（见 christoffel 注释）。

符号索引顺序（关键，弄混会得到差一个负号的结果）
--------------------------------------------------
Christoffel:  Gamma^k_{ij}，数组下标 [k][i][j]
Riemann:      R^i_{jkl} = d_k Gamma^i_{lj} - d_l Gamma^i_{kj}
                          + Gamma^i_{km} Gamma^m_{lj} - Gamma^i_{lm} Gamma^m_{kj}
              数组下标 [i][j][k][l]
Ricci:        R_{ij} = R^k_{ikj}
高斯曲率:      K = R_{theta phi theta phi} / det(g)   （二维专用，等价于里奇标量 / 2）

为什么需要两套平行移动实现
--------------------------------------------------
球坐标路线在 theta = 0 / pi（两极）处 cot(theta) 发散 —— 这是**坐标奇点**。
笛卡尔嵌入路线在极点照样成立 —— 因为极点处**几何**完全正常。
两套并存既是对拍手段，也是"坐标奇点 != 几何奇点"这个教学点的载体。

编号：本目录六个脚本按 fan-plan.md §3.1 的 **G1–G5** 编号，不再使用 T 系列。
"""

from __future__ import annotations

import os

import numpy as np

# ---------------------------------------------------------------- 路径与画图

# 图片统一落在 learning-todo 子仓（独立 git 仓，与 C++ 工程解耦）。
# 从本文件位置反推，避免写死绝对路径：
#   <repo>/graphics/demo/pyapps/sphere/_core.py
#   parents[0]=sphere [1]=pyapps [2]=demo [3]=graphics [4]=<repo>
_REPO = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
FIG_DIR = os.path.join(_REPO, "learning-todo", "imgs", "geo")


def setup_style():
    """配置 matplotlib：中文黑体 + 300 dpi + 负号正常显示。"""
    import matplotlib.pyplot as plt

    # PingFang SC 字重最全（含 600），排第一位可少触发 "找不到 bold" 警告
    plt.rcParams["font.sans-serif"] = ["PingFang SC", "Heiti SC",
                                       "Hiragino Sans GB", "Arial Unicode MS",
                                       "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False   # 否则负号渲染成方块
    plt.rcParams["figure.dpi"] = 100             # 屏幕预览
    plt.rcParams["savefig.dpi"] = 300            # 落盘 300 dpi
    plt.rcParams["savefig.bbox"] = "tight"
    plt.rcParams["font.size"] = 11
    return plt


def savefig(fig, name: str):
    """存图到 learning-todo/imgs/geo/<name>，回显绝对路径。"""
    os.makedirs(FIG_DIR, exist_ok=True)
    path = os.path.join(FIG_DIR, name)
    fig.savefig(path, dpi=300)
    return path


def draw_sphere_base(ax, R: float = 1.0, alpha: float = 0.09,
                     n_theta: int = 60, n_phi: int = 120, color: str = "gray"):
    """在 3D 轴上铺一层半透明球面作底（各脚本共用，免去重复五行）。"""
    X, Y, Z = sphere_surface(n_theta, n_phi, R)
    ax.plot_surface(X, Y, Z, alpha=alpha, color=color, linewidth=0)
    return ax


# ---------------------------------------------------------------- 坐标转换

def sph_to_cart(theta, phi, R: float = 1.0):
    """(theta, phi) -> 笛卡尔 (x, y, z)。逐元素运算，支持数组。"""
    return (R * np.sin(theta) * np.cos(phi),
            R * np.sin(theta) * np.sin(phi),
            R * np.cos(theta))


def cart_to_sph(p):
    """笛卡尔向量 -> (theta, phi)。p 可以是 (3,) 或 (N,3)。

    phi 返回 [0, 2pi)。**注意**：跨 phi = 0 的点列会得到 0/2pi 跳变，
    画平面展开图前必须用 unwrap_lon_deg 解卷（否则横跨整图的假跳线）。
    """
    p = np.atleast_2d(np.asarray(p, dtype=float))
    x, y, z = p[:, 0], p[:, 1], p[:, 2]
    r = np.linalg.norm(p, axis=1)
    theta = np.arccos(np.clip(z / np.where(r > 0, r, 1.0), -1.0, 1.0))
    phi = np.mod(np.arctan2(y, x), 2 * np.pi)
    if theta.size == 1:
        return float(theta[0]), float(phi[0])
    return theta, phi


def unwrap_lon_deg(lon_deg):
    """把 [0,360) 或折到 (-180,180] 的经度序列解卷成连续序列（度）。

    绘图专用。跨 180 度经线的路径若不解卷，会在图上横拉一条假跳线
    （北京->纽约的大圆就跨 180，实测跳变 358 度）。
    """
    lon = np.asarray(lon_deg, dtype=float)
    return np.degrees(np.unwrap(np.radians(lon)))


def sph_basis(theta, phi, R: float = 1.0):
    """坐标基向量 d_theta, d_phi 在笛卡尔系下的分量 (各为长度 3)。

    d_theta = R(cos th cos ph,  cos th sin ph,  -sin th)
    d_phi   = R sin th (-sin ph, cos ph, 0)
    注意 |d_phi| = R sin theta —— 在极点退化为 0，这就是坐标奇点的几何来源。
    """
    st, ct = np.sin(theta), np.cos(theta)
    sp, cp = np.sin(phi), np.cos(phi)
    e_theta = R * np.array([ct * cp, ct * sp, -st])
    e_phi = R * st * np.array([-sp, cp, 0.0])
    return e_theta, e_phi


def sphere_surface(n_theta: int = 60, n_phi: int = 120, R: float = 1.0):
    """生成球面网格 (X, Y, Z)，供 3D 画图做半透明底。"""
    th = np.linspace(0, np.pi, n_theta)
    ph = np.linspace(0, 2 * np.pi, n_phi)
    TH, PH = np.meshgrid(th, ph)
    X, Y, Z = sph_to_cart(TH, PH, R)
    return X, Y, Z


# ---------------------------------------------------------------- 度规与联络

def metric(theta, R: float = 1.0):
    """度规张量 g_ij (2x2)。ds^2 = g_ij dx^i dx^j。"""
    st = np.sin(theta)
    return np.array([[R * R, 0.0], [0.0, R * R * st * st]])


def metric_inv(theta, R: float = 1.0):
    """逆度规 g^{ij}（g 的矩阵逆）。用来升降指标。"""
    st = np.sin(theta)
    return np.array([[1.0 / (R * R), 0.0], [0.0, 1.0 / (R * R * st * st)]])


def christoffel(theta):
    """手写 Christoffel 符号 Gamma^k_{ij}，返回 (2,2,2)，下标 [k][i][j]。

    **不含 R**：Gamma 在度规整体乘常数时不变（R^2 缩放在
    (1/2) g^{kl}(...) 里上下相消）—— 所以球半径多大，联络都一样。
    单位球只有三个非零分量：
        Gamma^th_{ph ph} = -sin th cos th
        Gamma^ph_{th ph} = Gamma^ph_{ph th} = cot th
    几何直觉：Gamma 度量"坐标基向量沿着方向 i 移动时自己怎么转"。
    经度圈在赤道附近张得最开、在两极收拢，所以 cot(th) 在极点爆炸。
    """
    st, ct = np.sin(theta), np.cos(theta)
    G = np.zeros((2, 2, 2))
    G[0, 1, 1] = -st * ct                # Gamma^theta_{phi phi}
    G[1, 0, 1] = ct / st                 # Gamma^phi_{theta phi}
    G[1, 1, 0] = G[1, 0, 1]              # 无挠 => 下指标对称
    return G


def dchristoffel(theta):
    """d/dtheta 作用在 Christoffel 上（解析式）。d/dphi 恒为 0（绕轴对称）。同样不含 R。

        d_th Gamma^th_{ph ph} = d_th(-sin th cos th) = sin^2 th - cos^2 th
        d_th Gamma^ph_{th ph} = d_th(cot th)         = -1/sin^2 th
    用解析导数算曲率，精度能到机器级；有限差分版只用来做对拍。
    """
    st, ct = np.sin(theta), np.cos(theta)
    dG = np.zeros((2, 2, 2))
    dG[0, 1, 1] = st * st - ct * ct
    dG[1, 0, 1] = -1.0 / (st * st)
    dG[1, 1, 0] = dG[1, 0, 1]
    return dG


def christoffel_from_metric(x, g_func, h: float = 1e-5):
    """从度规数值反解 Christoffel（有限差分），用来验证手写的那个没写错。

    公式：Gamma^k_{ij} = (1/2) g^{kl} (d_i g_{jl} + d_j g_{il} - d_l g_{ij})
    只处理二维 (theta, phi) 且 g 与 phi 无关（球面绕轴对称）的情形 ⇒ d_phi g = 0，
    所以对 l = phi 那一支算出来恒为 0，是有意为之，不是漏算。
    """
    g = g_func(x[0])
    ginv = np.linalg.inv(g)
    # dg[l][i][j] = d_l g_{ij}；l=0 用中心差分，l=1 恒 0（g 不依赖 phi）
    dg = np.zeros((2, 2, 2))
    dg[0] = (g_func(x[0] + h) - g_func(x[0] - h)) / (2 * h)
    G = np.zeros((2, 2, 2))
    for k in range(2):
        for i in range(2):
            for j in range(2):
                acc = 0.0
                for l in range(2):
                    acc += ginv[k, l] * (dg[i][l, j] + dg[j][i, l] - dg[l][i, j])
                G[k, i, j] = 0.5 * acc
    return G


# ---------------------------------------------------------------- 曲率

def riemann(theta, h=None):
    """黎曼曲率张量 R^i_{jkl}，返回 (2,2,2,2)，下标 [i][j][k][l]。

    约定：R^i_{jkl} = d_k Gamma^i_{lj} - d_l Gamma^i_{kj}
                      + Gamma^i_{km} Gamma^m_{lj} - Gamma^i_{lm} Gamma^m_{kj}
    这个约定下球面给出 R^theta_{phi theta phi} = +sin^2 theta（单位球）。

    h=None  -> 用解析导数 dchristoffel（精度到机器级，默认）
    h=1e-4  -> 用中心差分（对拍用，误差 ~1e-8）
    """
    G = christoffel(theta)                # Gamma 与 R 无关（度规整体缩放不影响联络）

    # dG[k][i][a][b] = d_k Gamma^i_{ab}；球面绕轴对称 => 只有 k=theta(=0) 这一份非零
    dG = np.zeros((2, 2, 2, 2))
    if h is None:
        dG[0] = dchristoffel(theta)
    else:
        dG[0] = (christoffel(theta + h) - christoffel(theta - h)) / (2 * h)

    Rm = np.zeros((2, 2, 2, 2))
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    # R^i_{jkl} = d_k G^i_{lj} - d_l G^i_{kj}
                    #             + G^i_{km} G^m_{lj} - G^i_{lm} G^m_{kj}
                    val = dG[k][i][l][j] - dG[l][i][k][j]
                    for m in range(2):
                        val += G[i][k][m] * G[m][l][j]
                        val -= G[i][l][m] * G[m][k][j]
                    Rm[i, j, k, l] = val
    return Rm


def ricci(theta):
    """里奇张量 R_{ij} = R^k_{ikj}（对第一个上指标与第三个下指标缩并）。"""
    Rm = riemann(theta)
    out = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            out[i, j] = sum(Rm[k, i, k, j] for k in range(2))
    return out


def scalar_curvature(theta, R: float = 1.0):
    """里奇标量 R = g^{ij} R_{ij}。二维曲面上恒等于 2K = 2/R^2。"""
    ginv = metric_inv(theta, R)
    return float(np.sum(ginv * ricci(theta)))


def gaussian_curvature(theta, R: float = 1.0):
    """高斯曲率 K，走两条独立路径取平均并同时给出两者之差（自检）。

    路径 A：K = R_{th ph th ph} / det(g)              （全下标黎曼张量）
    路径 B：K = (里奇标量) / 2                          （二维恒等式 R = 2K）
    两者应当一致；若不一致说明 Riemann 的索引约定或 Christoffel 有一处写错。
    """
    Rm = riemann(theta)
    g = metric(theta, R)
    R_thphthph = np.einsum("i,i->", g[0], Rm[:, 1, 0, 1])   # g_{th i} R^i_{ph th ph}
    detg = np.linalg.det(g)
    K_A = R_thphthph / detg
    K_B = scalar_curvature(theta, R) / 2.0
    return K_A, K_B


# ---------------------------------------------------------------- 积分器

def rk4_step(f, s, y, h):
    """经典四阶龙格-库塔单步。y 与 f 的返回值同为 numpy 数组。"""
    k1 = np.asarray(f(s, y), dtype=float)
    k2 = np.asarray(f(s + h / 2, y + h / 2 * k1), dtype=float)
    k3 = np.asarray(f(s + h / 2, y + h / 2 * k2), dtype=float)
    k4 = np.asarray(f(s + h, y + h * k3), dtype=float)
    return y + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)


def integrate(f, y0, s_grid):
    """沿给定弧长网格 s_grid 用 RK4 积分，返回 (N, dim) 轨迹。"""
    y0 = np.asarray(y0, dtype=float)
    out = np.empty((len(s_grid), y0.size))
    out[0] = y0
    y = y0.copy()
    for idx in range(len(s_grid) - 1):
        h = s_grid[idx + 1] - s_grid[idx]
        y = rk4_step(f, s_grid[idx], y, h)
        out[idx + 1] = y
    return out


# ---------------------------------------------------------------- 大圆

def great_circle(u, n_hat, s):
    """大圆：gamma(s) = cos(s) u + sin(s) (n x u)，单位速度、弧长参数化。

    u     : 起点（单位向量）
    n_hat : 大圆所在平面的单位法向，n_hat = normalize(u x v)（v 为终点）
    s     : 弧长（单位球上 = 转过的圆心角）

    走的是 n_hat 决定的那一侧弧；要取短弧请用 great_circle_path（会自动择向）。
    千万别用"经纬度线性插值"代替 —— 那不是大圆（除非在赤道或同一经线上）。
    """
    n_hat = np.asarray(n_hat, dtype=float)
    n_hat = n_hat / np.linalg.norm(n_hat)
    u = np.asarray(u, dtype=float)
    t = np.cross(n_hat, u)               # 与 u 正交的切向，指向前进方向
    s = np.atleast_1d(np.asarray(s, dtype=float))
    return np.outer(np.cos(s), u) + np.outer(np.sin(s), t)


def great_circle_path(p1, p2, n: int = 500):
    """两点的**短弧**大圆点列（单位球），并返回中心角（弧度）。

    与 great_circle 的区别：这里自动挑方向保证取短弧 ——
    先把两点归一化，再比对 c.great_circle(u, +n_hat, L) 与 c.great_circle(u, -n_hat, L)
    哪个终点离 p2 更近。**必须用归一化后的点比对**：若拿单位球上的路径终点
    去比未归一化的 p2，距离恒 ≈ |p2|，判据永远为真，翻转分支形同虚设。
    """
    u = np.asarray(p1, dtype=float)
    u = u / np.linalg.norm(u)
    v = np.asarray(p2, dtype=float)
    v = v / np.linalg.norm(v)
    total = arc_length(u, v)
    s = np.linspace(0, total, n)
    n_hat = np.cross(u, v)
    n_hat = n_hat / np.linalg.norm(n_hat)
    path = great_circle(u, n_hat, s)
    if np.linalg.norm(path[-1] - v) > 1e-9:       # 走了长弧 ⇒ 翻向
        path = great_circle(u, -n_hat, s)
    return path, total


def arc_length(a, b):
    """单位球上两点的大圆距离（中心角）：arccos(a·b)，数值上做了裁剪。"""
    a = np.asarray(a, dtype=float) / np.linalg.norm(a)
    b = np.asarray(b, dtype=float) / np.linalg.norm(b)
    return float(np.arccos(np.clip(np.dot(a, b), -1.0, 1.0)))


def haversine(lat1, lon1, lat2, lon2, R: float = 6371.0):
    """haversine 大圆距离（km）。lat/lon 单位：度。R 默认地球平均半径。

    比 arccos 版在两点很近时数值稳定得多（arccos 在 1 附近丢有效位）。
    """
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp = p2 - p1
    dl = np.radians(lon2 - lon1)
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))


# ---------------------------------------------------------------- 平行移动

def transport_along(gamma_fn, s_grid, V0_cart):
    """沿**解析**曲线 gamma_fn(s) -> (位置, 单位切向) 平行移动，RK4。

    方程 dV/ds = -(V·gdot) gamma 的来历：
        球面法向 n_hat = gamma（单位球），nabla_gdot V = dV/ds 的切向投影；
        V·gamma = 0 求导得 dV/ds·gamma = -V·gdot，代入即得。
    几何读法：**向量唯一被迫做的改变，就是跟着曲面转的那部分。**

    RK4 的四个 stage 通过 gamma_fn(s0+s) 取到段内真实的中点曲线，
    所以整体四阶；s_grid 步长减半误差 /16。
    （教训：点列 + 两点间大圆插值 = 内接测地多边形，曲线本身被换掉了，
      误差 O(a^2) 锁死在 1e-6，步长再小也没用。）
    """
    s_grid = np.asarray(s_grid, dtype=float)
    N = len(s_grid)
    V = np.empty((N, 3))
    v = np.asarray(V0_cart, dtype=float)
    v = v / np.linalg.norm(v)
    V[0] = v
    for i in range(N - 1):
        s0 = s_grid[i]
        h = s_grid[i + 1] - s0

        def f(s, vv):
            g, gd = gamma_fn(s0 + s)           # gamma、gdot 随 stage 连续更新
            return -(vv @ gd) * g

        v = rk4_step(f, 0.0, v, h)
        g, _ = gamma_fn(s_grid[i + 1])
        v = v - (v @ g) * g                    # 投影回切平面（~1e-16 舍入漂移）
        v = v / np.linalg.norm(v)
        V[i + 1] = v
    return V


# ---------------------------------------------------------------- 报告辅佐

def angle_between(a, b, signed_about=None, eps: float = 1e-12):
    """两向量夹角。

    signed_about 给出法向时，返回带符号角（逆时针为正）；否则返回 [0, pi] 无符号角。
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    ca = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + eps)
    ca = np.clip(ca, -1.0, 1.0)
    ang = np.arccos(ca)
    if signed_about is not None:
        n = np.asarray(signed_about, dtype=float)
        if np.dot(np.cross(a, b), n) < 0:
            ang = -ang
    return float(ang)


def rotate_about(axis, angle):
    """罗德里格斯公式：绕单位轴 axis 转 angle 的 3x3 旋转矩阵。

    这是 lie 线（so(3) -> SO(3) 的指数映射）在本目录的唯一消费点：
    把"绕轴转 theta"和 expm(theta [u]_x) 当成同一件事。
    G2 用它把解析大圆从"cos/sin 线性组合"升级成"真的用旋转矩阵生成"，
    从而让"解微分方程 vs SO(3) 代数"这句两路对拍名副其实。
    """
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])     # [axis]_x
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)