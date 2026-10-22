"""
旋转表示公共数学层 —— 被同目录下各演示脚本共享。

三种旋转表示（本目录的主线）
--------------------------------------------------
    轴角 / 罗德里格斯    (u, theta)  ->  R = I + sinθ [u]× + (1−cosθ) [u]×²
    四元数               q = cos(θ/2) + sin(θ/2) u    （单位四元数，双覆盖 SO(3)）
    欧拉角               (α, β, γ) 配一个轴序列 + 内旋/外旋约定

同一个旋转有这三种写法，它们之间可以互相翻译 —— 这就是本目录要展示的东西。
**万向锁**是欧拉角这条路特有的缺陷，四元数与轴角都没有。

约定（全目录统一，改动会连锁影响所有脚本）
--------------------------------------------------
右手系。绕轴 u 按右手定则转正角。
[u]× 是叉积矩阵（反对称）：
        ⎡  0   −u_z   u_y ⎤
[u]× =  ⎢ u_z    0   −u_x ⎥
        ⎣−u_y   u_x    0  ⎦
外旋 (extrinsic)：绕世界固定轴依次转 —— R = R_末(a₃)·R_中(a₂)·R_首(a₁)
内旋 (intrinsic)：绕物体自身轴依次转 —— R = R_首(a₁)·R_中(a₂)·R_末(a₃)
（两者互为转置关系：同序列同角度的内旋与外旋，结果是彼此的逆。）

so(3) -> SO(3)
--------------------------------------------------
[u]× 属于李代数 so(3)（反对称 3×3），指数映射 exp(θ[u]×) 落到 SO(3)（旋转群）。
罗德里格斯公式就是这个指数映射的闭式。**[u]×³ = −[u]×**（u 为单位向量）
是它收敛成闭式的全部秘密。

编号：本目录六个脚本按 A5 切片的 pyapps 预演记作 **A5.1–A5.6**，
与球面目录（G1–G6）并列。
"""

from __future__ import annotations

import os

import numpy as np

# ---------------------------------------------------------------- 路径与画图

_REPO = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
FIG_DIR = os.path.join(_REPO, "learning-todo", "imgs", "rot")


def setup_style():
    """配置 matplotlib：中文黑体 + 300 dpi + 负号正常显示。"""
    import matplotlib.pyplot as plt

    plt.rcParams["font.sans-serif"] = ["PingFang SC", "Heiti SC",
                                       "Hiragino Sans GB", "Arial Unicode MS",
                                       "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["figure.dpi"] = 100
    plt.rcParams["savefig.dpi"] = 300
    plt.rcParams["savefig.bbox"] = "tight"
    plt.rcParams["font.size"] = 11
    return plt


def savefig(fig, name: str):
    os.makedirs(FIG_DIR, exist_ok=True)
    path = os.path.join(FIG_DIR, name)
    fig.savefig(path, dpi=300)
    return path


# ---------------------------------------------------------------- 基本工具

def hat(u):
    """叉积矩阵 [u]×（hat 映射 so(3) 的向量化逆：R³ -> so(3)）。"""
    u = np.asarray(u, dtype=float)
    return np.array([[0.0, -u[2], u[1]],
                     [u[2], 0.0, -u[0]],
                     [-u[1], u[0], 0.0]])


def vee(K):
    """[·]× 的逆：so(3) -> R³。"""
    K = np.asarray(K, dtype=float)
    return np.array([K[2, 1], K[0, 2], K[1, 0]])


def unit(v):
    """归一化；零向量抛异常（避免静默产生 NaN）。"""
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if n < 1e-15:
        raise ValueError("零向量不能归一化")
    return v / n


def rot_axis(axis, angle):
    """绕单位轴转 angle（弧度）的旋转矩阵 —— **罗德里格斯公式**。

        R = I + sinθ [u]× + (1 − cosθ) [u]×²

    来历：exp(θ[u]×) = Σ (θ[u]×)ⁿ/n!，用 [u]×³ = −[u]× 把奇偶次分别求和，
    得到 sinθ 与 (1−cosθ) 两项。
    """
    u = unit(axis)
    K = hat(u)
    return np.eye(3) + np.sin(angle) * K + (1.0 - np.cos(angle)) * (K @ K)


def rot_expm_series(axis, angle, n_terms=40):
    """同一件事的**级数版**：直接累加 exp(θ[u]×) 的 Taylor 级数。

    与闭式对拍用（两条实现路径互不参考）。级数对大 θ 收敛慢，
    本函数只用于「两条实现给出同一个矩阵」这个判据。
    """
    K = hat(unit(axis))
    term = np.eye(3)
    acc = np.eye(3)
    for n in range(1, n_terms + 1):
        term = term @ (angle * K) / n
        acc = acc + term
    return acc


def rot_xyz(axis_letter, angle):
    """绕坐标轴 X/Y/Z 转 angle。"""
    ax = {"X": [1, 0, 0], "Y": [0, 1, 0], "Z": [0, 0, 1]}[axis_letter.upper()]
    return rot_axis(ax, angle)


def is_rotation(R, tol=1e-10):
    """判 R 是否属于 SO(3)：Rᵀ R = I 且 det = +1。返回 (正交误差, det 误差)。"""
    R = np.asarray(R, dtype=float)
    orth = np.max(np.abs(R.T @ R - np.eye(3)))
    det = abs(np.linalg.det(R) - 1.0)
    ok = orth < tol and det < tol
    return ok, orth, det


# ---------------------------------------------------------------- 四元数

def quat_from_axis_angle(axis, angle):
    """(u, θ) -> 单位四元数 q = (w, x, y, z) = (cos θ/2, sin θ/2 · u)。

    **半角是四元数的全部秘密**：转 θ 需要 θ/2。倍数 q ↔ −q 给出同一个旋转。
    """
    u = unit(axis)
    h = angle / 2.0
    return np.array([np.cos(h), *(np.sin(h) * u)])


def quat_mul(a, b):
    """四元数乘法（Hamilton 约定：ij = k, ji = −k）。

    (w₁ + v₁)(w₂ + v₂) = (w₁w₂ − v₁·v₂) + (w₁v₂ + w₂v₁ + v₁×v₂)
    **不交换** —— v₁×v₂ 那一项就是非交换性的来源。
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    w1, v1 = a[0], a[1:]
    w2, v2 = b[0], b[1:]
    w = w1 * w2 - float(np.dot(v1, v2))
    v = w1 * v2 + w2 * v1 + np.cross(v1, v2)
    return np.array([w, *v])


def quat_conj(q):
    """共轭 q* = (w, −v)。单位四元数下 q* = q⁻¹。"""
    q = np.asarray(q, dtype=float)
    return np.array([q[0], -q[1], -q[2], -q[3]])


def quat_norm(q):
    return float(np.linalg.norm(np.asarray(q, dtype=float)))


def quat_to_matrix(q):
    """单位四元数 -> 旋转矩阵。"""
    q = unit(np.asarray(q, dtype=float))
    w, x, y, z = q
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)],
    ])


def matrix_to_quat(R):
    """旋转矩阵 -> 单位四元数（Shepperd 法：选最大对角元分支，避免除以小数）。"""
    R = np.asarray(R, dtype=float)
    t = np.trace(R)
    if t > 0:
        s = np.sqrt(t + 1.0) * 2
        w = 0.25 * s
        x = (R[2, 1] - R[1, 2]) / s
        y = (R[0, 2] - R[2, 0]) / s
        z = (R[1, 0] - R[0, 1]) / s
    elif R[0, 0] > R[1, 1] and R[0, 0] > R[2, 2]:
        s = np.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2]) * 2
        w = (R[2, 1] - R[1, 2]) / s
        x = 0.25 * s
        y = (R[0, 1] + R[1, 0]) / s
        z = (R[0, 2] + R[2, 0]) / s
    elif R[1, 1] > R[2, 2]:
        s = np.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2]) * 2
        w = (R[0, 2] - R[2, 0]) / s
        x = (R[0, 1] + R[1, 0]) / s
        y = 0.25 * s
        z = (R[1, 2] + R[2, 1]) / s
    else:
        s = np.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1]) * 2
        w = (R[1, 0] - R[0, 1]) / s
        x = (R[0, 2] + R[2, 0]) / s
        y = (R[1, 2] + R[2, 1]) / s
        z = 0.25 * s
    q = np.array([w, x, y, z])
    return q / np.linalg.norm(q)


def quat_rotate(q, v):
    """用四元数旋转向量：v' = q v q*。这里用标量+向量三元组直接算，避免构造 4×4。

        v' = v + 2 w (u × v) + 2 u × (u × v)     （q = (w, u)，单位）
    """
    q = unit(np.asarray(q, dtype=float))
    w, u = q[0], q[1:]
    v = np.asarray(v, dtype=float)
    return v + 2 * w * np.cross(u, v) + 2 * np.cross(u, np.cross(u, v))


def slerp(q0, q1, t):
    """球面线性插值。t 可为标量或数组。返回 (n,4) 或 (4,)。

    要点：先取**最短弧**（q0·q1 < 0 时翻 q1 的符号 ⇒ 走 q 的负代表元），
    再用 sin 权重而不是线性权重。
    """
    q0 = unit(np.asarray(q0, dtype=float))
    q1 = np.asarray(q1, dtype=float)
    d = float(np.dot(q0, q1))
    if d < 0.0:                       # 取短弧
        q1 = -q1
        d = -d
    d = min(d, 1.0)
    t = np.atleast_1d(np.asarray(t, dtype=float))
    if d > 0.9995:                    # 几乎同向：退化成 nlerp，避免除零
        out = q0[None, :] + t[:, None] * (q1 - q0)[None, :]
        return out / np.linalg.norm(out, axis=1, keepdims=True)
    th0 = np.arccos(d)
    th = th0 * t
    s0 = np.sin(th0 - th) / np.sin(th0)
    s1 = np.sin(th) / np.sin(th0)
    return s0[:, None] * q0[None, :] + s1[:, None] * q1[None, :]


def nlerp(q0, q1, t):
    """归一化线性插值：线性插值后归一化。比 slerp 快，但不是匀速。"""
    q0 = unit(np.asarray(q0, dtype=float))
    q1 = np.asarray(q1, dtype=float)
    if float(np.dot(q0, q1)) < 0.0:
        q1 = -q1
    t = np.atleast_1d(np.asarray(t, dtype=float))
    out = q0[None, :] + t[:, None] * (q1 - q0)[None, :]
    return out / np.linalg.norm(out, axis=1, keepdims=True)


def quat_angle(q):
    """单位四元数表示的旋转角（取 [0, 2π) 内的最短表示，即 θ ∈ [0, π]）。"""
    q = unit(np.asarray(q, dtype=float))
    if q[0] < 0:                       # 取 w ≥ 0 的代表元
        q = -q
    return 2.0 * np.arccos(np.clip(q[0], -1.0, 1.0))


# ---------------------------------------------------------------- 欧拉角

# 12 种约定 = 6 个轴序列 × {内旋, 外旋}。Tait-Bryan = 三个不同轴；proper = 首尾同轴。
EULER_SEQS = ["XYZ", "XZY", "YXZ", "YZX", "ZXY", "ZYX",
              "XYX", "XZX", "YXY", "YZY", "ZXZ", "ZYZ"]


def euler_to_matrix(seq, angles, intrinsic=False):
    """欧拉角 -> 旋转矩阵。

    seq      : 轴序列，如 "ZYX" 表示按 Z, Y, X 的次序施加
    angles   : 与 seq 对应的三个角（弧度）
    intrinsic: False = 外旋（绕世界固定轴），True = 内旋（绕物体当前轴）

    **内旋 = 把外旋的乘法次序倒过来**（实测 0.0 误差）：
        外旋 R = R_末·R_中·R_首 ；内旋 R = R_首·R_中·R_末

    三条常用的等价改写（全部实测核过，误差 ≤ 2.2e-16）：
        ① 内旋(seq, a) = 外旋(seq 倒序, a 倒序)
        ② 内旋(seq, a) = 外旋(seq, −a)ᵀ
        ③ 内旋(seq, a)ᵀ = 内旋(seq 倒序, −(a 倒序))

    **注意**：内旋(seq, a) ≠ 外旋(seq, a)ᵀ（实测差 1.4），别写成那个。
    """
    seq = seq.upper()
    if len(seq) != 3 or any(c not in "XYZ" for c in seq):
        raise ValueError(f"轴序列非法: {seq}")
    angles = np.asarray(angles, dtype=float)
    mats = [rot_xyz(c, a) for c, a in zip(seq, angles)]
    if intrinsic:
        R = mats[0] @ mats[1] @ mats[2]
    else:
        R = mats[2] @ mats[1] @ mats[0]
    return R


def matrix_to_euler_zyx(R):
    """把旋转矩阵分解成 ZYX **内旋**（= yaw-pitch-roll，航空航天惯例）。

        R = Rz(yaw) · Ry(pitch) · Rx(roll)
        pitch = asin(−R[2,0])
        roll  = atan2(R[2,1], R[2,2])
        yaw   = atan2(R[1,0], R[0,0])

    |pitch| = 90° 时 R[2,1] 与 R[2,2] 同时为 0 ⇒ roll/yaw 不可分离（**万向锁**）。
    此时约定 roll = 0，把全部自由给 yaw。

    锁死时 yaw 的公式（两支**同一个**，推导见下）：
        pitch = +90°:  R[0,1] = sin(roll−yaw),   R[1,1] = cos(yaw−roll)
                       ⇒ atan2(−R[0,1], R[1,1]) = yaw − roll
        pitch = −90°:  R[0,1] = −sin(yaw+roll),  R[1,1] = cos(yaw+roll)
                       ⇒ atan2(−R[0,1], R[1,1]) = yaw + roll
    两式都是 `atan2(−R[0,1], R[1,1])`，取 roll = 0 后即 yaw 本身。
    （曾在此处写成 `atan2(R[0,1], R[1,1])` 用于负支，符号错 —— 已修。）

    返回 (yaw, pitch, roll, locked)。
    """
    R = np.asarray(R, dtype=float)
    sp = np.clip(-R[2, 0], -1.0, 1.0)
    pitch = np.arcsin(sp)
    locked = abs(abs(sp) - 1.0) < 1e-10
    if locked:
        yaw = np.arctan2(-R[0, 1], R[1, 1])     # 两支同式，见 docstring
        roll = 0.0
    else:
        roll = np.arctan2(R[2, 1], R[2, 2])
        yaw = np.arctan2(R[1, 0], R[0, 0])
    return yaw, pitch, roll, locked


# ---------------------------------------------------------------- 画图辅佐

def draw_axes(ax, length=1.0, labels=("x", "y", "z"), alpha=0.85):
    """画三条坐标轴（带箭头）与轴标签。"""
    cols = ["#c0392b", "#27ae60", "#2980b9"]
    for i in range(3):
        v = np.zeros(3)
        v[i] = length
        ax.quiver(0, 0, 0, *v, color=cols[i], lw=1.6,
                  arrow_length_ratio=0.12, alpha=alpha)
    for i, lab in enumerate(labels):
        v = np.zeros(3)
        v[i] = length * 1.12
        ax.text(*v, lab, fontsize=11, color=cols[i])
    return ax


def draw_sphere_wire(ax, R=1.0, n_theta=26, n_phi=52, alpha=0.07, color="gray"):
    """半透明球面作底。"""
    th = np.linspace(0, np.pi, n_theta)
    ph = np.linspace(0, 2 * np.pi, n_phi)
    TH, PH = np.meshgrid(th, ph)
    X = R * np.sin(TH) * np.cos(PH)
    Y = R * np.sin(TH) * np.sin(PH)
    Z = R * np.cos(TH)
    ax.plot_surface(X, Y, Z, alpha=alpha, color=color, linewidth=0)
    return ax


def orthonormal_frame(R, scale=1.0):
    """刚体：把 R 的三列当成三条体轴，返回端点（从原点出发）。"""
    R = np.asarray(R, dtype=float)
    return [scale * R[:, i] for i in range(3)]


def fmt_mat(R, digits=6):
    """矩阵转可读字符串（用于终端打印）。"""
    R = np.asarray(R, dtype=float)
    return "\n".join("  [" + "  ".join(f"{x:+.{digits}f}" for x in row) + "]"
                     for row in R)


def angle_between(a, b):
    """两向量夹角（弧度）。"""
    a, b = unit(a), unit(b)
    return float(np.arccos(np.clip(np.dot(a, b), -1.0, 1.0)))


def rotation_angle_of(R):
    """旋转矩阵的转角 θ ∈ [0, π]。

    用 **atan2 良态口径**，不用 acos((tr−1)/2)：
        对任意旋转，tr = 1 + 2cosθ，且 ‖R − Rᵀ‖/2 = |sinθ|。
        于是  θ = atan2(‖vee(R−Rᵀ)‖/2, (tr−1)/2)

    为什么不用 acos：θ → 0 与 θ → π 两端 acos 的导数都发散
    （条件数 ~1/√(1−c²)），实测在 θ=π 处会给出 179.999998°（误差 2e-6°），
    进而让"反解轴+角重建矩阵"差 2e-8 —— 那是**判据病态**，不是实现错。
    atan2 口径在两端都是良态的（实测 θ=π 处误差回到机器精度）。
    """
    R = np.asarray(R, dtype=float)
    tr = float(np.trace(R))
    K = R - R.T
    s = np.linalg.norm(vee(K)) / 2.0            # = |sinθ|
    c = (tr - 1.0) / 2.0                        # = cosθ
    return float(np.arctan2(s, c))


def rotation_axis_of(R):
    """旋转矩阵的转轴：R 的特征值 1 对应的特征向量（叉积法，更稳）。"""
    R = np.asarray(R, dtype=float)
    # 对非对称矩阵用 (R − Rᵀ) 的 vee：R − Rᵀ = 2 sinθ [u]×
    K = R - R.T
    v = vee(K)
    n = np.linalg.norm(v)
    if n < 1e-12:                      # θ ≈ 0 或 π
        # θ ≈ π 时用 R + I 的列
        A = R + np.eye(3)
        idx = int(np.argmax(np.linalg.norm(A, axis=0)))
        return unit(A[:, idx])
    return v / n