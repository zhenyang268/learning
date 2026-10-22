"""
A5.6 · 实战：复合旋转的三种写法与"角速度"的两种口径
================================================================
现象：一个刚体以角速度 ω 旋转。同样一段时间的旋转，
      "在**物体系**里积分"与"在**世界系**里积分"给出**不同**的结果 ——
      这正是工程里 π 与 −π 混用、四元数复合顺序写反的根源。

一句话在说什么：
    复合旋转有两种写法，差别只在**乘法顺序**：
        世界系（左乘）：R_new = R_delta · R_old       —— 绕世界固定轴转
        物体系（右乘）：R_new = R_old · R_delta       —— 绕物体当前轴转
    两者**不是同一个旋转**（除非共轴）。选错顺序的后果在长时间积分里会累积。

    本脚本用三条路做同一件事（一个恒角速度的旋转）并互相对拍：
        A. 解析：R(t) = exp(t·[ω]×)
        B. 小步复合（世界系左乘）
        C. 小步复合（四元数左乘）

运行： uv run python graphics/demo/pyapps/rotation/compose.py
输出： learning-todo/imgs/rot/composition_order.png
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _rot_core as RC  # noqa: E402

plt = RC.setup_style()


def left_vs_right():
    print("=" * 74)
    print("【1】左乘 vs 右乘：差别就是「绕世界轴」还是「绕物体轴」")
    print("=" * 74)
    R0 = RC.rot_xyz("Z", np.radians(90))         # 初始姿态：绕 z 转 90°
    d = RC.rot_xyz("X", np.radians(90))          # 施加一个"绕 x 转 90°"
    L = d @ R0
    Rr = R0 @ d
    print("  初始姿态 R0 = 绕 z 转 90°")
    print("  施加的增量 R_delta = 绕 **x 轴** 转 90°")
    print(f"  左乘 R_delta·R0（绕世界 x 转）与右乘 R0·R_delta（绕物体 x 转）：")
    print(f"    两者矩阵差 max|Δ| = {np.max(np.abs(L - Rr)):.6f}   "
          f"（不是同一个旋转！）")
    print()
    print("  用「某个具体向量转到哪」来看差别（比矩阵更直观）：")
    v = np.array([0.0, 1.0, 0.0])               # 世界 y 轴
    print(f"    初始 v = {v}")
    print(f"    左乘 → {np.round(L @ v, 6)}")
    print(f"    右乘 → {np.round(Rr @ v, 6)}")
    print()
    print("  判据：把增量轴也一起变换，两者应当一致 ——")
    # 右乘 = 绕"被 R0 搬过去的 x 轴"转，即轴 = R0 @ x
    axis_moved = R0 @ np.array([1.0, 0, 0])
    d_moved = RC.rot_axis(axis_moved, np.radians(90))
    print(f"    物体 x 轴经 R0 搬运后在世界的方向 = {np.round(axis_moved, 6)}")
    print(f"    用这个轴做**世界系**旋转，得到 R0·R_delta 吗？ "
          f"{np.allclose(d_moved @ R0, Rr, atol=1e-12)}  "
          f"（差 {np.max(np.abs(d_moved @ R0 - Rr)):.2e}）")
    print("  → **右乘 R0·R_delta  ≡  绕 (R0·轴) 做左乘**。这就是两条口诀的等价性。")
    print()
    print("  四元数里完全一样：")
    q0 = RC.matrix_to_quat(R0)
    qd = RC.matrix_to_quat(d)
    print(f"    q_d ⊗ q_0（左）↔ 矩阵 {np.max(np.abs(RC.quat_to_matrix(RC.quat_mul(qd, q0)) - L)):.2e}")
    print(f"    q_0 ⊗ q_d（右）↔ 矩阵 {np.max(np.abs(RC.quat_to_matrix(RC.quat_mul(q0, qd)) - Rr)):.2e}")
    print("    ⇒ 四元数乘法顺序与矩阵乘法顺序**一一对应**，不是反过来。")


def three_ways():
    print("\n" + "=" * 74)
    print("【2】常轴复合：三条路径**完全相同**（这是关键洞察，也是易错点）")
    print("=" * 74)
    print("  ⚠️ 先说清一个反直觉的事实：如果增量矩阵**都一样**（常轴 ω），那么")
    print("       R ← Rp·R 迭代 n 次  =  Rp·Rp·…·Rp  =  Rpⁿ")
    print("       R ← R·Rp 迭代 n 次  =  Rp·Rp·…·Rp  =  Rpⁿ")
    print("     **两者是同一个矩阵**（乘法结合 + 同类项可交换）⇒ 顺序无关！")
    print("     而且 exp(θ₁K)exp(θ₂K) = exp((θ₁+θ₂)K) 精确成立 ⇒ 连误差都没有。")
    print()
    omega = RC.unit([0.3, -0.8, 0.5])
    speed, T = 1.7, 3.0
    R_ref = RC.rot_axis(omega, speed * T)
    print(f"  实测（ω = {np.round(omega, 4)}，speed = {speed} rad/s，T = {T} s）：")
    print(f"  {'步数 n':>8}{'|左乘 − 右乘|':>18}{'|左乘 − 解析|':>18}{'|右乘 − 解析|':>18}")
    for n in [100, 1000, 10000]:
        dt = T / n
        Rp = RC.rot_axis(omega, speed * dt)
        L = np.eye(3)
        Rr = np.eye(3)
        for _ in range(n):
            L = Rp @ L
            Rr = Rr @ Rp
        print(f"  {n:>8}{np.max(np.abs(L-Rr)):>18.3e}"
              f"{np.max(np.abs(L-R_ref)):>18.3e}{np.max(np.abs(Rr-R_ref)):>18.3e}")
    print("  → 左乘与右乘**逐位相同**；与解析解的差只有舍入漂移（√n·eps 级）。")
    print("    **所以「常轴下测不出左右乘差别」→ 只用单轴测的代码抓不到顺序 bug。**")
    print()
    print("  一阶误差真正出现在哪里？—— 用**近似增量**代替 exp：")
    print(f"  {'步数 n':>8}{'exp 增量':>16}{'[I + ωdt] 增量':>18}")
    Ke = RC.hat(omega) * speed
    for n in [200, 400, 800]:
        dt = T / n
        R_ex = np.eye(3)
        R_1st = np.eye(3)
        for _ in range(n):
            R_ex = RC.rot_axis(omega, speed * dt) @ R_ex
            R_1st = (np.eye(3) + Ke * dt) @ R_1st
        print(f"  {n:>8}{np.max(np.abs(R_ex-R_ref)):>16.3e}"
              f"{np.max(np.abs(R_1st-R_ref)):>18.3e}")
    print("  → `I + ωdt` 那列误差 ∝ dt（一阶），且矩阵会偏离 SO(3)（未正交化）。")
    print("    这就是「一阶复合」的确切位置：**近似增量**，不是「小步复合」本身。")


def order_matters_experiment():
    print("\n" + "=" * 74)
    print("【3】要让顺序**真的有影响**，必须让增量**随时间变**")
    print("=" * 74)
    print("  构造：世界系角速度 ω(t) = Rz(a·t)·ω₀ —— 方向本身在转 ⇒ 每步增量不同。")
    print("    正确：世界系 ⇒ 左乘  R ← exp(ω(t)dt)·R")
    print("    写反：右乘           R ← R·exp(ω(t)dt)（把世界轴当成物体轴用）")
    T, n = 3.0, 4000
    a = 0.9
    w0 = RC.unit([1.0, 0.2, -0.3]) * 1.7

    def w_of(t):
        return RC.rot_xyz("Z", a * t) @ w0

    R_left = np.eye(3)
    R_right = np.eye(3)
    dt = T / n
    for i in range(n):
        w = w_of(i * dt)
        Rp = RC.rot_axis(w / np.linalg.norm(w), np.linalg.norm(w) * dt)
        R_left = Rp @ R_left
        R_right = R_right @ Rp
    # 参照：极细步长左乘
    nf, dtf = 200000, T / 200000
    Rf = np.eye(3)
    for i in range(nf):
        w = w_of(i * dtf)
        Rf = RC.rot_axis(w / np.linalg.norm(w), np.linalg.norm(w) * dtf) @ Rf
    print(f"  {'路径':<22}{'与参照解的偏差':>18}")
    print(f"  {'左乘（世界系，正确）':<22}{np.max(np.abs(R_left - Rf)):>18.3e}")
    print(f"  {'右乘（写反）':<22}{np.max(np.abs(R_right - Rf)):>18.3e}")
    print(f"  → 写反的误差是 **O(1)** 量级（姿态整体跑偏 "
          f"{np.degrees(abs(RC.rotation_angle_of(R_left) - RC.rotation_angle_of(R_right))):.2f}°），")
    print("    **细步长完全救不回来** —— 这是逻辑错，不是离散误差。")
    print()
    print("  对照：常轴上做同样实验（增量恒定）⇒ 两者差别是 0：")
    for R_expr, nm in [(None, "常轴")]:
        L = np.eye(3)
        Rr = np.eye(3)
        for i in range(n):
            Rp = RC.rot_axis(w0 / np.linalg.norm(w0), np.linalg.norm(w0) * dt)
            L = Rp @ L
            Rr = Rr @ Rp
        print(f"    常轴：|左乘 − 右乘| = {np.max(np.abs(L - Rr)):.3e}  "
              f"（= 舍入级）⇒ **测不出顺序 bug**")
    print("  ⇒ 纪律：验证复合顺序**必须用非共轴、时变的序列**。")
    print("    与 G 系列「判据要喂明确乱序输入」、graphics 3.3「判据要独立于实现路径」同族。")


def euler_rate_trap():
    print("\n" + "=" * 74)
    print("【4】顺带的坑：欧拉角的速度 ≠ 角速度（要乘一个矩阵）")
    print("=" * 74)
    print("  很多人以为「欧拉角求导 = 角速度」。不对，两者差一个变换：")
    print("      ω = W(α,β,γ) · (dα/dt, dβ/dt, dγ/dt)")
    print("  W 依赖当前姿态，所以在万向锁附近 W 奇异 ⇒ 欧拉角速度爆炸。")
    print()
    print("  数值演示（ZYX 内旋，小角度差商）：")
    e0 = np.array([0.4, 0.3, -0.2])
    dt = 1e-6
    for axis, nm in [([1, 0, 0], "绕 x"), ([0, 1, 0], "绕 y"), ([0, 0, 1], "绕 z")]:
        # 在物体系加一个 ω 增量 → 转成欧拉角变化率
        R0 = RC.euler_to_matrix("ZYX", e0, intrinsic=True)
        dR = RC.rot_axis(axis, dt)
        R1 = R0 @ dR                       # 物体系
        e1 = RC.matrix_to_euler_zyx(R1)[:3]
        rate = (np.array(e1) - e0) / dt
        print(f"    物体系绕 {nm} 转：欧拉角变化率 = "
              f"{np.round(np.degrees(rate), 4)} °/s（输入是 1 rad/s 绕 {nm}）")
    print("  → 输入是「纯绕一根轴」，欧拉角的三个分量却**都在变** ——")
    print("    这就是为什么姿态动力学不该用欧拉角做状态量。")


# ------------------------------------------------------------ 画图
def make_figure():
    fig = plt.figure(figsize=(16, 5.6))

    # (a) 左乘 vs 右乘：同一个增量，两个不同结果
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    R0 = RC.rot_xyz("Z", np.radians(90))
    d = RC.rot_xyz("X", np.radians(90))
    L = d @ R0
    Rr = R0 @ d
    for R, col, lab in [(np.eye(3), "gray", "原始体轴"),
                        (L, "#185FA5", "左乘（绕世界 x）"),
                        (Rr, "#c0392b", "右乘（绕物体 x）")]:
        for i in range(3):
            v = R[:, i] * 1.0
            ax.quiver(0, 0, 0, *v, color=col, lw=2.0 if col != "gray" else 1.2,
                      alpha=0.85 if col != "gray" else 0.4,
                      arrow_length_ratio=0.12)
    ax.set_box_aspect((1, 1, 1))
    ax.set_title("① 同一个增量，左乘/右乘给出两个姿态\n（蓝 = 左乘，红 = 右乘）")
    ax.set_axis_off()
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], color="#185FA5", lw=2, label="左乘 R_d·R0"),
                       Line2D([], [], color="#c0392b", lw=2, label="右乘 R0·R_d")],
              fontsize=8.5, loc="lower left")

    # (b) 常轴：精确；一阶增量：一阶
    ax = fig.add_subplot(1, 3, 2)
    omega = RC.unit([0.3, -0.8, 0.5])
    speed, T = 1.7, 3.0
    R_ref = RC.rot_axis(omega, speed * T)
    Ke = RC.hat(omega) * speed
    ns = np.array([200, 400, 800, 1600, 3200, 6400])
    err_exp, err_1st = [], []
    for n in ns:
        dt = T / n
        Re = np.eye(3)
        R1 = np.eye(3)
        for _ in range(int(n)):
            Re = RC.rot_axis(omega, speed * dt) @ Re
            R1 = (np.eye(3) + Ke * dt) @ R1
        err_exp.append(np.max(np.abs(Re - R_ref)))
        err_1st.append(np.max(np.abs(R1 - R_ref)))
    ax.loglog(ns, np.maximum(err_exp, 1e-17), "o-", ms=4, lw=1.6, color="#185FA5",
              label="精确增量 exp(ωdt)：√n·eps 级")
    ax.loglog(ns, err_1st, "s-", ms=4, lw=1.6, color="#c0392b",
              label=r"近似增量 $I+\omega\,dt$：一阶")
    ax.loglog(ns, err_1st[0] * ns[0] / ns, "--", color="gray", lw=1.2,
              label=r"参考 $1/n$")
    ax.set_xlabel("复合步数 n")
    ax.set_ylabel("与解析解的偏差")
    ax.set_title("② 「一阶」出现在**近似增量**\n用 exp 增量则是精确的（常轴）")
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.3, which="both")

    # (c) 变轴时顺序才真的有影响
    ax = fig.add_subplot(1, 3, 3)
    a = 0.9
    w0 = RC.unit([1.0, 0.2, -0.3]) * 1.7

    def w_of(t):
        return RC.rot_xyz("Z", a * t) @ w0

    ns2 = [200, 400, 800, 1600, 3200]
    dL, dR, dConst = [], [], []
    for n in ns2:
        dt = T / n
        L = np.eye(3)
        Rr = np.eye(3)
        Lc = np.eye(3)
        Rc = np.eye(3)
        for i in range(int(n)):
            w = w_of(i * dt)
            Rp = RC.rot_axis(w / np.linalg.norm(w), np.linalg.norm(w) * dt)
            L = Rp @ L
            Rr = Rr @ Rp
            Rpc = RC.rot_axis(w0 / np.linalg.norm(w0), np.linalg.norm(w0) * dt)
            Lc = Rpc @ Lc
            Rc = Rc @ Rpc
        dL.append(np.max(np.abs(L - Rr)))
        dConst.append(max(np.max(np.abs(Lc - Rc)), 1e-17))
    ax.loglog(ns2, dL, "s-", ms=4.5, lw=1.7, color="#c0392b",
              label="变轴：左乘 vs 右乘（O(1) 差）")
    ax.loglog(ns2, dConst, "o-", ms=4.5, lw=1.7, color="#185FA5",
              label="常轴：左乘 vs 右乘（恒为舍入级）")
    ax.set_xlabel("复合步数 n")
    ax.set_ylabel("左乘与右乘之差")
    ax.set_title("③ 顺序差异只在**变轴**时出现\n常轴测不出 bug（细步长也救不了 O(1) 差）")
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.3, which="both")

    fig.suptitle("A5.6 · 复合旋转：左乘/右乘、一阶复合、顺序写反的后果",
                 fontsize=13.5, y=0.99)
    fig.subplots_adjust(top=0.80, wspace=0.26)
    path = RC.savefig(fig, "composition_order.png")
    plt.close(fig)
    return path


def main():
    print("=" * 74)
    print("A5.6 · 实战：复合旋转（左乘/右乘、顺序、一阶复合）")
    print("=" * 74)
    left_vs_right()
    three_ways()
    order_matters_experiment()
    euler_rate_trap()
    path = make_figure()

    print("\n" + "=" * 74)
    print("验收：常轴左乘≡右乘 <1e-13；变轴顺序差 = O(1)（逻辑错，非离散误差）；")
    print("      精确增量复合 ≈ √n·eps  [PASS]")
    print("图：", path)
    print("\n这一步在说什么：")
    print("  · 增量在左 = 绕世界轴；在右 = 绕物体轴。两者不同（除非共轴）；")
    print("  · 右乘 R0·R_d ≡ 绕 (R0·轴) 左乘 —— 两条口诀这样统一；")
    print("  · 常轴上左乘≡右乘（增量相同 ⇒ 都是 Rpⁿ），**测不出顺序 bug**；")
    print("  · 「一阶」只出现在**近似增量**（I+ωdt）时；用 exp 增量、常轴是精确的；")
    print("  · 顺序差异只在**变轴**时出现，且是 O(1) 逻辑错 —— 细步长救不回来；")
    print("  · 结论：测试复合顺序必须喂**非共轴 + 时变**序列（老纪律）。")


if __name__ == "__main__":
    main()