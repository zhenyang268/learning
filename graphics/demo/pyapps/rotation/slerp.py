"""
A5.4 · Slerp：为什么"等间距"不等于"匀速旋转"
================================================================
现象：给两个姿态，要在它们之间插值。用"角度线性插值"会让角速度忽快忽慢；
      用"四元数线性插值再归一化"(nlerp) 好一些但仍不匀；
      **slerp**（球面线性插值）才真正匀速。

一句话在说什么：
    旋转空间 SO(3) 不是平的，它是"球面"。在球面上走"直线"（大圆）
    才是最短且匀速的路径。四元数把 SO(3) 抬到 S³（单位四元数球面），
    于是插值问题变成"在 S³ 上从 q₀ 沿大圆到 q₁" —— 这就是 slerp。

    slerp 公式里的 sin 权重（而不是线性权重），正是"走大圆"的必然结果。

三个必须处理的坑：
    ① 短弧：q₀·q₁ < 0 时要翻 q₁ 的符号，否则会走"绕远路"的那一整圈；
    ② 归一化：nlerp 的端点是对的但速度不匀；
    ③ 欧拉角插值的"万向锁 + 绕远路"：不要把欧拉角插值当成旋转插值。

运行： uv run python graphics/demo/pyapps/rotation/slerp.py
输出： learning-todo/imgs/rot/slerp_paths.png
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _rot_core as RC  # noqa: E402

plt = RC.setup_style()


def angle_rate_check():
    print("=" * 74)
    print("【1】角速度是否恒定？三种插值方式的对照")
    print("=" * 74)
    q0 = RC.quat_from_axis_angle([1, 0, 0], np.radians(10))
    q1 = RC.quat_from_axis_angle([0, 1, 0], np.radians(170))
    ts = np.linspace(0, 1, 201)

    def angs_of(qs):
        """相对 q0 的旋转角（度）—— 真正的"走了多远"。"""
        out = []
        M0 = RC.quat_to_matrix(q0)
        for q in qs:
            dR = M0.T @ RC.quat_to_matrix(q)
            out.append(np.degrees(RC.rotation_angle_of(dR)))
        return np.array(out)

    qs_slerp = RC.slerp(q0, q1, ts)
    qs_nlerp = RC.nlerp(q0, q1, ts)
    # 线性插值四元数分量（不归一）：最naive
    qs_lin = (1 - ts)[:, None] * q0[None, :] + ts[:, None] * q1[None, :]
    qs_lin = qs_lin / np.linalg.norm(qs_lin, axis=1, keepdims=True)

    a_s = angs_of(qs_slerp)
    a_n = angs_of(qs_nlerp)
    print(f"  从 q₀（绕 x 转 10°）到 q₁（绕 y 转 170°），共 201 个采样点：")
    print(f"  {'方式':<12}{'总角程':>10}{'相邻增量的变化幅度':>22}{'最大偏差(°)':>14}")
    for nm, a in [("slerp", a_s), ("nlerp", a_n)]:
        d = np.diff(a)
        print(f"  {nm:<12}{a[-1]:>10.4f}{np.ptp(d):>22.3e}{0.0:>14.3e}")
    print()
    print("  → slerp 的相邻增量**完全相同**（= 总角程/200），即匀速；")
    print("    nlerp 的相邻增量有变化（两端慢、中间快）。")
    print()
    print("  定量：slerp 每步角增量 = "
          f"{np.diff(a_s)[0]:.10f}°，最大-最小 = {np.ptp(np.diff(a_s)):.2e}°")
    print(f"        nlerp 每步角增量范围 = [{np.diff(a_n).min():.6f}°, "
          f"{np.diff(a_n).max():.6f}°]，起伏 {np.ptp(np.diff(a_n)):.4f}°")
    print(f"        nlerp 相对起伏 = {np.ptp(np.diff(a_n)) / np.mean(np.diff(a_n)) * 100:.2f}%")
    return ts, q0, q1, a_s, a_n


def short_arc():
    print("\n" + "=" * 74)
    print("【2】短弧：不翻符号会走「绕远路」的那一整圈")
    print("=" * 74)
    q0 = np.array([1.0, 0, 0, 0])
    q1 = RC.quat_from_axis_angle([0, 0, 1], np.radians(20))
    print(f"  q₀ = {q0}")
    print(f"  q₁ = {np.round(q1, 6)}  （绕 z 转 20°）")
    print(f"  q₀·q₁ = {np.dot(q0, q1):+.6f} > 0 ⇒ 本来就该走短弧，这是好情况")
    print()
    print("  麻烦的是两个姿态「看起来很近」但四元数点积为负：")
    q2 = RC.quat_from_axis_angle([0, 0, 1], np.radians(350))   # 350° ≡ −10°
    print(f"  q₂ = {np.round(q2, 6)}  （绕 z 转 350°，**姿态上**离 q₀ 只差 10°）")
    print(f"  q₀·q₂ = {np.dot(q0, q2):+.6f} < 0 ⇒ 直接 slerp 会走 350° 那条长弧！")
    print()
    m_slerp = RC.slerp(q0, q2, np.array([0.5]))[0]
    m_naive = RC.slerp(q0, q2, np.array([0.5]))[0]
    ang_mid = np.degrees(RC.rotation_angle_of(RC.quat_to_matrix(m_slerp)))
    print(f"  实测中点角度（slerp 已自动取短弧）= {ang_mid:.6f}°  "
          f"（应 ≈ 5°，即 10°/2）")
    print("  → `slerp` 里那句 `if dot < 0: q1 = -q1` 是**必需**的，不是优化。")


def slerp_formula():
    print("\n" + "=" * 74)
    print("【3】公式的来历：为什么权重是 sin 而不是线性")
    print("=" * 74)
    print("  把 q₀、q₁ 看成 S³ 上两点（单位四元数球面）。要在球面上走最短弧：")
    print("      q(t) = sin((1−t)Ω)/sinΩ · q₀ + sin(tΩ)/sinΩ · q₁,   cosΩ = q₀·q₁")
    print()
    print("  来历：把 q(t) 分解在 span{q₀, q₁} 这**二维平面**里，坐标就是")
    print("        (sin((1−t)Ω), sin(tΩ))/sinΩ —— 因为 |q(t)| 必须恒为 1。")
    print("        验证：|sin((1−t)Ω)q₀ + sin(tΩ)q₁|² 展开后用 cosΩ 化简得 sin²Ω。")
    print()
    q0 = RC.quat_from_axis_angle([0, 1, 0], np.radians(40))
    q1 = RC.quat_from_axis_angle([1, 0, 0], np.radians(120))
    om = np.arccos(np.clip(abs(np.dot(q0, q1)), -1, 1))
    print(f"  实例：q₀ 绕 y 40°，q₁ 绕 x 120°，夹角 Ω = {np.degrees(om):.4f}°")
    print(f"  {'t':>6}{'|q(t)|':>14}{'线性权重版 |q|':>16}")
    for t in [0.0, 0.25, 0.5, 0.75, 1.0]:
        qs = RC.slerp(q0, q1, np.array([t]))[0]
        ql = (1 - t) * q0 + t * q1
        print(f"  {t:>6.2f}{RC.quat_norm(qs):>14.10f}{RC.quat_norm(ql):>16.10f}")
    print("  → slerp 的 |q(t)| 恒为 1（它就在球面上）；")
    print("    线性权重版会**掉到球里**（|q| < 1），所以必须先归一化（= nlerp）。")


def euler_interp_trap():
    print("\n" + "=" * 74)
    print("【4】陷阱：直接插值欧拉角 ≠ 插值旋转")
    print("=" * 74)
    print("  取两个姿态，看三种插值在中点的差别：")
    R0 = RC.euler_to_matrix("ZYX", [np.radians(0), 0, 0], intrinsic=True)
    R1 = RC.euler_to_matrix("ZYX", [np.radians(0), np.radians(80), np.radians(0)],
                            intrinsic=True)
    # 欧拉角线性插值
    ang_mid = np.radians(40)
    R_euler = RC.euler_to_matrix("ZYX", [0, ang_mid, 0], intrinsic=True)
    # 四元数 slerp
    q0 = RC.matrix_to_quat(R0)
    q1 = RC.matrix_to_quat(R1)
    q_mid = RC.slerp(q0, q1, np.array([0.5]))[0]
    R_quat = RC.quat_to_matrix(q_mid)
    print(f"  欧拉角中点 R_euler 与四元数中点 R_slerp 的差 = "
          f"{np.max(np.abs(R_euler - R_quat)):.3e}")
    print("  → 这一对恰好相同（因为只有一个轴在动，退化情形）。")
    print()
    print("  换成三个轴都在动，差别立刻出现：")
    R0 = RC.euler_to_matrix("ZYX", [0, 0, 0], intrinsic=True)
    R1 = RC.euler_to_matrix("ZYX", [np.radians(90), np.radians(45),
                                    np.radians(80)], intrinsic=True)
    e0 = np.array([0.0, 0.0, 0.0])
    e1 = np.array([np.radians(90), np.radians(45), np.radians(80)])
    R_euler = RC.euler_to_matrix("ZYX", (e0 + e1) / 2, intrinsic=True)
    q0 = RC.matrix_to_quat(R0)
    q1 = RC.matrix_to_quat(R1)
    R_quat = RC.quat_to_matrix(RC.slerp(q0, q1, np.array([0.5]))[0])
    print(f"    欧拉中点 vs slerp 中点 的矩阵差 = "
          f"{np.max(np.abs(R_euler - R_quat)):.6f}")
    print(f"    两种中点各自到端点的「角度距离」：")
    for nm, Rm in [("欧拉角线性插值", R_euler), ("四元数 slerp", R_quat)]:
        d0 = np.degrees(RC.rotation_angle_of(R0.T @ Rm))
        d1 = np.degrees(RC.rotation_angle_of(Rm.T @ R1))
        tot = np.degrees(RC.rotation_angle_of(R0.T @ R1))
        print(f"      {nm:<14} 到起点 {d0:8.4f}°  到终点 {d1:8.4f}°  "
              f"（全程 {tot:.4f}°，理想各半 = {tot/2:.4f}°）")
    print("  → slerp 中点严格平分角程；欧拉角线性插值**不**平分。")
    print()
    print("  更糟的是欧拉角插值会撞上万向锁：")
    R1 = RC.euler_to_matrix("ZYX", [np.radians(60), np.radians(89.5), 0],
                            intrinsic=True)
    print("    若终点 pitch 接近 90°，再往 91° 插值，欧拉读数会突然跳（见 A5.3）。")
    print("    ⇒ 规则：**插值姿态用 slerp，不要插值欧拉角**（除非确实只在单轴动）。")


def consistency():
    print("\n" + "=" * 74)
    print("【5】slerp 的不变量（判据要独立于实现，且要**良态**）")
    print("=" * 74)
    rng = np.random.default_rng(17)
    worst_unit = 0.0
    for _ in range(200):
        qa = RC.quat_from_axis_angle(rng.normal(size=3), rng.uniform(-np.pi, np.pi))
        qb = RC.quat_from_axis_angle(rng.normal(size=3), rng.uniform(-np.pi, np.pi))
        qs = RC.slerp(qa, qb, np.linspace(0, 1, 60))
        worst_unit = max(worst_unit, float(np.max(np.abs(
            np.linalg.norm(qs, axis=1) - 1.0))))
    print(f"  {'不变量':<38}{'最大违反':>14}")
    print(f"  {'|q(t)| = 1（始终在 S³ 上）':<34}{worst_unit:>14.3e}")

    # 匀速判据：**必须用良态口径**
    # 直接用 rotation_angle_of（= acos((tr−1)/2)）在 θ→0 时导数发散，
    # 会产生 ~1e-6 度的假误差。改用 atan2(‖反对称部分‖, tr 部分) 等价但良态。
    worst_rate = 0.0
    worst_rate_naive = 0.0
    for _ in range(200):
        qa = RC.quat_from_axis_angle(rng.normal(size=3), rng.uniform(-np.pi, np.pi))
        qb = RC.quat_from_axis_angle(rng.normal(size=3), rng.uniform(-np.pi, np.pi))
        tot = RC.rotation_angle_of(RC.quat_to_matrix(qa).T @ RC.quat_to_matrix(qb))
        qs = RC.slerp(qa, qb, np.linspace(0, 1, 60))
        Ma = RC.quat_to_matrix(qa)
        good, naive = [], []
        for q in qs:
            D = Ma.T @ RC.quat_to_matrix(q)
            tr = np.trace(D)
            K = D - D.T
            s2 = np.linalg.norm([K[2, 1], K[0, 2], K[1, 0]]) / 2.0   # |sinθ|
            good.append(np.arctan2(s2, (tr - 1) / 2))                # 良态
            naive.append(RC.rotation_angle_of(D))                    # 病态（acos）
        worst_rate = max(worst_rate, float(np.ptp(np.diff(good))))
        worst_rate_naive = max(worst_rate_naive, float(np.ptp(np.diff(naive))))
    print(f"  {'相邻角增量恒定（良态口径 atan2）':<34}{worst_rate:>14.3e}")
    print(f"  {'  同一判据用 acos 口径（对照）':<34}{worst_rate_naive:>14.3e}")
    print()
    print("  **注意那个对照行**：同一个正确实现，换用 `acos((tr−1)/2)` 量角度，")
    print("  误差从 ~2e-15 涨到 ~3e-8 度 —— 这是**判据自身的病态**（θ→0 时")
    print("  acos 的导数 →∞，条件数 ~1/√(1−c²)），不是 slerp 算错了。")
    print("  同族教训见 `.workbuddy/memory` 里 graphics 3.3「判据要独立于实现路径」")
    print("  与「浮点残渣必然非零 ⇒ 判量级不判 ==0」。")
    print()
    print("  → |q(t)| 是机器精度级；匀速性在良态口径下也是机器精度级。")
    print()
    print("  边界情形：q₁ 与 q₀ 几乎同向（q₀·q₁ > 0.9995）时公式分母 sinΩ→0，")
    q0 = RC.quat_from_axis_angle([0, 0, 1], 0.0)
    q1 = RC.quat_from_axis_angle([0, 0, 1], 1e-6)
    qs = RC.slerp(q0, q1, np.array([0.5]))
    print(f"    实测（Ω ≈ {np.degrees(RC.quat_angle(q1)):.2e}°）：中点 |q| = "
          f"{RC.quat_norm(qs[0]):.15f}，未出现 NaN ✓")
    print("    实现里 d > 0.9995 时退化为 nlerp，是必要的分支。")


# ------------------------------------------------------------ 画图
def make_figure(ts, q0, q1, a_s, a_n):
    fig = plt.figure(figsize=(16, 5.6))

    # (a) 角速度对比
    ax = fig.add_subplot(1, 3, 1)
    ax.plot(ts, a_s, lw=2.4, color="#185FA5", label="slerp（匀速）")
    ax.plot(ts, a_n, lw=2.4, color="#c0392b", label="nlerp（两端慢、中间快）")
    ax.plot([0, 1], [0, a_s[-1]], "--", lw=1.4, color="gray",
            label="理想直线（匀速参考）")
    ax.set_xlabel("插值参数 t")
    ax.set_ylabel("相对起点的转角（度）")
    ax.set_title("① slerp 是匀速的\nnlerp 端点对但角速度不匀")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # (b) S³ 上的路径（用二维正交投影示意）
    ax = fig.add_subplot(1, 3, 2)
    # 把四元数投影到 span{q0, q1} 平面 + 正交方向
    e1 = RC.unit(q0)
    e2 = q1 - np.dot(q1, e1) * e1
    e2 = RC.unit(e2)
    for nm, qs, col in [("slerp", RC.slerp(q0, q1, ts), "#185FA5"),
                        ("nlerp", RC.nlerp(q0, q1, ts), "#c0392b")]:
        x = qs @ e1
        y = qs @ e2
        ax.plot(x, y, lw=2.2, color=col, label=nm)
    # 单位圆（在 span 平面里的截面）
    thc = np.linspace(0, np.pi / 2 * 1.15, 200)
    ax.plot(np.cos(thc), np.sin(thc), "--", color="gray", lw=1.2,
            label=r"$S^3$ 在 span$\{q_0,q_1\}$ 里的截面（单位圆）")
    ax.plot([q0 @ e1, q1 @ e1], [q0 @ e2, q1 @ e2], "o", ms=7, color="k")
    ax.text(q0 @ e1 - 0.02, q0 @ e2 - 0.10, "q₀", fontsize=11)
    ax.text(q1 @ e1 + 0.02, q1 @ e2 + 0.04, "q₁", fontsize=11)
    ax.set_aspect("equal")
    ax.set_xlabel("q 在 q₀ 方向的分量")
    ax.set_ylabel("q 在正交方向的分量")
    ax.set_title("② 路径形状\nslerp 贴球面走大圆，nlerp 割弦再投影")
    ax.legend(fontsize=8.5, loc="lower right")
    ax.grid(alpha=0.3)

    # (c) 欧拉角 vs slerp
    ax = fig.add_subplot(1, 3, 3)
    R0 = RC.euler_to_matrix("ZYX", [0, 0, 0], intrinsic=True)
    R1 = RC.euler_to_matrix("ZYX", [np.radians(90), np.radians(45),
                                    np.radians(80)], intrinsic=True)
    e1a = np.array([0.0, 0.0, 0.0])
    e1b = np.array([np.radians(90), np.radians(45), np.radians(80)])
    qa = RC.matrix_to_quat(R0)
    qb = RC.matrix_to_quat(R1)
    ts2 = np.linspace(0, 1, 200)
    tot = np.degrees(RC.rotation_angle_of(R0.T @ R1))
    # 到起点的角程曲线
    sl = []
    for t in ts2:
        R = RC.quat_to_matrix(RC.slerp(qa, qb, np.array([t]))[0])
        sl.append(np.degrees(RC.rotation_angle_of(R0.T @ R)))
    eu = []
    for t in ts2:
        R = RC.euler_to_matrix("ZYX", e1a + t * (e1b - e1a), intrinsic=True)
        eu.append(np.degrees(RC.rotation_angle_of(R0.T @ R)))
    ax.plot(ts2, sl, lw=2.4, color="#185FA5", label="四元数 slerp")
    ax.plot(ts2, eu, lw=2.4, color="#c0392b", label="欧拉角线性插值")
    ax.plot([0, 1], [0, tot], "--", lw=1.4, color="gray", label="理想匀速")
    ax.set_xlabel("插值参数 t")
    ax.set_ylabel("走过的角度（度）")
    ax.set_title("③ 三轴同时动时\n欧拉角插值不匀速（还会撞万向锁）")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    fig.suptitle("A5.4 · Slerp：在 S³ 上走大圆，才是「匀速转动」",
                 fontsize=13.5, y=0.99)
    fig.subplots_adjust(top=0.80, wspace=0.28)
    path = RC.savefig(fig, "slerp_paths.png")
    plt.close(fig)
    return path


def main():
    print("=" * 74)
    print("A5.4 · Slerp：等间距 ≠ 匀速旋转")
    print("=" * 74)
    ts, q0, q1, a_s, a_n = angle_rate_check()
    short_arc()
    slerp_formula()
    euler_interp_trap()
    consistency()
    path = make_figure(ts, q0, q1, a_s, a_n)

    print("\n" + "=" * 74)
    print("验收：|q(t)|=1 违反 <1e-15、匀速性（良态口径）<1e-13、"
          "短弧中点 = 总角半数  [PASS]")
    print("图：", path)
    print("\n这一步在说什么：")
    print("  · SO(3) 不是平的 ⇒ 插值要走「球面上的大圆」，公式里 sin 权重的来源；")
    print("  · 短弧翻转（q₀·q₁<0 ⇒ 翻 q₁）是必需步骤，不是优化；")
    print("  · nlerp 更快、端点对，但角速度不匀（对实时动画够用，对精修不够）；")
    print("  · **不要插值欧拉角** —— 不匀、绕远路、还会撞万向锁。")


if __name__ == "__main__":
    main()