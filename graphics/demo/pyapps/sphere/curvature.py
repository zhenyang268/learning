"""
G5 · 曲率张量：度规 -> Christoffel -> 黎曼 -> 里奇 -> 高斯曲率
=========================================================================
现象：球面的 Christoffel 符号在两极发散（cot theta -> ∞），
      但由它组合出来的高斯曲率处处恒等于 1/R^2。

一句话在说什么：
    Christoffel（联络系数）**不是张量** —— 它随坐标而变，可以把一个好好的点
    写成无穷大；曲率才是张量 —— 坐标怎么换，它的值不动。这就是"张量"这个概念
    被逼出来的原始动机：**我们需要一种不随坐标说谎的量**。

Γ 为什么会炸？不是球有问题，是**坐标基退化了**：
    坐标基 |d_phi| = sin(theta)，在极点缩成 0。用两个长度趋于 0 的箭头
    去丈量一个正常的切向量，分量自然要趋于无穷 —— 分量无穷大，
    向量本身有界。这就是"分量不是几何"的最朴素例子。

运行： uv run python graphics/demo/pyapps/sphere/curvature.py
输出： learning-todo/imgs/geo/curvature_basics.png
"""

from __future__ import annotations

import os
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _core as C  # noqa: E402

plt = C.setup_style()


# ------------------------------------------------------------ 1. 符号推导
def symbolic_derivation():
    """用 sympy 从度规出发完整推一遍 g -> Γ -> R^i_{jkl} -> Ric -> K。

    返回值里的 K 应当化简为 1/R^2（半径 R 的球面高斯曲率）。
    """
    th, ph = sp.symbols("theta phi", positive=True)
    R = sp.Symbol("R", positive=True)
    x = [th, ph]
    names = ["theta", "phi"]

    # 度规：ds^2 = R^2 dθ^2 + R^2 sin^2 θ dφ^2
    g = sp.Matrix([[R ** 2, 0],
                   [0, R ** 2 * sp.sin(th) ** 2]])
    ginv = g.inv()

    # Γ^k_{ij} = (1/2) g^{kl} (∂_i g_{lj} + ∂_j g_{li} − ∂_l g_{ij})
    G = [[[sp.simplify(sum(
        ginv[k, l] * (sp.diff(g[l, j], x[i]) + sp.diff(g[i, l], x[j])
                      - sp.diff(g[i, j], x[l])) / 2
        for l in range(2))) for j in range(2)] for i in range(2)]
        for k in range(2)]

    # R^i_{jkl} = ∂_k Γ^i_{lj} − ∂_l Γ^i_{kj} + Γ^i_{km} Γ^m_{lj} − Γ^i_{lm} Γ^m_{kj}
    def Riem(i, j, k, l):
        v = (sp.diff(G[i][l][j], x[k]) - sp.diff(G[i][k][j], x[l])
             + sum(G[i][k][m] * G[m][l][j] for m in range(2))
             - sum(G[i][l][m] * G[m][k][j] for m in range(2)))
        return sp.simplify(v)

    Rm = [[[[Riem(i, j, k, l) for l in range(2)] for k in range(2)]
           for j in range(2)] for i in range(2)]

    # Ric_{ij} = R^k_{ikj}
    Ric = sp.Matrix(2, 2, lambda i, j: sp.simplify(
        sum(Rm[k][i][k][j] for k in range(2))))
    # 里奇标量 R = g^{ij} Ric_{ij}
    Rs = sp.simplify(sum(ginv[i, j] * Ric[i, j] for i in range(2)
                         for j in range(2)))

    # 全下标：R_{ijkl} = g_{im} R^m_{jkl}
    def Rlow(i, j, k, l):
        return sp.simplify(sum(g[i, m] * Rm[m][j][k][l] for m in range(2)))

    # 球面应当满足 R_{ijkl} = K (g_{ik} g_{jl} − g_{il} g_{jk})
    K = sp.simplify(Rs / 2)

    print("=" * 66)
    print("【1. 符号推导】从度规一路推到 K（sympy 全程符号运算）")
    print("=" * 66)
    print("  度规 g_ij =")
    sp.pprint(g)
    print("\n  非零 Christoffel 分量 Γ^k_ij（注意：**不含 R**，只含 θ）：")
    for k in range(2):
        for i in range(2):
            for j in range(2):
                if G[k][i][j] != 0:
                    print(f"    Γ^{names[k]}_{names[i]}{names[j]} = "
                          f"{sp.simplify(G[k][i][j])}")
    print("\n  里奇张量 Ric_ij =")
    sp.pprint(Ric)
    print(f"\n  里奇标量  R  = {Rs}")
    print(f"  高斯曲率  K  = R/2 = {K}")
    print("  → g 含 R²，Γ 不含 R，K = 1/R²：**R² 在曲率里才重新出现**。")

    # 全下标恒等式：逐分量验，不只看一个
    worst = 0
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    d = sp.simplify(Rlow(i, j, k, l)
                                    - K * (g[i, k] * g[j, l] - g[i, l] * g[j, k]))
                    worst = max(worst, abs(complex(d.evalf(subs={th: 1.234, R: 2.0}))))
    print(f"  恒等式 R_ijkl = K(g_ik g_jl − g_il g_jk)：16 个分量最大偏离 = "
          f"{worst:.3e}  (应为 0)")
    print("  → 二维曲面上曲率只有**一个**自由度 K，所有分量被它锁死。")
    return th, R, K


# ------------------------------------------------------------ 2. 数值对拍
def numeric_checks():
    """四条独立路径对照，任一不一致就说明公式或代码有一处写错。"""
    print("\n" + "=" * 66)
    print("【2. 数值对拍】手写 Γ  vs  从度规有限差分反解 Γ")
    print("=" * 66)
    worst = 0.0
    for th in [0.15, 0.5, 1.0, np.pi / 2, 2.0, 2.9]:
        d = np.max(np.abs(C.christoffel(th)
                          - C.christoffel_from_metric([th, 0.0], C.metric)))
        worst = max(worst, d)
        print(f"    θ={th:5.2f}   max|ΔΓ| = {d:.3e}")
    print(f"    → 全区间最大偏差 {worst:.3e}（有限差分的理论误差量级）")

    print("\n" + "=" * 66)
    print("【3. 数值对拍】解析导数 Riemann  vs  有限差分 Riemann")
    print("=" * 66)
    for th in [0.4, 1.0, 2.0]:
        d = np.max(np.abs(C.riemann(th) - C.riemann(th, h=1e-4)))
        print(f"    θ={th:5.2f}   max|ΔR| = {d:.3e}")

    print("\n" + "=" * 66)
    print("【4. 数值对拍】K 的两条独立算法（单位球应当全等于 1）")
    print("=" * 66)
    print(f"    {'θ':>6}  {'K=R_θφθφ/det g':>16}  {'K=R_scalar/2':>14}  {'两者之差':>11}")
    worst = 0.0
    for th in np.linspace(0.05, np.pi - 0.05, 9):
        KA, KB = C.gaussian_curvature(th)
        worst = max(worst, abs(KA - KB), abs(KA - 1.0))
        print(f"    {th:6.3f}  {KA:16.12f}  {KB:14.12f}  {abs(KA - KB):11.2e}")
    print(f"    → 偏离 1 的最大值 = {worst:.3e}")
    print("    注：第二条路 = 里奇标量/2，二维里 R = 2K；它是**独立于**第一条路的")
    print("    张量运算路径，不是把同一个数再算一遍。")

    print("\n" + "=" * 66)
    print("【5. 缩放检验】半径 R 的球，K 应当 = 1/R²（Γ 不变、K 变）")
    print("=" * 66)
    for R in [0.5, 1.0, 2.0, 6371.0]:
        KA, _ = C.gaussian_curvature(0.9, R=R)
        print(f"    R={R:8.1f}   K={KA:.6e}   1/R²={1 / R ** 2:.6e}")
    print("    → 半径放大 2 倍，Γ 一动不动，K 缩到 1/4：几何在变，坐标没变。")
    return worst


# ------------------------------------------------------------ 3. 画图
def make_figure():
    """四连图：Γ 发散 → 坐标基退化（原因）→ K 恒定（结论）→ Ric = K·g。"""
    th = np.linspace(0.012, np.pi - 0.012, 900)
    deg = np.degrees(th)
    Gth_phph = np.array([C.christoffel(t)[0, 1, 1] for t in th])
    Gph_thph = np.array([C.christoffel(t)[1, 0, 1] for t in th])
    ephi_len = np.sin(th)                      # |∂_φ| = sin θ
    Ks = np.array([C.gaussian_curvature(t)[0] for t in th])
    Ric00 = np.array([C.ricci(t)[0, 0] for t in th])
    Ric11 = np.array([C.ricci(t)[1, 1] for t in th])
    g11 = np.array([C.metric(t)[1, 1] for t in th])

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.6))
    fig.subplots_adjust(top=0.86, hspace=0.34, wspace=0.24)
    (ax0, ax1), (ax2, ax3) = axes

    # ---- ① Γ 分量发散 ----
    ax0.plot(deg, Gth_phph, lw=2.2, label=r"$\Gamma^\theta_{\varphi\varphi}=-\sin\theta\cos\theta$")
    ax0.plot(deg, np.clip(Gph_thph, -60, 60), lw=2.2,
             label=r"$\Gamma^\varphi_{\theta\varphi}=\cot\theta$（纵轴截断到 ±60）")
    ax0.axhline(0, color="k", lw=0.6)
    ax0.set_xlabel("极角 θ（度）")
    ax0.set_ylabel("分量值")
    ax0.set_title("① 联络 Γ 在两极 → ∞\n（分量炸了）", color="#7f1d1d")
    ax0.legend(fontsize=9, loc="lower center")
    ax0.grid(alpha=0.3)
    ax0.set_ylim(-62, 62)
    for x0, ha in [(0, "left"), (180, "right")]:
        ax0.annotate("→ +∞", xy=(x0, 55), xytext=(x0 + (8 if x0 == 0 else -8), 40),
                     ha=ha, fontsize=13, color="#7f1d1d", weight="bold")

    # ---- ② 原因：坐标基退化 ----
    ax1.plot(deg, ephi_len, lw=2.5, color="#2980b9")
    ax1.fill_between(deg, 0, ephi_len, color="#2980b9", alpha=0.12)
    ax1.axhline(1.0, ls="--", color="gray", lw=1)
    ax1.set_xlabel("极角 θ（度）")
    ax1.set_ylabel(r"坐标基长度 $|\partial_\varphi|=\sin\theta$")
    ax1.set_title("② 为什么炸：坐标基在极点缩成 0\n（不是球坏了，是尺子坏了）",
                  color="#1e3a8a")
    ax1.grid(alpha=0.3)
    ax1.set_ylim(0, 1.12)
    for x0, ha in [(0, "left"), (180, "right")]:
        ax1.annotate("|∂φ| = 0\n（极点）", xy=(x0, 0.02),
                     xytext=(x0 + (30 if x0 == 0 else -30), 0.35), ha=ha,
                     fontsize=9, color="#1e3a8a",
                     arrowprops=dict(arrowstyle="->", lw=0.9, color="#1e3a8a"))
    ax1.annotate("赤道：坐标基最长", xy=(90, 1.0), xytext=(96, 0.82), fontsize=9,
                 arrowprops=dict(arrowstyle="->", lw=0.9))

    # ---- ③ K 恒定 ----
    ax2.plot(deg, Ks, lw=2.6, color="#c0392b")
    ax2.axhline(1.0, ls="--", color="gray", lw=1, label="理论值 K = 1")
    ax2.set_xlabel("极角 θ（度）")
    ax2.set_ylabel("高斯曲率 K")
    ax2.set_title("③ 曲率 K ≡ 1，一点没炸\n（曲率是张量：坐标换不动它）",
                  color="#7f1d1d")
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.3)
    ax2.set_ylim(0.9, 1.1)

    # ---- ④ Ric = K·g ----
    ax3.plot(deg, Ric00, lw=2.2, label=r"$R_{\theta\theta}$")
    ax3.plot(deg, Ric11, lw=2.2, label=r"$R_{\varphi\varphi}$")
    ax3.plot(deg, g11, ls="--", lw=1.8, color="green",
             label=r"$g_{\varphi\varphi}=\sin^2\theta$")
    ax3.set_xlabel("极角 θ（度）")
    ax3.set_ylabel("分量值")
    ax3.set_title("④ 里奇张量 = K · 度规\n（分量随坐标变，但比例恒为 K）")
    ax3.legend(fontsize=9)
    ax3.grid(alpha=0.3)
    ax3.annotate("两条曲线处处重合\n⇒ $R_{ij}=1\\cdot g_{ij}$", xy=(90, 1.0),
                 xytext=(104, 0.55), fontsize=9,
                 arrowprops=dict(arrowstyle="->", lw=0.9))

    fig.suptitle("G5 · 联络会发散，曲率不会 —— 这就是“张量”被逼出来的原因",
                 fontsize=13.5, y=0.965)
    path = C.savefig(fig, "curvature_basics.png")
    plt.close(fig)
    return path


def main():
    print("=" * 66)
    print("G5 · 曲率张量：ds² = R²dθ² + R²sin²θ dφ²  →  K = 1/R²")
    print("=" * 66)
    symbolic_derivation()
    worst = numeric_checks()
    path = make_figure()

    print("\n" + "=" * 66)
    print("验收：K 偏离 1 的最大值 = "
          f"{worst:.3e}  [{'PASS' if worst < 1e-10 else 'FAIL'}]")
    print("图：", path)
    print("\n这一步在说什么：")
    print("  度规给定了「怎么量长度」；联络 Γ 是从度规算出来的「基向量怎么转」；")
    print("  曲率 R 是从 Γ 算出来的「转一圈回不到原处多少」。")
    print("  Γ 在极点发散只是坐标的锅（端点处坐标基 ∂_φ 长度 = sinθ → 0），")
    print("  几何本身完好 —— 判断依据就是：由 Γ 组装出来的 K 处处是 1，一点没炸。")
    print("  **分量不是几何，张量才是。**")


if __name__ == "__main__":
    main()