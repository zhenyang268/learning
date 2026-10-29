"""rk4_step / integrate 的数值判据 —— 这是 C++ 侧 math/ode.hpp 的对拍基准。

为什么先测 Python 这一侧: `NOW.md` 里 ode 的"算做完"写着要与
`graphics/demo/pyapps/sphere/_core.py` 的**手写 RK4 逐点对拍**。
对拍要求 Python 侧给出的是**可被程序调用的数值**，而不是一张人眼看的图 ——
所以先把这层数值核钉住，C++ 实现才有真基准可对。

契约 (读 `_core.py` 得到，不是猜的):
  rk4_step(f, s, y, h) -> 新状态;  y 与 f 的返回值都是 numpy 数组
  integrate(f, y0, s_grid) -> (N, dim) 轨迹;  y0 标量时 dim=1, 取 `[-1, 0]`

判据边界表
----------
段 | 阈值                  | 实测值 (2026-10-30)        | 专抓什么错
---|-----------------------|----------------------------|--------------------------
A  | 相对误差 <= 1e-11     | h=1e-2 时 8.26e-13         | 四阶权重写错 (漏 /6, 2/6<->1/6)
B  | 每级观测阶 >= 3.9     | 3.98 / 3.99 / 4.00         | 退化成欧拉或 RK2 (阶掉到 1/2)
C  | rtol <= 1e-5          | h=0.1 走 10 步 9.06e-7     | 只对标量正确, 数组广播/索引错位

阈值为什么留这些余量:
  · A 的 1e-11 对实测 8.26e-13 留 ~12 倍余量。卡到 1e-12 只剩 1.2 倍,
    换台机器换 BLAS 就会假红。
  · C 的 1e-5 对实测 9.06e-7 留 ~11 倍。h=0.1 的 RK4 全局误差本就 ~h^4,
    要求 1e-10 是在测舍入, 不是在测算法。
"""

from __future__ import annotations

import numpy as np

from _core import integrate, rk4_step


def test_rk4_step_matches_analytic():
    """A 段: 单步对 y' = y 的解析解 e^h，相对误差 <= 1e-11。"""
    f = lambda s, y: y
    h = 1e-2
    got = float(rk4_step(f, 0.0, 1.0, h))
    want = np.exp(h)
    assert abs(got - want) / want <= 1e-11


def test_rk4_order_is_four():
    """B 段: h 减半，误差降约 16 倍 ⇒ 每级观测阶 >= 3.9。"""
    f = lambda s, y: y
    t_end = 1.0
    errs = []
    for n in (40, 80, 160, 320):
        traj = integrate(f, 1.0, np.linspace(0.0, t_end, n + 1))
        errs.append(abs(float(traj[-1, 0]) - np.exp(t_end)) / np.exp(t_end))

    orders = [float(np.log2(errs[i] / errs[i + 1])) for i in range(len(errs) - 1)]
    assert min(orders) >= 3.9, f"观测阶 {orders} 掉出四阶，实测误差序列 {errs}"


def test_rk4_array_state():
    """C 段: 状态是数组时，两个独立分量必须各自走对 (抓广播/索引错位)。"""
    # y0' = -y0, y1' = y1 —— 两个分量方向相反，广播错位会立刻可见
    f = lambda s, y: np.array([-y[0], y[1]])
    y = np.array([1.0, 1.0])
    for _ in range(10):
        y = rk4_step(f, 0.0, y, 0.1)

    want = np.array([np.exp(-1.0), np.exp(1.0)])
    rel = np.abs((y - want) / want)
    assert np.all(rel <= 1e-5), f"各分量相对误差 {rel}"
