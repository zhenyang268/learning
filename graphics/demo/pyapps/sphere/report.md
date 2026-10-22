# 球面几何与曲率 · Python 数值验证报告

> 位置：`graphics/demo/pyapps/sphere/`（纯 Python 现象演示区，不进 ctest、不进 CMake）
> 运行方式：`uv run python graphics/demo/pyapps/sphere/<名字>.py`，每个文件独立可跑
> 图片输出：`learning-todo/imgs/geo/`
> 日期：2026-10-17（修订）｜ 依赖：numpy + matplotlib + sympy（未用 scipy —— ODE 全部手写 RK4）
> 编号依据：`learning-todo/fan-plan.md` §3.1 的 **G 系列**（早期 T 系列编号已废弃）

---

## 0. 主线：五步走完「张量入门」

**度规定长度 → 联络定搬运 → 测地线定直线 → 平行移动定转角 → 曲率是转角的来源**

| 步 | 文件 | 一句话 | 验收数字 |
|---|---|---|---|
| G5 | `curvature.py` | Γ 会发散，K 不会 —— 张量被逼出来的原因 | K 偏离 1 < 4.4e-14 |
| G1 | `great_circle.py` | 投影会撒谎，度规不会 | 北京→纽约 10 991 km，恒向线多飞 30.9% |
| G2 | `geodesic_eq.py` | Γ 当发动机，跑出来的是大圆 | vs 解析大圆 5.1e-10；SO(3) 与 cos/sin 差 1.1e-16 |
| G3 | `parallel_transport.py` | 方向没有绝对意义，搬一圈的差角就是曲率 | \|V\| 漂移 3.2e-15，holonomy 对拍 ≤1.7e-11° |
| **G4★** | `octant.py` | 角超 = 面积 = holonomy = 90° | 三者偏离 0.00e+00 |
| G6 | `deviation.py` | 曲率怎么变成力：ξ'' + Kξ = 0 | R² = 1.00000000，偏差理论式逐位吻合 |

> **编号说明**：fan-plan §3.1 把球面系列定为 G1–G5，其中 G1=大圆、G2=测地线方程、
> G3=平行移动、G4=八分仪/角超、G5=曲率张量；`deviation.py`（测地偏离）是 G3 与 G4 之间的
> 延伸件，此处记作 **G6** 并已在 fan-plan 登记。**旧的 T1–T4 编号已作废**，
> 早期报告里的「T2★」即现在的 G4★。

**建议阅读顺序**（由直觉到抽象）：
G1（最短路径）→ G3（平行移动，最直观的曲率）→ G4★（三重验证收口）→ G2（Γ 方程）
→ G5（张量动机）→ G6（偏离方程 → 广义相对论预告）。

---

## 1. G5 曲率张量（`curvature.py`）

单位球面度规：$ds^2 = d\theta^2 + \sin^2\theta\, d\varphi^2$（半径 R 时整体乘 $R^2$）。

手写 Christoffel（从度规算，同时让 sympy 独立推一遍，两路一致）：

$$\Gamma^\theta_{\varphi\varphi} = -\sin\theta\cos\theta,\qquad
\Gamma^\varphi_{\theta\varphi} = \Gamma^\varphi_{\varphi\theta} = \cot\theta,\qquad \text{其余} = 0$$

**三个关键事实**（这是原版没讲透的地方）：

1. **g 含 $R^2$，Γ 不含 R** —— 度规整体乘常数时 $(1/2)g^{kl}(\cdots)$ 上下相消，
   所以球半径多大联络都一样。
2. **$\Gamma$ 在极点发散不是球的错，是坐标基退化**：$|\partial_\varphi| = \sin\theta$，
   在极点缩成 0。用长度趋零的箭头去丈量正常的切向量，分量自然趋无穷。
3. **$R^2$ 在曲率里重新出现**：$K = 1/R^2$。

曲率张量（约定 $R^i{}_{jkl} = \partial_k\Gamma^i_{lj} - \partial_l\Gamma^i_{kj} + \Gamma^i_{km}\Gamma^m_{lj} - \Gamma^i_{lm}\Gamma^m_{kj}$）：

$$R_{\theta\varphi\theta\varphi} = \sin^2\theta,\qquad
R_{ijkl} = K\,(g_{ik}g_{jl} - g_{il}g_{jk}),\qquad
K = \frac{1}{R^2},\qquad R_{ij} = K\,g_{ij}$$

**验收**：K 两条独立算法（$R_{\theta\varphi\theta\varphi}/\det g$、里奇标量/2）+ sympy 符号推导，
全等于 1，最大偏离 **4.396e-14**；16 个分量恒等式 $R_{ijkl}=K(g_{ik}g_{jl}-g_{il}g_{jk})$
全过；缩放检验 $K=1/R^2$ 对 R = 0.5/1/2/6371 全对。

**图（curvature_basics.png，四联）**：① Γ 两极发散 → ② **原因**：$|\partial_\varphi|=\sin\theta$ 在极点归零
→ ③ K 恒为 1 → ④ $R_{ij}=K g_{ij}$。①②的因果对照是本项目最想让人带走的一张图：
**联络不是张量（随坐标说谎），曲率是张量（坐标换不动它）**。

## 2. G1 大圆 vs 等角航线（`great_circle.py`）

北京 (39.9°N, 116.4°E) → 纽约 (40.7°N, 74.0°W)：

| 航线 | 长度 | 说明 |
|---|---|---|
| 大圆 | **10 991 km** | 几何法与 haversine 对拍，差 5.5e-12 km |
| 等角航线 | 14 383 km | 墨卡托空间直线，多飞 **30.9%** |

大圆最高爬到 **83.9°N** —— 出发地才 40°N，为了抄近路飞机往北极绕。

**图（great_circle_vs_rhumb.png，三联）**：① 球面视角 ② 经纬度平面（短的弯、直的贵）
③ **墨卡托平面**：恒向线在这里是**直线**。

**新增的量化判据**：曲线是直线 ⟺ 斜率处处相同。
经纬度平面 $d\lambda/d\varphi$ 相对波动 1.2e-2；墨卡托平面 $d\lambda/d\psi$ 相对波动 **2.0e-11**（舍入级）。
—— 墨卡托的非线性拉伸 $z=\psi(\varphi)$ 恰好把恒向线拉成直线，代价是长度被放大。

## 3. G2 测地线方程（`geodesic_eq.py`）

$$\ddot\theta = \sin\theta\cos\theta\,\dot\varphi^2,\qquad
\ddot\varphi = -2\cot\theta\,\dot\theta\dot\varphi$$

**白话**：坐标加速度不为零不是"有力"，是基向量自己在转；两者相抵就是"直线"。

**三条路会师同一条曲线**（判据独立于实现路径）：

| 路线 | 做法 | 与数值解偏差 |
|---|---|---|
| A | RK4 积分测地线方程 | —— |
| B | SO(3)：$\gamma(s) = R(s\hat n)p_0$（罗德里格斯公式） | 5.1e-10 |
| C | cos/sin 线性组合 | 5.1e-10 |

且 **B 与 C 逐位相同（差 1.1e-16）** —— "绕轴转一个角"与"线性组合"是同一件事的两种写法。
（原版只在注释里声称"用旋转矩阵生成"，代码里并没有真的用；现在真的用了 `rotate_about`。）

- Clairaut 守恒量 $\sin^2\theta\,\dot\varphi = \mathrm{const}$ 漂移 **2.1e-10** —— 测地线自带角动量守恒；
- 收敛比 15.9~16.0（= 2⁴，四阶确认）。

## 4. G3 平行移动（`parallel_transport.py`）

球坐标方程（由 $\frac{dV^k}{ds} + \Gamma^k_{ij}V^i\dot\gamma^j = 0$ 展开）：

$$\frac{dV^\theta}{ds} - \sin\theta\cos\theta\, V^\varphi\dot\varphi = 0,\qquad
\frac{dV^\varphi}{ds} + \cot\theta\,(V^\theta\dot\varphi + V^\varphi\dot\theta) = 0$$

**第二套实现**（笛卡尔嵌入，极点无奇点）：$\dfrac{dV}{ds} = -(V\cdot\dot\gamma)\,\gamma$。
**白话**：平行移动 = 向量唯一被迫做的改变，就是跟着曲面转的那部分。

**holonomy = $-2\pi\cos\theta_0$**（在正交标架下连续展开），9 个纬度扫描最大偏差 **1.7e-11°**：

| $\theta_0$ | 10° | 60° | 90° | 170° |
|---|---|---|---|---|
| 转角 | −354.5° | **−180°** | 0° | +354.5° |

- 赤道（大圆 = 测地线）转角 0°；贴北极的小圈 ≈ −360° ≡ 0°（局部近似平面）；
- 脚本同时给出「直接比首末向量」的 mod 360° 代表元，并逐点校验两者同余（全 ✓）；
- |V| 全程漂移 **3.2e-15**（阈值 1e-8）；RK4 步长减半误差 ÷16.0 稳定。

**图（transport_latitude.png）**：左 = 纬线上的向量演化；右 = holonomy 随纬度 vs 理论曲线。

## 5. G4★ 八分仪三重验证（`octant.py`）

顶点 A=(1,0,0)、B=(0,1,0)、C=(0,0,1)，三条 π/2 大圆弧。

| 量 | 数值 | 来源 |
|---|---|---|
| 内角和 − π | **90.000000°** | 每顶点两边切向量夹角（全是 90°） |
| 面积 | **π/2** | 解析 4π/8 与 ∫∫sinθ dθdφ（2000 格中点规则）互证，差 4.0e-8 |
| holonomy | **90.000000°** | 首末向量夹角（两者都在 A 的切平面内） |

三者最大偏离 **0.00e+00** —— 球面角超公式 $(\alpha+\beta+\gamma) - \pi = \text{面积}/R^2$ 的实例化，
即 Gauss–Bonnet 测地多边形版本：**角超 = 面积 = holonomy，是同一个东西的三种量法**。

**必须说清的两点**（原版含糊或写错）：

1. 面积本身单位是 sr，没有"度"。汇总里换算成角当量只是为了与另两个角同轴对比 ——
   用的正是"面积 = 角超"这条恒等式，所以那个 90° **不是独立测量**，这一点在输出里已注明。
2. 「转角 = 90°」用的是**首末向量夹角**口径。若改用「沿途切向标架下连续展开」的口径，
   得到的是 **−180°**。两口径**不 mod 360 同余**（参照系不同：一个是向量间夹角，
   一个是相对局部标架的落后量）。原版注释写的「−270°」是凭记忆填的，代码里从没算过，
   已删除。

- 动画 `octant_transport.gif`：长矛沿八分仪走一圈，末帧标注转角。

## 6. G6 测地偏离（`deviation.py`）

赤道上经度差 Δλ 的两点同时向北，分离距离精确解 $\sin(\xi/2) = \sin\theta\,\sin(\Delta\lambda/2)$，
线性化 $\xi(s) = \Delta\lambda\cos s$，满足**雅可比方程**

$$\xi'' + K\xi = 0 \quad (K=1)\qquad \text{—— 简谐振动方程}$$

**白话**：正曲率 = 会聚透镜（经线在北极=共轭点聚焦），负曲率 = 发散透镜。
把测地线换成自由落体，这就是广义相对论的潮汐力方程。

**原版写错了两个理论式，现已修正并坐实**：

| 量 | 正确理论式 | 实测 vs 理论 |
|---|---|---|
| 精确解 − 线性化 | $-\dfrac{\Delta\lambda^3}{24}\sin\theta\cos^2\theta$ | max 2.005e-06，逐点差 **2.2e-10** |
| 精确解代雅可比残差 | $-\dfrac{\Delta\lambda^3}{12}\cos 3s$，max $= \Delta\lambda^3/12$ | 1.042e-05 **完全吻合** |

（原版注释写「偏差 ~ Δλ³/24·Δλ = Δλ⁴/24」= 2.6e-07，与实测 2.0e-06 差 8 倍，是错的。）

- 数值（两条 RK4 测地线实测距离）vs 精确球面距离：**4.2e-15**；
- 偏差按 **Δλ³** 缩放：相邻比恒 125 = 5³，末行累计比 15625 = 25³；
- 拟合 R² = **1.00000000**（> 0.999 ✓）；共轭点 s = π/2；
- 半径缩放 R = 1/2/5 全对。

**图（deviation_meridians.png）**：左 = 经线束向北汇聚；右 = ξ(s) 数值/精确/线性化三条线叠合。

---

## 7. 方法学备忘（本次踩过的真坑）

1. **点列 + 两点间大圆插值 = 内接测地多边形**：曲线本身被换掉了，误差 O(a²) 锁死在 2e-6，
   步长再小没用。解法：解析曲线 $\gamma(s)$ 接口，让 RK4 的 stage 取到真实中点 → 四阶恢复（5e-13）。
   教训与 graphics 3.2「判据要独立于实现路径」同族。
2. **坐标奇点是教学素材，但必须换路线**：八分仪顶点 C 恰在极点，球坐标 cotθ 爆炸 ——
   平行移动换笛卡尔嵌入路线；测地偏离只积分到 θ=0.02。
3. **理论对比要用精确公式而不是小量近似**：$\xi \approx \Delta\lambda\sin\theta$ 与精确
   $\sin(\xi/2)=\sin\theta\sin(\Delta\lambda/2)$ 差 O(Δλ³)。**两者要分清是"哪个"三阶项**：
   偏差式是 $\sin\theta\cos^2\theta$，雅可比残差是 $\cos 3s$ —— 系数不同，别混用。
4. arccos 型距离公式在小距离丢位，换 arcsin（haversine 型）稳。
5. **跨 180° 经线的路径画平面图必须解卷**：北京→纽约大圆跨 180，不解卷会拉出一条
   横跨整图的假跳线（实测跳变 358°）。已封装为 `_core.unwrap_lon_deg`。
6. **"短弧判定"要拿归一化后的点比**：若用单位球上的路径终点去比未归一化的 p2，
   距离恒 ≈ |p2| ≈ 6371，判据永远为真 —— 翻转分支形同虚设（原 `great_circle_path` 的 bug）。

## 8. 待做（下一轮）

- G7 三种几何对照（平面 K=0 / 球面 K=+1 / 庞加莱圆盘 K=−1，同边长比内角和）
- 锥面缺失角（绕锥顶转角 = α，锥面除顶点 K=0 —— 曲率可以集中在测度零的点上）
- 交互式：可拖动的球面三角形实时显示内角和与 holonomy

## 9. 运行一览

```bash
cd /Users/apple/develop/learning
uv run python graphics/demo/pyapps/sphere/great_circle.py        # G1
uv run python graphics/demo/pyapps/sphere/geodesic_eq.py         # G2
uv run python graphics/demo/pyapps/sphere/parallel_transport.py  # G3
uv run python graphics/demo/pyapps/sphere/octant.py              # G4★ + GIF
uv run python graphics/demo/pyapps/sphere/curvature.py           # G5
uv run python graphics/demo/pyapps/sphere/deviation.py           # G6
```

公共数学层 `_core.py`：坐标转换（含经度解卷）、度规/Γ/黎曼（解析导数 + 差分对拍开关）、
RK4、大圆（含短弧择向）、haversine、平行移动积分器、SO(3) 旋转（罗德里格斯）、
3D 球面底图辅助。