# plan.md — 数学 × 图形学 学习实施计划

> **读者**：agent（执行工程任务）+ 人 ｜ **最后更新**：2026-09-22
>
> **真源声明**
> - 本文只管**任务定义与验收阈值**，**不记录任何状态**（不勾选、不写「进行中」）。
> - 进度（现在在哪 / 下一步）→ `dev-progress.md`（唯一真源）。
> - 矩阵与线代的**设计理由与红线** → `learning-todo/matrix.md`，本文不重述，只引用。
> - 目录分层与硬规则 → `AGENTS.md`。
>
> **验收编号约定**：`V{阶段}.{序号}`（如 `V1.3` = 阶段 1 第 3 条）。
> 旧文中的裸字母编号 `A1`–`A7`、`B1`–`B7`、`C1`–`C4`、`D1`–`D5`、`E1`–`E5`、`F1`–`F4` **已废弃**——
> 它们与笔记代号（`B1`、`D1`、`E1`…）同名，是 agent 误读的主要来源。对照表见文末。
>
> **符号约定**：本文属工程文档，矩阵算子一律写 ASCII —— `A^T`（转置）、`A^H`（共轭转置）。
>
> 主线原则：以 graphics 项目为骨架，每实现一个模块前精读 notebook 对应章节；
> 每个 demo 完成后回到 notebook 对应 `.tex` 追加「实现笔记 + 踩坑」一节，
> 让 AI 原稿笔记变成亲手验证过的版本。

---

## 0. 设计红线（摘要，理由见 `matrix.md`）

违反以下任一条视为 bug；**要改先改 `matrix.md` 并写清理由，再动代码**。

1. C++20，`Scalar` concept + `ScalarTraits<T>`；**实复差异只存在于 traits 特化，算法层零分支**。
2. Matrix 是**哑容器**：扁平行主序存储、`operator[]` 即时派生行指针（不存 `vector<T*>`）、
   逐元素四则 + Transpose/Adjoint。**其余运算全是 linalg/math 自由函数，不进类**。无缓存层。
3. **函数与算子不进 `Matrix<T>`**：`Scalar` 额外要求乘法交换（`CommutativeScalar`）；
   微分算子走 `LinearOperator{repr, basis, truncation}`（算子 = 基 + 边界条件 + 截断矩阵）。
4. 高维数据靠**展平 / 切片**归约为二维线性代数，**不造 Tensor 类**；固定尺寸 `Matrix<T,M,N>` 与 SIMD 推迟。
5. 复数自实现 `Complex<T>`，不用 `std::complex`；求逆/特征值/SVD 每阶段只实现一种通用简单方法。

---

## 阶段 1：地基 —— 泛型标量 + 矩阵容器 + 最小 linalg

**对应笔记**：B1（矩阵→群论）、§26–27（SO(2)/SO(3)/SU(2) 部分）、G3（数值线代）

**任务**

1. `core/scalar.hpp`：Scalar 概念 + float/double traits。
2. `complex/complex.hpp`：Complex<T> 四则/模/共轭/exp（欧拉公式），特化 traits。
3. 重写 `core/matrix.hpp`：清 P0（`#pragma once`、删行指针层、构造歧义改静态工厂
   `Zero/Identity/Diagonal`、friend 不匹配、删死代码、const 正确性、缓存机制整体删除）。
4. 修 `autodiff/dual.hpp`：isZero/比较走 traits；`SumLine/MultLine/MultDiagnal`
   迁出到 `autodiff/jacobian.hpp`，Matrix 中零残留。
5. `linalg/solve.hpp`：LU 部分选主元 + `Solve(A, b)`。
6. `linalg/inverse.hpp`：`Inverse(A)` = LU + 逐列回代（唯一方法）。
7. `linalg/eigen.hpp`：Jacobi 旋转法，仅对称/Hermitian。
8. `linalg/svd.hpp`：经 `A^H·A` 特征分解组装 SVD。

**验收节点（全部通过才进阶段 2）**

- **V1.1** `Matrix<double>` 与 `Matrix<Complex<double>>` 使用**同一份** linalg 源码
  零修改编译通过（traits 设计正确性的硬验收）。
- **V1.2** 随机矩阵求逆残差 ‖A·A^-1 − I‖∞ < 1e-8（double）；复矩阵同样通过。
- **V1.3** 对称阵 Jacobi 重构 max‖A − Q·Λ·Q^H‖ < 1e-8；Hermitian 复矩阵通过同一测试。
- **V1.4** SVD 重构 ‖U·diag(s)·V^H − A‖ < 1e-8；奇异值非负有序。
- **V1.5** `Dual<Complex<double>>` 作为 Matrix 元素编译通过，2×2 复矩阵 LU 的
  value 与直接用 `Complex<double>` 计算一致（为波包演化铺路）。
- **V1.6** 对易关系数值验证：so(3) 生成元 exp(Ât) 作用向量 = 旋转；`[Jx, Jy] = Jz`。
- **V1.7** 复数乘法交互 demo：平面上单位复数乘法 = 旋转（纯演示，不作硬验收）。

**不过时先查什么**

- V1.1 编译不过 → 大概率是某处**直接比较元素**或**漏了 `conj`**，而不是 concept 写错了。
- V1.2/V1.3/V1.4 残差偏大 → 先查主元选取与回代下标（差一错误），再查
  `abs2` 是否被误写成开方、容差是否用了全局 EPSILON 而非 `numeric_limits<real_t>::epsilon()` 的倍数。

**测试需完善**

- demo 注册进 ctest；公共断言工具 `demo/common.hpp`：`IsApprox(A,B,eps)`、
  随机矩阵生成（实/复/Dual 三种元素）、行列式/残差范数打印。
- LU 主元交换正确性：与 naive 高斯消元在小矩阵上对拍。
- 边界用例：1×1、奇异矩阵抛异常、非方阵（Transpose/Adjoint/SVD 输入）。

---

## 阶段 2：线代核心 + 旋转（李群实操）

**对应笔记**：B1 表示论基础、§32.1 四元数与 SU(2)、§22.1–22.2（泰勒↔拉普拉斯算子、矩阵指数）、G4（Padé）

**任务**

1. `linalg/qr.hpp`：Householder QR。
2. `linalg/cholesky.hpp`：`LL^H` 分解（Hermitian 正定）。
3. `linalg/expm.hpp`：scaling-and-squaring + Padé；先用 Taylor 级数版做对照实现。
4. `linalg/classify.hpp`：`Classify(A)`（Diagonal / Orthogonal(酉) / RigidTransform / SPD）
   + `Inverse(A, hint)` 快速路径：正交→Adjoint、对角→倒数、刚体→分块公式。
5. `quat/quaternion.hpp`：`Quaternion<T>`、与旋转矩阵互转、Slerp。
6. `linalg/operator.hpp`：`LinearOperator` + 基/边界条件工厂
   （`DifferentiationMatrix` / `LaplacianMatrix` / `MultiplicationMatrix`）。
   配套手推：`learning.md` 第二阶段 2.2 全部四题。
7. （可选，了解层）一般非对称特征值：Hessenberg + 位移 QR —— 不阻塞主线。

**验收节点**

- **V2.1** QR 重构 ‖QR−A‖ 小、Q 列正交性 ‖Q^H·Q−I‖ 小（实/复双类型）。
- **V2.2** Cholesky 与 LU 在 SPD 矩阵上解同一方程组，解之差 < 1e-8；
  对比两者 flops 计时，验证「快约 1/3」（对应 `matrix.md` 的求逆方法表）。
- **V2.3** `Classify` 正确识别构造出的各类特殊矩阵；快速求逆与 LU 求逆结果一致。
- **V2.4** `expm(A)` 与 Taylor 级数版在小矩阵上差 < 1e-8；
  `expm(t·Jx)` 作用向量 = 绕 x 轴旋转 t（李代数→李群指数映射数值版）。
- **V2.5** 四元数↔矩阵往返一致；`Slerp(q0→q1)` 端点正确、插值轨道为匀速大圆弧
  （逐帧 ‖Δ角‖ 恒定）。
- **V2.6** 万向锁 demo：绕 Y 90° 后 X/Z 轴重合，欧拉角→矩阵可视化演示。
- **V2.7** `LinearOperator` 三项验证：多项式基下 D 矩阵 = 移位阵、
  `expm(t·D)` 作用在多项式系数向量 = Taylor 移位（与解析对拍）；
  Dirichlet 正弦基下 Laplacian 为对角阵 `diag(-n^2)`；
  对易子 `Compose(D,X) − Compose(X,D) → I` 的残差随截断 N 收敛到 0（记录收敛率）。

**不过时先查什么**

- Q 不正交 → 先查 Householder 反射的符号约定与 `v` 的归一化，再查是否漏了在复数情形取相位。
- V2.4 偏差大 → 先查 Padé 阶数与 scaling 的平衡指数取整方式（`ceil` / `floor` 差一即失阶）。
- V2.7 对易子残差不收敛 → 先查基的截断是**两胞**还是**单胞**（非交换性来自边界项）。

**测试需完善**

- 旋转群性质回归：det(R)=1、`R·R^T = I` 对所有生成的旋转成立。
- Slerp 退化为 Lerp 的小角度误差界测试。

---

## 阶段 3：ODE = 动力学的可视化 + 渲染窗口骨架

**对应笔记**：de1–de4（第一/三层）、G1.1（RK4/辛积分器）、D1–D3（拉格朗日/Noether/哈密顿）、M7（变分积分器）、cmath/pendulum.c 升级

**任务**

1. `gfx/` 最小骨架：SDL2 窗口 + ImGui 面板 + 帧缓冲画点线（只到「能画曲线」，不做光栅化管线）。
2. `math/ode.hpp`：Euler、RK4、Verlet/辛欧拉，状态用 `Matrix<T>` 列向量表示。
3. 交互 demo：单摆参数滑条（长度/重力/初角），实时相图 (θ, ω)。
4. N 体重力模拟 demo（粒子渲染 + 轨迹拖尾）。
5. `x' = Ax` 线性系统：数值解 vs `expm(A t)` 解析解对拍。

**验收节点**

- **V3.1** 同一单摆，RK4 长时间能量漂移 vs 辛积分器能量有界振荡的曲线对比
  （matplotlib 或 ImGui 曲线绘制均可）。
- **V3.2** `x' = Ax` 数值解与 `expm` 解析解差 < 局部截断误差理论阶（RK4 为 O(h^4)：
  h 减半误差 ÷16 近似成立）。
- **V3.3** 两体问题出圆形/椭圆轨道，周期满足开普勒第三定律（误差 < 1%）。
- **V3.4** Noether 演示：旋转对称势场中角动量数值守恒（辛积分器下漂移 < 1e-6/千步）。

**不过时先查什么**

- 能量单调漂移 → 先查是不是用了非辛格式；有界振荡是**正常**的，不要「修」掉。
- V3.2 阶数不足 → 先查是否在积分过程中复用了同一 RHS 求值（RK4 要 4 次独立求值）。

**测试需完善**

- 积分器收敛阶自动化测试（步长序列误差斜率拟合）。
- 刚性方程（van der Pol 大 μ）展示显式格式失稳 —— 为阶段 4 隐式方法埋伏笔。

---

## 阶段 4：PDE 数值 —— 场和波 + 变分视角

**对应笔记**：de 第五~八层、G2（FDM/迭代法）、E1–E2（微分形式统一 grad/curl/div，
在 §21.5 弱解处回头精读）、I1/I2/I9（流体了解层）

**任务**

1. `math/pde.hpp`：一维热/波动方程显式 FDM；二维 Poisson 的 Jacobi/Gauss-Seidel/CG。
2. 可视化：热方程 = 图像模糊（`e^{tΔ}` 半群视角，帧缓冲直接显示）；
   二维波动方程 = 鼓膜振动。
3. CG 求解过程可视化：每次迭代画能量泛函 Φ(x) = (1/2)·x^T·A·x − b^T·x 下降 ——
   把「变分法求弱解」从公式变成动画。
4. （了解层，可选）Stam stable fluids 简版：平流 + 投影，作为阶段 1–4 的期末考试。

**验收节点**

- **V4.1** CFL 失稳实测：超临界步长立即爆炸，亚临界收敛，验证稳定条件与理论一致。
- **V4.2** 解析解（分离变量构造）对拍：热方程收敛阶符合格式理论。
- **V4.3** CG 在 SPD 矩阵上迭代次数 ≤ n 精确收敛；能量 Φ 单调下降；
  对比 Jacobi / Gauss-Seidel 收敛速度曲线。
- **V4.4** 离散余弦/正弦变换解 Poisson 与 CG 解一致（谱方法惊鸿一瞥，衔接 H5）。
- **V4.5**（若做流体）散度 ‖∇·u‖∞ 投影后 < 1e-6（对应笔记「散度自由条件的保持」）。

**不过时先查什么**

- 立即爆炸 → 先查 CFL 步长与**边界条件组装**（Dirichlet 行是否被正确消去）。
- Φ 不单调 → 先查 A 是否真 SPD（用 Cholesky 试跑一次即可判定），再查步长公式里的 `<Ap,p>` 下标。

**测试需完善**

- 迭代求解器统一测试脚手架：矩阵条件数 vs 迭代次数曲线。
- 边界条件（Dirichlet/Neumann）组装矩阵的单元对拍（小网格手工验证）。

---

## 阶段 5：采样与渲染 —— 概率测度论落地

**对应笔记**：C1/C5（概率测度论口径重读）、H2（渲染方程）/H4（蒙特卡洛）/H6（采样）/H7（辐射度）、J1/J5、M2（黑体辐射）

**任务**

1. `gfx` 补帧缓冲与光栅化基本件（三角形填充），打通 `demo/a_triangle`。
2. `math/mc.hpp`：LCG/xorshift RNG、Sobol/Halton 低差异序列、求积（中点/高斯）。
3. 蒙特卡洛积分 π：伪随机 vs 低差异 vs 分层，收敛率曲线 `1/√N` vs `N^-1` 边界。
4. CPU 光线追踪器：球求交、Lambert/Phong、镜面递归、progressive 累积显示。
5. 黑体辐射 Planck 曲线 → 颜色映射渲染成渐变条（M2 可视化）。

**验收节点**

- **V5.1** π 估计三种采样器收敛率 log-log 斜率与理论吻合。
- **V5.2** 半球余弦采样对 Lambert 积分方差为 0（重要性采样教科书验证）。
- **V5.3** 光线追踪球图与解析参考（镜面球对场景的镜像）像素级一致。
- **V5.4** progressive 渲染下图像噪声随 `1/√N` 消退（误差曲线验证）。
- **V5.5** `a_triangle` 上屏（渲染管线里程碑，纯工程不作数学验收）。

**不过时先查什么**

- 收敛率斜率不对 → 先查随机数是否在多次估计间被**重复使用同一种子**，再查低差异序列的维度顺序。
- V5.2 方差不严格为 0 → 先查 pdf 的归一化因子（余弦分布是 `cosθ/π` 而非 `cosθ`）。

**测试需完善**

- 求交/采样的解析对拍单测（射线-球、球面点均匀性的统计检验 χ^2）。
- RNG 周期与低差异序列偏差（discrepancy）数值测试。

---

## 阶段 6：交叉应用层（球谐 / 复数矩阵物理 / 离散几何）

**对应笔记**：H5（球谐）、H1（离散微分几何）、§32.3（SH = SO(3) 不可约表示基）、F6、K5（路径积分↔路径追踪类比）、M6（DEC）

**任务**

1. `math/fourier.hpp`：一维/二维 DFT（radix-2 + 朴素版对拍）。
2. 球谐函数：实 SH 基计算 + 辐照度系数投影/重建，渲染成光照贴图。
3. Schrödinger 波包演化（split-step 或显式 + 隐式）：`Complex<T>` 与
   `Dual<Complex<T>>` 的实战收官（自动求导验证波函数对参数的梯度）。
4. 网格曲率（高斯/平均）离散化 + 着色可视化（了解层，配合 H1）。

**验收节点**

- **V6.1** DFT vs 朴素 O(n^2) 变换一致；Parseval 定理数值验证。
- **V6.2** SH 重建低频环境光与直接蒙特卡洛半球积分结果吻合（低频段误差 < 5%）。
- **V6.3** 波包自由演化色散符合解析高斯波包公式；概率守恒 ∫|ψ|^2 = 1 漂移 < 1e-6。
- **V6.4** 球面离散高斯曲率积分 ≈ 4π（Gauss-Bonnet 数值验证，拓扑不变量的震撼时刻）。

**不过时先查什么**

- Parseval 不成立 → 先查归一化约定（`1/N` 放在正变换还是逆变换）与 `2π` 因子。
- ∫|ψ|^2 漂移 → 先查 split-step 是否用了非幺正的算子分裂顺序（半-全-半 而非 全-全）。

**测试需完善**

- SH 正交性数值检验（离散球面积分 ⟨Yl,Ym⟩ ≈ δ）。
- 复数特征分解（Hermitian）在量子 demo 中回归测试。

---

## 暂缓清单（不做承诺，遇到回读）

| 主题 | 处置 |
|---|---|
| A1 测度论/勒贝格积分、A2 泛函分析 | 阶段 4「CG = 能量下降」、阶段 5「1/√N」处间接学，不单独推进 |
| B2 流形/纤维丛、微分形式严格化 | 仅阶段 4/6 用其计算面（外微分 = 矩阵） |
| E4 规范场、F 量子公理化、L 费曼图 | 了解层，读笔记不写代码 |
| I5–I8 湍流/激波多相流 | 阶段 4 的 V4.5 为限 |
| `Matrix<T,M,N>` 固定尺寸特化、SIMD、GPU | 阶段 5 渲染 profiling 出现真实瓶颈后再立项 |
| `Tensor<T,rank>` | 仅当出现展平无法处理的场景 |

---

## 验收编号新旧对照（历史引用用）

| 旧编号 | 新编号 | 阶段 |
|---|---|---|
| A1–A7 | V1.1–V1.7 | 阶段 1 |
| B1–B7 | V2.1–V2.7 | 阶段 2 |
| C1–C4 | V3.1–V3.4 | 阶段 3 |
| D1–D5 | V4.1–V4.5 | 阶段 4 |
| E1–E5 | V5.1–V5.5 | 阶段 5 |
| F1–F4 | V6.1–V6.4 | 阶段 6 |

---

## 待办 / 已知问题（只追加，不删除）

- [ ] `graphics/CMakeLists.txt` 仍是 C++17，未升 C++20（阶段 1 任务 0）。
- [ ] `src/core/perspective.hpp`、`rotate.hpp` 为旧接口，阶段 1 重写 matrix 时一并处置。
- [ ] `demo/autodiff/compute_graph.cpp` 依赖旧 Matrix friend 接口，矩阵重写后需适配。
- [x] 仓库根目录的 `requirements`（无扩展名，内容为 pip 源说明）**已于 2026-09-22 删除** ——
  它与 `pyproject.toml` + `uv.lock` 重复，且会误导 agent 走 pip 路径。依赖一律以 uv 为准。
- [x] `notebook/` 文件名两处拼写问题**已于 2026-09-22 修正**：
  `B_algbra.tex` → `B_algebra.tex`、`M_addtion.tex` → `M_addition.tex`（`main.tex` 的 `\input` 已同步）。
- [ ] `§NN.N` 是编译后全局编号，**插入新 `.tex` 就会整体位移** —— 文档引用优先用文件内编号（`B1`、`G1.1`）。
