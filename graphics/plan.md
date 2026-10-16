# plan.md — 数学 × 图形学 学习实施计划

> **读者**：agent（执行工程任务）+ 人 ｜ **最后更新**：2026-10-16
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
  配套用例（对应 `matrix.md` §7.2 复数求导图景）：
  - **z² 单种子**（全纯）：derive 种子 = {1,0}，传播后 derive = 2z 的复数值；
  - **|z|² 双种子**（非全纯，z̄·z）：两个种子独立种（{1,0}、{0,1}），
    两次传播的组合出 2z̄，与解析值对拍；
  - 阈值：value 一致 1e-12、derive 一致 1e-10（double）。
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

- [x] `graphics/CMakeLists.txt` 仍是 C++17，未升 C++20（阶段 1 任务 0）。
      → **2026-10-13 核实：已升**（`set(CMAKE_CXX_STANDARD 20)` + `CXX_STANDARD_REQUIRED ON`），本条作废。
- [ ] `src/core/perspective.hpp`、`rotate.hpp` 为旧接口，阶段 1 重写 matrix 时一并处置。
- [x] `demo/autodiff/compute_graph.cpp` 依赖旧 Matrix friend 接口，矩阵重写后需适配。
      → **2026-10-13 完成**：真实待适配的不是 friend 接口，而是① `GradOp` 转 `enum class` 后
      `grad_graph.hpp` 的 35 处 `case` 标签缺作用域前缀；② 类内 CTAD `Matrix(1,1)` 须写 `Matrix<T>(1,1)`；
      ③ 缺 `Output`/`Mult` 分支（编过也输出恒 0）。`ctest --test-dir graphics/build` 1/1 通过。
- [x] 仓库根目录的 `requirements`（无扩展名，内容为 pip 源说明）**已于 2026-09-22 删除** ——
  它与 `pyproject.toml` + `uv.lock` 重复，且会误导 agent 走 pip 路径。依赖一律以 uv 为准。
- [x] **标量契约与设计文档漂移**（2026-10-13 修）：`scalar.hpp` 里写的是文档中不存在的
      `const_n(内层 T)`，而 §7.1 第 2 条要的是 `fromReal(real_t)` —— 后果是
      `Scalar<Complex<double>>` / `Scalar<Dual<float>>` 双双不成立，任何
      `requires CommutativeScalar<T>` 的 linalg 函数都拒绝它们（**V1.1 的直接拦路石**）。
      已改名并按设计补进 concept；`ScalarTraits<Dual<T>>::conj/abs2` 的成员名
      （`x.real`/`x.dual` → `x.value`/`x.derive`）一并修；`jacobian.hpp` 的 `Derive`
      改为按值收参（原先接不了 `Apply(...)` 的临时量）。
      新增 ctest 用例 **`scalar_traits`**（`demo/core/`）：静态断言钉契约 + 交换律单测
      + Complex/Dual 原函数对拍，44 项检查。<br>→ 注意：**这不等于 V1.1 通过**，
      V1.1 还要求 linalg（inverse/eigen/svd）真的落地并跑出残差。
- [ ] **V1 收尾体检（2026-10-14 探针实测）：V1.1–V1.7 一条都还没通过，不能收尾。**
      根因不是"还差几个函数"，而是 **linalg 从未被任何编译单元 include** —— `demo/` 只有
      `scalar_traits` 与 `compute_graph` 两个 ctest，528 行 `linalg/*.hpp` 一次都没进过编译器，
      所以"看着写完了"。`-fsyntax-only` 实例化探针实测（两轮共 32 个编译错误）：
  - [x] `solve.hpp`（**2026-10-14 修复，S1**）：`SwapRow` 定义 2 参而调用处传 3 参；函数体里用
    未声明的 `max_row`；`Tol` 未定义（改为随规模缩放的 eps 倍数 → `DefaultTol<T>(n)`）；
    `piv[j]` 应为 `result.piv[j]`；主元比较写裸 `>`（违反 D7 —— Complex/Dual 没有 `operator>`）。
    另有两处**编过也全错**的真 bug：① 消元整段被 `if (max_row != i)` 包住 —— 主元行恰好落在 i
    上时整段跳过，U 少消一行、L 少一个乘子，静默出错；② **换行时 L 的乘子列没有跟着换**
    （`r1`/`r2` 已带乘子时 `P·A != L·U`）。顺带定下运算符口径：Complex/Dual 无复合赋值，
    泛型 linalg 只能写 `a = a - b`。
    验收：新增 ctest 用例 **`solve`**（`demo/linalg/solve.cpp`，34 项检查）+ 测试基建
    **`demo/common.hpp`**（`MaxAbs2`/`NormInf`/`IsApprox`/实-复-Dual 随机矩阵/`Checker`）。
    **注入反例实测**：去掉乘子换行 → 15 条 FAIL；恢复 `if (max_row != i)` → 26 条 FAIL。
  - [x] `inverse.hpp`（**2026-10-14 修复，S2**）：`scalarTraits` 拼写（应为 `ScalarTraits`）；
    `Identity<T>(n)` 应为 `Matrix<T>::Identity(n)`；`diagonal_inverse` 里 `tol` 未声明；
    `Cholesky(A)`/`Solve` 均未 include；`Inverse` 只有声明、无定义。另有两处**编过会算错**的
    探测 bug：① `is_symmetric` 把比较对象写成 `AH(j,i)`（= conj(A(i,j))），实际判的是
    "A(i,j) 的虚部是否为 0"，实数矩阵恒真 → **所有实方阵都被判成 Symmetric**；
    ② `is_orthogonal` 只查对角线，`A^T·A` 对角全 1 但非对角非零的矩阵被误判成正交。
    另 `is_positive_definite`/`cholesky_inverse`/`qr_inverse` 三处依赖未接入的 cholesky/qr，
    改为**显式 throw**（不留静默 no-op），`Classify` 的 SPD 档暂时关闭。
    验收：新增 ctest 用例 **`inverse`**（`demo/linalg/inverse.cpp`，48 项检查）—— V1.2 硬指标
    ‖A·A⁻¹ − I‖∞ < 1e-8（double n=2..6 与 Complex n=2..4 共用同一份源码），实测最差
    double 1.18e-13 / complex 9.16e-16。
    **注入反例实测**：`is_symmetric` 还原成 `AH(j,i)` → 2 条 FAIL；`is_orthogonal` 还原成
    "只查对角" → 1 条 FAIL。
  - `eigen.hpp`：`eig` 调 `is_symmetric(A)` 缺 `tol` 实参；`Identity(n, n)` 只收 1 参；
    `ScalarTraits()::one()/zero()` 语法错；`jacobi` 定义在使用点之后（非 ADL 不可见）；
    **且实现是 2x2 单次旋转，不是 Jacobi 迭代**（无收敛循环、无非对角元判据）。
  - `qr.hpp`（不在阶段 1 清单，但已建文件）—— **2026-10-14 编译层 + 算法层全部修复（S6b）**：
    编译层 9 处：`real_t` 未定义、`aplha`/`Cows` 拼写、`QRComposition` 大小写、
    `ScalarTraits<T>::Zero/One` 大小写、`Identity` 多传一参、`atan2` 单参、`-=` 复合赋值、缺 `return`。
    **算法层（编译过也会错）**：① `QR_Householder` 分母用了**改前**的 ‖x‖²（应为改后 vᴴv，
    否则 H 不再对合、消不掉 x）、`+sqrt` 无符号选择（x₁<0 时灾难性相消）、`H.Transpose()` 复数域应为 `Adjoint`；
    ② `QR_givens` 循环写成"外层行 i、内层列 j>i"——**是正确结构的镜像**（应外层列 j、内层行 i>j），
    且 j 上界用 `A.Rows()` 而 R 只有 n 列 → **堆越界**（ASan 实证）；`atan2` 取 `real()` 丢虚部，复数域无解；
    ③ `QR_Shimit` 把投影系数塞进 Q（Q 应是正交基矩阵）、拿 R 的列当"已正交基"、无归一化、R 不更新 ——
    与 Gram-Schmidt 无关。**三个函数全部重写**，并统一 thin 口径（Q: m×n，R: n×n）。
    验收：新增 ctest **`qr`**（`demo/linalg/qr.cpp`，189 项）—— 见下文 S6b。
  - `operator.hpp` —— **2026-10-14 修复（S6b 前置）**：`Dot`/`Abs2` 是成员函数残骸
    （自由函数里用 `this->rows/cols/data`，`Abs2` 还写了大写成员 `matrix.Data`）。
    改为两参 `Dot(a,b) = Σ conj(aᵢ)·bᵢ`（**补上共轭** —— 原版没有，是第 4 个「double 免疫」点）、
    `Abs2` 返回 `real_t`。
  - `cholesky.hpp` 是**唯一**能通过实例化的新文件。
  - 空文件：`svd.hpp`（V1.4 必需）、`expm.hpp`（V1.6 必需）、`householder.hpp`；
    `dot.hpp` 已删（无引用残留）。
  - 本阶段要求的 `demo/common.hpp`（`IsApprox` / 实-复-Dual 随机矩阵生成）**不存在**。
  - 收尾步骤 **S1–S7**（2026-10-14 定稿：一步一文件一用例一提交，不并发开两个文件）：
    - **S1 `solve.hpp`**（LU + Solve）+ `demo/common.hpp` 测试基建 + `demo/linalg/solve.cpp`
      —— ✅ 已完成
    - **S2 `inverse.hpp`**（补 `Inverse` 分发 + 修 `scalarTraits`/`Identity<T>`/`tol`/`Cholesky`）
      + `demo/linalg/inverse.cpp` → **V1.2 硬验收** —— ✅ 已完成（2026-10-14）
      （SPD / QR 分支本刀显式 throw，待 S2a / S6 接入）
    - **S2a `cholesky.hpp` 分解侧** —— ✅ 已完成（2026-10-14，**方案甲**：删递归版、只留元素公式版）：
      修 `L(i,j) -=`/`L(i,j) /=`（Complex/Dual 无复合赋值，D7）、**两处 `L(·,·)*L(·,·)` 漏 `conj`**、
      `L(i,i) <= 0`（复数无 `<=`，改走新加 `ScalarTraits<T>::real`）、`sqrt` 未限定命名空间、
      错误码空矩阵改 throw（D12）；删半成品 `cholesky_r` / `cholesky_recursive`
      （递归分块版由用户后续自行补）；新增 `cholesky_solve`（V2.2 要"解同一方程组"）。
      配套结构性整理：`ScalarTraits<T>::real`（4 个特化 + concept 加一条契约）与
      `DefaultTol<T>(n)` 上移到 `core/scalar.hpp` 标量层（linalg 各文件不再为取一个容差互相 include）。
      验收：新增 ctest **`cholesky`**（`demo/linalg/cholesky.cpp`，98 项）—— V2.2 硬指标
      ‖X_chol − X_lu‖∞ < 1e-8 且各自残差 < 1e-10，重构 ‖L·Lᴴ − A‖∞ < 1e-8，实/复共用同一份源码。
      **注入反例**：两处漏 `conj` → **23 条 FAIL 全在 complex，double 毫发无损**（印证"实数是这个
      bug 的盲区"）；去掉正定性判据 → 7 条 FAIL 全在异常/SPD 判据。
    - **S2a-2 `inverse.hpp` SPD 档** —— ✅ 已完成（2026-10-14）：`is_positive_definite` 改走
      "cholesky 不选主元走通即真"（只吃 `runtime_error`，别的异常照常外传；非方阵直接答否）；
      `cholesky_inverse` 走 `cholesky_solve(A, I)`（前代 + 回代，复用已验收的分解路径，
      不另写一遍三角求逆）；`Classify` 的 SPD 档**嵌在 `is_symmetric` 分支内部** ——
      平铺成同级 if 会让它永远不可达（SPD ⊂ Symmetric，会被上一行先 return 掉）。
      验收：`inverse` 用例 48 → **78 项**（新增 G 段 30 项）全过；V1.2 硬指标最差
      double 1.18e-13 / complex 9.16e-16 不变；随机 SPD(n=2..6) 与 LU 路径逐元差 ≤ 1.11e-16，
      自动探测入口与 hint 结果一致。**注入反例**：① Classify 平铺（制造顺序陷阱）→ **8 条 FAIL**
      （全为 SPD 判据）；② 正定判据恒 true → **5 条 FAIL**（全在非正定一侧）。
    - **S3 `eigen.hpp`** 重写为**真 Jacobi 迭代**（原为 2×2 单次旋转）+ `demo/linalg/eigen.cpp`
      → **V1.3** —— ✅ 已完成（2026-10-14）：cyclic Jacobi，一轮轮扫全部 i<j，直到
      off(A)² = Σ_{i≠j}|A(i,j)|² <= tol²（tol = 8n·eps·‖A‖_F；酉相似保 Frobenius 范数，标尺只算一次）；
      每步只更新第 p、q 两行两列（O(n)，不是整体旋转的 O(n³)）；Hermitian 用酉阵
      U = [[c, -s·ph], [s·conj(ph), c]]（ph = z/|z| —— 复数域比实对称只多这一个相位，
      实对称下 ph = ±1 退化成普通实旋转）；对角元走解析式 (a+b)/2 ± hypot(a−b, 2|z|)/2，
      与 θ = 0.5·atan2(2|z|, a−b) 严格配套（取 θ+π/2 也能消零但会对调特征值）；
      100 轮不收敛则 throw（D12，不返不准的结果）。删掉 `rotate()` 整体旋转与硬编码 double 的
      `jacobi()`；`eig(A)` 补 tol 实参、非对称改显式 throw（原先调不存在的 `eig_general`）。
      验收：新增 ctest **`eigen`**（`demo/linalg/eigen.cpp`，119 项）—— V1.3 硬指标最差重构
      **double 1.78e-15 / complex 1.33e-15**（判据 1e-8）；Q 正交/酉 ≤ 1.4e-15；对角化 ≤ 1.1e-15；
      逐列特征对 ‖A·q_j − λ_j·q_j‖∞ ≤ 1.6e-15；迹不变 |Σλ − tr(A)| ≤ 2.2e-16。
      **注入反例**：① 只扫一轮 → **16 条 FAIL**（**n=2 免疫、n ≥ 3 全中** —— 正是"必须循环"
      的硬证据：2×2 只有一对非对角元，一次旋转就彻底清零）；② 漏 `conj` → **24 条 FAIL
      全在 complex、实对称 0 条**，连 Q 的酉性都一起坏（漏 conj 直接破坏 U 的酉性）。
    - S4 `svd.hpp`（经 A^H·A 特征分解组装）+ `demo/linalg/svd.cpp` → **V1.4** —— ✅ 已完成（2026-10-14）：
      原文件是 19 行空壳（只有 `SVDcomposition` 声明，无实现、无引用点）。`B = A^H·A` 走 `eig_symmetric`，
      `σ_i = √λ_i`（非负根 + 降序）、`u_i = A·v_i/σ_i`；经济型 U(m×r) / Σ(r×r) / V(n×r)，r = min(m,n)。
      σ ≈ 0 的列（秩亏）用 Gram-Schmidt 扫标准基补正交补（`linalg_detail::FillOrthonormalColumn`）。
      **σ 的精度不拖累重构**：`U·Σ·V^H = A·(V·V^H)`，σ 与 1/σ 精确相消 —— 第③条注入反例就是它的数值证据。
      判据分层：U 列正交判 1e-10（方向精度 eps·κ(A)，实测最差 2.456e-13）、V 判 1e-12（实测 ~1e-16）。
      验收：新增 ctest **`svd`**（`demo/linalg/svd.cpp`，266 项）—— V1.4 最差重构
      **double 1.332e-15 / complex 7.772e-16**（判据 1e-8）；`Σσ² = ‖A‖_F²` 相对差 ≤ 2.9e-16；
      `σ_i²` 与 `eig(AᴴA)` 谱逐项吻合 ≤ 1.9e-16；覆盖实/复 × 方阵/长方阵 + 秩亏/零矩阵/1×1/1×2 + 解析对拍。
      **注入反例**：① `Adjoint` → `Transpose` → **34 条 FAIL 全在 complex、实数域 0 条**（实矩阵下二者等价、
      天然免疫 —— 与 S2a/S3 漏 conj 第三次同构）；② 不排序 → 5 条 FAIL 且**重构一条不挂**，
      **并暴露判据覆盖漏洞**（原对角用例 `diag(3,−2)` 本就降序）→ 已补「升序输入」用例
      （`diag(1,5,3)` 实 / `diag(1,4,2)` 复）；③ **σ 不开方 → 102 条 FAIL 但「重构」0 条挂** ——
      证明 V1.4 的字面判据「重构 < 1e-8」单独抓不到 σ 算错，必须有独立判据兜底。
    - S5 `expm.hpp`（0 字节 → 落地）+ `demo/linalg/expm.cpp` → **V1.6**（so(3) 对易）
      **（上层刀序见 `learning-todo/lie-plan.md` 的八刀路线：定义 → 李括号 → 闭式 → 非单射 →
      SU(2) 双覆盖 → BCH → 极分解 → 幺正演化。S5 = 刀 1，本文只保留本刀的 V1.6 阈值；
      刀 2–8 是扇面加餐，不新增 V 节点。）**
    - S6 `operator.hpp` / `qr.hpp` / `householder.hpp` 派支线处置决定 + `demo/linalg/smoke.cpp`
      （每个公开入口在 double 与 Complex 下各实例化一次 —— 同类问题防复发闸门）
      → **✅ 2026-10-16 全部完成**。`smoke.cpp` 见下 V1.1 条；三条派支线定案 ——
      `Dot` 由 smoke 收编（不再零调用方）、`MatrixType` **用户定案：挪到应用层 `demo/common.hpp`**
      （它按用途划分，「怎么用是上层的责任」；core 只留按**代数结构**分的 `MatrixKind`）、
      `qr_inverse` 补回代（= S-App6 ✅）。
    - **S6b `qr.hpp` 三函数修复（2026-10-14，用户选「三个都修」）** —— ✅ 已完成：
      `operator.hpp` 的 `Dot`/`Abs2` 修正 + `qr.hpp` 编译层 9 处 + **三个函数算法层重写**
      （householder：分母/相位/`Adjoint`；givens：循环结构反了 + 堆越界；shimit：重写为 MGS）。
      统一 thin QR 口径（Q: m×n 列正交，R: n×n 上三角，`A = Q·R`；`m < n` 三个入口统一 throw）。
      验收：ctest **8/8**，`qr` 用例 **189 项 / FAIL 0**；V2.1 最差重构 **9.93e-16**、
      Q 列正交 **8.88e-16**（判据分别 1e-8 / 1e-12）。判据三条：重构 ‖QR−A‖、列正交 ‖QᴴQ−I‖、
      R 严格下三角 < 1e-14（浮点残渣必然非零，判"== 0"是判错了）。
      **注入反例 4 条**：① householder 分母退回改前 ‖x‖² → **35 条 FAIL 全在 householder**
      （givens/mgs 0 条）；② givens `Adjoint`→`Transpose` → **3 条全在 complex、实数域 0 条**
      （第 4 次同构的「double 免疫」）；③ MGS 内积漏 `conj` → **3 条全在 complex**（Q 酉性崩）；
      ④ 去掉相位/符号选择 → 5 条，**首轮实数域 0 条 → 暴露判据覆盖漏洞**（实数用例 x₁ 全为正）
      → 补「病态 x₁<0」用例（`[[-1,0],[1e-6,1]]` 实 / 复）后重跑，**实数域 0 → 1 条**，覆盖补上。
    - **V1.5 `Dual<Complex>` 验收：`demo/linalg/dual_complex.cpp`（2026-10-16）** —— ✅ 已完成：
      `demo/common.hpp` 补 `RandFill<Dual<Complex<U>>>` 特化 —— 原来复用 `RandFill<Dual<U>>`
      写 `static_cast<U>(a)`，`U = Complex<U'>` 时 `Complex` 无单参构造/无默认参数 ⇒ 编译失败；
      偏序上 `Dual<Complex<U>>` 更特化，两条并存即可。新增 ctest **`dual_complex`**（22 项 / FAIL 0）。
      **V1.5 三条**：① `Matrix<Dual<Complex<double>>>` 零修改可用（含矩阵层 `d(A²) = A·A' + A'·A`）；
      ② derive 置零时 `Solve` 的 value 与纯 Complex **精确 0**（n=2/4/6）；③ z² 单种子
      （两种子 `derive/s` 一致 = 2z₀）/ |z|² 双种子（两方向 `derive/s` 不同、重组出 2z̄，残差 0）。
      另加 E 段「实用形态」：解 `A(θ)x = b` 时 `dx/dθ` 由 Dual 自动传播，与解析式 −A⁻¹A₁x
      （5.2e-18）、中心差分（6.9e-12）三方对拍 —— 这是 V1.5「为波包演化铺路」的真实形状。
      **注入 3 条**：① `ScalarTraits<Dual<T>>::conj` 漏 `conj(derive)` → 1 条 FAIL（残差 2.0，
      与手算一致，其余 8 用例不受影响）；② `Dual::operator*` 漏 `derive*other.value` → 4 条 FAIL
      **全在复数域**（同时 `scalar_traits`/`compute_graph` 挂 —— 那两个测的是实数域乘法，本刀补的
      正是复数域覆盖）；③ 把 D 段的「非全纯」换成 `z*z` → 2 条 FAIL，**暴露我自己的判据漏洞**：
      原判据比的是**原始** `derive`，而全纯的两种子原始值也恰好不同（2z₀ vs 2z₀·i）⇒ 漏判；
      改为比归一化后的 `derive/s` 后 2 → **3 条 FAIL**，漏洞补上。验收：ctest **9/9**。
    - **V1.1 显式闸门：`demo/linalg/smoke.cpp`（2026-10-16）** —— ✅ 已完成：
      V1.1 要的不是「某个入口算得对」（那是 V1.2–V1.6 各自的事），而是
      「**每个公开入口都在 double 与 Complex 下各实例化过一次**」。两段：
      **[A] 实例化闸门**对 **40 个公开入口**（27 自由函数 + 13 个 `Matrix<T>` 成员，
      `Inverse` 两个重载须 `static_cast` 消歧）**取函数指针**（odr-use ⇒ 强制实例化）、
      刻意不调用，配 `kEntryCount = 40` 计数（漏抄一行即报警）；**[B] 最小数值闸门**
      每个入口跑一次 + 一条轻量判据（阈值 1e-10，专项用例才是权威）。
      验收：ctest **`smoke`**，**54 项 / FAIL 0**，`-c -w` 从零 + `-Wall -Wextra` **0 告警**。
      **注入反例两条**（本条最值钱的证据）：① `Dot` 去掉 `conj` → **全项目只有 smoke 挂**
      （1 条，complex，残差 2.048e-01），其余 10 个用例全绿（`Dot` 此前零调用方，旧体系看不见）；
      ② 往 `Dot` 体内写只在实例化时非法的表达式 → **只 include 的探针编得过**、
      **其余 10 个 target 照常构建**，只有 `smoke.cpp` 编不过 —— 两面钉死「include ≠ 实例化」。
      副产物：`Dot` 从此有了第一个调用方（不再零调用方）；`qr_inverse` 当日随后补上回代（S-App6 ✅），
      闸门里「如实断言它 throw」的临时判据随之换成数值判据 [B6d] + 奇异抛 [B6e]。
    - S7 文档收尾：V1.1–V1.7 勾选 + dev-progress §3/§4/§6 + 记忆
- [ ] **应用扇面（S-App 系列，2026-10-14 用户定调「先铺开再收束」）**：
      主线顺序不变（S5 → V1.5/V1.7 → V1 收口 → 阶段 2/3），**不插入 S 序列**；
      扇面件与主线**并行**、每件半天到一天、共用已验收的 linalg 零件。
      规则：扇面只加 demo 与薄封装（≤ ~60 行），**不新建库文件、不新增 V 验收节点**
      （plan 红线：验收阈值只跟主线；扇面的「验收」= demo 自身的打印/对拍断言）。
      - [ ] **S-App1 PCA demo**（`demo/linalg/pca.cpp`）：生成二维正态点云（Box–Muller）→
            协方差矩阵 → `eig_symmetric` 对角化 → 主轴画 ASCII 散点或导 CSV 给 matplotlib。
            概念收获：多元正态 N(μ,Σ) 等高线=椭圆、主轴=特征向量、半轴=√λ；中心化先行。
            数学正课仍在 C6–C10（CLT/特征函数），demo 只是几何直觉。
      - [ ] **S-App2 牛顿法 demo**（`demo/linalg/newton.cpp`）：`jacobian` + `Solve` 拼装，
            一个二维非线性方程组 + 收敛阶打印（探针已验：二次收敛 1.8→…→1.1e-15）。
            概念收获：多元泰勒一阶截断 = 线性化；J 奇异时的行为。
            正课在阶段 3 隐式积分器（van der Pol），此处见「它是什么」。
      - [ ] **S-App3 伪逆 + 左逆右逆**（`svd.hpp` 薄封装 `Pseudoinverse(A, tol)`，~15 行）：
            列满秩 A⁺A=I（左逆）、行满秩 AA⁺=I（右逆）、秩亏取截断 Σ⁺。
            **数值秩判据用 `√eps·σ_max` 或 `1e-6·σ_max`，不用 `DefaultTol`**
            （2026-10-14 探针实测：8n·eps 口径对 7×5 真秩 4 矩阵只判对 73/100，
            √eps 口径 100/100 —— `DefaultTol` 是唯一性判据，不是数值秩判据）。
      - [ ] **S-App4 最小二乘 `lstsq`**（`qr.hpp` 薄封装，QR 直解 R⁻¹Qᴴb，~20 行）：
            直线/二次多项式拟合 + 残差打印。**默认 QR 不走 SVD**（本库实测：κ=1e8 时
            正规方程爆 1e9、QR 1.3e-8、SVD 4.2e-1 —— A^H·A 路线条件数平方）。
            κ 对比一节正是「正规方程为什么被弃用」的教科书一课。
      - [ ] **S-App5 旋转扇面**（等 S5 expm 落地后做，天然衔接 V1.6/V2.4）：
            Rodrigues vs `expm(t·[u]×)` 对拍（探针已验 3.3e-16）+ 欧拉角万向锁复现 +
            四元数预演（阶段 2 正课）。概念收获：李代数=切空间线性化、exp=切空间→流形。
      - [x] **S-App6 `qr_inverse` 补回代** —— ✅ 2026-10-16 完成（`inverse.hpp`，~25 行含注释）：
            `A = Q·R ⇒ A⁻¹ = R⁻¹Qᴴ`；R 已是上三角，**直接回代**求 R⁻¹（不复用 LU、不显式求逆），
            奇异判据 `|R_ii| ≤ DefaultTol` 抛。验收：`inverse` 用例 **[H] 段**（+8 项 → 85 项 FAIL 0）
            + `smoke` **[B6d/B6e]**。实测残差 double ≤ **1.11e-15** / complex **2.48e-16**，
            与 `lu_inverse` 逐元差 ≤ **6.2e-15**。注入「`Qᴴ`→`Qᵀ`」→ **只挂复数域**（第 5 次 double 免疫）。
      收束条件（何时回主线全时）：V1.5/V1.7 + 上述扇面全部提交后 → S7 文档收尾 →
      **扇面清零、主线独占**，进阶段 2。扇面不做完不阻塞阶段 2（App5 依赖 S5，其余随时可插）。
- [x] **✅ 反向模式（`ComputationGraph::Backward` / `GradNode::Backward`）已落地**（2026-10-17）：
      上一版（2026-10-13）把它标成「待立项功能、非 bug」是**判断错了** —— 实测发现那 130 行不是
      「草稿未验证」，而是**有 4 类确定性 bug**（22 个 GradOp 里只有 5 个伴随公式正确）：
      ① `case Mult` / `case Output` **缺失** ⇒ 静默 no-op（Forward 里 Mult ≡ CWiseMult，
      但 Backward 只写了 CWiseMult）；② `case Transpose` 用 `node1->matrix` 而非上游伴随、且用
      `=` 而非 `+=`；③ `case Divide` / `Function` / `MatrixFunction` 把**逐元素**写成**矩阵乘**；
      ④ `case MultAll` 维度不匹配 + 偏导藏在 Dual 的 `derive` 分量里没提取。
      更根本的是 `adjoint` **从不按形状初始化** ⇒ 第一个节点就抛 `Matrix dimension mismatch`
      —— 所以 `ComputationGraph::Backward` 的 `throw` 是「知情拒绝」，但既然知道了就该初始化。
      **现修**：`ComputationGraph::Backward()` 做「形状初始化 → 输出播种 → 拓扑逆序回传」，
      并补 5 个 case（`-Wswitch` 归零）＋ 修上述 4 处公式 ＋ `jacobian.hpp` 加
      `localFunc`/`LocalDeriv`（局部导数**与种子解耦**：`Derive(Apply(...))` 得到的是
      f'(x)·种子，不是 f'(x)）。**新增 ctest `backward`（29 项，含最小闭环三方对拍）**。
      **仍未动**：8 个 `throw` 占位（CWisePow/Inverse/Sigmoid/Relu/Tanh/Softmax/CrossEntropy/
      LnMultAll）与 4 个激活函数的**空 body** Backward（`[E3]` 已如实断言它仍静默 0）。
      **关于 tape**：§7.5 第 4 条的「独立文件」仍是终态方向，但**不是现在必须做** ——
      「顺手存局部导数」是工程惯例（省一次 O(n) 计算），不是正确性要求；
      本架构 `shared_ptr` 持有 `matrix` ⇒ 中间量一直在 ⇒ 永远可重算。
      详见 `dev-progress.md` §6 的 2026-10-17 条。
- [x] `notebook/` 文件名两处拼写问题**已于 2026-09-22 修正**：
  `B_algbra.tex` → `B_algebra.tex`、`M_addtion.tex` → `M_addition.tex`（`main.tex` 的 `\input` 已同步）。
- [ ] `§NN.N` 是编译后全局编号，**插入新 `.tex` 就会整体位移** —— 文档引用优先用文件内编号（`B1`、`G1.1`）。
- [x] **✅ autodiff：`Apply` 的域契约已写进签名**（2026-10-16 落地，用户令「你改吧」）：
      按域拆成 **4 个重载** —— Dual 专用 ×2（签名 `Matrix<Dual<T>>`，走 `dualFunc`）
      + 标量专用 ×2（签名 `Matrix<T>`，走新增的 `elemFunc`）。两者共存靠
      [temp.deduct.partial] **部分排序**分流（`Matrix<Dual<U>>` 比 `Matrix<T>` 更特化，
      标准行为、非编译器扩展）；标量版体内 `static_assert(!IsDual<T>::value, ...)` 作守卫
      （域判别 trait `IsDual` 落在 `dual.hpp`）。原报错
      `error: use of undeclared identifier 'dualFunc'` 的根因是 `dualFunc` 为 `Dual` 的
      **in-class friend**（ADL-only、无法前置声明），而签名贪心写成了任意 `Matrix<T>`。
      **新增 ctest `apply`**（`demo/autodiff/apply.cpp`，17 项 FAIL 0）：4 重载 × 4 域
      （double / Complex / Dual / Dual&lt;Complex&gt;）+ `GradOp::Function` / `MatrixFunction`
      节点真实路径（此前**从未跑过**）+ 嵌套对偶二阶导 + 维度不匹配 throw。
      **顺带修** `dual.hpp` 的 `ln`/`log` ADL 断链（`ln` 体内用 `log`、而友元叫 `ln` ⇒
      嵌套对偶 = 二阶导时 `no matching function for call to 'log'`）。
      四条注入对照（含偏序证据）见 `dev-progress.md` §6 的 2026-10-16 条。
      **验证**：ctest **12/12**、`-c -w`（从零 + `-Wall -Wextra`）0 告警。**已了结。**
