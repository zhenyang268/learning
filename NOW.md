# NOW.md — 近期要做什么（唯一滚动文件）

> **本文件是全仓唯一写「接下来」的文件。**
> **规则**：1–3 周一版；过期就**整篇重写**，不追加、不归档、不留历史。
> 历史只进 `docs/logs/<YYYY-MM>.md`（append-only，最新在上）。
> **不写阶段、不写月份、不写"长期目标"** —— 那些一律视为噪音。
>
> **本文件不重复别处的规范**：验收阈值 → 各 `demo/*.cpp` 头部「判据边界表」；
> 判据方法论 → `graphics/demo/README.md` §4；数值契约 → `graphics/numeric-contract.md`。
> 这里只写**做什么、什么算做完、不做什么、待定**（一条一句）。**不给先后顺序**（2026-10-30 用户定）。

**本期起点**：2026-10-29（第四次重写。原六条中 quat 已完成并入日志，故本期为五条 + 一条并行线）
**2026-10-30 修订**：用户定「内部先后由我自己并行安排」⇒ **本文不再给顺序意见**，条目编号只是标识。

---

## 本期五条（做完就重写本文件）

**编号不代表先后。** 这五条彼此无依赖（`ode` 吃 `Matrix` 列向量，不吃 `Tensor`），
先做哪个、并行几条由你自己安排。本文只写每条的「做什么 / 算做完 / 不做 / 待定」。
第 2–4 条落在**新建的** `graphics/src/math/`；第 5 条是纯数学，不阻塞任何代码。

---

### 1. 张量 `Tensor` / `TensorView` 实现

**前置更正（2026-10-28 实测，仍未变）**：`/tmp/tensor.full.{hpp,cpp}`、`.wip-tensor/`、
`graphics/demo/core/tensor.full.bak/` **均已空** —— 887 行完整实现已于 2026-10-25 按用户指示删除。
**只能从零实现，不要再找备份。**

- **做什么**：`graphics/src/core/tensor.hpp`（框架：`throw` 语句 **35 处**、`待实现` 63 处，2026-10-30 实测，**只数语句不数注释**）按
  ① 几何工具 → ② 构造与索引 → ③ 视图 → ④ 深拷贝 → ⑤ 逐元素 顺序落地；
  `demo/core/tensor.cpp` 的 `[P1]`–`[P24]` 占位判据**逐条换真判据**
  （真判据规格写在同一文件末尾「升级路径」一节，照抄编号即可）。
- **算做完**：ctest `tensor` 从「断言 throw」变真判据全绿、`[G]` 契约段打开、
  从零构建 **0 告警**；且**两条注入对照**做过了（方案写在 `graphics/demo/core/tensor.cpp` 末尾，本文件不重复）。
- **开工前必须定三个待定**（不定就返工；`graphics/src/core/tensor.hpp` 头注释留的尾巴）：
  1. **`strides_` 是否收敛为单一真源**
     (a) `Tensor` 不存 strides，`Strides()` 每次现算；
     (b) 保留成员，但只留**一条写入路径**（构造时 `ComputeStrides(shape_)`，不给 setter）。
     约束：`Strides()` 现在的签名返回 `const vector&`；且 `operator()` 必须与 view 共用
     `TensorOffset`（索引算术只准一处实现）。
  2. **广播** (a) 显式自由函数 `BroadcastTo`（stride=0 视图，须**只读**）；
     (b) 进 `operator+`。
     约束：`graphics/src/core/tensor.hpp` 已把「形状不一致必须 throw，不静默广播」写成契约。
  3. **`shape_` / `strides_` 是否「内联容量 4 + 溢出才堆」** (a) 做；(b) 不做。
     约束：头注释把「单线程、最简单版本」定为**设计目标，不是待补的欠债**。
- **不做**：指标收缩 / 乘法 / einsum **一律不进这个类**（是 `src/math/` 的自由函数）；
  不做编译期 rank（会把「不同 rank 的 view 互转」堵死）；不做并发 / 引用计数 / COW / 表达式模板。
- **注意**：设计口径的真源是 `graphics/src/core/tensor.hpp` 头注释 + `graphics/acceptance.md` §0
  （原 `docs/matrix.md` 已于 2026-10-28 删除，不要按旧路径找）。

---

### 2. `math/ode.hpp` —— Euler / RK4 / Verlet·辛欧拉

**为什么不能延后**（四份文档唯一交集）：
① **Gate A 的最后一块**（`dev-progress.md` §2）；
② 变分段（段 5）的**曲线出口**（`docs/route-math.md` §7）；
③ **CAE 60 #3 + HPC 60 #4** 的共同要求（`docs/eng-competency.md`）；
④ 数学物理工程 60 #1（`docs/baseline-levels.md` §2.2）。
**没有任何其它件能做到「一份时间喂四条线」。**

**起手成本（2026-10-30 实测更正）**：`graphics` 是 header-only，`src/` 只是 include 路径 ——
新建 `graphics/src/math/` 不需要动根 `CMakeLists.txt`，只加 `demo/math/CMakeLists.txt` +
`demo/CMakeLists.txt` 一行 `add_subdirectory`，代价约等于两个小文件。
（本文不据此给顺序意见 —— 先后由你自己安排。）

- **做什么**：本件是 `graphics/src/math/` 的**第一件**。
  状态用 `Matrix<T>` 列向量（D19）；`f` 签名 `Matrix<T> f(real_t, const Matrix<T>&)`。
  三入口：`Euler`（一阶对照）、`RK4`（主力）、`Verlet`/辛欧拉（保结构）。
- **算做完**：收敛阶测出来了（`h` 减半时 RK4 误差降约 16 倍、Verlet 约 4 倍）；
  与 `graphics/demo/pyapps/sphere/_core.py` 的手写 RK4 **逐点**对拍；
  保结构判据（辛格式能量有界振荡 vs 非辛单调漂移）两条曲线同图。
- **不做**：自适应步长（RK45）、隐式方法（BDF/Radau）、事件检测 —— 都不在段 5 消费点上。

---

### 3. `math/fourier.hpp` —— DFT / FFT / 卷积 / 傅里叶级数

**为什么做**：① 它是数学**段 2（L^2 与傅里叶，枢纽）的工程出口** ——
`docs/route-math.md` §7 直接点名 `math/fourier.hpp`（卷积定理 `F(f*g) = F(f)F(g)`、Parseval）；
② `docs/eng-competency.md` 现状表里「FFT / 卷积」一栏是 `无`，这是明确缺口；
③ 它和 `ode.hpp` 一起把 `src/math/` 的目录能力搭起来。

- **做什么**（四入口）：
  1. 朴素 `DFT`（O(N^2)，**只作对照基线**，不做优化）
  2. `FFT`（radix-2，迭代式位反转）
  3. `Convolve` —— **三条路径**：时域直接算 / 频域乘法（补零到 >= N+M-1）/
     **循环矩阵乘**（复用 `linalg`，这是唯一跨模块的独立路径，比 DFT 对拍 FFT 强得多）
  4. 傅里叶级数：`FourierCoefficients`（系数投影）+ `Reconstruct`（截断重构）
- **共享原语（不许各写一份）**：逐点复乘 + 归约求和。`Convolve` = 逐点乘后接逆变换；
  Parseval / 内积 = 逐点乘后接归约。**内积不在本文件实现**（`linalg` 的 `Dot`/`Adjoint` 已覆盖），
  只复用 —— 同一算术只准一处实现。
- **算做完**：N = 1024 下 FFT vs 朴素 DFT **复数逐点**相对误差 <= 1e-12；
  卷积**三条路径**逐点一致 + Parseval（**带 `1/N`**）<= 1e-12；
  方波重构的 Gibbs 过冲约 8.95% 且**不随 N 减小**（只向跳变点收缩）；
  N = 2^20 的计时比报出 O(N log N) vs O(N^2) 的量级差。
- **不做**：任意长度（混合基 / Bluestein）、多维 FFT、实数 rfft、窗函数、频谱估计。
- **注意**：必须能吃 `Complex<T>`（库已有 `core/complex.hpp`）；实数 `T` 域也要能编过。
  阈值均为**初判据，实现后按实测收紧**。
- **符号与归一化已定案**：`graphics/numeric-contract.md` §11（**D20**）。
  归一化跟 `numpy.fft`（正向不归一化、逆带 `1/N`），**不选酉变换**；
  高斯不动点用 `exp(-pi*x^2)`，其阈值必须实测后填（连续恒等式，离散只是近似）。
- **卷积不进 `Einsum`**：它带平移核，不是纯下标收缩；要进得先 Toeplitz 展开成 N×N，
  那是另一个东西。这条也补上了「einsum 薄层无消费者」的理由。

---

### 4. `math/` —— erf / 正态分布 / Box–Muller 采样

**为什么做**：段 1「概率三件」的工程出口 —— `docs/route-math.md` §7 点名
`math/mc.hpp`（采样 / CLT）；`docs/eng-competency.md` §4 把「随机采样 / 蒙特卡洛」
列为 HPC 赛道检查点。

- **做什么**（三件）：
  1. `Erf(x)` / `Erfc(x)` + `NormalCdf(x)`（= `(1 + erf(x/sqrt2))/2`）
  2. `NormalQuantile(p)`（Φ 的逆）—— 二分反解或 Acklam 近似，二选一
  3. `BoxMuller` 采样（**自带可复现种子**）
- **算做完**：与 Python 标准库 `math.erf` 逐点对拍（**不需要 scipy**），最大绝对误差进日志；
  `Phi(NormalQuantile(p)) - p` <= 1e-12（p 取 [1e-9, 1-1e-9]）；
  N = 10^6 的样本均值/方差偏差落在 3/sqrt(N) 量级内。
- **不做**：多元正态、Cholesky 相关采样、拒绝采样族、KS 检验（只做矩对照）。
- **落点已定（2026-10-30）**：按「**复值 vs 实值**」切，不按主题切 ——
  `math/special.hpp` 放 `Erf`/`Erfc`/`NormalCdf`/`NormalQuantile`（**实值纯函数**，
  不许依赖 `Complex`）；`math/prob.hpp` 放**正态特征函数** `phi_X(t)` 与 `BoxMuller` 采样
  （复值，属傅里叶域，符号见 D20）。
- **不依赖 `autodiff`**（2026-10-30 定）：矩的交叉验证 `phi^{(k)}(0) = i^k * E[X^k]`
  **放在 testcase 里**做 —— 用例同时 include `autodiff/dual.hpp` 与 `prob.hpp`，
  用 `Dual` 求 `phi` 在 0 点的导数，与解析矩对拍。**库际零依赖，验证力一点不减。**
- **跨文件三角判据**：两个高斯密度的 `Convolve` 结果，与 `phi_X * phi_Y` 的逆变换对拍 ——
  同时压 `fourier` 与 `prob` 两侧，任一边符号错都会红。这也是「一份时间喂两条线」的落点。

---

### 5. 测度第 3 关收尾（数学，无代码）

- **做什么**：Carathéodory 可测性判据 → Vitali 不可测集构造 → Lebesgue 积分四步构造。
- **算做完**：能讲清「外测度 → 限制成测度」为什么需要 Carathéodory 判据；
  Vitali 集为什么必须用选择公理。允许跳证明，**不允许跳定义**（`docs/route-math.md` §1）。
- **失败模式**：把 Carathéodory 判据当技术条件背下来 —— 讲不出「外测度不满足可数可加性，
  所以必须先限制到可测集族」；或说 Vitali 集「取代表元没有构造性」却不提**选择公理**。
- **位置**：明细记 `docs/math-log.md`。
- **排期性质**：不阻塞任何代码，也不进任何赛道；工程线卡住时它随时可做，反之不然。

---

## 并行线：测试框架迁移（已完成，2026-10-30）

基础设施：依赖脚本 `graphics/scripts/fetch_deps.sh`、gtest 接线、
`./graphics/run_tests.sh --coverage`、pytest 架子（细节见 `docs/logs/2026-10.md` 2026-10-30 条）。

**17 个 ctest 用例里 16 个已是 gtest**，逐个做过注入对照（明细同上日志）。
唯一还是 legacy 的是 **`core/tensor`** —— 它**迁移与实现合并做**：
24 条 `[P]` 是「断言仍 throw」的占位，先迁到 gtest 再实现等于写两遍。

迁移中新发现的**契约缺口**（要不要补由你定，均已登记 `graphics/test-map.md` §4.3）：
`expm_pade` 的非方阵/空矩阵 throw、`inverse` 的 `diagonal_inverse` 近奇异与三个
`logic_error` 分支、`eigen` 的 `eig_symmetric` 守卫与 Jacobi 迭代上限、`svd` 的内部错误、
`jacobian` 3 参数版的两个守卫、`grad_graph` 里 `Cos`/`Log` 的局部导数分支。

---

## 并行线：微基准仓（`learning/learning-benchmark/`）

> **定位（2026-10-28 实测）**：该仓**就在父仓内**（不是 `~/develop/learning-benchmark/`），
> 独立 git 仓、父仓 `.gitignore` 第 71 行已登记、父仓不跟踪。
> 真源 = 该仓 `docs/bench-plan.md`（排期）+ `docs/bench-workflow.md`（两机分工）。
> **本仓不重复记它的进度**，只记「它在哪、卡在哪」。

**它的下一步**（该仓 `bench-plan.md` §4）：M1 的 **T1.1–T1.6** ——
在 **WSL2**（Remote SSH `192.168.50.15 -p 2222`）写：
`harness/env` · `harness/timing` · `cases/gemm_ref_f64` + `tolerance` ·
`cases/gemm_generic` · `harness/main` · `tests/meta_schema`。
**前置阻塞已全部清除。**

**两条已定案的实测结论**（勿再排查）：

| 项 | 结论 |
|---|---|
| **硬件 PMU** | 不可用（`unsupported p6 CPU model 183 no PMU driver`）⇒ 归因**只走唯一变量对照**，`attribution.method` 永不为 `"perf"` |
| **绑核** | 不绑核（WSL2 拓扑是合成的：`lscpu` 报 10 核×2 线程，与真实 6P+8E 不符）⇒ `bind.kind` 恒 `none`，多线程数据 `reliable=false` + `caveats` 必填 |

**与 `docs/baseline-levels.md` 的联动**：因不绑核，**HPC 赛道 60 #3（性能分析）**只能靠
**单线程数据 + 对照实验**满足；**多线程加速比不许写进简历**。

---

## 本期不做（明确挂起）

| 不做 | 理由 |
|---|---|
| 稀疏 CSR / SpMV / CG | CAE 入口，等 `ode.hpp` 与 `Tensor` 之后 |
| Gauss–Legendre 求积 | **零行代码**。注意库里三个「Gauss」不是一回事（`docs/eng-competency.md` §3 第 2 条） |
| 渲染（BRDF/路径追踪/BVH）、几何（NURBS/半边） | 属其它赛道 |
| GR 计算、G 系 C++ 化 | 依赖张量实现 |
| `einsum` 薄层 | 无消费者（张量侧的收缩产品还没定） |

---

## 下一期候选（无序，按需挑，不承诺）

| 候选 | 触发条件 | 一次投入喂几条线 |
|---|---|---|
| **稀疏 CSR + SpMV + CG + 1D Poisson FEM + Gauss–Legendre** | 本期五条完成后 | CAE 60 #1/#2 + HPC 60 #1/#2 |
| **`math/pde.hpp` + `LaplacianMatrix` + CG** | 想开 Gate B（PDE 可视化） | CAE + 数学段 4–5 |
| **微基准仓 M2**（变体铺开） | 该仓 G1 闸门过后 | HPC 求职主线 |
| **李群刀 4–6**（exp 非单射 / SU(2) 双覆盖 / BCH 测阶） | 纯线代，成本最低 | 只在基线里，不喂赛道 |

---

## 待你确认（我查得到、但不能替你决定）

1. **岗位抽样**（`docs/eng-competency.md` §3 第 3 条）：取 20 个真实 JD，
   只记三列 ——「学历要求 / 是否明写稀疏·预条件·MPI·GMRES·BRDF·BVH / 是否要求专业对口」。
   **这是市场事实，agent 判断替代不了。**

2. **本期还剩三个契约空白**（未定前对应条目不开工；论述在该条目内，这里只做清单）：
   - `Tensor` 三个待定：`strides_` 单一真源 / 广播是否进 `operator+` / 内联容量 4 —— 见「张量」条

   > 已于 2026-10-30 结清：erf / ppf 落点（`special.hpp` + `prob.hpp`）、
   > 傅里叶符号与归一化（`numeric-contract.md` D20）、`prob.hpp` 是否依赖 `autodiff`（不依赖）。

---

## 「做完」的判据模板（每个交付物都套）

> **条件 + 数字 + 失败模式**，三件缺一不算完。
> 来源：`docs/eng-competency.md` §0.2 —— **这是 60 分与「我跑通了」的唯一分界。**

- **条件**：机器、编译器与优化级、规模、容差、基线是什么
- **数字**：误差 / 收敛阶 / 迭代数 / 加速比，要能和基线对照
- **失败模式**：什么输入会让它错、错成什么样、你怎么知道的

反例：「CG 快 3.2×」→ 噪音。
正例：「N=10^5、float64、-O3、容差 1e-8，对比 Jacobi 预条件，迭代数 812 → 96，预条件耗时占比 7%」。