# NOW.md — 近期要做什么（唯一滚动文件）

> **本文件是全仓唯一写「接下来」的文件。**
> **规则**：1–3 周一版；过期就**整篇重写**，不追加、不归档、不留历史。
> 历史只进 `docs/logs/<YYYY-MM>.md`（append-only，最新在上）。
> **不写阶段、不写月份、不写"长期目标"** —— 那些一律视为噪音。

**本期起点**：2026-10-28（第三次重写。本期由三条扩为六条，加入 quat / fourier / erf；
前两条为骨架填空，后四条含两件新建目录）

---

## 本期六条（做完就重写本文件）

顺序：quat → Tensor → ode → fourier → erf → 测度。
前两条**骨架与判据都已写全**，开工即填空；第 3–5 条落在新建的 `graphics/src/math/`；
第 6 条是纯数学，沉底是因为它不阻塞任何代码。

### 1. `quat/quaternion.hpp` —— 四元数与三维旋转（10 处 `throw` 清零）

**为什么排第一**：本期唯一一个「骨架 + 占位判据 + Python 对拍源全齐」的件 ——
骨架 `src/quat/quaternion.hpp`，占位判据 `[P2]`–`[P10]` 在 `demo/quat/quaternion.cpp`，
对拍源 `demo/pyapps/rotation/`（A5.1–A5.6 + `report.md`）已跑过并留下数字。
**开工即填空，不需要先做设计决策**，所以放在最前面当热手。
赛道出处：`docs/eng-competency.md` CG 60 #1（变换与空间：矩阵 / 四元数 / SLERP）。

- **做什么**：10 处占位按依赖序清零，四步，每步单独跑一遍：
  1. `Inverse` / `Normalize`（`[P2][P3]`）
  2. `FromAxisAngle` / `ToMatrix` / `Rotate`（`[P4][P5][P7]`）
  3. `FromMatrix`（`[P6]`）—— **必须四分支 Shepperd**：按 `tr(R)` 与三个对角元中最大者选分支
  4. `AngleOf` / `Slerp` / `Nlerp`（`[P8][P9][P10]`）
  每步：实现 → 对应 `[P]` 变红 → **立刻**换成真判据 → 绿 → 再进下一步。
- **验收**（条件：AppleClang 17 / Debug、`T = double` 与 `float` **各跑一遍**、容差除注明外 1e-12；
  数字全部取自 `rotation/report.md`，不是新造）：
  - **同态**：`M(a*b) = M(a)M(b)` 误差 <= 1e-12；`Rotate(v)` 与 `ToMatrix()*v` 逐分量一致 <= 1e-12
    （两条独立路径，能抓 `ToMatrix` 的转置错）
  - **往返**：`FromMatrix(ToMatrix(q))` <= 1e-12（double）；四表示互转往返 < 1e-12
  - **双覆盖**：`M(q) = M(-q)` 精确 0（1000 组）
  - **非交换**：`q_x q_z != q_z q_x`，差 O(1)（实测 1.0）
  - **Slerp**：中点 `|q| = 1`（<= 1e-12）；`q0 ≈ q1` 退化 nlerp 不出 NaN
  - **失败模式**：① 叉积符号写反 —— **往返测试抓不到**（正反抵消），
    只有 `q_x q_z != q_z q_x` 与 `M(a*b) = M(a)M(b)` 能抓；
    ② `FromMatrix` 用单分支 `w = sqrt(1+tr)/2` 在 `tr -> -1`（θ→π）处精度崩；
    ③ `AngleOf` 用 `2*arccos(w)` 在 θ→0 导数发散（必须 `2*atan2(|v|, w)`）。
- **不做**：对偶四元数、欧拉角互转、squad 插值。
- **注意**：`Quaternion<T>` 有 `static_assert(is_floating_point_v<T>)` ——
  **只吃浮点，不吃 `Complex` / `Dual`，不参与全库泛型铁律**。
  D16：四元数封死在 `quat/`，不做矩阵元素、不需要 `ScalarTraits<Quaternion>`。
- **契约空白**：`Normalize(零四元数)` 抛还是回单位元 —— 见「待你确认」第 3 条。

---

### 2. 张量 `Tensor` / `TensorView` 实现

**注意：前置更正（2026-10-28 实测）**：`/tmp/tensor.full.{hpp,cpp}`、`.wip-tensor/`、
`graphics/demo/core/tensor.full.bak/` **均已空** —— 887 行完整实现已于 **2026-10-25 按用户指示删除**。
**只能从零实现，不要再找备份。**

- **做什么**：`graphics/src/core/tensor.hpp`（507 行框架，实测 **51 处 `throw`**）按
  ① 几何工具 → ② 构造与索引 → ③ 视图 → ④ 深拷贝 → ⑤ 逐元素 顺序落地；
  `demo/core/tensor.cpp` 的 `[P1]`–`[P24]` 占位判据**逐条换真判据**
  （真判据规格写在同一文件末尾「升级路径」一节，照抄编号即可）。
- **验收**（条件：AppleClang 17、Debug 构建、`-Wall -Wextra`、`T = double` 与 `float` 各跑一遍）：
  ctest `tensor` 从「断言 throw」变真判据全绿；`[A]`–`[F]` 六组有具体数字；
  `[G1]`–`[G8]` 契约段打开；`-c -w`（从零 + 全告警）**0 告警**。
  **失败模式**：① `Slice` 若偷偷复制数据 → 只有 `[C1][C2][C3]` 会红，其余全绿（所以要注入对照）；
  ② 非连续 view 上 `Reshape` 若「需要时静默拷一份」→ `[D3]` 绿而别名语义已错；
  ③ 秩 0 的 0 参 `operator()` 若用零长数组写法 → `-Wpedantic` 告警，且 `t()` 越界读。
- **两条注入对照**（不做就不知道闸门有没有牙）：
  ① 把 `TensorOffset` 的和式写成 `idxs[d]*strides[n-1-d]` ⇒ **应只有 `[B1][C6][C9]` 变红**；
  ② 把 `View()` 改成「拷一份再给 view 指过去」⇒ **应只有 `[C1][C2][C3]` 变红，其余全绿**。
- **开工前必须先定三个待定**（不定就返工；`tensor.hpp` 头注释留的尾巴）：
  1. **`strides_` 是否收敛为单一真源**
     (a) `Tensor` 不存 strides，`Strides()` 每次现算；
     (b) 保留成员，但只留**一条写入路径**（构造时 `ComputeStrides(shape_)`，不给 setter）。
     约束：`Strides()` 现在的签名返回 `const vector&`；且 `operator()` 必须与 view 共用
     `TensorOffset`（索引算术只准一处实现）。
  2. **广播** (a) 显式自由函数 `BroadcastTo`（stride=0 视图，须**只读**）；
     (b) 进 `operator+`。
     约束：`tensor.hpp` 已把「形状不一致必须 throw，不静默广播」写成契约。
  3. **`shape_` / `strides_` 是否「内联容量 4 + 溢出才堆」** (a) 做；(b) 不做。
     约束：头注释把「单线程、最简单版本」定为**设计目标，不是待补的欠债**。
- **不做**：指标收缩 / 乘法 / einsum **一律不进这个类**（是 `src/math/` 的自由函数）；
  不做编译期 rank（会把「不同 rank 的 view 互转」堵死）；不做并发 / 引用计数 / COW / 表达式模板。
- **注意**：设计口径的真源是 `tensor.hpp` 头注释 + `graphics/acceptance.md` §0
  （原 `docs/matrix.md` 已于 2026-10-28 删除，不要按旧路径找）。

---

### 3. `math/ode.hpp` —— Euler / RK4 / Verlet·辛欧拉

**为什么不能延后**（四份文档唯一交集）：
① **Gate A 的最后一块**（`dev-progress.md` §2）；
② 变分段（段 5）的**曲线出口**（`docs/route-math.md` §7）；
③ **CAE 60 #3 + HPC 60 #4** 的共同要求（`docs/eng-competency.md`）；
④ 数学物理工程 60 #1（`docs/baseline-60.md` §2.2）。
**没有任何其它件能做到「一份时间喂四条线」。**

**为什么排第 3 而不是第 1**：它是**新建目录**（`graphics/src/math/` 目前不存在），
没有骨架、没有占位判据可对拍，起手成本高于前两条。

- **做什么**：本件是 `graphics/src/math/` 的**第一件**。
  状态用 `Matrix<T>` 列向量（D19）；`f` 签名 `Matrix<T> f(real_t, const Matrix<T>&)`。
  三入口：`Euler`（一阶对照）、`RK4`（主力）、`Verlet`/辛欧拉（保结构）。
- **验收**（条件：`T = double`、与解析解对拍（谐振子 + 单摆小角）、步长 `h = 1/2^k`、k = 4..10）：
  - **收敛阶**：`h` 减半，RK4 误差降 **≈16 倍**（斜率 ≈ −4）；Verlet ≈ 4 倍（≈ −2）
  - **对拍**：与 `graphics/demo/pyapps/sphere/_core.py` 的手写 RK4 **逐点**对拍
  - **保结构**：辛积分器能量**有界振荡**；非辛**单调漂移**（两条曲线同图）
  - **失败模式**：步长过大 ⇒ 辛格式**也会爆**（不是无条件稳定）；刚性问题显式格式失稳
- **不做**：自适应步长（RK45）、隐式方法（BDF/Radau）、事件检测 —— 都不在段 5 消费点上。
- **注意**：`T` 必须泛型（`double` / `Complex` / `Dual` 都要编过）—— 全库铁律。

---

### 4. `math/fourier.hpp` —— DFT / FFT / 卷积 / 傅里叶级数

**为什么做**：① 它是数学**段 2（L^2 与傅里叶，枢纽）的工程出口** ——
`docs/route-math.md` §7 直接点名 `math/fourier.hpp`（卷积定理 `F(f*g) = F(f)F(g)`、Parseval）；
② `docs/eng-competency.md` 现状表里「FFT / 卷积」一栏是 `无`，这是明确缺口；
③ 它和 `ode.hpp` 一起把 `src/math/` 的目录能力搭起来。

- **做什么**（四入口）：
  1. 朴素 `DFT`（O(N^2)，**只作对照基线**，不做优化）
  2. `FFT`（radix-2，迭代式位反转）
  3. `Convolve` —— 时域直接算 + 频域乘法**两条路径**（补零到 >= N+M-1）
  4. 傅里叶级数：`FourierCoefficients`（系数投影）+ `Reconstruct`（截断重构）
- **验收**：
  - **对拍一**：N = 1024、`Complex<double>`，FFT vs 朴素 DFT **复数逐点**相对误差 <= 1e-12
  - **对拍二**：卷积两条路径逐点 <= 1e-12；Parseval 能量守恒 <= 1e-12
  - **级数**：方波重构的 Gibbs 过冲 **≈ 8.95%（不随 N 减小）**，随 N 增大过冲位置向跳变点收缩
  - **规模**：N = 2^20 时 FFT / DFT 计时比，报 O(N log N) vs O(N^2) 的量级差异
  - **失败模式**：① N 非 2 的幂若静默截断 ⇒ 结果错而不报错（必须 throw）；
    ② 位反转索引写错只错**相位**，幅度判据 `|X_k|` 抓不到，必须比复数逐点；
    ③ 补零不足 `N+M-1` ⇒ 循环卷积卷回头部。
- **不做**：任意长度（混合基 / Bluestein）、多维 FFT、实数 rfft、窗函数、频谱估计。
- **注意**：必须能吃 `Complex<T>`（库已有 `core/complex.hpp`）；实数 `T` 域也要能编过。
  新条目阈值均为**初判据，实现后按实测收紧**。

---

### 5. `math/` —— erf / 正态分布 / Box–Muller 采样

**为什么做**：段 1「概率三件」的工程出口 —— `docs/route-math.md` §7 点名
`math/mc.hpp`（采样 / CLT）；`docs/eng-competency.md` §4 把「随机采样 / 蒙特卡洛」
列为 HPC 赛道检查点。

- **做什么**（三件）：
  1. `Erf(x)` / `Erfc(x)` + `NormalCdf(x)`（= `(1 + erf(x/sqrt2))/2`）
  2. `NormalQuantile(p)`（Φ 的逆）—— 二分反解或 Acklam 近似，二选一
  3. `BoxMuller` 采样（**自带可复现种子**）
- **验收**（条件：`T = double`、x 取 [-6, 6] 上 1001 个等距点、采样 N = 1e6、固定种子）：
  - **对拍**：与 Python 标准库 `math.erf` 逐点对拍（**不需要 scipy**），最大绝对误差记进日志
  - **互逆**：`Phi(NormalQuantile(p)) - p` <= 1e-12，p 取 [1e-9, 1-1e-9]
  - **样本矩**：N = 1e6，样本均值与 0 的偏差、样本方差与 1 的偏差都落在 **3/sqrt(N)** 量级内
  - **失败模式**：① 大 x 用 `1 - erf(x)` 相消 ⇒ 尾部全变 0（必须走 `erfc`）；
    ② `NormalQuantile` 牛顿法在 `p -> 0 / 1` 处溢出；
    ③ `Box–Muller` 里 `log(u1)` 遇 `u1 = 0` 得 `-inf`；
    ④ 用 `std::rand` ⇒ 不可复现，矩验收无法重跑。
- **不做**：多元正态、Cholesky 相关采样、拒绝采样族、KS 检验（只做矩对照）。
- **落点待定**：erf / ppf 是纯函数，可另立 `math/special.hpp`；Box–Muller 按
  `route-math.md` 归 `math/mc.hpp`。边界见「待你确认」第 3 条。
  阈值均为**初判据，实现后按实测收紧**。

---

### 6. 测度第 3 关收尾（数学，无代码）

- **做什么**：Carathéodory 可测性判据 → Vitali 不可测集构造 → Lebesgue 积分四步构造。
- **验收**：能讲清「外测度 → 限制成测度」为什么需要 Carathéodory 判据；
  Vitali 集为什么必须用选择公理。
- **判据**：允许跳证明，不允许跳定义（`docs/route-math.md` §1）。
  **本条无数字验收**（纯数学条目），失败模式就是典型误解：把 Carathéodory 判据当技术条件背下来 ——
  讲不出「外测度不满足可数可加性，所以必须先限制到可测集族」；或说 Vitali 集「取代表元没有构造性」
  却不提**选择公理**。
- **位置**：明细记 `docs/math-log.md`；本件是段 0 的收口。
- **为什么排末位**：它不阻塞任何代码，也不进任何赛道；工程线卡住时它随时可做，反之不然。

---

## 并行线：微基准仓（`learning/learning-benchmark/`）

> **注意：定位更正（2026-10-28 实测）**：该仓**就在父仓内**（不是 `~/develop/learning-benchmark/`），
> 独立 git 仓、父仓 `.gitignore` 第 71 行已登记、父仓不跟踪。
> **10 次提交，设计全部闭合，代码零。**
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

**注意：与 `docs/baseline-60.md` 的联动**：因不绑核，**HPC 赛道 60 #3（性能分析）** 只能靠
**单线程数据 + 对照实验**满足；**多线程加速比不许写进简历**。

---

## 本期不做（明确挂起）

| 不做 | 理由 |
|---|---|
| 稀疏 CSR / SpMV / CG | CAE 入口，等 `ode.hpp` 与 `Tensor` 之后 |
| Gauss–Legendre 求积 | **零行代码**。注意库里三个「Gauss」不是一回事（`eng-competency.md` §3 第 2 条） |
| 渲染（BRDF/路径追踪/BVH）、几何（NURBS/半边） | 属其它赛道 |
| GR 计算、G 系 C++ 化 | 依赖张量实现 |
| `einsum` 薄层 | 无消费者（张量侧的收缩产品还没定） |

---

## 下一期候选（无序，按需挑，不承诺）

| 候选 | 触发条件 | 一次投入喂几条线 |
|---|---|---|
| **稀疏 CSR + SpMV + CG + 1D Poisson FEM + Gauss–Legendre** | 本期六条完成后 | CAE 60 #1/#2 + HPC 60 #1/#2 |
| **`math/pde.hpp` + `LaplacianMatrix` + CG** | 想开 Gate B（PDE 可视化） | CAE + 数学段 4–5 |
| **微基准仓 M2**（变体铺开） | 该仓 G1 闸门过后 | HPC 求职主线 |
| **李群刀 4–6**（exp 非单射 / SU(2) 双覆盖 / BCH 测阶） | 纯线代，成本最低 | 只在基线里，不喂赛道 |

---

## 待你确认（我查得到、但不能替你决定）

1. **岗位抽样**（`docs/eng-competency.md` §3 第 3 条）：取 20 个真实 JD，
   只记三列 ——「学历要求 / 是否明写稀疏·预条件·MPI·GMRES·BRDF·BVH / 是否要求专业对口」。
   **这是市场事实，agent 判断替代不了。**

2. **父仓 33 项未提交**（13 M / 8 D / 12 ??；数字以 `git status` 为准），含 `D graphics/plan.md`、
   `D src/main.cpp`、`D src/view/*` 等删除，以及 `?? NOW.md`、`?? docs/baseline-60.md`、`?? docs/eng-competency.md`。
   **文档重构成果只在工作区** —— 建议提交一次。
   注意：`NOW.md` 是未跟踪文件，本版覆盖不可恢复。

3. **本期五个契约空白**（未定前对应条目不开工；论述在各条目内，这里只做清单）：
   - `Tensor` 三个待定：`strides_` 单一真源 / 广播是否进 `operator+` / 内联容量 4 —— 见本期第 2 条
   - `Normalize(零四元数)` 抛还是回单位元 —— 见本期第 1 条
   - erf / ppf 的落点：另立 `math/special.hpp` 还是并进 `math/mc.hpp` —— 见本期第 5 条

---

## 「做完」的判据模板（每个交付物都套）

> **条件 + 数字 + 失败模式**，三件缺一不算完。
> 来源：`docs/eng-competency.md` §0.2 —— **这是 60 分与「我跑通了」的唯一分界。**

- **条件**：机器、编译器与优化级、规模、容差、基线是什么
- **数字**：误差 / 收敛阶 / 迭代数 / 加速比，要能和基线对照
- **失败模式**：什么输入会让它错、错成什么样、你怎么知道的

反例：「CG 快 3.2×」→ 噪音。
正例：「N=10^5、float64、-O3、容差 1e-8，对比 Jacobi 预条件，迭代数 812 → 96，预条件耗时占比 7%」。