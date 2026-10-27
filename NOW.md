# NOW.md — 近期要做什么（唯一滚动文件）

> **本文件是全仓唯一写「接下来」的文件。**
> **规则**：1–3 周一版；过期就**整篇重写**，不追加、不归档、不留历史。
> 历史只进 `docs/logs/<YYYY-MM>.md`（append-only，最新在上）。
> **不写阶段、不写月份、不写"长期目标"** —— 那些一律视为噪音。

**本期起点**：2026-10-28（第二次重写。上一版同日首版有两处事实错误，本轮已修）

---

## 本期三条（做完就重写本文件）

### 1. `math/ode.hpp` —— Euler / RK4 / Verlet·辛欧拉

**为什么排第一**（四份文档唯一交集）：
① **Gate A 的最后一块**（`dev-progress.md` §2）；
② 变分段（段 5）的**曲线出口**（`docs/route-math.md` §7）；
③ **CAE 60 #3 + HPC 60 #4** 的共同要求（`docs/eng-competency.md`）；
④ 数学物理工程 60 #1（`docs/baseline-60.md` §2.2）。

**没有任何其它件能做到「一份时间喂四条线」。**

- **做什么**：`graphics/src/math/` 目录**不存在**，本件是第一件。
  状态用 `Matrix<T>` 列向量（D19）；`f` 签名 `Matrix<T> f(real_t, const Matrix<T>&)`。
  三入口：`Euler`（一阶对照）、`RK4`（主力）、`Verlet`/辛欧拉（保结构）。
- **验收**（条件 + 数字 + 失败模式，三件缺一不算完）：
  - **收敛阶**：`h` 减半，RK4 误差降 **≈16 倍**（斜率 ≈ −4）；Verlet ≈ 4 倍（≈ −2）
  - **对拍**：与 `graphics/demo/pyapps/sphere/_core.py` 的手写 RK4 **逐点**对拍
  - **保结构**：辛积分器能量**有界振荡**；非辛**单调漂移**（两条曲线同图）
  - **失败模式**：步长过大 ⇒ 辛格式**也会爆**（不是无条件稳定）；刚性问题显式格式失稳
- **不做**：自适应步长（RK45）、隐式方法（BDF/Radau）、事件检测 —— 都不在段 5 消费点上。
- **⚠️** `T` 必须泛型（`double` / `Complex` / `Dual` 都要编过）—— 全库铁律。

---

### 2. 张量 `Tensor` / `TensorView` 实现

**⚠️ 前置更正（2026-10-28 实测）**：`/tmp/tensor.full.{hpp,cpp}` **与 `.wip-tensor/` 均已空** ——
887 行完整实现已于 **2026-10-25 按用户指示删除**。**只能从零实现，不要再找备份。**

- **做什么**：`graphics/src/core/tensor.hpp`（507 行框架，**实测 51 处 `throw`**）
  按 ①形状工具 → ②构造与索引 → ③视图 → ④深拷贝 顺序落地；
  `demo/core/tensor.cpp` 的 `[P1]`–`[P24]` 占位判据**逐条换真判据**。
- **验收**：ctest `tensor` 从「断言 throw」变真判据全绿；`[A]`–`[F]` 六组有具体数字；
  `-c -w`（从零 + `-Wall -Wextra`）0 告警。
- **开工前必须先定三个待定**（`docs/matrix.md` §6 第三次订正留的尾巴；不定就返工）：
  1. `strides_` 是否收敛为**单一真源**（`ComputeStrides(shape_)` 派生，只有 view 允许任意 strides）
  2. 广播做显式 `BroadcastTo`（stride=0 只读视图）还是进 `operator+`
  3. `shape_`/`strides_` 是否「内联容量 4 + 溢出才堆」
- **不做**：指标收缩/乘法/einsum **一律不进这个类**（是 `src/math/` 自由函数）；
  不做编译期 rank；不做并发/引用计数/COW/表达式模板。

---

### 3. 测度第 3 关收尾（数学，无代码）

- **做什么**：Carathéodory 可测性判据 → Vitali 不可测集构造 → Lebesgue 积分四步构造。
- **验收**：能讲清「外测度 → 限制成测度」为什么需要 Carathéodory 判据；
  Vitali 集为什么必须用选择公理。
- **判据**：允许跳证明，不允许跳定义（`docs/route-math.md` §1）。
- **位置**：明细记 `docs/math-log.md`；本件是段 0 收口。

---

## 并行线：微基准仓（`learning/learning-benchmark/`）

> **⚠️ 定位更正（2026-10-28 实测）**：该仓**就在父仓内**（不是 `~/develop/learning-benchmark/`），
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
| **硬件 PMU** | ❌ **不可用**（`unsupported p6 CPU model 183 no PMU driver`）⇒ 归因**只走唯一变量对照**，`attribution.method` 永不为 `"perf"` |
| **绑核** | ❌ **不绑核**（WSL2 拓扑是合成的：`lscpu` 报 10 核×2 线程，与真实 6P+8E 不符）⇒ `bind.kind` 恒 `none`，多线程数据 `reliable=false` + `caveats` 必填 |

**⚠️ 与 `docs/baseline-60.md` 的联动**：因不绑核，**HPC 赛道 60 #3（性能分析）** 只能靠
**单线程数据 + 对照实验**满足；**多线程加速比不许写进简历**。

---

## 本期不做（明确挂起）

| 不做 | 理由 |
|---|---|
| **`quat` 剩余 10 处 `throw`** | 骨架/判据/对拍源全齐，是**填空型** —— 直通 V2.5/V2.6，但**不喂任何赛道**，等三条主线做完再当填空 |
| 稀疏 CSR / SpMV / CG | CAE 入口，等 `ode.hpp` 之后 |
| Gauss–Legendre 求积 | **零行代码**。注意库里三个「Gauss」不是一回事（`eng-competency.md` §3 第 2 条） |
| FFT / 卷积 / `fourier.hpp` | 段 1–2 枢纽，但工作量最大（4 入口 + 两套对拍） |
| erf / ppf / Box–Muller | 段 1 概率三件数值侧 |
| 渲染（BRDF/路径追踪/BVH）、几何（NURBS/半边） | 属其它赛道 |
| GR 计算、G 系 C++ 化 | 依赖张量实现 |
| `einsum` 薄层 | 无消费者 |

---

## 下一期候选（不承诺）

| 候选 | 触发条件 | 一次投入喂几条线 |
|---|---|---|
| **稀疏 CSR + SpMV + CG + 1D Poisson FEM + Gauss–Legendre** | 上面三条完成后 | **CAE 60 #1/#2 + HPC 60 #1/#2**（重叠最高） |
| **`fourier.hpp`**（DFT/FFT/卷积） | 想推数学段 1–2 | 段 1–2 枢纽 |
| **微基准仓 M2**（变体铺开） | 该仓 G1 闸门过后 | HPC 求职主线 |

---

## 待你确认（我查得到、但不能替你决定）

1. **岗位抽样**（`docs/eng-competency.md` §3 第 3 条）：取 20 个真实 JD，
   只记三列 ——「学历要求 / 是否明写稀疏·预条件·MPI·GMRES·BRDF·BVH / 是否要求专业对口」。
   **这是市场事实，agent 判断替代不了。**

2. **父仓 33 项未提交**（13 M / 8 D / 12 ??），含 `D graphics/plan.md`、`D src/main.cpp`、
   `D src/view/*` 等删除，以及 `?? NOW.md`、`?? docs/baseline-60.md`、`?? docs/eng-competency.md`。
   **文档重构成果只在工作区** —— 建议提交一次。

---

## 「做完」的判据模板（每个交付物都套）

> **条件 + 数字 + 失败模式**，三件缺一不算完。
> 来源：`docs/eng-competency.md` §0.2 —— **这是 60 分与「我跑通了」的唯一分界。**

- **条件**：机器、编译器与优化级、规模、容差、基线是什么
- **数字**：误差 / 收敛阶 / 迭代数 / 加速比，要能和基线对照
- **失败模式**：什么输入会让它错、错成什么样、你怎么发现的

反例：「CG 快 3.2×」→ 噪音。
正例：「N=10⁵、float64、-O3、容差 1e-8，对比 Jacobi 预条件，迭代数 812 → 96，预条件耗时占比 7%」。