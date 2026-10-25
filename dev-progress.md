# dev-progress.md — 进度与排期

> **读者**：人 + agent ｜ **最后更新**：2026-10-23 ｜ **本文是「进度」的唯一真源**
>
> **真源声明**
> - 本文是**总纲式进度真源**：回答「整体现在在哪、下一步做什么」。
> - `graphics/plan.md` 只负责任务定义与验收阈值（验收项记作 `V{阶段}.{序号}`），**不记录状态**。
> - **数学线的关卡明细，真源是 `learning-todo/ns-progress.md`**（逐次会话日志 + 关卡状态）。
>   本文只保留一行指针，**不重复记录**，以免两处漂移。
> - 冲突处理：工程线以本文为准；数学线以 `ns-progress.md` 为准。
>
> 分工（2026-10-23 用户明确）：**AI 只负责整个工程与文档** —— 排期 / 定门槛 / 验收阈值 /
> 设计契约 / 代码审查；**不写代码**（用户原话：「你永远只负责整个工程，以及文档，
> 不负责代码，代码我会在其它 agent 里自己写」）。优先学习，工程增量推进。

---

## 0. 当前指针（每次会话先看这里）

> **数学线在 L1（测度与积分）中段；工程线在 `graphics` 阶段 2 进行期 —— V1.1–V1.6 全验收 + V2.4（Padé 一半）✅ + pyapps 两系列十二脚本落地（球面 G1–G6 + 旋转 R1–R6，图 13 张），ctest **16/16**（autodiff 新增 `jacobian` 18 项；`quat` 59 项；`so3` 51 项）；库侧剩 quat 实函数 / LinearOperator / classify 三件（见 §4）。**

| 线 | 当前位置 | 下一个动作 |
|---|---|---|
| **数学（优先）** | **段 0（测度与积分）进行中** —— 路线真源改为 `learning-todo/math-plan.md`（十段表；原 L0–L5 作废）。当前位置明细见 `learning-todo/ns-progress.md` 的「📍 当前位置」 | 段 0 收口（Carathéodory / Vitali / 积分四步）→ **段 1 概率三件**（条件期望=正交投影、卷积=前推）→ 段 2 L²/傅里叶 |
| **工程** | `graphics` HEAD `fb71bec`（**V2.4 Padé 一半 ✅**，expm_pade = [13/13] scaling），ctest **17/17**（**新增 `tensor` = Test #2，24 条占位判据全绿**，2026-10-23 实测；`quat` 59 项，9 个入口仍占位 throw；**新增 `so3` 51 项**）+ **autodiff 布局翻案**（全库列向量 + 分子布局，新增 `jacobian` 18 项）。**2026-10-20 上午新增 `src/linalg/so3.hpp`**（`Cross` / `Hat` / `Vee`，只依赖 `core/matrix.hpp`）—— 四元数 `operator*` 与 `Rotate` 的叉积抽层落点，**`quat/` 故意不 include 它**（D16）。**2026-10-23 新增 `src/core/tensor.hpp`**（507 行**框架**：`Tensor<T>` 哑容器 + `TensorView<T>` 非占有视图，契约全、实现体全 `throw`）+ `demo/core/tensor.cpp`（347 行，`[P1]`–`[P24]` 断言"仍然 throw"；实现掉一条就变红）；契约真源 = `learning-todo/matrix.md` §6 第三次订正 + 决策 #21 + 红线 #4。**einsum 仍未立项，`src/math/` 目录不存在**。**2026-10-17 上午 pyapps 两系列落地（未提交，整目录 untracked）**：`demo/pyapps/sphere/` **G1–G6**（`_core.py` = 度规/Γ/黎曼/RK4/平行移动公共层，六脚本全过：K 偏离 1 仅 4.4e-14、holonomy 对拍 1.7e-11°）+ `demo/pyapps/rotation/` **R1–R6**（五表示互转/万向锁/Slerp，`_rot_core.py` 15 组自检；**R2 已含四元数全套数学**）；图 13 张 → `learning-todo/imgs/{geo,rot}/`（rot 6 张未提交）。数学侧两次口径修正已入 fan-plan（G3 极角/纬度、R 系判据病态 atan2） | ① **lie 刀 4/5/6**（expm.cpp 续写；**R2 已把刀 4 的双覆盖 Python 前哨做掉**——`M(q)=M(−q)`）② **N2 `math/ode.hpp`**（Gate A 最后一块积木；**G 系列 `_core.py` 已手写 RK4，即现成对拍源**）③ **B1 quat 实函数**（骨架已立，`operator*` ✔，9 个入口待实现）④ **B2 LinearOperator**（弦/鼓膜/热模糊本体） |
| **小项目** | 未解锁 | 等 Gate A（= matrix ✔ + linalg ✔ + **一个积分器** ← 就差 `ode.hpp`） |
| **扇面（并行）** | 真源 `learning-todo/fan-plan.md`（§3.1 G / §3.2 R 已含实测状态与口径警告）；**A5「四元数」pyapps 预演完成**（R1–R6），A1 C++ 版已出（string_modes） | **扇面收束计分（§7）**：A 面 **3 条 ✔**（G+R+string，超额）；B/C/D 面 0 条 → **D2 CLT 动画是下一优先**（D 面第一条 + 特征函数困惑的正解）；C 面 C1 共形映射可顺延 R/G 经验 |
| **李群线（扇面 A4 的展开）** | 刀 1/2/3 已落地（V1.6 ✅）；刀 4 的双覆盖 Python 前哨已在 R2 完成一半（`M(q)=M(−q)` + 漂移 √n·eps） | **刀 4/5/6**（expm.cpp 续写）：θ 扫 [0,4π] 非单射 + **SU(2) 双覆盖 C++ 版**（θ=2π 处 `−I`）+ BCH 测阶（缩半降 8 倍） |
| **笔记 web 化** | ✅ **已完成**（父仓此前漏记，2026-10-17 补正实况）：`learning-web` 独立仓 **139 次提交**、**json 词条 318/318 全入库**、前端四层实现齐全（`shell`/`renderers`/`modules`/`viz` 27+ 档1 模块）+ `dist/` 产物（`index.html` / `knowledge.html` / `graph/` / `json/`）已出，**2026-10-04 收工**（HEAD `3ec63a6`）| 无（已收束）；后续按需加餐，不占主线 |
| **微基准/归因框架（求职线）** | ✅ **已独立成仓**（2026-10-25）：`~/develop/learning-benchmark/`（独立 git 仓，init 提交 `a3234b2`；父仓 `.gitignore` 已登记，不跟踪）。**设计层闭合、实现层为零**。真源 = 该仓 `docs/bench-scope.md`（设计 608 行）+ `docs/bench-handoff.md`（执行清单 238 行）。目标岗位 = **算法岗配套的测试开发**（非通用 QA）。口径 = **单机**（14600K + 3080 / WSL2；M3 降为 P3 可选） | 见该仓 `docs/bench-handoff.md` §2 的 P0 清单；**唯一阻塞** = WSL2 里 `dmesg \| grep -i PMU` 实测。**代码不归我写** —— 我负责审五步依赖表 + 代码审查 + 验收 |

> **2026-09-22 修正**：本表原先写「L1 第 3 关 —— 外测度已学、勒贝格积分刚学」，
> 与 `ns-progress.md` 的「外测度 / Carathéodory **未接触**」**互相矛盾**。
> 已改为**只留指针不做重复记录** —— 数学线的关卡状态一律以 `ns-progress.md` 为准。
> 同日用户确认：**外测度确已学**，滞后的是 `ns-progress.md` 一侧，该文件已一并更新。

**更新方式**：每次会话用户报「学到哪 / 写到哪」→ 改上表 + 在第 6 节日志追加一条。
上表只保留**当前**状态（覆盖式），历史一律进日志（追加式，不删改）。

---

## 1. 协作模式

- 数学线（慢，优先）与工程线（快，够用即做）并行。
- 每次会话：用户报「学到哪 / 写到哪」 → AI 更新本文件并给下一步。
- **matrix 改造**与**自动微分**由用户自行抽时间完成，AI 只给验收点。

---

## 2. 小项目启动门槛（Gate）

> **⚠️ 2026-10-23 口径更新**：数学路线已改为 `learning-todo/math-plan.md` 的**十段表**
> （原 L0–L5 作废）。**真正的关卡从「L1 完成」挪到「段 5 变分两出口」** ——
> 那里是「能读教材」与「能读场方程」的分界（math-plan.md §5 冲突 ③）。
> 下表按新段号重定义；旧 Gate A/B/C 的数学条件**作废**。

| 门槛 | 数学（按 math-plan 段号） | 工程 | 解锁的小项目 |
|---|---|---|---|
| **Gate A** | 段 0（测度积分）+ ODE 数值解直觉（Euler/RK4/Verlet） | matrix 改造 ✔ + 最小 linalg（LU/solve）+ **一个积分器（`ode.hpp`）** | 动力系统可视化：单摆 / 行星轨道 / 双摆 / N 体 |
| **Gate B** | 段 1–2（概率三件 + L²/傅里叶）+ 段 5 前半（δ 规则、EL） | Gate A + FDM 基础（`pde.hpp`）+ `fourier.hpp` | PDE 可视化：热方程（图像模糊）、波动（鼓膜） |
| **Gate C** ★真关卡 | **段 5 全程**（两出口 + 边界项 + Noether）+ 段 3–4（泛函/Sobolev） | 阶段 4 完整（`pde.hpp` + 边界条件组装 + CG） | 弱解类（通往 NS / 场方程） |
| **Gate D** | 段 6–7（流形与张量场、联络与曲率） | `Tensor`/`TensorView` + `Einsum` 薄层 | GR 侧：测地线、曲率、场方程数值 |

> **第一个小项目定在 Gate A**；数学线继续按 `math-plan.md` 推进，不互相阻塞。
>
> **共用的「主干五件」**（`ode.hpp` / `pde.hpp` / `fourier.hpp` / `Tensor+Einsum` /
> `LaplacianMatrix`）是 **CG / GR / NS 三条支线全部共用**的（math-plan.md §7.1）。
> **不是「选一条支线」，而是「把主干修宽」** —— 主干修好，三个终点都能接着长，
> 任何一支都不浪费。工程排序的唯一判据：**这件东西在不在主干上。**

---

## 3. 现状明细（截至 2026-10-15）

- **数学线**：**只以 `learning-todo/ns-progress.md` 为准**，本文不重复关卡明细。
- **工程线**（`graphics`，分支 `graphics_0`；2026-10-19 实地盘点到 ctest 15/15）：

  **linalg 已全部落地、且每个公开入口都真的进过编译器** —— 2026-10-14 那次「528 行从未被任何 TU
  include」的体检问题，已随 S1–S6 全部修复；当时 32 个编译错误的逐条明细保留在 §6 的
  2026-10-14 条，本节只保留结果：

  | 文件 | 行数 | 验收用例 | 状态 |
  |---|---|---|---|
  | `operator.hpp` | 40 | `smoke`（`Dot` 的第一/唯一调用方） | ✔ S6 收编 |
  | `expm.hpp` | 54 | `expm` 16 项 | ✔ S5 → **V1.6** |
  | `cholesky.hpp` | 116 | `cholesky` 98 项 | ✔ S2a |
  | `solve.hpp` | 129 | `solve` 34 项 | ✔ S1 |
  | `svd.hpp` | 137 | `svd` 266 项 | ✔ S4 → **V1.4** |
  | `inverse.hpp` | 188 | `inverse` **85 项**（含新增 `qr_inverse` [H] 段）+ `smoke` B6d/B6e | ✔ S2 / S2a-2 → **V1.2**；**S-App6 ✔ 2026-10-16** |
  | `eigen.hpp` | 167 | `eigen` 119 项 | ✔ S3 → **V1.3** |
  | `qr.hpp` | 258 | `qr` 189 项 | ✔ S6b（V2.1） |

  合计 **1089 行**（`qr_inverse` 补回代 +27 行）。`householder.hpp` **已删除并随 `6999ca2` 提交**（内容由 `qr.hpp` 自带的
  Householder 接手）。**全局实例化闸门** `demo/linalg/smoke.cpp`（2026-10-16 新增）：
  `EntryTable<T>` 对 40 个公开入口取函数指针（odr-use ⇒ 强制实例化）+ 54 项最小数值判据，
  double 与 Complex 各跑一遍 —— **这是 V1.1 的显式证据**。

  - **测试入口**：`./graphics/run_tests.sh` → ctest **16/16**、构建 0 告警 0 错误
    （含 `-c -w` 从零 + `-Wall -Wextra` 复验）；附加区 `./graphics/run_apps.sh`
    （`demo/apps/`，**不注册 ctest**，首个成员 `string_modes`，即 fan-plan 扇面 A1 的 C++ 版）。
  - **V1 验收落点**（`graphics/plan.md` 的 V1.1–V1.7，逐条）：
    - **V1.1 ✔** 2026-10-16 通过（ctest **`smoke`**，54 项 / FAIL 0，构建 0 告警）：
      40 个公开入口 × {`double`, `Complex<double>`} 各实例化一次（同一份源码零修改），
      数值侧最差 Solve 残差 5.8e-17、Inverse 2.2e-16、QR 重构 3.4e-16、expm 自洽 6.7e-16。
      **注入反例两条**（见 §6 同日条）：漏 `conj` → 全项目**只有 smoke 挂**（1 条）；
      写入只在实例化时暴露的编译错误 → 「只 include」的探针**编得过**、其余 10 个 ctest
      **全部照常构建**，只有 smoke 编不过。
    - **V1.2 ✔** 求逆残差最差 double **1.18e-13** / complex **9.16e-16**（判据 1e-8）
    - **V1.3 ✔** Jacobi 重构最差 double **1.78e-15** / complex **1.33e-15**
    - **V1.4 ✔** SVD 重构最差 double **1.332e-15** / complex **7.772e-16**；σ 非负降序
    - **V1.5 ✔** 2026-10-16 通过（ctest **`dual_complex`**，22 项 / FAIL 0）：
      `Matrix<Dual<Complex<double>>>` 零修改走通（含 `Solve`，value 与纯 Complex **精确 0**，n=2/4/6）；
      z² 单种子（两种子 `derive/s` 一致 = 2z₀）/ |z|² 双种子（重组出 2z̄）均**精确 0**。
      另加「A(θ)x=b 求 dx/dθ」自动传播，与解析式 5.2e-18 / 中心差分 6.9e-12 三方对拍。
    - **V1.6 ✔** 2026-10-16 通过（ctest **`expm`**，16 项 / FAIL 0）：`expm.hpp`（S5 刀 1）+
      刀 2 `[Ji,Jj]=ε_ijk Jk` 判**精确 0** + 刀 3 **三线对拍**（级数 / 谱路线 / Rodrigues）
      在 θ∈{0,0.3,1,2,3,−0.9,2π} 上 ≤1.2e-15 + 几何判据（`R·u=u`、`tr R=1+2cosθ`、`det R=1`、
      `RᵀR=I`）与超强判据 `det(exp X)=exp(tr X)` ≤1.0e-15。
    - **V1.7 ❌ 未做** —— 纯演示（单位复数乘法 = 平面旋转，**不作硬验收**），可落 `demo/apps/`；
      **用户 2026-10-16 明确「先不用管」**，故不阻塞 V1 收口。
  - **其余工程事实**：`compute_graph` 基线（forward-mode）成立；**反向模式明确未实现**
    （抛 `logic_error`，按 `matrix.md` §7.5 属**待立项功能**而非 bug）；
    `scalar_traits` 44 项 + `compute_graph` 36 项；`a_triangle` 空文件、无 CMakeLists、
    未注册（等阶段 3 渲染骨架）。
  - **autodiff `Apply` 域分裂（2026-10-16 落地；§3 此前「编不过」的记录随之作废）**：
    `Apply<T>` 原先只有一个 `Matrix<T>` 版本，体内却假定 `T = Dual<U>` ⇒ 标量域必然报
    `error: use of undeclared identifier 'dualFunc'` —— `dualFunc` 是 `Dual` 的 **in-class
    friend**（**ADL-only**，无法前置声明），`T = Dual<U>` 时实例化点 ADL 命中 ✅，标量域 ❌。
    **现按域拆成 4 个重载**：`{FunctionType, Matrix<FunctionType>} × {Dual 专用, 标量专用}`。
    - Dual 版签名写 `Matrix<Dual<T>>`（走 `dualFunc`）；标量版 `Matrix<T>`（走新增的
      `elemFunc` —— 标量域逐元素分派，`using std::sin;` + ADL）。
    - **分流靠 [temp.deduct.partial] 部分排序**（`Matrix<Dual<U>>` 比 `Matrix<T>` 更特化），
      不是任何 hack；标量版体内的 `static_assert(!IsDual<T>::value, ...)` 是**守卫**
      （新增 trait `IsDual`，落在 `dual.hpp`）。
    - **顺带修掉一个真 bug**：`dual.hpp` 的 `Dual<T>::ln` / `pow` 体内写的是 `log(x.value)`，
      而 `Dual` 的友元叫 **`ln` 不是 `log`** ⇒ 对 `T = Dual<U>`（嵌套对偶 = **二阶导**）
      ADL 断链（`no matching function for call to 'log'`）。补一个 `log` 转发即通
      （实测 `ln(2)` 二阶导 −0.25 = −1/x²、`x³` 二阶导 12 = 6x，全中）。
    - **实测**：4 域（`double` / `Complex` / `Dual<double>` / `Dual<Complex<double>>`）
      × 4 重载全部实例化且数值对；derive 对拍解析导数 50 组**误差 0**；中心差分对拍最差
      **3.666e-09**。
  - **新增 ctest `apply`（17 项 FAIL 0）** —— 补上了此前最大的覆盖缺口：`GradOp::Function` /
    `MatrixFunction` 两个节点**从未跑过** ⇒ `Apply` 数值零覆盖（连前向都没验）。现覆盖
    4 重载 × 4 域 + Function/MatrixFunction 节点真实路径 + 嵌套对偶二阶导 + 维度不匹配 throw。
  - **未动的一条**：`GradNode<T>` 仍然只支持 `T = Dual<U>`（`SumAll` / `MultAll` 等签名
    写死 `Matrix<Dual<T>>`，`GradNode<double>` 实例化报 4 条错）—— 属另一件事。
  - **工作区未提交改动**（`HEAD = 6999ca2`，**未 push**）：`demo/linalg/smoke.cpp`（新，V1.1 闸门）、
    `demo/linalg/CMakeLists.txt`（注册 smoke）、`src/autodiff/{dual.hpp,jacobian.hpp}`（域分裂 +
    log 修复）、`demo/autodiff/{apply.cpp,CMakeLists.txt}`（新用例）、`plan.md`（S 清单）、本文件。

  > 2026-10-13 的「基线接口层 bug 第二批修复」明细（标量契约 `const_n → fromReal`、
  > `ScalarTraits<Dual<T>>::conj/abs2` 成员名、`Derive` 收临时量、`Dual::operator=` 弃用告警、
  > arity / `matrixF` 越界保护、`SumAll` 占位形状）**已移入 §6 的 2026-10-13 条**，本节不再重复。

---

## 4. 下一步

**数学线（优先）**

> 关卡明细与逐条状态**只以 `learning-todo/ns-progress.md` 为准**，本文不重复记录。
> （2026-09-23 移除本节原有的 `- [ ]` 勾选清单：它长期挂着「外测度未学」，
> 与 §0 指针及 `ns-progress.md` 直接矛盾 —— 属「两处记录」漂移，按真源约定删除。）

**工程线**（用户自定节奏；任务定义与阈值见 `graphics/plan.md` 阶段 1）

- [x] 修 `grad_graph.hpp` / `compute_graph.cpp` 编译错误 → **2026-10-13 完成**（真因见 §6 2026-10-13 条：
      `enum class` 作用域前缀 35 处 + `Matrix(1,1)` CTAD 2 处 + 缺 `Output`/`Mult` 分支）；
      `ctest --test-dir graphics/build` **2/2** 通过（当时的用例数，现为 8/8）
- [x] 修「接口层」bug（同一会话第二批，见 §6 2026-10-13 条）：标量契约 `const_n → fromReal` + concept 补条目、
      `ScalarTraits<Dual<T>>::conj/abs2` 成员名、`Derive` 收临时量、`Dual::operator=` 弃用告警、
      arity/`matrixF` 越界保护、`SumAll` 占位形状；反向模式明确拒绝执行。
      新增 ctest `scalar_traits`（44 项）+ `compute_graph`（36 项）
- [x] 落「基线」提交 → **2026-10-13/14 完成**：父仓 `graphics_0` 三次提交 ——
      `375f618`（基线修复）、`c82c7ed`（新增 `graphics/run_tests.sh` + AGENTS.md 登记）、
      `1f33b7a`（脚本补 `--build-only` / `--test-only`）；子仓 `learning-todo` `e733363`（`matrix.md`）。
      **均未 push**（推送由用户自行决定）。

> **S1–S6b 已全部提交**（`solve` `inverse` `cholesky` `eigen` `svd` `qr` 六个入口 + `demo/common.hpp`
> 测试基建），提交序 `17ea31b`(S1) → `c0ab458`(S2) → `daefb3e`(S2a) → `f73790b`(S2a-2) →
> `5751b05`(S3) → `edfc572`(S4) → `bae45df`(S6b) → **`6999ca2`（V1 收尾一刀：expm.hpp=刀1 +
> `dual_complex`（V1.5）+ `expm.cpp`（V1.6）+ apps 附加区与 `run_apps.sh`，ctest 10/10）**；
> `HEAD = 6999ca2`，**均未 push**。过程明细见 §6；**当前缺口只剩下面两条。**

**V1 剩余两步（V1.1–V1.6 已全部验收，ctest 16/16）**：

- [x] **① V1.5 用例 `demo/linalg/dual_complex.cpp`** —— ✅ 已完成（2026-10-16）：`demo/common.hpp`
      补 `RandFill<Dual<Complex<U>>>` 特化（原 `RandFill<Dual<U>>` 的 `static_cast<U>(a)` 在
      `U = Complex` 时编译失败 —— `Complex` 无单参构造），新增 ctest `dual_complex`（22 项 / FAIL 0）。
      三条注入验证过（含一条**暴露并修掉我自己判据漏洞**的）。明细见 §6 同日条（第 1 条）。
- [x] **① V1.6 demo `demo/linalg/expm.cpp`** —— ✅ 已完成（2026-10-16，与 `lie-plan.md` 刀 2/3
      合并做）：16 项判据（刀 1 expm 本体 / 刀 2 `[Ji,Jj]=ε_ijk Jk` 判精确 0 / 刀 3 **三线对拍** /
      几何 + `det(exp X)=exp(tr X)` 超强判据），已随 `6999ca2` 提交。
- [x] **② V1.1 显式闸门 `demo/linalg/smoke.cpp` + S6 收口** —— ✅ 已完成（2026-10-16，ctest 11/11 → 12/12）：
      40 个公开入口取函数指针强制实例化（`kEntryCount` 计数防漏抄）+ 54 项最小数值判据，
      double 与 Complex 各一遍。**注入两条**证明它真的补上了覆盖（漏 `conj` → 全项目只有它挂；
      只在实例化时报的编译错误 → 只 include 的探针编得过、其余 10 个 ctest 照常构建，只有它编不过）。
      **支线**：`Dot` 由本闸门收编（不再零调用方）；`qr_inverse` 当日随后补上回代（**S-App6 ✔**，
      闸门 B6d/B6e 由「断言 throw」换成数值判据）；`MatrixType` 由用户定案挪到 `demo/common.hpp`。
- [ ] **③ V1.7 + S7 收尾**：V1.7（单位复数乘法 = 平面旋转，纯演示，落 `demo/apps/`；**用户已定
      「先不用管」**）；S7 = V1.1–V1.7 在 `graphics/plan.md` 逐条勾选 + 本文 §0/§3/§4/§6 + 记忆同步。
      **注：V1 的实质内容（V1.1–V1.6）已全绿，V1.7 是唯一不作硬验收的演示件 —— 阶段 2 不必等它。**

- [x] **S6 支线定案 —— ✅ 2026-10-16 三条全结算**：① **`Dot`** 由 `smoke` 收编（不再零调用方）；
      ② **`MatrixType`** —— **用户定案：挪到应用层 `demo/common.hpp`**（它按**用途**划分，代码原注释
      自陈「不能按用途去分，怎么用是上层的责任」）；core 层只留按**代数结构**分的 `MatrixKind`
      （带 tol 判定）—— 与 `matrix.md` §3.7 实测结论一致（不给稠密 Matrix 存运行期结构标签）；
      ③ **`S-App6` `qr_inverse`** 补 `A⁻¹ = R⁻¹Qᴴ`（回代）—— 已实现，见 §6 今日条。
- [x] **✅ autodiff：`Apply` 的域契约已写进签名**（2026-10-16 落地，用户令「你改吧」）：
      按域拆成 **4 个重载** —— Dual 专用 ×2（签名 `Matrix<Dual<T>>`，走 `dualFunc`）
      + 标量专用 ×2（签名 `Matrix<T>`，走新增的 `elemFunc`）；两者共存靠
      [temp.deduct.partial] **部分排序**分流（`Matrix<Dual<U>>` 更特化），标量版体内
      `static_assert(!IsDual<T>::value, ...)` 作守卫。**新增 ctest `apply`（17 项 FAIL 0）**。
      顺带修掉 `dual.hpp` 的 `ln`/`log` ADL 断链（嵌套对偶 = 二阶导）。四条注入对照见 §6 今日条。
      **已了结。**
- [ ] **遗留（不阻塞 V1 收口）**：`src/core/perspective.hpp` / `rotate.hpp` 旧接口未处置
      （阶段 1 重写 matrix 时本应一并处理）。

### 下一版基线（2026-10-17 03:00 定，以「可视化 demo」为目标）

**基线原点 = `fb0e78b` + ctest 13/13**（含本日 03:00 修复：`expm.cpp` 的 `E_pade_vs_taylor`
误贴两份导致编译炸，删重复块后全绿；备份 `/tmp/expm.cpp.wipbak`）。**收基线纪律：动过
demo/ 或 src/ 的会话结束前必跑 `run_tests.sh`，全绿才作数。**

阶段 2 收窄为「可视化 demo 期」——qr/cholesky/expm 三个阶段 2 库件已提前完成，剩下的
库侧缺口按「服务哪个可视化 demo」倒排：

**库侧三件（要加进底层库的全部）**：

- [ ] **B1 `src/quat/quaternion.hpp`**（~150 行）：四元数四则、↔旋转矩阵互转、Slerp。
      为谁服务：**Slerp 插值 demo**（匀速大圆）+ **万向锁 demo**（V2.5/V2.6，姿态演示的
      唯一缺口）。注意 D16：四元数封死 quat/ 目录、不做矩阵元素。验收 V2.5/V2.6 归
      `graphics/plan.md`（阈值不动）。
      - **2026-10-17 09:34 骨架已立**：`src/quat/quaternion.hpp`（新目录）+
        `demo/quat/quaternion.cpp`（ctest `quaternion`，32 项 FAIL 0；ctest 13/13 → **14/14**，
        `-c -w` 0 告警）。**已实现**：① 的 `+ - 一元- 标量乘` + ② 的
        `Conjugate / SquaredNorm / Norm` + `Identity()`。**全部占位（`throw std::logic_error`）**：
        `operator*`（Hamilton 积 = 第 ① 步的魂）、`Inverse` / `Normalize`、
        `FromAxisAngle` / `ToMatrix` / `FromMatrix`（Shepperd 四分支）/ `Rotate`、
        `Slerp` / `Nlerp` / `AngleOf` —— 共 **10 个入口**，用例用 [P1]–[P10]
        **如实钉住"仍不可用"**（照 `smoke.cpp` 对 `qr_inverse` 的规矩：实现掉一条，对应 [P] 变红，
        红就是提醒换成真判据）。四步划分、每步判据、易错点写在头文件注释里。
        域守卫：类内 `static_assert(is_floating_point_v<T>)`，实测 `Quaternion<int>` 与
        `Quaternion<Complex<double>>` 均被拦下并给出中文原因。
- [ ] **B2 `src/linalg/operator.hpp` 扩正体**（现在只有 Dot/Abs2 暂驻）：
      `LinearOperator` + `DifferentiationMatrix`（多项式基，移位阵）/ `LaplacianMatrix`
      （1D 三对角 / 2D 五点差分，Dirichlet）工厂。为谁服务：**弦模态 Python 版的对拍源**、
      鼓膜 Chladni（C4）、热模糊（C3）—— 全部阶段 4 demo 的矩阵本体；G 系列 Python
      手写的离散 Laplacian 将来要跟它对拍。验收 V2.7 归 plan.md。
- [ ] **B3 `src/linalg/classify.hpp` + 收编**：`Classify(A)` + `Inverse(A, hint)` 快速路径，
      把 `inverse.hpp` 里散住的 `is_diagonal/is_orthogonal/is_positive_definite` 搬进来
      收编成结构探测层（V2.3）；`qr_inverse` 回代已落地（S-App6 ✔），此刀只做搬家不动数值。
      为谁服务：让「结构 → 快速路径」的语义有个正式的家（教学价值 > demo 价值，顺位可后）。

**demo 侧（刀序）**：

- [ ] **D1 lie 刀 4/5/6**（`expm.cpp` 续写，半天）：θ 扫 [0,4π] 非单射 + **SU(2) 双覆盖**
      （θ=2π 处 `−I` vs `I`，lie-plan 收束判据）+ BCH 测阶（缩半降 8 倍）。
- [ ] **D2 `demo/pyapps/clt.py`**（D 扇面第一条，纯 numpy）——图已出过概念版，本刀是
      「有图 + 有判据」的正式版（φⁿ(t) 逐点趋 e^{−t²/6}）。
- [ ] **D3 G 系列 .py 源入库 `demo/pyapps/`**：图已出（`learning-todo/imgs/geo/` 7 件，
      子仓 `5580227`）但 **.py 源没进版本库**——按 pyapps 纪律（文件名=演示名、
      头部三行注释）收编 great_circle/geodesic_eq/parallel_transport 三件。
- [ ] **D4 quat demo**（等 B1）：Slerp 轨迹 + 万向锁复现（apps 或 ctest 按 V2.5/V2.6 定）。
- [ ] **D5 V1.7 rotate demo**（apps，纯演示，用户此前「先不用管」，随时可做）。

**不进本基线**（明确排除，防范围蠕变）：pybind11 绑定（等第一个联动 demo 触发）、
SDL2/gfx 骨架（等 Gate A）、lstsq/pinv 薄封装（S-App3/4，等 PCA/拟合 demo 真要）、
Hessenberg 非对称特征值（plan.md 明确缓做）。

**扇面线（并行，不阻塞主线）**：`graphics/plan.md` 的 **S-App1–5**（PCA / 牛顿法 / 伪逆 / lstsq /
旋转扇面）+ `learning-todo/fan-plan.md` 的四扇面。收束条件见 fan-plan §7 —— 扇面件必须带对拍断言，
否则就退化成抄 demo。

---

## 5. 里程碑

- [ ] Gate A 达成 → 启动第一个小项目（动力系统可视化）
- [x] ~~第一个小项目完成 → 拉 `viz` / `web` 多仓基建~~
      → **已提前至 2026-09-22**：不等 Gate A，直接以「笔记 web 化」启动两仓（理由见第 6 节日志）
- [x] ~~笔记 web 化 A1 垂直切片走完七步（`learning-web/AGENTS.md` §12 / `DESIGN.md` §5）~~
      → **✅ 已远超目标完成（2026-10-17 补记）**：不止 A1，**318/318 全部单元**入库、
      前端四层 + `dist/` 产物齐全，独立仓 139 次提交、2026-10-04 收工。
      （父仓此前指针停在「实现未开工」，是 2026-09-24 骨架期写下的过期状态，已订正。）

---

## 6. 进度日志（append-only：只追加，不修改历史条目）

### 2026-10-23（夜 · 三条定调：终点不止一个 / 数学跳着学 / 角色边界）

用户对上一轮三个风险给出答复，三条都改变了文档口径：

- **① 终点不止一个 —— 我上一轮「终点改成 GR」的判断作废。**
  用户原话：最终落点不止一个，目的是把**数学 / 数学物理 / CG / 编程**都学成，
  数学路线相当于**费曼学习法的输出端**。
  ⇒ `math-plan.md` 新增 **§0 定位：终点不止一个，路线是「共用的主干」**，
  并**重写 §7.1**：CG / GR / NS 是主干之上的**三条支线**，不是二选一；
  工程排序判据从「服务哪条支线」改为「**在不在共用主干上**」。
  §8 终点验收也加了说明：四句话只是 **GR 支线的样板**，CG/NS 各有验收。
- **② 数学是跳着学的、状态没持续记录，缺口自己慢慢补。**
  ⇒ `math-plan.md` §6 加前置说明：`ns-progress.md` 的「当前位置」**只能当「大概在哪」**，
  不是精确进度条；**段号是路线图不是进度闸门**，允许非连续跳跃（可同时摸段 0 与段 6）；
  **agent 不据此拦截或纠偏**。禁区清单与终点验收仍有效，但**只在用户主动要求盘点时使用**。
- **③ 角色边界（硬约束，已提升为最高优先级）。**
  用户原话：「**你永远只负责整个工程，以及文档，不负责代码，代码我会在其它 agent 里自己写。**」
  ⇒ 已写入 `MEMORY.md`**新的第零节**（最高优先级）+ 用户级 `MEMORY.md` 的「工作方式」首条。
  **具体含义**：不主动写 `src/` 实现（`.hpp` / `.cpp` / 演示 `.py`）；
  可做——设计契约（骨架/规格注释）、验收规格与阈值、代码审查、指出 bug、文档、排期。
  **上一轮我提的「我可以直接开始写 `ode.hpp`」作废** —— 那不是我的活。
  **与「要动手不要只给建议」的关系**：后者指**文档与工程决策要落地**，不是写代码的授权。
- **对「风险 3（段 3 三失效是否进库）」的答复**：用户「你自己权衡，我不懂」⇒
  **我的裁决：不进库，纯 pyapps 演示**（三条各一个脚本，成本近零）。
  理由：三失效的价值在**反例的直观**（$\{\sin nx\}$ 有界无强收敛子列、微分算子无界、
  谱离散 vs 连续），用 20 行 numpy 画出来比封一层 C++ API 有效得多；
  且 `Tensor`/`linalg` 侧没有消费者 —— 按「主干判据」它不在主线上。

**本轮文档改动**：`learning-todo/math-plan.md`（§0 新增、§6 前置说明、§7.1 重写、
§8 验收说明、章节编号重排 0–10）；`dev-progress.md`（分工行、§2 主干口径、本日志）；
`MEMORY.md`（第零节角色边界）；用户级 `MEMORY.md`（工作方式首条）。
**无代码改动**（符合新的角色边界）。

### 2026-10-23（晚 · 数学路线重排：十段表取代 L0–L5，工程排期随之改）

- **事件**：用户带来与另一 agent 讨论的数学路线（判据 / 主轴 / 十段表 / 不做清单 / 终点验收）。
  骨架完整采用，**新建真源 `learning-todo/math-plan.md`**（约 230 行）。
- **路线要点**：一条判据（**看懂允许跳证明、不允许跳定义**）+ 双线并行
  （分析线 A：测度→概率→L²→泛函→Sobolev→变分；几何线 B：流形→张量场→联络→曲率），
  **只在 L² 与变分两处交汇**，且 **A 必须先到变分，B 才能接场方程**。
- **三处冲突已落地订正**：
  ① **概率整块缺失** → 新增段 1（条件期望=正交投影、卷积=群前推）；
  ② **「Lee 前 11 章」是错的读法** → ISM 该读 **1–3、8–9、12（张量）、13（黎曼度量）**，
     曲率在另一本书（*Riemannian Manifolds*）；
  ③ **L2 泛函不该当关卡** → 真关卡挪到**段 5「变分两出口」**（§2 Gate 表已按此重定义，
     旧 Gate A/B/C 的数学条件作废，新增 Gate D = 张量/曲率）。
- **对本文件 §0 的改动**：数学行改为「段 0 进行中」并指向 `math-plan.md`；
  §2 Gate 表按段号重定义。
- **工程侧影响（math-plan §6 消费点映射，这是本次最值钱的产出）**：
  - 原 `plan.md` 的阶段顺序与数学段**几乎正交**：阶段 5/6（渲染、球谐、波包）在新路线里
    **没有消费点**；而段 1/2 最需要的 **`fourier.hpp` 原排在最后（阶段 6）** ⇒ **应大幅提前**。
  - **共用件只有五件**：`ode.hpp` / `pde.hpp` / `fourier.hpp` / `Tensor+Einsum` /
    `LaplacianMatrix` 是 **CG 与 GR 两条终点线共用**；分叉在「CG：渲染光追球谐」vs
    「GR：张量场联络测地线」。**建议先做齐共用五件，再按兴趣选分叉**。
  - 新增一件原计划没有的：**段 1 的 `lstsq`（条件期望=正交投影）** —— 原为 S-App4 扇面件，
    现在升为段 1 消费点。
- **未动**：`graphics/plan.md` 的阶段顺序（只给建议，实际调动由用户定；阈值一字不改）。
- 纯文档改动，无代码变更。

### 2026-10-23（Tensor 容器骨架落地 + 张量契约口径统一）

- **落地的代码（未提交）**：`graphics/src/core/tensor.hpp`（507 行）+ `demo/core/tensor.cpp`（347 行）+
  `demo/core/CMakeLists.txt`（注册 `tensor`，ctest 变 **17/17**，`tensor` = Test #2）。
  **契约 100%，实现 0%**：函数体全是 `throw std::logic_error("… 待实现 …")`；demo 侧 24 条
  `[P1]–[P24]` 专门断言"仍然 throw"，**实现掉一条就变红** —— 这是"故意绿"，别当成功能已完成。
- **契约翻案（`matrix.md` §6 第三次订正）**：从「不造 Tensor 类」改为「不造**胖** Tensor 类」。
  翻案逻辑是用户点破的：10-20 那条的理由是「存储层免费、运算层值钱」，**这个理由恰恰支持放行容器**——
  既然存储免费，「不值得写」不成立，只是「不值得写成**完整容器类**」。
  **放行的是容器那一半；运算那一半红线一字未动**：收缩 / 乘法 / 求逆 / einsum 一律自由函数，
  einsum 不进 `core/`。同步落 `plan.md` 红线 #4 新文本 + 决策 #21（修订 #9/#10）。
- **本次文档修补（6 处口径漂移 + 2 处悬空引用）**：
  - `.workbuddy/memory/MEMORY.md` §六整段重写 —— 原标题仍是「不造 Tensor 类」，与 10-23 定案**直接矛盾**，
    会让后续 agent 反过来劝阻立项。
  - `matrix.md` 宪法（68–70 行）「高维靠展平**不造张量**」→「高维靠展平（**Tensor 只存不算**）」；
    决策 #9/#10 补修订标记（此前只有 #21 反向括注"修订 #9/#10"，**被修订方无标记**）。
  - **悬空引用**：`matrix.md` 三处引 `§3.7`（该节**不存在**，§3 只到 §3.4）→ 改锚 `§2.3 黑名单 + §3.1`；
    两处引 `§6.1`（**不存在**，§6 在 476 行结束）→ 改指该节「具体形态」三条。
  - `plan.md` 暂缓清单：`Tensor<T,rank>` 那行改为「**编译期 rank** 仍暂缓，哑容器已立项」。
  - `plan.md` B2 那行补限定：此处"暂缓"**只约束工程验收**，理论线不暂缓（`math-notebook/partN.tex:171`
    把 B2 列为第一层承重墙）。两处此前看着矛盾，实为两个不同问题（代码要做到哪一步 / 书要不要读）。
  - `so3.hpp` 头注释：去掉「与 §六『不造 Tensor 类』同一口味」，并澄清决策 #21 立项的是**运行期长度**的
    `Vector<T>`，**不是**固定 3 维的 `Vector3`（否则与"不造 Vector3"看着打架）。
  - `demo/pyapps/tensor-convention.md` 头部：旧口径（"不进 C++ 库，因为不造 Tensor 类"）换新口径 ——
    仍不进库，但理由换成「它管的是 numpy 层的轴序/上下标/升降规约，C++ 记号另见 `tensor.hpp`」。
- **实测**：`ctest --test-dir graphics/build` → **100% tests passed, 17/17**。
- **遗留（未做，等用户口令 —— 张量设计）**：实现前三个待定 ① `strides_` 是否收敛为**单一真源**
  （`Tensor` 由 `ComputeStrides(shape_)` 派生、不给外部设置入口；只有 view 允许任意 strides）；
  ② 广播是否做成显式 `BroadcastTo` 自由函数（返回 stride = 0 的**只读**视图），而非进 `operator+`；
  ③ `shape_`/`strides_` 是否改「内联容量 4 + 溢出才堆」（动态 rank 不变，只是去掉 view 上的堆分配）。
- **附（结构问题，只记不搬）**：本文件 §2026-09-29 小节末尾还压着两条**误置**记录 ——
  2026-10-20（`so3.hpp`）、2026-10-23（`tensor-convention.md` 初版），没有自己的日期标题。
  按日期检索会漏，此处留指针。

### 2026-10-19（autodiff 布局翻案：分子格式 + 列向量输入输出）

- **用户令**：「重新搞我的 autodiff，换成分子格式和列向量输入输出」→「好的全改，代码文档全部同步，
  还要检查测试用例是否有问题以及是否符合规定」。这是 **D19 第二次翻案**（10-01 定、10-07 第一次修订）。
- **诊断阶段两条实测（决定工作量分布）**：
  - **`MatrixMult` 的伴随公式本来就与布局无关** —— 未改一行库代码，直接用列向量
    `W(2×3)·x(3×1)` 跑通：前向 `[6,15]`、`dW` 2×3 全 1、`dx=[5,7,9]`。原因：y=AB 的
    Ā=āBᵀ、B̄=Aᵀā 在两种约定下形式相同，翻案只换"谁是 A 谁是 B"。
  - **`Apply(A, f)` 的形状检查在列向量下失效**：`A.Cols()=1 == f.Cols()=1` 恒真 ⇒ 长度不符
    不抛错、按 `A.Size()` 读 `f.Data()[i]` **越界**。这是潜伏 bug，不是风格问题。
- **库侧改动（2 行 + 1 个新函数）**：
  - `jacobian.hpp`：形状检查换成 `CheckFuncMatrixShape`（1×1 广播 / 否则**逐维**同形），
    报错信息带两边形状。
  - `jacobian.hpp` 新增 **`Jacobian(f, x[, m])`**：分子布局 `J_ij = ∂y_i/∂x_j`，形状 m×n，
    一次 sweep 填**第 j 列**；x 非列向量 / 输出非列向量 / m 不符均 throw。
    —— 此前 `jacobian.hpp` 里**根本没有 `Jacobian()` 函数**，分子布局长期是句没有代码承载的口号。
  - `grad_graph.hpp`：**代码零改动**，只改注释（布局说明块、MatrixMult 的 Forward/Backward 推导）；
    顺带把 `Inverse` Forward 过期的「待 linalg::Inverse 落地」改成「本节点未接入（linalg::Inverse 已可用）」。
- **测试全部转列向量**：`apply.cpp`（11 处非 1×1 构造 + [F] 形状判据改 Rows）、
  `backward.cpp`（helper `MV` + 13 处维度实参 + Transpose 形状 + **MatrixMult 换成非方阵
  W(2×3)·x(3×1)**，两侧都查）、`compute_graph.cpp`（**整文件重写**：3 处构造 +
  56 处 `[0][k]`→`[k][0]` + `Cols()`→`Rows()`，顺带 **float→double**，
  差分步长 h=1e-2 → 1e-6、eps 1e-3 → 1e-8）。
- **新增 ctest `jacobian`（18 项）**：形状 m×n + 逐元解析 + 中心差分 + **链式法则自然序
  `J_{f∘g}=J_f·J_g`**（分子/分母的判别式）+ **图反向逐元对拍**（seed e_i ⇒ x.adjoint 是 J 第 i 行）
  + 退化（m=1 的 1×n 行 / n=1）+ 非法输入 throw；double 与 float 双域。
- **审计发现并修掉的用例问题（3 条）**：
  1. **`[P1] operator* 仍是占位` 已红** —— 用户 10-17~10-19 之间实现了 Hamilton 积。
     这正是占位判据的设计意图（实现即变红）。已按设计**升级为真判据 [C1]–[C14]**：
     9 条基元乘法表 + 反交换 + 单位元/结合律/分配律 + **共轭反自同构 (ab)*=b*·a*** +
     **范数乘性 |ab|²=|a|²|b|²** + 同轴相加 q(θ1)q(θ2)=q(θ1+θ2) + **双覆盖 q(θ+2π)=−q**。
     理由：叉积符号写反时往返测试抓不到（正反抵消），只有这组能钉死约定。`quat` 32 → 59 项。
  2. **`backward.cpp` 只跑 `Dual<double>` 一个域** —— 违反 graphics 3.3 第 1 条（分域探针）。
     10-17 对 float / `Dual<Complex<double>>` 的验证只存在于 `/tmp` 探针、**没进 ctest**。
     已补 `[F]` 段：同一份 grad_graph 换域再实例化，两域各跑最小闭环。29 → 31 项。
  3. **`apply.cpp` 的 [F] 段只测"不匹配会抛"，没有反向对照**。已补 [F3]：1×1 函数矩阵
     必须**广播不抛** —— 否则"见到非同形就抛"也能让 [F1] 通过，判据不独立。
- **验证**：ctest **15/15**；`-c -w` 从零 + `-Wall -Wextra` **0 告警 0 错误**。
- **文档同步**：`matrix.md` §7.4 全节重写（保留 10-07 历史 + 新增 10-19 第二次修订六条）+
  决策表 #19 + 头部增量 + §6 雅可比表 + §6 Field 展平 + §7.5 第 2 条接口公式；
  `dev-progress.md` §0 头部/工程行/§3/§6；`MEMORY.md` 新增 §3.10 布局约定。
- **未提交**，`HEAD = aac7b27`。

### 2026-10-17（晚 · 状态真源订正：learning-web 早已完成）

- **订正事件**：用户指出「learning-web 不是完成了吗」——核查属实，**是父仓真源漏记**。
  - **实况（`learning-web` 独立仓核验）**：HEAD `3ec63a6`（2026-10-04 14:33），
    **提交总数 139**；`json/` 词条 **318/318 全入库**（与 `HANDOFF-json.md` 的「85/318」
    快照已隔 233 单元）；前端四层全实现（`src/shell` / `renderers` / `modules` / `viz`
    含 27+ 档1 模块）；`dist/` 产物齐全（`index.html` / `knowledge.html` / `graph/` / `json/`）；
    `git status` **干净**（无未提交）。
  - **根因**：本文件 §0 指针表「笔记 web 化」行是 **2026-09-24 骨架期**写下的
    「实现未开工」，此后 web 仓在 **9/24–10/04 密集推进并收工**（9/25 单日 93 次提交），
    但**从未回写父仓** —— §5 里程碑的 `[ ] A1 垂直切片走完七步` 也一直空挂着。
    本次按实况补齐 §0 该行 + §5 该条，并同步 ctest 计数 13/13 → **14/14**（`quat` 骨架已立）。
  - **教训**：独立子仓的完成事件**必须由该仓会话主动回写父仓**，否则父仓真源会长期停在
    启动态；引用子仓状态时**以 `git log` + 目录实况为准，不照抄父仓摘要**（本次即栽在这）。
  - 纯文档改动，无代码/数据变更。

### 2026-10-17（上午 · R/G 两系列落地后的文档同步）

- **实测盘点（09:23）**：`git log` 显示 `aac7b27`（expm_pade，**V2.4 Padé 一半 ✅**，
  [E] 段 16→25 项：相对差 5.3e-16 ~ 1.5e-13、幂零逐位精确）已提交；**pyapps 两系列整目录
  untracked**（sphere G1–G6 + rotation R1–R6，共 14 个 .py + 2 份 report.md）；
  `learning-todo` 子仓 `imgs/rot/` 6 张图未提交、`imgs/{deriv,diff,matrix}/` 三目录 untracked。
  fan-plan.md §3.1/§3.2 已含 G/R 实测状态与两条口径警告（G3 极角/纬度混用、
  R 系 acos 病态→atan2），本会话此前未同步 §0 —— 现已对齐。
- **文档更新**：§0 头部（阶段 2 进行期 + V2.4 + pyapps 十二脚本）、工程行（含
  「R2 已把 lie 刀 4 双覆盖的 Python 前哨做掉」「`_core.py` 是 ode.hpp 的现成对拍源」
  两个新事实）、扇面行（A 面 3 条 ✔ 超额，B/C/D 0 条 → D2 CLT 下一优先）、
  李群行（刀 4 前哨状态）。
- **给用户的提交建议（攒批分组）**：① 父仓 `graphics/demo/pyapps/`（G+R 十四脚本两报告，
  一刀「pyapps: 球面 G1–G6 + 旋转 R1–R6」）② 子仓 imgs（rot 6 张 + deriv/diff/matrix
  三目录 + fan-plan.md + matrix.md 的 real() 增量——matrix.md 的改动是 10-14 就有的工作区
  遗留，与 fan-plan 可分开两刀）③ `dev-progress.md`/`plan.md` 的文档改动随手带。
- **下一步排序（见 §0 工程行）**：lie 刀 4/5/6 → ode.hpp（Gate A 钥匙，`_core.py` RK4
  即对拍源）→ LinearOperator → quat C++ 化。

### 2026-10-17（凌晨 · 下一版基线定版）

- **扇面分层调整 + G 系列立项（文档层，无代码改动）** —— 用户两项决定落地：
  ① **viz 层暂缓启用**：纯 Python 数学演示改挂 **`graphics/demo/pyapps/`**（与 `demo/apps/`
  同构：不进 ctest、无 CMake、文件名 = 演示名）；「Python×C++ 联动才去 viz 层」，
  触发条件 = 第一个联动 demo（届时按既定方案建 `graphics/bind/linalg.cpp`：
  pybind11 绑 7 算法 × 双域、Matrix 走 numpy↔Matrix 壳函数转换、`run_bind.sh`）。
  ② **G 系列（球面测地线）立项**：用户选题「三维球体（地球）测地线入门旋转、张量、曲率」——
  G1 大圆 vs 纬线 / G2 测地线方程（RK4 预演）/ G3 平行移动 holonomy（转角=2π·cosφ）⭐ /
  G4 Gauss–Bonnet（球面三角形内角和>π + 网格离散版，即 V6.4 的 pyapps 预演）/
  G5 曲率张量 2×2 实例。依据：`learning.md` 204 行手推题 + 607–608 Gauss–Bonnet 数值验证
  + 阶段一 1.2 清单「曲率张量、测地线」；数学正课仍在 `learning.md` 阶段一（Lee 前 11 章）。
  - 改动：`learning-todo/fan-plan.md` §1 分层表加 pyapps 行 + viz 暂缓说明与触发条件、
    §3 加「↳ G 系列」行 + **§3.1 G 系列五刀表**、§9 开局三刀按最新状态修订
    （A1 C++ 版已出 / S5 ✅ 转 lie 刀 4/5/6 / D2 挂 pyapps）；本文件 §0 扇面行同步。
  - fan-plan.md 在 learning-todo 子仓（与 matrix.md 改动同处未提交区），提交须在该子仓进行。
  - **git 状态提醒（父仓，23 项未提交）**：V1.5/V1.1 两刀（dual_complex/smoke/apply 等
    ctest 11/11 的成果）仍在工作区（HEAD `6999ca2` 之后无提交）；另有 11408-notebook
    新增 4 ipynb + modules/408 两文件 + learning-408 整目录 untracked。攒批提交时按刀分组。

- **autodiff 反向模式落地：`grad_graph.hpp` 修 4 类 bug + 新增 ctest `backward`（12/12 → 13/13）**
  —— 用户令「落实吧，最小化修改我的逻辑，并且把我的 bug 要说出来」。承接本日诊断（见下条），
  把反向模式从「能跑但基本是错的」修到「12 个可用算子全对」。
  - **`GradNode::Forward` 一行未动**（只修 `Backward` 侧 + `ComputationGraph::Backward`）。
  - **bug 清单（全部实测确认；每条都写进了 `demo/autodiff/backward.cpp` 的头注释）**：

    | # | 类别 | 位置 | 症状 |
    |---|---|---|---|
    | 1 | 架构 | `ComputationGraph::Backward()` | 直接 `throw` 拒绝执行 |
    | 2 | 架构 | `GradNode::adjoint` | **从不按形状初始化**（默认 0×0），而所有分支都写 `+=` ⇒ 第一个节点就抛 |
    | 3 | 静默 | `case Mult` **缺失** | Forward 里 Mult ≡ CWiseMult，Backward 只写了 CWiseMult ⇒ **乘法节点反向静默 no-op** |
    | 4 | 静默 | `case Output` **缺失** | 透传节点不回传伴随 |
    | 5 | 静默 | `case Transpose` | `node1->adjoint = node1->matrix.Transpose()` —— 用了**原值**而非上游伴随，且 `=` 而非 `+=` |
    | 6 | 抛异常 | `case Divide` 的 node1 | 逐元素除法误写成矩阵乘 `adjoint * node2->matrix.Transpose()` |
    | 7 | 抛异常 | `case MultAll` | 1×n 与 1×1 做 `cwiseProduct`；且偏导在 Dual 的 `derive` 分量，需先 `Derive` 提取 |
    | 8 | 抛异常 | `case Function`/`MatrixFunction` | `adjoint * Derive(Apply(...))` 用了矩阵乘；逐元素函数应逐元素乘 |
    | 9 | 语义 | 同上 | `Derive(Apply(A,f))` = **f'(x)·A.derive**（依赖种子），不是 f'(x) ⇒ 改用 `LocalDeriv` |
    | 10 | 静默 | `case Sigmoid/Relu/Tanh/Softmax` | case 写了但 **body 是空的**，`-Wswitch` 抓不到（**本次未动**，Forward 侧本来就是占位） |
    | 11 | 注释过期 | `case Inverse` 的 Forward | throw「待 linalg::Inverse 落地 (阶段 1)」，而它 **V1.2 就落地验收了** |

  - **改动（`grad_graph.hpp` +71/−17；`jacobian.hpp` 加 `localFunc`/`LocalDeriv`）**：
    ① `ComputationGraph::Backward()` 换真实现（**形状初始化 → 输出播种 → 拓扑逆序回传**），
    并加 `Backward(const Matrix<T>& seed)` 重载；② 补 5 个 case（`Input`/`Constant`/`Output`/
    `Mult`/`CWisePow`），**顺带 `-Wswitch` 归零**；③ 修上表 §5–§8 四处的公式；
    ④ `jacobian.hpp` 加 `localFunc`（局部导数标量分派，7 行）+ `LocalDeriv`（与种子解耦）。
  - **`LocalDeriv` 为什么不用 `Derive(Apply(...))`**：实测 `x.derive` = 0 / 1 / 3 会让
    `Derive(Apply(x, sin))` 得 **`[0,0]` / `cos(x)` / `3·cos(x)`** —— 局部导数**必须与种子无关**。
  - **新增 ctest `backward`（29 项 FAIL 0）**：`[A]` 最小闭环 `y = Σ(w ⊙ x) + b` 三方对拍
    （解析 / 中心差分 / **前向模式**）+ `[B]` 12 个可用算子逐条 + `[C]` adjoint 自动初始化 +
    `[D]` 幂等 + `[E]` **如实断言占位算子仍不可用**（照抄 `smoke.cpp` 的规矩，不假装它可用）。
  - **两条注入对照**：① 拿掉 `case Mult` → **只有 `backward` 挂，且只挂 [B3]**；
    ② 删掉形状初始化 → `backward` **退出码 134（SIGABRT）**，正是原版「一跑就抛」的活体复现。
    恢复用 `cp` + `shasum`（`3fd3f50a…` 一致）、`grep src/` 0 残留、复跑 13/13。
  - **验证**：`./run_tests.sh` **13/13**；`-c -w`（从零 + `-Wall -Wextra`）**告警 0、错误 0、13/13**。
  - **未动**：8 个 `throw` 占位（CWisePow / Inverse / Sigmoid / Relu / Tanh / Softmax /
    CrossEntropy / LnMultAll）按用户令「占位的先不管」保持原样 —— 但 `[E3]` 已**如实断言**
    `Sigmoid::Backward` 仍是**空 body 静默 0**（将来实现后该判据会 FAIL，提醒更新用例）。
- **口径修正：「要不要 tape」（承接同日讨论）** —— 上一轮把 §7.5 的「顺手存局部导数」讲成了
  硬约束，**是我说错了**。用户指出：反向时重算 `f'(x.value)` 与「前向时存」**等价**（纯函数），
  「存」只剩「省一次 O(n) 计算」；且本架构 `shared_ptr` 持有 `matrix` ⇒ 中间量一直在 ⇒
  **永远可重算**。tape 真正的独占优势是**省内存**（不保留全部中间量）与动态控制流。
  §7.5 三条重新分层：**① 局部导数与种子解耦 = 硬（正确性）**；② 每算子第二份代码 = 事实但免不掉；
  ③ 伴随独立存储 = 工程取舍。用户的取舍（`adjoint` 仍用 `Matrix<Dual<U>>` 以保 `GradNode<T>`
  单模板参数）**成立**，代价只有 `MultAll` 需 `Derive(...)` 提取一次。

- **V2.4 完成：`expm_pade`（[13/13] Padé + scaling）落地（代码改动 + 用例扩展；ctest 仍 13/13）**
  —— 承接用户令「接着往下研究 expm 的拓展」。
  - `src/linalg/expm.hpp`（54 → 125 行）：新增 `expm_pade` + `#include "linalg/solve.hpp"`。
    与原 `expm`（Taylor + scaling 到 0.5）**并列**、互为对照基线 —— 符合本库「多条路径并列、
    互不复用实现」的风格（同 `lu_inverse` / `cholesky_inverse` / `qr_inverse`）。
  - **关键差异是 scaling 阈值**：Taylor 折半到 ‖As‖ ≤ 0.5；Padé [13/13] 折半到 θ₁₃ ≈ 5.372
    ⇒ **s 少约 3 个数**（平方少 3 次）。代价是每步要解一个线性系统（`Solve`，不显式求逆）。
  - `demo/linalg/expm.cpp` 新增 [E] 段 7 条（16 → **25 项**）：相对差 5.34e-16 / 2.52e-15 /
    **1.54e-13（5×5, |a|≤100）** / 1.14e-15（复 4×4）；Padé→Rodrigues 1.11e-16；幂零**精确 0**。
  - **踩坑（又一次老坑）**：第一版判据用了**绝对差** ⇒ `|a|≤100` 的 5×5 报 `1.1e+29 FAIL`，
    那是**假警报**（元素量级是 e^‖A‖）。改相对判据后同一组只有 1.54e-13。
    **`matrix.md` §3.6 早写过「expm 判据必须相对」，这次是没照做。**
  - **操作事故（已修）**：`expm.cpp` 的 [E] 段被插入**两次**（同一脚本对该文件执行了两遍，
    「读—改—写」不幂等）⇒ `redefinition of 'E_pade_vs_taylor'`、构建失败。
    修法：`git show HEAD:<path>` 取回干净基线 → 重插一次。
    **教训：「读—改—写」类脚本必须幂等**（插入前 `assert 目标串不存在`），否则重跑即毁。
  - **验证**：`./run_tests.sh` **13/13**；`-c -w` **告警 0、错误 0**。

### 2026-10-16

- **V1.5 验收通过：`demo/linalg/dual_complex.cpp`（代码改动；ctest 8/8 → 9/9）** —— 用户令「按照计划来」，
  并先问「V1.5 再解决什么问题」。
  - **在解决什么（三层，不只是"复数求导"）**：① **工程** —— 泛型域覆盖的最后一格（`double` /
    `Complex` / `Dual<double>` 都验过，`Dual<Complex<double>>` **从未进过编译器**）;
    ② **数学** —— `Dual` 的规则 `(fg)' = f'g + fg'` 隐含「导数唯一」；实数域可导即唯一，
    **复数域只有全纯函数才唯一**，非全纯上 Dual 给的是「沿某个方向的数」（随方向变）⇒ 是错答案；
    ③ **用法** —— 非全纯场合的补救 = **双种子**（两个线性无关方向各传播一次再重组）。
  - **代码改动两处**：① `demo/common.hpp` 补 `RandFill<Dual<Complex<U>>>` 特化 —— 原来复用
    `RandFill<Dual<U>>`（其内写 `static_cast<U>(a)`），`U = Complex<U'>` 时 **`Complex` 无单参构造、
    无默认参数** ⇒ 实测 `no matching conversion for static_cast from 'double' to 'Complex<double>'`。
    偏序上 `Dual<Complex<U>>` 比 `Dual<U>` 更特化，两条并存即可。② 新增 `demo/linalg/dual_complex.cpp`
    + `demo/linalg/CMakeLists.txt` 注册四行。
  - **判据 6 段 / 22 项全过**：[A] 泛型域（含**矩阵层**链式法则 `d(A²) = A·A' + A'·A`）；
    [B] 跨域一致（derive 置零时 `Solve` 的 value 与纯 Complex **精确 0**，n=2/4/6）；
    [C] 全纯 z²（两种子 `derive/s` 一致 = 2z₀）；[D] 非全纯 |z|²（双种子重组出 2z̄，残差 0）；
    [E] **实用形态**（解 `A(θ)x = b` 的 `dx/dθ` 自动传播：解析 −A⁻¹A₁x 5.2e-18 / 中心差分 6.9e-12）；
    [F] 边界（1×1 / 非方阵 throw / 主元近奇异 throw）。
  - **注入 3 条**：① `ScalarTraits<Dual<T>>::conj` 漏 `conj(derive)` → **1 条 FAIL（残差 2.0，与手算一致）**，
    其余 8 用例不受影响；② `Dual::operator*` 漏 `derive*other.value` → **4 条 FAIL 全在复数域**
    （同时 `scalar_traits`/`compute_graph` 挂 —— 那两个测的是实数域乘法，本刀补的正是**复数域覆盖**）；
    ③ 把 D 段的「非全纯」换成 `z*z` → 首轮只挂 2 条，**暴露我自己的判据漏洞**：原判据比的是**原始**
    `derive`，而全纯的两种子原始值也恰好不同（2z₀ vs 2z₀·i）⇒ 漏判；改为比归一化后的 `derive/s`
    后 2 → **3 条 FAIL**，漏洞补上。注入全部恢复，两个库文件 `shasum` 与备份一致、`grep 注入` 0 残留。
  - **文档同步**：`graphics/plan.md` S 清单插入 V1.5 完成条（S6b 后 / S7 前）；本文件 §0 头部 +
    §3 V1.5 行 + §4「剩余四步 → 三步」（① 划掉，②③④ 重编号 ①②③）。
  - **未提交**：改动仍在工作区（`HEAD` 仍 `bae45df`）。

- **V1.1 显式闸门落地：`demo/linalg/smoke.cpp`（代码改动 + 文档同步；ctest 10/10 → 11/11）** ——
  用户令「v1.7 那个 rotate demo 先不用管，看接下来要搞什么」，故选 **V1 最后一个缺硬证据的节点**
  （V1.1），即 §4 的第 ② 步。**V1.7 是唯一不作硬验收的演示件，不阻塞收口。**
  - **在解决什么**：V1.2–V1.6 各自只覆盖**自己用到的**入口（solve/inverse/cholesky/eigen/svd/qr/expm）。
    V1.1 要的证据不是「某个入口算得对」，而是「**每一个公开入口都在 double 与 Complex 下各实例化过一次**」
    —— 与 2026-10-14 那次「528 行从未被任何 TU include」是同一件事的收尾：模板 bug 只在实例化时暴露。
  - **两段设计（不合并）**：**[A] 实例化闸门** —— `EntryTable<T>` 对 40 个公开入口（27 自由函数 +
    13 个 `Matrix<T>` 成员，含 `Inverse` 两个重载的显式消歧 `static_cast`）**取函数指针**（odr-use
    ⇒ 强制实例化），刻意**不调用**；配 `kEntryCount = 40` 计数，漏抄一行会被抓。
    **[B] 最小数值闸门** —— 每个入口跑一次 + 一条轻量判据（阈值 1e-10，比专项用例松）。
    合计 **54 项 / FAIL 0**：double 侧 Solve 5.8e-17、Inverse 2.2e-16、QR 重构 3.4e-16；
    complex 侧 Solve 5.8e-17、Inverse 4.4e-16、expm 自洽 6.7e-16。
  - **顺带结算的支线**：① **`Dot` 零调用方** —— 本闸门成为它的第一个调用方，判据 `<v,v> = Σ|v_i|²`
    把它与 `Abs2` 钉在一起（顺带覆盖「复数域必须共轭第一个参数」这个被测过 4 次的「double 免疫」点）；
    ② **`qr_inverse`** —— 闸门**如实断言它 throw `logic_error`**（不假装它可用），补回代留给 S-App6；
    ③ **`MatrixType` 孤儿** —— 不属 V1.1，仍在 §4 挂着等定案（建议删，结构类归 `MatrixKind`）。
  - **注入反例两条（本次最值钱的证据）**：
    - **① 行为型**：`operator.hpp::Dot` 去掉 `conj` → 全项目 **只有 smoke 挂**（1 条，[B12] complex，
      残差 2.048e-01），其余 **10 个 ctest 全绿** —— 因为 `Dot` 此前零调用方，
      **这个 bug 在旧测试体系下完全不可见**。
    - **② 编译型**：往 `Dot` 体内写一句只在实例化时非法的 `a.OnlyAtInstantiation()` →
      只 `#include "linalg/operator.hpp"` 的探针**编得过**（证明「include 过」≠「编译器看过」），
      而 `run_tests.sh` 在 **smoke.cpp 处编不过**、**其余 10 个 target 全部照常 `Built`** ——
      正好把「实例化才暴露」这件事两面都钉死。
      恢复用 `cp` 备份 + `shasum` 比对（`7e9a1092…` 前后一致），`grep 注入` 0 残留，复跑 11/11。
  - **验证**：`./run_tests.sh` **11/11**；另跑 **`-c -w`**（从零 + `-Wall -Wextra`）→ **告警 0、错误 0、11/11**。
  - **🆕 顺带发现（未修，不属本刀）**：`autodiff/jacobian.hpp` 的 `Apply<T>` 在 **double 域探针**下
    报 `use of undeclared identifier 'dualFunc'`（`:22`/`:42`）—— 当时据此记为「两个重载都编不过」。
    ⚠️ **该结论已于本日第二次探针更正**（见本日「autodiff 域契约勘定」条）：`Apply` 在
    **`Dual` 域是好的**，只是标量域越界；`Derive`/`SumAll`/`MultAll`/`MultAllBackward` 在 Dual 域
    实例化干净。autodiff 按 §1 归用户自行推进，只记录不改动。（探针：`/tmp/probe_jacobian.cpp`）
  - **文档同步**：`graphics/plan.md` S 清单补 V1.1 条 + 待办追加 jacobian 发现；本文件 §0 头部与
    工程行、§3 表格与 V1 落点、§4「三步 → 两步」。**未提交**（`HEAD` 仍 `6999ca2`）。

- **S6 三条派支线结算（代码改动 + 文档同步；ctest 仍 11/11）** —— 用户令「MatrixType 我放到 demo 的
  common 里了，Dot 是啥，还有 qr_inverse 那个你也搞一下」，三件一并收口。
  - **① `MatrixType` → 应用层（用户定案）**：从 `src/core/common.hpp` 挪到 `demo/common.hpp`
    （用户改的，core 侧已删除），core 层只留按**代数结构**分的 `MatrixKind`。原注释「不能按用途去分，
    怎么用是上层的责任」即理由 —— 与 `matrix.md` §3.7 的实测结论同向（不给稠密 `Matrix` 存运行期
    结构标签）。代码侧零引用（`grep -rn MatrixType src/ demo/` 只剩定义处 + smoke 的说明注释），
    对编译无影响。
  - **② `qr_inverse` 补回代（S-App6）** —— `inverse.hpp` 的 `throw` 占位换成真实现（+27 行含注释）：
    `A = Q·R ⇒ A⁻¹ = R⁻¹Qᴴ`。用 `QR_Householder`（新增 `#include "linalg/qr.hpp"`；qr.hpp 只
    include core/* + operator.hpp，**不反向依赖 inverse.hpp ⇒ 无循环**），对上三角 R **直接回代**
    求 R⁻¹（右端项 B = I 原地被解成 R⁻¹；**不复用通用 `Solve`**——R 已是上三角，再走一遍 LU 等于
    白花一遍 O(n³)；也不显式求逆再乘），末尾一次 `X * Q.Adjoint()`。奇异判据 `|R_ii| <= DefaultTol`
    即抛（与 `lu_inverse` / `diagonal_inverse` 同口径 D12）。
    它**不接 `Inverse` 分发**：`MatrixKind` 里没有、也不该有「QR 档」（QR 是算法不是代数结构），
    故与 `lu_inverse` / `cholesky_inverse` 并列，是一条显式点名的独立路径。
  - **验收**：`demo/linalg/inverse.cpp` 新增 **[H] 段**（+8 项 → `inverse` **85 项 / FAIL 0**）：
    随机 n=2..5 残差 + **与 `lu_inverse` 逐元对拍**（只判残差不钉住"算的是同一个逆"；两条独立实现
    互为判据）+ 复 2×2（D7 老坑：`Qᴴ` 极易写成 `Qᵀ`）+ 奇异 / 非方阵 / 0×0 三条 throw 路径。
    实测残差 double ≤ **1.11e-15**、complex **2.48e-16**，与 `lu_inverse` 逐元差 ≤ **6.22e-15**；
    汇总行加了 `qr_inverse` 一列。`smoke.cpp` 的 [B6d]（原「如实断言 throw」）改为数值判据，
    并新增 [B6e]（零矩阵 → R 对角奇异 → throw），`smoke` **54 → 56 项**。
  - **注入反例两条（本次最值钱的证据）**：① `Q.Adjoint()` → `Q.Transpose()`（D7 经典坑）→
    **只有复数域挂**（`inverse` [H]② 残差 8.889e-01、`smoke` [B6d] complex 残差 4.088e-02），
    double 侧全绿、11 个 ctest 里只红 `inverse`/`smoke` —— **第 5 次「double 免疫」现场**；
    ② 删掉奇异判据 → **恰好只挂那 3 条 throw 判据**（`inverse` [H] ×1、`smoke` [B6e] ×2），
    其余全绿。恢复用 `cp` + `shasum`（`f03aea15…` 前后一致）、`grep -r 注入 src/` **0 残留**。
  - **验证**：`./run_tests.sh` **11/11**；另跑 **`-c -w`**（从零 + `-Wall -Wextra`）→ **告警 0、错误 0**。
  - **文档同步**：`graphics/plan.md` S6 条、S-App6 待办（勾选）、V1.1 条副产物；本文件 §0 工程行、
    §3 表格（`inverse.hpp` 161 → **188** 行、`inverse` 78 → **85** 项；linalg 合计 1062 → **1089** 行）
    与 V1 落点、§4 S6 支线条。**未提交**（`HEAD` 仍 `6999ca2`；工作区里的 `demo/common.hpp` /
    `src/core/common.hpp` 是用户自己改的 MatrixType 迁移）。

- **autodiff 域契约勘定（纯诊断，零代码改动；ctest 仍 11/11）** —— 用户令「先把 autodiff 调通吧，
  不改代码，告诉我要怎么搞」。**结论：autodiff 没坏 —— 是域契约没写进签名 + 覆盖全缺。**
  - **推翻了本日早先的一条记录**：`Apply<T>` **不是「两个重载都编不过」**，而是
    **只在 `T = Dual<U>` 域编得过**。早先结论出自 double 域探针 —— 与 §3.3 第 1 条「分域探针」
    同源、**方向相反**：那次是「double 探针 0 错掩盖真 bug」，这次是「double 探针报错」被误读成
    「整个模板坏了」。**同一个坑的第二次踩法。**
  - **机制（实测）**：`dualFunc` 是 `Dual<T>` 的 **in-class friend** ⇒ **ADL-only**、无法前置声明。
    `Apply` 体内的 `dualFunc(A.Data()[i], f)` 是**依赖调用**，名字查找 = 定义点普通查找（找不到）
    + **实例化点 ADL**。`T = Dual<U>` ⇒ ADL 命中 `Dual<U>` 的友元 ✅；`T = double`/`Complex`
    ⇒ 无关联类含此友元 ⇒ `use of undeclared identifier 'dualFunc'` ❌。
  - **实测矩阵**（同一份头文件，只换域；`clang++ -std=c++20 -I src -fsyntax-only`）：
    | 探针 | 域 | 结果 |
    |---|---|---|
    | `Apply(A, FunctionType::Sin)` | `double` / `Complex<double>` | ❌ `undeclared identifier 'dualFunc'` |
    | 同上 | `Dual<double>` / `Dual<float>` / `Dual<Complex<double>>` | ✅ 全过（退出码 0） |
    | `Apply(A, Matrix<FunctionType>)` | `Dual<double>` | ✅ 过 |
    | `GradNode<double>::Forward()` | `double` | ❌ 4 条（`SumAll`/`MultAll` no matching + `Apply<double>` ×2） |
    | `GradNode<Dual<float>>::Forward()` | `Dual<float>` | ✅ 过 |
  - **数值已验**：`Apply<Dual<double>>` 五个 `FunctionType` 的 `derive` 与解析导数**逐条误差 0**
    （seed=1：Sin 0.764842…/Cos −0.295520…/Exp 3.320116…/Log 0.5/Sqrt 0.267261…）——
    所以 `Apply` 在它自己的域里是**正确**的。
  - **真正缺口是覆盖，不是正确性**：① `compute_graph` 只跑 Input/Mult/Add/Output/SumAll 五种
    GradOp —— **`Function`/`MatrixFunction` 从未跑过**，`Apply` 的数值行为**零覆盖**（连前向都没验）；
    ② `GradNode<T>` 整体**只支持 `Dual<U>`**（其余签名早写死 `Matrix<Dual<T>>`），
    但**没有任何地方声明这一点**，越界时报错分散在四个位置、信息晦涩。
  - **三条处置路线（用户选）**：① 给 `Apply` 加域约束（`requires`/`static_assert` 把 `Dual` 写进
    签名）——错误信息从「undeclared identifier」变「只支持 Dual 域」，**最小改动**；
    ② 给标量域补一层逐元素分派（走 `std::sin` 等，复数用 ADL）——**扩能力**；
    ③ 只补闸门不动实现。**三条都缺同一件东西**：autodiff 实例化烟测（照抄 `EntryTable`，
    但**按域铺开**——同一模板在每个适用域各取一次函数指针）+ 一个 `Function` 节点数值用例。
  - **文档同步**：`dev-progress.md` §0 工程行 / §3 新发现段 / §4 待办条 / 本 §6 条；
    `graphics/plan.md` 待办条 —— 均把「两个重载都编不过」更正为域契约表述、附**两条对照探针**
    （`/tmp/pa.cpp` double 域报错 vs `/tmp/pb.cpp` Dual 域 0 错）。**零代码改动**，`HEAD = 6999ca2`。

- **autodiff `Apply` 域分裂落地 + `dual.hpp` log 断链修复（代码改动 + 新 ctest；11/11 → 12/12）**
  —— 用户令「你改吧」，在上一条给出的三选一里走**路线 ②**（标量域补逐元素分派）：它是唯一
  「既把域契约钉死、又不丢能力」的一条；路线 ①（只加约束）会让标量域彻底不可用。
  - **`src/autodiff/dual.hpp`（152 → 169 行）**：① 新增域判别 trait `IsDual`（前置声明 + 特化，
    置于 `Dual` 类之前），供守卫使用；② 新增 `friend Dual log(const Dual&)` 转发 `ln` ——
    **修一个真 bug**：`ln` 体内写的是 `log(x.value)`，而 `Dual` 的友元叫 `ln` 不叫 `log`，
    对 `T = Dual<U>`（嵌套对偶 = **二阶导**）ADL 断链（`no matching function for call to 'log'`）；
    `pow` 体内同样用 `log`，一并受益。
  - **`src/autodiff/jacobian.hpp`（100 → 162 行）**：新增 `elemFunc`（标量域逐元素分派，
    `using std::sin;` + ADL）、补 `#include <cmath>`/`<type_traits>`；`Apply` 拆成 **4 个重载**
    —— `{FunctionType, Matrix<FunctionType>} × {Matrix<Dual<T>> 版, Matrix<T> 版}`。
    Dual 版走 `dualFunc`、标量版走 `elemFunc`；标量版体内 `static_assert(!IsDual<T>::value, ...)`
    作守卫。**分流机制**：`Matrix<Dual<U>>` 比 `Matrix<T>` 更特化，[temp.deduct.partial] 部分排序
    必然把 Dual 域判给 Dual 版 —— 标准行为，非编译器扩展。
  - **新增 ctest `apply`（`demo/autodiff/apply.cpp`，17 项 FAIL 0）**：4 重载 × 4 域 +
    Function/MatrixFunction 节点真实路径 + 嵌套对偶二阶导 + 维度不匹配 throw。
    **关键数字**：derive 对拍解析导数 **50 组误差全 0**；中心差分对拍最差 **3.666e-09**（阈值 1e-8）；
    `ln(2)` 二阶导 **−0.25**（= −1/x²）、`x³` 二阶导 **12**（= 6x）；Complex 域对拍 `std::complex` ≤ 1e-15。
  - **四条注入对照（本次最值钱的证据）**：
    | # | 注入 | 结果 |
    |---|---|---|
    | ① | 行为型：`elemFunc` 的 `Log` 错映射成 `exp` | **只有 `apply` 挂，其余 11 个全绿** —— 证明标量域此前**零覆盖** |
    | ② | 编译型：标量版体内写 `A.OnlyAtInstantiation()` | 只 include 的探针**编过**；`compute_graph` 照常 Built；**只有 `apply` 编不过**，报 `Matrix<Complex<double>>` 无此成员 ⇒ 标量版的 double 与 Complex **两个域都被实例化** |
    | ③ | 偏序：`#if 0` 掉 Dual 专用版 | Dual 域**确实落到标量版**，被守卫拦住并给出 `static assertion failed … : Apply: Dual 域应由 Dual 专用重载接管`（`Dual<double>` + `Dual<Complex<double>>` 各一次）⇒ 守卫**必要**且**有效** |
    | ④ | 枚举扩展：临时给 `FunctionType` 加 `Tanh` | `elemFunc`（`jacobian.hpp`）+ apply 的 4 个 switch **全报 `-Wswitch`**，而 `dual.hpp` 的 `dualFunc` **命中 0** ⇒ 坐实它的 `default: return x;` 会**静默吞掉**新枚举 |
    恢复一律 `cp` + `shasum` 比对（`c614c460…` 一致）、`grep -r 注入 src/` **0 残留**、复跑全量。
  - **验证**：`./run_tests.sh` **12/12**；`-c -w`（从零 + `-Wall -Wextra`）→ **告警 0、错误 0、12/12**。
    （中途出过一条 `-Wunused-function`：`apply.cpp` 里写了 `FName` 却没用上 —— 已改为在 [B] 段
    逐函数打印最差误差，顺带让输出更有信息量。）
  - **未动**：`GradNode<T>` 仍只支持 `T = Dual<U>`；`GradNode::Backward` 的 `-Wswitch`
    （5 个枚举未处理，含 `Mult` **静默 no-op**）—— 都属 `grad_graph.hpp`，另议。
  - **文档同步**：本文件 §0 头部与工程行、§3 段（改为"已落地"）、§4 待办条（→ `[x]`）、§6 本条；
    `graphics/plan.md` 待办条。**未提交**（`HEAD` 仍 `6999ca2`）。

### 2026-10-15

- **`learning-todo/lie-plan.md` 落地（文档层，无代码改动）** —— 新建「用 `expm` 学李群李代数」
  的八刀路线。用户诉求：「我打算用这一个 `expm.hpp` 学完李群李代数」（并附一份自写的
  《李群/李代数/矩阵指数》总结），随后令「落成 plan」。
  - **定位**：`fan-plan.md` **A4 切片（expm + so(3)）的展开件**。扇面层说「A4 是什么」，
    lie-plan 说「怎么一刀刀走完」。**刀 1 = S5**（主线停点，V1.6 阈值仍在 `graphics/plan.md`），
    刀 2–8 是扇面加餐，**不新增 V 节点、不新建库文件**（`lie.hpp` 标记「跑通后再决定」）。
  - **两处必修表述（照原文写会写出错的用例）**：
    ① **极分解那行不成立** —— 「`X = S+A` ⟹ `exp(X) = R·P`」是错的；`exp(对称)` 是正定对称、
    `exp(反对称)` 是旋转（各自都真），但 `G = exp(X)` 的极分解因子**一般** ≠ `exp(A)`/`exp(S)`，
    除非 `[S,A] = 0`（差的正是 BCH 首项那一阶）。
    ② **「exp 是局部一一对应」要补全局性** —— `d·exp₀ = id` 只保证局部微分同胚；
    `so(3)→SO(3)` 满射但非单射（θ 与 θ+2π 同像），一般李群上 exp 连满射都不是（`SL(2,ℝ)` 有反例）。
  - **能/不能数值落地的分界表**（本计划最有价值的一段）：`exp` 级数、结构常数、单参数子群、
    BCH 首项（可测**阶**）、exp 非单射、**SU(2) 双覆盖**、谱路线、Jacobi 公式 → 能做；
    流形结构 / 表示论 / 根格 / `π₁` 的**证明** → 回 `li1–li3.tex` + `Q_topology.tex`。
  - **八刀**：① expm 本体（= S5）② so(3) 三基 + 对易子 ③ `exp(θ[u]×)` vs Rodrigues（**三线对拍**：
    级数 / 谱 / 闭式）④ θ 扫 [0,4π] 证非单射 ⑤ **SU(2) 双覆盖**（θ=2π 处回 `−I`，SO(3) 回 `I`；
    `π₁(SO(3))=ℤ₂` 的数值影子）⑥ BCH 首项（**测阶**：缩半降 8 倍，不是设绝对阈值）⑦ 极分解 / 谱路线
    对拍 ⑧ `exp(−iHt)` 时间演化（复数域）。
  - **两个超强判据**：`det(exp X) = exp(tr X)`（Jacobi，任意方阵成立，so(3) 上自动排除反射，
    比「与 Rodrigues 对拍」更不依赖专门知识）+ 三线对拍。
  - **收束判据（三样各用一次）**：Jacobi 公式 ≥2 个域、三线对拍在 so(3) 闭合、SU(2) 双覆盖跑出来。
  - **指针已接**：`fan-plan.md` §3 A4 行后加「↳ 展开件」行、§8 加一条、§9 开局三刀第 2 条补注；
    `graphics/plan.md` S5 条目补注（只保留刀 1 的 V1.6 阈值）；本文件 §0 表加「李群线」一行。
  - **代码未动**（写这条时的状态）：`expm.hpp` 仍是 0 字节，**按用户指令由用户自己写**（S5）；
    工作区另有两项 10-14 遗留未提交：`graphics/src/linalg/qr.hpp` 的 `<q_k,q_k>≡1` 注释、
    `graphics/dev-progress.md` 的 S6b 哈希回填；`householder.hpp` 在工作区显示为 `D`（已删，非本会话所为）。
    → **同日稍后 `expm.hpp` 已落地**（Taylor + scaling-and-squaring，54 行），见本日第 3 条。
- **`demo/apps` 附加区 + `run_apps.sh` 落地（工程结构改动）**：用户问「demo 是不是该分成普通 demo 和
  testcase」→ 结论「该分，但分的是**有没有判据**」。落地方式：**不搬已有 8 个 ctest 用例**，只加附加区 ——
  新建 `demo/apps/`（每个 app 只 `add_executable(... EXCLUDE_FROM_ALL ...)`、**不注册 ctest**）、
  新建 `run_apps.sh`（扫 `demo/apps/*.cpp` 列清单 / 单独构建 target / 运行 / `--` 后参数透传；
  `-c -w -j -b -h` 与 `run_tests.sh` 对齐），`demo/CMakeLists.txt` 加 `add_subdirectory(apps)`。
  判定规则（写进注释）：**要不要出现在 ctest 的绿灯里** → 要 `demo/<模块>/`，不要 `demo/apps/`。
  同步三处文档：`demo/CMakeLists.txt:1` 口径、`AGENTS.md` §八目录用途（原文把 `test/` 与 ctest 混为一谈，
  实际 `graphics/test/` 是**早期 C 实验与 notebook**）、`run_tests.sh` 头注释。
  第一个成员 `demo/apps/string_modes.cpp`（fan-plan 扇面 **A1** 的 C++ 版）：离散 Laplacian
  `tridiag(-1,2,-1)/h²` 特征值 → `(kπ)²`，N=8→128 比值 0.989887→**0.999951**，第一模态认出 `sin(πx)`。
  **验证 4 条**：`run_tests.sh` 仍 **8/8**、`ctest -N` 总数 **8**（无 app）、`run_apps.sh string_modes` 跑通、
  **注入 `#error` 到 `demo/apps/zz_inject_broken.cpp` 后 8/8 仍绿，而单独 `--target` 构建该文件确实报错**
  （`EXCLUDE_FROM_ALL` 两面都成立），注入已清理、`grep zz_inject` 无残留。**未提交**。

- **S5 刀 1：`expm.hpp` 落地**（工程线，代码改动；`expm.hpp` **54 行**；**未提交**）。
  用户先自写了第一版骨架（41 行），本刀按 6 处编译问题 + 2 个边界修干净 ——
  **算法逻辑一行未动**（scaling + Taylor + squaring 的骨架本来就对），改的全是「名字/API 记错」
  与「数学等价的改写」：`from_real`→`fromReal`、`ldexp(1,-s)`→`ldexp(real_t(1),-s)`、
  `DefaultTol(A.rows())`→`DefaultTol<T>(n)`、`term*As/fromReal(k)`→`term*As*fromReal(1/k)`
  （**`Matrix / T` 这个运算符本库不存在** —— `matrix.hpp` 只有 `operator*(const T&)`）、
  `Abs2(sum)`→`Abs2(result)`（`sum` 根本没声明）。**真正改行为的只有两处**：
  `while (n2 >= 0.25)` → `>`（只影响 2 的幂边界上多做/少做一次缩放+平方，**数学结果不变**）、
  **补空矩阵 throw**（原版 `Matrix(0,0)` 因 `Rows()==Cols()` 会静默返回空矩阵）。
  - **三域实例化**：`-Wall -Wextra` 0 告警，double / `Complex<double>` / `Dual<double>` 一次过。
  - **基础段（探针 `/tmp/probe_expm_check/probe.cpp`，全 PASS）**：`expm(0)=I` **精确 0**；
    `expm(diag)` 对角 6.4e-14 / **非对角精确 0**（级数自带逐元退化，**不需要为对角写特化**）；
    幂零 `A²=0 → I+A` **精确 0**；`det(exp X)=exp(tr X)` rel 1.1e-15。
  - **用例构造教训（我第一版 A3 假 FAIL，残差 7.0）**：**「严格上三角」≠ `A²=0`** ——
    `[[0,a,b],[0,0,c],[0,0,0]]` 的 `A²(0,2) = a·c ≠ 0`，幂零指数是 **3**。
    判「逐位精确」的用例必须先确认输入在浮点下真满足代数前提（取 2×2，或只留单个超对角元），
    否则是**判据**的错、不是实现的错。
  - **注入①（s 恒 0 = 不做 scaling）对照**（6×6；实对称对谱路线参考 / 随机非正规阵自洽残差）：

    | ‖A‖_F | 正确版·对称 | 注入版·对称 | 正确版·非正规自洽 | 注入版 |
    |---|---|---|---|---|
    | 1 | 1.6e-15 | 9.3e-16 | 6.7e-16 | 4.4e-16 |
    | 5 | 7.2e-15 | 1.0e-15 | 2.8e-14 | 2.0e-13 |
    | 10 | 1.4e-14 | **4.7e-9** | 1.5e-11 | 5.2e-9 |
    | 20 | 2.8e-14 | **2.4e-3** | 1.8e-7 | **1.0** |
    | 40 | **5.6e-14** | **7.8e-1** | **1.0** | 1.0 |

    **两件事必须分开读**：① **算法侧没问题** —— 正确版在对称（良条件）下到 `‖A‖=40` 仍 5.6e-14；
    注入版 ‖A‖=10 就烂到 4.7e-9 ⇒ **scaling 就是 ‖A‖≳5 能算对的原因**。
    ② **判据侧病态** —— 非正规那两列是 `exp(A)exp(−A)` 自身的条件数 e^{2L} ＋ 平方放大
    （非正规阵每次平方按 `2^s` 放大舍入），`‖A‖=40` 时**正确版也读到 1.0**，**不是算法错**。
    ⇒ 重申：expm 判据一律归一化（`‖Xs−Xp‖∞/‖Xs‖∞`）或改用标量恒等式，
    **大范数下别拿自洽式当唯一判据**。
  - **验证**：`./run_tests.sh` 仍 **8/8**（`expm.hpp` 目前**无任何 TU include** ⇒ 不进 ctest，符合预期）；
    注入用 `cp` 备份 + `shasum` 比对还原（`c5d5576a…` 前后一致），`grep 注入` 0 残留。
  - **未做**：`demo/linalg/expm.cpp` 与 ctest 注册 —— 即 §4 的第 ② 步（V1.6），按用户节奏留待下一刀。

- **文档基线同步（本会话收尾）**：`dev-progress.md` §0 指针（工程行 + 李群线行）、§3 现状明细、
  §4 下一步全部按 2026-10-15 实际状态覆盖重写 —— 原文还停在「linalg 528 行从未被任何 TU include」
  「ctest 2/2」「V1.1–V1.7 一条都没通过」，已与实况（linalg 1062 行全部落地、ctest 8/8、
  V1.2/V1.3/V1.4 通过）严重脱节。被覆盖的 2026-10-13/14 明细在 §6 同日条中已完整保留，
  未丢历史。

### 2026-10-14

- **V1 收尾体检（agent 会话，纯探针，无代码改动）**：用户问「v1.0 可以收尾了吗」→ **答案：不能**。
  数据已同步到 `graphics/plan.md` 同日条与本文 §0/§3/§4。
  - **根因：`linalg/` 528 行从未被编译器看过** —— `demo/` 只有 `scalar_traits`、`compute_graph`
    两个 ctest，没有任何 TU include 过 `linalg/*.hpp`。所以先前记的"已落地"实际含义是"文件存在"。
  - 两个 `-fsyntax-only` 实例化探针共 **32 个编译错误**（`solve`+`inverse`+`eigen` 一轮 16 个；
    `cholesky`+`qr`+`operator` 一轮 16 个）。**`cholesky.hpp` 是唯一能通过实例化的**。
    空文件：`svd.hpp`（V1.4 必需）、`expm.hpp`（V1.6 必需）、`householder.hpp`；
    `dot.hpp` 已删（`grep` 确认无引用残留）。
  - **编译错误之外还有算法错误**（编过也会算错）：`solve.hpp` 消元写成 `U(i,k) -= L(i,j)*U(i,k)`
    （应为 `U(j,k)`），且整段被 `if (max_row != i)` 包住 —— 早先记的"`if (max_row != i)` 可疑"
    被实测证实；主元比较写裸 `>`（违反 D7，Complex/Dual 没有 `operator>`）；
    `eigen.hpp` 是 2x2 单次旋转，**不是 Jacobi 迭代**（无收敛循环、无非对角元判据）。
  - V1.1–V1.7 逐条全未通过：V1.4 缺 `svd.hpp`、V1.6 缺 `expm`、V1.1/V1.2/V1.5 卡在 linalg 编不过、
    V1.3 卡在 eigen 非迭代且编不过、V1.7 无 demo；plan 明确要求的 `demo/common.hpp` 不存在。
    8 步收尾最小路径见 §4。
- **新增 `graphics/run_tests.sh`（10-13 傍晚会话，已提交）**：配置 + 构建 + ctest 一条链，
  选项 `-c` 从零 / `-d` 明细 / `-w` 告警 / `-b` 指定目录 / `-j` 并行度 / `-l` 列表 /
  `--build-only` 只构建 / `--test-only` 只测试 / `-h`；退出码 0 全绿 / 1 构建失败 / 2 用例失败 / 64 参数错。
  用例清单由 `ctest -N -V` 自省（新用例自动带出），cwd 无关，非 tty 关色。
  已登记进 `AGENTS.md` 第三节「命令」。
- **提交（10-13 傍晚，均未 push）**：父仓 `graphics_0` —— `375f618`（基线修复，12 文件）、
  `c82c7ed`（`run_tests.sh` + AGENTS.md 登记，含先前已写未提交的 learning-408 条款）、
  `1f33b7a`（脚本补 `--build-only`/`--test-only`）；子仓 `learning-todo` —— `e733363`（`matrix.md`）。
- **S1 `solve.hpp` 数值层修复完成（代码改动；`ctest` 3/3，`solve` 用例 34/34）**。
  节奏按用户定的「一步一条 → 只解释清楚 → 再动手 → 单独编译/跑用例」：
  - **③ 容差**：删全局 `Tol`（**float** 的 eps 1.19e-7，用于 `Matrix<double>` 差 9 个数量级）
    → `DefaultTol<T>(n) = 8n·eps`，判据 `abs2(p) <= tol*tol`。老写法还藏了第二条缺陷：
    左边是平方量、右边不平方 → 等效阈值 √Tol = 3.45e-4（探针 `/tmp/probe_tol.cpp` 实测：
    条件数 1e12 的良定 2×2 被老判据误判奇异）。依据 D12 + §10.8 ① + §2.4。
  - **④ 边界检查**：`LU` 非方阵 / `Solve` 的 B 行数不符 → `throw std::invalid_argument`（D12 四分类）。
    **没给 `matrix.hpp` 的 `operator()` 加检查** —— D12：越界属编程 bug → assert，`at()` 才是检查版，
    维度校验放 API 边界即可。此前 `LU(Matrix(2,3))` 越界读、`Solve(I₃,B₂ₓ₁)` 越界写 = 静默 UB（用例 SIGABRT）。
  - **① 删 `if (max_row != i)`**（**用户自己判断出来的**）：`max_row` 是换行**前**的行号，换行后主元行
    已坐在 `j` 上，"不消主元行"由 `i > j` 天然保证；该判定实际**跳错了行**（跳掉的是被换进来的旧第 j 行，
    那行此刻正带着非零 `U(i,j)`）。只在 `max_row > j` 且发生换行时触发 → **n=2 永不触发，n≥3 才炸**。
  - **② 乘子随换行**：新增 `SwapMultipliers(L, r1, r2, upto)`，换 U 行时同步换 L 的 `[0, j)` 乘子列；
    **不能整行换**（单位下三角的对角 1 与行号绑定，整行换会把 1 搬到上三角）。
  - 修复路径：OK **3 → 7**（③④）→ **19**（①）→ **34**（②）。最差相对误差
    `P·A−L·U` 1.187e0 → **1.135e-16**，`A·X−B` 2.996e0 → 1.317e-15，复数 3.926e0 → **6.218e-16**。
  - **注入反例实测**：去掉 `SwapMultipliers` 调用 → **15 条 FAIL**（恰是 ② 的下游全部）。
  - 期间的替代思路讨论（用户提「先把 A 按主元编号排好再减」）：方向对，但主元是**自适应**选的
    （第 j 列主元取决于前 j−1 列消元后的数值），排列无法一次排完，只能与消元交替；
    该思路的干净实现是**原地 LU**（L/U 同一矩阵、单位对角隐式、整行换），探针 `/tmp/probe_inplace_lu.cpp`
    实测 `max|L·U − P·A|` ≤ 4.4e-16（手算 3×3 + 随机 4–8 阶），得到的 L 与方案 (a) 完全一致。
    **用户最终选 (a)**（`SwapMultipliers`，不动 `LUResult` 接口）。
- **S2 `inverse.hpp` 落地完成（代码改动；`ctest` 4/4，`inverse` 用例 48/48）—— V1.2 硬验收通过**。
  按用户定的节奏「先搞 `inverse.hpp`，QR/cholesky 那些慢慢引入」，本刀自包含（不碰 cholesky/qr）：
  - **编译层 5 处**：打开 `#include "linalg/solve.hpp"`；`scalarTraits` → `ScalarTraits`；
    `Identity<T>(n)` → `Matrix<T>::Identity(n)`（它是静态成员，不是自由函数）；
    `diagonal_inverse` 补 `tol` 形参；`Cholesky(A)` 调用移除。
  - **两条"编过会算错"的探测 bug**（比编译错误更该先看）：
    ① `is_symmetric` 把比较对象写成 `AH(j,i)` = conj(A(i,j))，实际判的是"A(i,j) 的虚部是否为 0"
    → **实数矩阵恒真，所有实方阵都被判成 Symmetric**（正确判据是比 `AH(i,j)`）；
    ② `is_orthogonal` 只查 `A^H·A` 的对角线 → `A^T·A` 对角全 1、非对角非零的矩阵被误判成正交。
    两处统一改成全元素 `abs2 <= tol*tol` 口径（与 `solve.hpp` 一致）。
  - **分发**：`Inverse(A)` = `Classify(A, DefaultTol<T>(n))` → `Inverse(A, hint)`；
    `hint` 分支 Diagonal→`diagonal_inverse`、Orthogonal→`A.Adjoint()`、Symmetric/Generic→
    `lu_inverse`（LU）。**SPD / QR 分支显式 `throw logic_error`**（未接入 cholesky/qr），
    `is_positive_definite` / `cholesky_inverse` / `qr_inverse` 三处同样 throw —— 不留静默 no-op；
    `Classify` 的 SPD 档暂时关闭，对称矩阵落到 Generic 走 LU（结果依然正确，只是没用上 Cholesky 优化）。
  - **验收**：新增 ctest `inverse`（`demo/linalg/inverse.cpp`，48 项）—— V1.2 硬指标
    ‖A·A⁻¹ − I‖∞ < 1e-8，double n=2..6×4 seed 与 Complex n=2..4×3 seed **共用同一份源码**；
    实测最差 double **1.18e-13** / complex **9.16e-16**。另含结构分支解析对拍、Classify 语义、
    hint 一致性、异常路径（奇异 / 非方阵 / 0×0 / SPD / QR）。
  - **注入反例实测**：`is_symmetric` 还原成 `AH(j,i)` → 2 条 FAIL；`is_orthogonal` 还原成
    "只查对角" → 1 条 FAIL。闸门有效。
  - **工具坑（记录）**：注入脚本用 `shutil.copy` 恢复文件后，CMake 可能因 **mtime 粒度**判
    "up-to-date" 而不重建，跑出的仍是注入版二进制（我一度误读为"恢复后仍 1 FAIL"）——
    恢复后必须 `touch <源文件>` 或删构建产物再跑。
  - 新增 **S2a**（`cholesky.hpp` 修复 + 打开 SPD 档）已登记进 `graphics/plan.md` 的 S 清单。

- **S2a `cholesky.hpp` 分解侧完成**（工程线，代码改动；`ctest` **5/5**，`cholesky` 用例 **98/98**）。
  用户定的是**方案甲**：删递归版、只留元素公式版，先把 V2.2 拿到手（递归分块版用户后续自行补）。
  - **结构性整理（本刀前置，单独构建验证无回归）**：
    ① `ScalarTraits<T>::real(const T&) -> real_t` 新增：float/double 恒等、`Complex<T>` 递归取
    `.real`、`Dual<T>` 只看 `value`（与 `abs2` 同口径），concept 同步加一条契约 —— 于是
    `Complex<Dual<double>>` / `Dual<Complex<double>>` 这类嵌套也成立（探针 6/6 实测）。
    ② `DefaultTol<T>(n)` 从 `linalg/solve.hpp` 上移到 `core/scalar.hpp` 标量层：它是 D12 的容差策略，
    solve / inverse / cholesky / eigen / svd 共用，不该让 cholesky 为了拿一个容差去 include LU。
  - **修掉的问题**（编译层 6 类 + 数学层 2 处）：
    `-=` / `/=`（Complex/Dual 无复合赋值，D7）；**两处 `L(·,·)*L(·,·)` 漏 `conj`**（实数下免费、
    复数下必错）；`L(i,i) <= 0`（复数无 `<=`，改走新加 `real`，且写成 `!(real(d) > tol)` 让 NaN 也落 throw）；
    `sqrt` 未限定命名空间（两步法 —— Complex/Dual 的 sqrt 是类内 friend，ADL 可见）；
    返回 `Matrix<T>()` 当错误码 → 改 throw（D12）；删半成品 `cholesky_r` / `cholesky_recursive`
    （后者自己写着 TODO「只处理了对角线」，且被定义在它之后）。
  - **新增** `cholesky_solve`（前代 + 回代；V2.2 要求「解同一方程组」）。
  - **验收**：新增 ctest `cholesky`（`demo/linalg/cholesky.cpp`，98 项）。实测最差重构
    **double 1.78e-15 / complex 1.78e-15**，解之差 **double 5.55e-17 / complex 2.78e-17**；
    double n=2..6×4 seed、Complex n=2..4×3 seed **共用同一份源码**。含手算 2×2 / 3×3 / Hermitian 2×2、
    结构判据（严格上三角精确 0；对角元为正实数）、秩亏边界（A = M·Mᵀ，M 3×2 —— 最后一个主元才塌）、
    接口契约（只读下三角）、SPD 判据语义。
  - **注入反例**：① 两处漏 `conj` → **23 条 FAIL，全部落在 complex 组；double 重构仍是 1.78e-15
    毫发无损** —— 直接印证「实数测试是这个 bug 的盲区」；② 去掉正定性判据 → 7 条 FAIL 全在
    异常 / SPD 判据上；③ 消元那步除法写成乘法 → 71 条 FAIL（并在 E 段之后 abort，因为 E 段没包 try）。
  - **两个新踩的坑（都记下）**：
    ① **clang 默认 `-ffp-contract=fast` 会把 `a - b*c` 融成 FMA**，于是「虚部」判据若写成
    `abs2(d) - re*re`，算出来的是 `fl(d²) − d²`（即 d² 的舍入残渣 ~1e-16），在 double 上恒大于 `tol²`
    —— 实测加 `-ffp-contract=off` 后残渣精确为 0，坐实成因。正确写法是 `abs2(d - fromReal(re))`
    （减法精确抵消，不受 contract 影响；复数下正好是 `Im(d)²`）。**这是判据写法的坑，不是实现的错**。
    ② 判据函数不包 try 的话，实现一旦抛异常会让**整个用例 abort**，其余判据一条都跑不到 ——
    已给 `CheckReconstruct` / `CheckAgainstLU` 加 try，并把异常记成一条 FAIL。
  - **未动**：`inverse.hpp` —— SPD 档的打开、`cholesky_inverse` / `is_positive_definite` 的接入
    留作 **S2a-2**（当前三处仍显式 throw），按「一刀一文件」纪律分开做。
  - **提交**：S2a 由用户放行后单独一刀 `daefb3e`（2026-10-14，未 push）。

- **S2a-2 `inverse.hpp` SPD 档完成**（工程线，代码改动；`ctest` **5/5**，`inverse` **48 → 78 项**全过）。
  用户指令「可以开启 spd 吧」；只在 `inverse.hpp` + `demo/linalg/inverse.cpp` 动刀，`cholesky.hpp` 未碰。
  - **三处改动**：
    ① `is_positive_definite(A)` 由 `throw logic_error` 改为"cholesky 走通即真"（matrix.md §3.2 的
    规格就这一句）。**只 catch `std::runtime_error`**（cholesky 判出的非正定/近奇异），别的异常照常
    外传不吞；非方阵/空矩阵**提前返回 false**（不问 cholesky，免得把入参错解读成"非正定"）。
    ② `cholesky_inverse(A)`：`A⁻¹ = (L⁻¹)ᴴ·L⁻¹`，实现为 `cholesky_solve(A, Identity(n))` ——
    展开正是"L·Y = I 前代 + Lᴴ·X = Y 回代"，同一条路径，直接复用那份压过 98 项验收的代码。
    ③ `Classify`：SPD 档**嵌进 `is_symmetric` 分支内部**，不是平铺的同级 if。
  - **本刀的核心坑（单独记）**：`inverse.hpp` 原注释写 `Diagonal -> Orthogonal -> Symmetric ->
    (SPD) -> Generic`，按字面**平铺**成同级 if，SPD 档永远不可达 —— 因为 **SPD ⊂ Symmetric**，
    上一行 `is_symmetric` 先 return 就把它吃掉了。正确写法是在对称分支内部筛子类。
    附带结论：`Diagonal` 排在 SPD 之前是有意的（单位阵落 Diagonal，逐元倒数比 Cholesky 更省）。
  - **取舍（已认可）**：改完后每个对称矩阵都会多跑一次 Cholesky（非正定则 throw 后 catch，当判据用），
    白付一次 O(n³/3)。可接受 —— `Classify` 本质就是"付一次 O(n³) 探测换后续用上最优算法"，
    且 `is_orthogonal` 那步本来就要算 `Aᴴ·A`，总体量级不变。
  - **验收**：`inverse` 用例 **48 → 78 项**（新增 G 段 30 项）。V1.2 硬指标最差
    **double 1.18e-13 / complex 9.16e-16** 与 S2 持平；随机 SPD（n=2..6，S = M·Mᵀ + n·I）走 cholesky
    与走 LU 的逐元差 **≤ 1.11e-16**、残差 ≤ 4.44e-16；`is_positive_definite` 真/假两侧都测
    （正定 / 不定 / 全负 / 半正定奇异 / 单位阵 / 非方阵）；复 Hermitian 正定走同一份源码。
  - **注入反例**：① Classify 的 SPD 档平铺成同级 if（制造顺序陷阱）→ **8 条 FAIL**，全部是 SPD 判据；
    ② 正定判据恒 true（吞掉 cholesky 的否定结论）→ **5 条 FAIL**，全部落在"非正定"一侧。
    两次注入均恢复后复跑 5/5，`grep 注入` 无残留。
  - **提交**：S2a-2 单独一刀 `f73790b`（2026-10-14，未 push）。

- **S3 `eigen.hpp` 真 Jacobi 迭代完成**（工程线，代码改动；`ctest` **6/6**，`eigen` **119 项**全过；**V1.3 通过**）。
  用户开工前自己把理解纠正过来了（原以为"把非对角元改成 0 就完事"，意识到"还是得套个 while 循环"）——
  这个判断是对的，下面第①条注入反例正是它的硬证据。
  - **原实现三层问题**（40 行的 `eig_symmetric`）：① **只扫一遍**（无收敛循环、无非对角判据）；
    ② `rotate()` 做整体 `R·A·R^T`（每次 O(n³)，且用 `Transpose` 而非 `Adjoint`，复数域方向就错）；
    ③ `jacobi()` 硬编码 `double`，`ScalarTraits()::one()` 语法错，`Identity(Rows, Cols)` 多传一参，
    `eig(A)` 调 `is_symmetric(A)` 缺 tol、且调了根本不存在的 `eig_general(A)`。
  - **重写结构**（cyclic Jacobi）：
    · 外层 `for sweep` 一轮轮扫，内层遍历全部 i<j；收敛判据 `off(A)² = Σ_{i≠j}|A(i,j)|² <= tol²`。
    · **容差标尺只算一次**：酉相似变换保持 Frobenius 范数，故 `tol = DefaultTol<T>(n) * ‖A‖_F`（D12）。
    · **单元素跳过阈值** `tol²/nd`（nd = 非对角元个数）—— 即使剩余元素全被跳过，总和也不超过 tol²。
    · 每步**只更新第 p、q 两行两列**（O(n)），不再整体旋转。
    · **Hermitian 的相位**：酉阵 `U = [[c, -s·ph], [s·conj(ph), c]]`，`ph = z/|z| = e^{iφ}`；
      实对称下 `z` 为实、`ph = ±1`，同一份代码自动退化 —— 这就是 V1.1「同一份源码」的落点。
    · **对角元走解析式** `(a+b)/2 ± hypot(a−b, 2|z|)/2`，与 `θ = 0.5·atan2(2|z|, a−b)` **必须配套**：
      取 θ+π/2 同样能消零，但会把两个特征值对调（这一条在稿纸上先算错过一次，靠 2×2 数值验证找回）。
    · 100 轮仍不收敛 → throw（D12：不返不准的结果）。
  - **踩到的坑（S2a 教训复用）**：第一次注入（只扫一轮）时用例直接 **SIGABRT** —— `eigen.cpp` 主函数
    没包 try，实现抛出的 `runtime_error` 把整个进程 terminate 掉，一条判据都看不到。
    已加 `RunEig()` helper 把异常记成一条 FAIL（B/C 批量段），之后注入才表现为 16 条 FAIL 计数。
    **又一次印证：判据函数与批量循环必须包 try。**
  - **验收**：`eigen` **119 项**全过。V1.3 最差重构 **double 1.78e-15 / complex 1.33e-15**（判据 1e-8）；
    Q 正交/酉 ≤ 1.4e-15；`‖QᴴAQ − Λ‖∞` ≤ 1.1e-15；逐列特征对 `max‖A·q_j − λ_j·q_j‖∞` ≤ 1.6e-15；
    迹不变 `|Σλ − tr(A)|` ≤ 2.2e-16（与特征值取值无关，是独立校验）。
    含 2×2 闭式对拍（(5±√5)/2、2±√2）、已对角/单位阵/零矩阵/1×1、非对称与非方阵 throw。
  - **注入反例**：
    ① **只扫一轮** → **16 条 FAIL**，分布极有说明力：**n=2 全过、n ≥ 3 全中** ——
    2×2 只有一对非对角元，一次旋转就彻底清零；从 3×3 起，旋转 (0,1) 会在 (0,2)/(1,2) 上
    **重新灌入非零**（fill-in），所以必须一轮轮扫。这直接回答了"为什么会完"的疑问。
    ② **漏 `conj`**（`s·conj(ph)` → `s·ph`）→ **24 条 FAIL 全落在 complex、实对称 0 条**，
    且连 `‖QᴴQ − I‖∞` 都坏到 7.4e-3（漏 conj 破坏了 U 的酉性）——
    与 S2a 那次漏 conj 注入（23 条全在 complex）完全同构。
    两次注入恢复后复跑 5/5 → 6/6，`grep 注入` 无残留。
  - **提交**：S3 单独一刀 `5751b05`（2026-10-14，未 push）。

- **S4 `svd.hpp` 组装完成**（工程线，代码改动；`ctest` **7/7**，`svd` **266 项**全过；**V1.4 通过**）。
  原文件是 **19 行空壳**（只有 `SVDcomposition` 结构体声明，无任何实现、无引用点）。
  - **路线**（`plan.md` §3.3 指定）：经 `A^H·A` 的 Jacobi 特征分解组装 —— `B = A^H·A`（n×n Hermitian 半正定）
    走刚验完的 `eig_symmetric`，`σ_i = √λ_i`（非负根 + 降序），`u_i = A·v_i/σ_i`。
    经济型：U 为 m×r、Σ 为 r×r、V 为 n×r，`r = min(m,n)`。
  - **为什么 σ 的精度不拖累重构**（已写进头文件注释）：`U·Σ·V^H = Σ_i (A·v_i/σ_i)·σ_i·v_i^H = A·(V·V^H)`，
    σ 与 1/σ **精确相消** —— 重构只取决于 V 的列正交归一、且张成 A 的行空间。
    这条不只是理论好看：第③条注入反例就是它的数值证据（σ 算错 102 条挂，重构一条不挂）。
  - **σ ≈ 0 的列**（秩亏时）：`u = A·v/σ` 数值不可用，而零奇异值对应的左奇异向量本就在
    正交补里任取 —— 用 Gram-Schmidt 扫标准基补一个确定性结果（`linalg_detail::FillOrthonormalColumn`）。
  - **判据分层（U 判 1e-10、V 判 1e-12）**：U 的列方向精度是 `eps·κ(A)` —— σ 走 `A^H·A` 只有 √eps 级，
    除以 σ_i 后放大，实测最差 **2.456e-13**；V 直接出自 `eig_symmetric`，实测 ~1e-16。
    两者都判 1e-12 会留下一个只有 4 倍余量的脆阈值，故按各自路线的固有精度分层。
  - **验收**：`svd` **266 项**全过。V1.4 最差重构 **double 1.332e-15 / complex 7.772e-16**（判据 1e-8）；
    `‖UᴴU − I‖∞` ≤ 2.456e-13、`‖VᴴV − I‖∞` ≤ 5e-16；`Σσ² = ‖A‖_F²` 相对差 ≤ 2.9e-16；
    `σ_i²` 与 `eig(AᴴA)` 谱的逐项相对差 ≤ 1.9e-16。
    覆盖：实/复 × 方阵/长方阵（3×5、5×3、1×4、4×1、2×7、6×2）、秩亏（全 1 阵、列 2 = 2×列 1）、
    零矩阵、1×1、1×2，以及解析对拍（σ = 3/1、1/1、4/1、5、重排 5/3/1 与 4/2/1）。
  - **注入反例（三条，均恢复后复跑 7/7，`grep 注入` 无残留）**：
    ① **`Adjoint` → `Transpose`**（`B = Aᵀ·A`）→ **34 条 FAIL 全在 complex、实数域 0 条**
    （实矩阵下 `Transpose ≡ Adjoint`，天然免疫）—— 与 S2a/S3 的漏 conj 注入**第三次同构**。
    ② **不排序**（`ord` 保持恒等）→ **5 条 FAIL，且「重构」判据一条不挂**（σ 与配对都没错，只是顺序不同）。
    ⚠️ **这条暴露了判据覆盖漏洞**：原有用例的对角矩阵 `diag(3,−2)` 本就降序，`eig_symmetric` 对随机矩阵
    的输出顺序又常恰好降序 —— 于是"不排序"几乎抓不到。已补「**升序**输入」用例
    （`diag(1,5,3)` 实、`diag(1,4,2)` 复，要求输出重排为降序），补后该注入 FAIL 数 3 → 5，两条新用例全中。
    ③ **σ 不开方**（`σ_i = λ_i`）→ **102 条 FAIL，但「重构」判据 0 条挂**。原因正是上面那条相消性质：
    u 用错 σ 算、再乘错误的 σ，`U·Σ·V^H` 依旧精确等于 A。**结论：V1.4 的字面判据「重构 < 1e-8」
    单独抓不到 σ 算错**，必须靠 `Σσ² = ‖A‖_F²`、`σ_i² = λ_i(AᴴA)`、解析对拍这些独立判据兜底 ——
    这是本刀最有价值的一条反例。
  - **提交**：S4 单独一刀 `edfc572`（2026-10-14，未 push）。

- **S6b `qr.hpp` 三函数修复**（用户要求提前做 qr、跳过 S5 顺序；`operator.hpp` + `qr.hpp` +
  新增 `demo/linalg/qr.cpp` 189 项；ctest **8/8**）：
  - `operator.hpp`：`Dot`/`Abs2` 是**成员函数残骸**（自由函数里用 `this->rows/cols/data`，
    `Abs2` 还写了大写 `matrix.Data`）。改为两参 `Dot(a,b) = Σ conj(aᵢ)·bᵢ`（**补上共轭**，
    原版没有 —— 第 4 个「double 免疫」点）、`Abs2` 返回 `real_t`。
  - `qr.hpp` 编译层 9 处（`real_t` 裸用 / `aplha`·`Cows`·`QRComposition` 拼写 / `One`·`Zero` 大小写 /
    `Identity` 多传一参 / `atan2` 单参 / `-=` / 缺 `return`）。
  - **算法层三处（编译过也会错，这才是值钱的）**：
    ① `QR_Householder` 分母 `Abs2(X)` 算的是**改前**的 ‖x‖²，而 `H = I − 2vvᴴ/(vᴴv)` 的分母必须是
    **改后**的 vᴴv（= 2‖x‖(‖x‖+Re x₁)）—— 量用错，H 不再对合、消不掉 x；`+sqrt(alpha)` 恒取正号、
    无符号选择（x₁<0 时灾难性相消）；`H.Transpose()` 复数域应为 `Adjoint`。
    ② `QR_givens` 内层循环 `for (j=i+1; j<i; …)` 条件永不成立（死码），用户改成 `j < A.Rows()` 后
    **立刻堆越界**（ASan: heap-buffer-overflow @ qr.hpp）—— 因为 R 只有 n 列。**根因是循环结构本身是镜像**：
    正确结构是外层列 j、内层行 i>j（消下三角），原版是外层行 i、内层列 j>i（去消上三角）。
    ③ `QR_Shimit` 把投影系数塞进 `Q(k,j)`（Q 应是正交基矩阵）、拿 `res.R` 的第 k 列当"已正交基"
    （该取 Q 的列）、无归一化、`res.Q = Q * res.Q` 语义完全不对 —— 与 Gram-Schmidt 无关。
  - **修复**：三个函数重写。householder 用 `v = x + phase(x₁)·‖x‖·e₁`（分母取改后 vᴴv、`Adjoint`）；
    givens 改成外层列内层行、复相位 `s = (f/|f|)·conj(g)/r`；shimit 重写为**修正 Gram-Schmidt (MGS)**
    （逐次用当前 v_j 算系数、列归一化落 R 对角、秩亏走正交补补齐）。
    **统一 thin QR 口径**：Q: m×n 列正交，R: n×n 上三角，`A = Q·R`；`m < n` 三个入口统一 throw。
  - **验收**：`qr` 189 项 / FAIL 0；V2.1 最差重构 **9.93e-16**、Q 列正交 **8.88e-16**（判据 1e-8 / 1e-12）；
    覆盖实/复 × 方阵/长方阵/列向量 + 解析对拍（对角/单位/2×2 闭式）+ 秩亏/零矩阵/1×1 + 异常（空/`m<n`）。
  - **注入反例 4 条**：① 分母退回改前 ‖x‖² → **35 条 FAIL 全在 householder**（givens/mgs 0 条）；
    ② givens `Adjoint`→`Transpose` → **3 条全在 complex、实数域 0 条**（第 4 次「double 免疫」）；
    ③ MGS 内积漏 `conj` → **3 条全在 complex**（Q 酉性崩到 4.3e-1）；
    ④ 去掉相位/符号选择 → 5 条，**首轮实数域 0 条 → 暴露判据覆盖漏洞**（实数用例 x₁ 全为正）
    → 补「病态 x₁<0」用例（`[[-1,0],[1e-6,1]]`）后重跑，**实数域 0 → 1 条**（相消残渣 8.89e-11）。
  - **提交**：S6b 单独一刀 `bae45df`（2026-10-14，未 push）。
- **目标口径修正（用户定调）**：用户明确「NS 方程只是打开整个学习路线的钥匙，不代表这些基础
  就不重要了」——即 **NS 是选题与排序的钥匙，基础本身也是目标**。此前多轮咨询（V1 数学基础、
  高数/线代/概率回补、最小二乘/PCA/左逆右逆/牛顿/叉积落点）里 agent 反复强调「窄而深、
  概念在消费点按需到货」，该口径**只对排序有效，不构成降级基础的理由**；推迟 = 排到消费点
  再学，≠ 可以不学。已写进 `AGENTS.md` §二（agent 规约层，防未来 session 再误读）、
  `README.md` 总纲句。数学线的概念正课（多元正态/CLT/特征函数 → C6–C10；左逆右逆的
  无穷维版本 → L2 泛函）不受此影响，仍是必经路径。
  当日提交序（均未 push）：`f73790b`（S2a-2）→ `5751b05`（S3）→ `edfc572`（S4）→ `bae45df`（S6b）。
- **学习模式定调：「先铺开再收束」（用户提出，agent 落地）**。用户想最大程度从应用入手
  （PCA、旋转、李代数→李群……），先把广度铺开获得「每写一个东西就是巨大收获」的体验，
  再回主线逐深；并点名了麦克斯韦链（多元函数 → 复变 → 场）与概率→测度论收束。
  落地形态（与「NS=钥匙、基础=目标」的新口径兼容）：
  - **工程侧 = 应用扇面 S-App1–6**（已写进 `graphics/plan.md` 待办）：PCA/牛顿/伪逆/lstsq/
    旋转扇面/qr_inverse，每件 0.5–1 天，全部是**已验收 linalg 零件上的薄封装 + demo**，
  不新建库文件、不新增 V 节点；主线 S5→V1 收口→阶段 2 顺序不变、与之并行。
  收束条件：V1.5/V1.7 + 扇面提交后清零回主线。
  - **数学侧 = 两条「视角线」走数学线既有材料**（不新建文件）：
  ① 麦克斯韦链 = E1–E2（微分形式）+ Q/R 笔记的读线 —— 高数(多元)→复变(R)→场(E)→麦克斯韦
  是一条真实的物理收束线，E 笔记本来就是为此写的；特征函数/ODE/PDE 的疑惑在
  `learning.md` 12 阶段题库里都有对应手推题，届时按题进入。
  ② 概率→测度论收束 = C6–C10 的概率深化（CLT/特征函数/均匀可积）正是测度论
  （L1 正在学的收敛定理）在概率上的投影 —— 这条收束**恰好就是 L1 的动机**，无需新安排。
  - 风险条款（agent 提醒，用户知情）：扇面每件都「甜」，铺开期的纪律是**扇面件必须见验收
  断言（对拍/打印），不然就退化成抄 demo**；收束点定死在 V1 收口，不许滑动。
  - **升级（同日稍后，合并另一 session 的「四扇面 + 枢纽」设计）**：用户贴来另一份扇面计划
  （四扇面 A/B/C/D + 谱枢纽：本征函数 × 傅里叶 × 李群在 L2 谱定理汇流），确认合并。
  落地：**新建 `learning-todo/fan-plan.md`**（扇面唯一真源：切片清单 + 分层规则 + 收束判据）；
  S-App 清单保留为 C++ 薄封装执行清单（对应关系：S-App1→A2、S-App5→A4 等，见 fan-plan §8）。
  **引入 `learning-viz` Python 层当「实验室」**（用户判断「现在引入 viz 层的 python 也可以了」）：
  扇面切片 Python 优先（matplotlib 分钟级出图，图存 `learning-todo/imgs/`），C++ 只收需要
  V 验收的本体；`viz_modules/` 目录即为此就位。Python 纪律：只 `uv run`、每切片一图、
  有图才算完、数值进主线须走 V 节点不反向。开局三刀：**A1 弦的模态**（Python）→
  **S5 expm**（C++ 主线停点）→ **D2 CLT 动画**（特征函数困惑 + L1 动机）。
  §0 指针表已加「扇面」行。fan-plan.md 在 learning-todo 子仓，**提交须在该子仓内进行**。

### 2026-10-13

- **工程线实地盘点**（agent 会话，无代码改动）：发现实际代码远超 2026-09-29 日志记录——
  C++20、scalar/complex/matrix 重写、dual 修复、LU+Solve 均已提交在 `graphics_0` 分支
  （`a5ff34e`~`2156bb1`，用户自行推进，未走会话）。同时发现三个卡点：
  ① `compute_graph` demo 编译不过（`grad_graph.hpp` 引用未声明的 `LnMultAll`/`Function`）；
  ② `inverse.hpp` 有 `scalarTraits<T>` 拼写错误并引用未定义的 `Cholesky(A)`；
  ③ `Solve` 的消元循环 `if (max_row != i)` 逻辑可疑（与 piv 交换语义冲突，且 `SwapRow(result.U, j, max_row)` 只换 U 不换 L 已换行）。
  V1.1–V1.7 验收全部未跑。§0 指针、§3 现状、§4 工程待办已按实际状态同步。

- **compute_graph 基线落地**（工程线，代码改动；`cmake --build graphics/build` + `ctest` 1/1 通过，
  4 组用例 36 项检查全 OK）。**根因更正**：不是旧记的「`LnMultAll`/`Function` 未声明」——
  这两个名字是 `GradOp` 的合法枚举值；真因是 `enum class` 迁移留下的三类适配问题。
  改动清单：① `grad_graph.hpp` 两处 switch 的 35 个 `case X:` 加 `GradOp::` 前缀
  （其中 `SumAll`/`MultAll` 还会撞 `jacobian.hpp` 同名函数，报错信息极具误导性）；
  ② 类内 CTAD 失败 2 处：`Matrix(1, 1)` → `Matrix<T>(1, 1)`；
  ③ 补两个**缺失分支**：`Output`（透传 —— 缺它导致所有用例输出恒 0）与 `Mult`
  （逐元素乘，与 `CWiseMult` 同义；不补则 `y=w*x` 这类用例算出全零却"编译通过"）；
  ④ 未实现算子（`CWisePow`/`LnMultAll`/`CrossEntropy`/`Sigmoid`/`Relu`/`Tanh`/`Softmax`）
  由静默 `break` 改为 `throw std::logic_error("...: 未实现")`，杜绝静默零值；
  ⑤ demo 有限差分步长 `h` 由 1e-3 提到 1e-2（float 中央差分舍入 ~ eps·|y|/h ≈ 3.7e-3，
  原值在 `∂y/∂x2[0]` 一项假报 FAIL：1.99890 vs 2.0；四个用例对被扰动变量都是线性的，
  截断误差为 0，故加大 h 纯赚）；⑥ demo 头部注释的 `g++ -std=c++17` 单编说明改为 cmake/ctest。
  新增待改项（见 §4）：`scalar.hpp` 的 `ScalarTraits<Dual<T>>::conj/abs2` 成员名写成
  `x.real`/`x.dual`（Dual 实为 `value`/`derive`）；`jacobian.hpp` 的 `Derive` 收不了临时量。
  **本次未提交**（留用户核验后自行 commit/push）。

- **基线"接口层"bug 清完**（工程线，代码改动；`ctest` **2/2** 通过）。把上一条列的"已知边界"
  里属真 bug 的部分修掉，逐条都有实测证据：

  1. **标量契约与设计文档漂移（最重的一条）**：`scalar.hpp` 里是文档中不存在的
     `const_n(内层 T)`，而 `matrix.md` §7.1 第 2 条要的是 `fromReal(real_t)`（且
     `complex.hpp` 的 `sqrt` 已经在调 `fromReal` —— 一实例化就编译错，只是没人调过）。
     后果用 `static_assert` 实测：`Scalar<Complex<double>>`、`Scalar<Dual<float>>`
     **双双不成立**，于是任何 `requires CommutativeScalar<T>` 的 linalg 函数都拒绝
     Complex/Dual —— **V1.1 与 V1.5 直接卡死在这一条上**。已按设计改名、补进 concept
     （顺带堵住 §7.1 点名的"泛型代码裸写字面量"假设），并同步 `matrix.md` §1.2 的
     concept 片段与特化表（learning-todo 子仓，另提交）。
  2. `ScalarTraits<Dual<T>>::conj/abs2` 里的成员名 `x.real`/`x.dual` → `x.value`/`x.derive`
     （`Matrix<Dual<float>>::Adjoint()` 实测：修前报 `no member named 'real' in 'Dual<float>'`，
     修后 `Adjoint[0][0] = (1.0, 0.1)` 对）。
  3. `jacobian.hpp` 的 `Derive(Matrix<T>&)` → 按值收参返回（`Backward` 实例化报
     `no matching function for call to 'Derive'`，两处；修后编译通过）。
  4. `dual.hpp` 手写 `operator=` → 删除。**触发 `-Wdeprecated-copy` 的条件是
     "有 user-provided 拷贝赋值 + 无 user-declared 拷贝构造"，把参数改成 `const&` 不解除**
     （先按 const& 改、告警仍在，实测）；删掉后拷贝/移动四个特化成员全由编译器合成，
     `d = 0.5f` 语义不变（探针验证），`-Wall -Wextra` 下四个 TU 告警数全 0。
  5. `grad_graph.hpp` 补两处"明确失败"：二元算子 arity ≠ 2、`Function`/`MatrixFunction`
     的 `matrixF` 为空（原来分别是 nullptr 解引用段错误、`matrixF[0][0]` 越界读）；
     实测分别抛出 `GradNode::Forward: 二元算子需要 2 个输入, 实际 1` 与
     `Function: matrixF 未设置 (空 0x0)`。
  6. `demo` 里 `SumAll` 节点的占位形状 `1x2` → `1x1`：建图期与运行期形状不符的隐患，
     实测 `y = sum(a)*k` 修后得 50（修前会在 `cwiseProduct` 抛尺寸不匹配）。
  7. **反向模式明确标记未实现**：`ComputationGraph::Backward` 改为抛 `logic_error`
     —— 实测跑它必然在第一节点抛 `Matrix dimension mismatch`（节点 adjoint 从未按形状初始化），
     这个报错会把人误导到 Matrix 上去。按 §7.5 属**待立项功能**，不是 bug；
     `GradNode::Backward` 保留为未验证草稿。
  8. **新增 ctest 用例 `scalar_traits`**（`demo/core/`，44 项检查）：静态断言钉契约
     （含 `Dual<Complex<double>>`）、`matrix.md` §1.3 铁律要求的 `a*b == b*a` 交换律单测
     （Complex / Dual / Dual<Complex> 三类型）、Complex 与 Dual 原函数对拍。
     这个用例本来是缺的 —— 有它，第 1 条那种漂移在编译期就会被抓住。
  **本次仍未提交**（留用户核验后自行 commit/push；父仓 4 文件 + 新增 2 文件，
  `matrix.md` 属 learning-todo 子仓需单独提交）。

### 2026-09-29

- **设计落档（工程线，文档层）**：`learning-todo/matrix.md` 并入第二轮泛型数值库评审定稿
  **D11–D18**（learning-todo 子仓 commit `a6f089c`，+133/−13）：
  - 新增 §1.4 标量层细则（D11 Smith 除法 / D12 除零两层 / D13 禁 is_xxx 探测 /
    D15 混合运算 / D16 四元数封死 / D18 命名收敛）、§2.5 BlockView 视图（D17）、
    §3.4 特征值 tag dispatch（D14）；§0 宪法补「分块靠视图不嵌套」；§8 决策表扩到 18 条。
  - §7 重写为 Dual 六条细化清单（对着 `dual.hpp` 现码逐条锚定：删 operator/ 的 ==0
    吞错分支、字面量走 traits、std::sin 改 ADL 两步法、sumLine/multLine 迁
    jacobian.hpp 并修 multLine 零因子除零 bug、DualFunc 签名改 T(const T&)、
    新增 ScalarTraits<Dual<T>> 逐分量共轭）。
- `graphics/plan.md` V1.5 补复数求导用例（父仓 commit `a7277d4`）：
  z² 单种子（全纯）、|z|² 双种子（非全纯）与解析值对拍，阈值 value 1e-12 / derive 1e-10。
- `matrix.md` 新增 **§10 C++ 特性注解**（learning-todo 子仓 commit `e7ec45d`）：
  用户尚未掌握文中涉及的 C++ 特性，按「先讲解决什么问题、再讲最小语法、标注本库落点」
  写了 9 小节（模板全特化/偏特化/typename/concept+requires/if constexpr/ADL 两步法/
  tag dispatch/numeric_limits·enum class·static_assert/两条写码纪律），正文 7 处加指路。
- 工程代码仍未动（`scalar.hpp` / `complex/` / `linalg/` 均不存在）；
  下一步按 `matrix.md` §9 顺序动代码：CMake 升 C++20 → scalar.hpp → complex.hpp →
  matrix.hpp 重写 → dual.hpp 修 → linalg 四件套 → V1.1–V1.7。

- 文档治理：统一四层结构（入口 `AGENTS.md` / 真源 `README.md`+本文 / 展开 `plan.md`+`learning.md`+`matrix.md` / 产物 `learning-todo/`）。
- 本次改动不动学习进度，仅治理文档；`plan.md` 的验收项编号由 `A1`–`F4` 改为 `V{阶段}.{序号}`，
  以消除与笔记代号（`B1`、`E1` 等）的撞名。
- 删除仓库根目录的 `requirements`（与 uv 工作流重复，且内容是 pip 源说明，会误导 agent 走 pip 路径）。
- 清理工作区：删除多余的 worktree `graphics-0-2fe0063d` 及其分支 `workbuddy/graphics-0-2fe0063d`
  （确认无独有提交，分支 `12dfde4` 已合并进 `graphics_0`）。现在只剩主仓与 `graphics-0-0c774982`。
  被删 worktree 的 `.workbuddy/memory/` 已备份到
  `graphics-0-0c774982/.workbuddy/memory/archive/2fe0063d-2026-09-22/`（含 md5 校验）。
- `notebook/` 笔记体检与修补（commit `b899d33`，18 文件）：
  - 修 bug：`diffrential/de4.tex` 重复的空标题「知识网络」（原 135/136 行），后续小节重编号为 25.2–25.5。
  - 内容缺口：`A_analysis.tex` 补 σ-代数「可数并而非任意并」的理由；「可数可加性」由三行定义
    扩为完整论证（有限可加太弱的三反例表 + 失效定理清单 + 不可数可加太强 + 集合列连续性 + Vitali 三难）。
  - 消冗余：新增 11 个 `\label` 与 6 处指路句，打通 F→K、G→de4、H→I、B→lie1/lie2、A→de3 的重复链路。
  - 新增 `partN.tex` §15.5「章节与教材对照索引」（教材定位表 + 精读/略读分级）。
  - `main.tex` 给两个专题加 `\part` 分界（专题一 p.41 / 专题二 p.61），与 A.–O. 字母编号拉开层次。
  - 顺手修两处既有缺陷：全篇 38 处 ASCII 直引号规范为中文引号；9 处被静默丢字的字符
    （7×`↔`、2×`✓` → `\leftrightarrow` / `\checkmark`，Latin Modern 无此字形）。
  - 校验：`latexmk -xelatex` 通过，76 页，0 缺字 / 0 Overfull / 0 未定义引用 / 0 重复标签；
    另用 pdfpages+sips 导出关键页目视复核（p.11–12 新论证、p.39 索引表、p.41 `\part`、p.48/53 符号修复）。
  - `.gitignore` 补 `*.xdv`（XeLaTeX 中间产物）。
- `notebook/` 内容补缺（第二次会话）：用户要求把 `learning-todo/chats/` 里的「计算必修清单」落成正式章节，
  并盘查其余缺口。新增/扩写五块，`main.pdf` 由 76 页增至 **112 页**（A–R 字母连续，两个专题顺延到 §19 起）：
  - **P 计算必修**（`P_computation.tex`，P1–P10）：由黑板 `ses_f5e6…ipynb` md-24～md-34 转写为 LaTeX。
    组号 **G1–G9 → P1–P9**（`G1`–`G5` 已被 `G_numeric.tex` 占用，撞名已消）；自检题按约定留在聊天/黑板，不入笔记。
  - **C6–C10 概率深化**（扩写 `C_partial.tex`，不重编号 C1–C5，避免 A_analysis 里「→ C4」「C2 整节失效」的引用失效）：
    收敛模式与层级 + 反例博物馆 + 均匀可积性 + Egorov/Scheffé/Vitali 收敛/单调类定理 + Borel–Cantelli + 两种大数定律、
    特征函数与 CLT 证明、马尔可夫链与 MCMC/MLT、蒙特卡洛与方差缩减、Itô 微积分与 Fokker–Planck。
  - **A4 Sobolev 空间与弱解**（扩写 `A_analysis.tex`）：弱导数、Sobolev 嵌入与紧嵌入、Poincaré、Lax–Milgram、正则性提升。
  - **Q 拓扑学**（新建 `Q_topology.tex`）：点集拓扑、紧性（含与弱紧性/Arzelà–Ascoli/Rellich–Kondrachov 的对接）、
    商空间、同伦与 $\pi_1$、单纯同调与 Euler 特征、de Rham 上同调；与 `lie2` 的「拓扑结构」明确分工。
  - **R 复分析**（新建 `R_complex.tex`）：Cauchy 理论、留数（补上 `de3` 缺的 Laplace 反演一环）、共形映射与共形参数化、
    解析延拓（接 K4 的 Wick 转动）、复势与二维流体。
  - 顺带：表格定宽列统一为 `L{}`（`\raggedright`，定义在两个 `main.tex`），消除 CJK 窄列拉伸；
    `partN.tex` 教材对照索引补 P/Q/A4/C6–C10 五行；`AGENTS.md` 代号速查补 P/Q/R。
  - 校验：0 缺字 / 0 Overfull / 0 Underfull / 0 未定义引用 / 0 重复标签；子文档 `diffrential/main.tex` 也能独立编译（0 报错）。
- 仍缺、**尚未动手**的内容（下次可从这些里挑）：抽象代数（群环域、模、Galois）、图论与组合、
  凸优化与变分法系统化、数值 PDE 理论（Lax 等价定理、von Neumann 稳定性分析）、有限元理论、连续介质力学/弹性、
  复几何与辛几何、测度论进阶（Radon–Nikodym 的完整证明、乘积测度与 Fubini 的构造）、偏微分方程弱解的存在唯一（Evans Ch.6 级）。
- **多仓基建提前启动**（本日第三次会话）：`learning-viz` 与 `learning-web` 各落一个 init 提交，
  不再等「第一个小项目完成」。起因是拿到一份「笔记 web 化」设计（另一 session 产出，已存为
  `learning-viz/python/notes_pipeline/DESIGN.md` 并作为该模块的设计真源）：
  - `learning-viz`（commit `f782f4c`）：`README.md`（定位/红线/跨仓路径坑）、`.gitignore`（`dist/` 不入库）、
    `python/` 骨架（`notes_pipeline/vendored`、`app`、`viz_modules`、`contracts`、`tests`）、
    `python/notes_pipeline/DESIGN.md`。
  - `learning-web`（commit `3e44024`）：`README.md`、`.gitignore`、`src/` 骨架（`shell`、`renderers`、`viz`、`modules`）+ `public/`。
  - **不做** CDN、不做前端计算；代号（`A1`、`G1.2`）作路由 ID，不用 `§NN.N`。
  - 设计经评审后**就地校正 10 处**（文中标 `[校正 Cn]`），要点：`pandoc --katex` 实为**客户端**渲染须改机制
    （改走构建期 Node 预渲染 KaTeX）；filter 正则 `[A-O]` → `[A-R]`；插槽代号 `A1.1.3`/`A1.1.4` →
    `A1.3`/`A1.4`（不存在三级代号）；`\providecommand{\vizslot}` 需加父/子**两个** `main.tex`；
    R1「Python 3.14 缺 wheel」经 PyPI 核实**已解除**（`pypandoc-binary 1.17` 有 `py3-none-macosx_11_0_arm64.whl`、
    `flask 3.1.3` 为 `py3-none-any`）；笔记规模由「21 文件 / 72 页 / 33 节」更正为
    **24 内容文件 / 112 页 / 202 个 `\subsection`**。
  - 本次只做骨架与文档治理，未装任何依赖、未改 `notebook/`；`\vizslot` 兜底与步骤 1 尚未开始。
  - 父仓文档同步：`README.md` 两仓状态「空仓」→「开发中」；`AGENTS.md` 文档地图新增 `DESIGN.md` 一行、
    第八节目录用途改写；本节即第 5 节里程碑提前的说明。
- **上一条的修正：笔记转换链路改归 `learning-web`**（本日第四次会话）。用户带来一份新的
  `learning-web` 入口规约（比上一版细得多，含四问分级判据、掌握度 L0–L5、四类图边），
  并拍板了它 §10 的两个未决问题：
  - **Q1 → 本仓**：知识图谱语义边真源放 `learning-web/graph/<chapter>.yaml`（原方案放父仓 `notebook/graph/`）。
    脱节风险由「yaml 头记 `notebook_commit` + 构建时与 `provenance.json` 比对、不一致告警」兜住。
  - **Q2 → `build/`**：转换脚本全部放 `learning-web/build/`。**连带口径收窄**：
    `learning-viz` 只做**计算模块**的后端（`/api/viz/*` + 本地渲染），**不再有 `/api/notes/*`**。
    → 上一轮放在 `learning-viz/python/notes_pipeline/` 的 `DESIGN.md` 与目录**已删除**（commit `4d093e3`）。
  - `learning-web`（commit `2a28d9f`）：新增 `AGENTS.md`（入口规约）、`DESIGN.md`（决策 D1–D8 + 契约样例 +
    A1 切片七步验收）、`build/`（含 `vendored/`）、`graph/`；`README.md` 改为短版人读总纲。
  - 新规约同样经评审**就地校正 9 处**（`AGENTS.md` §0.1，标 `[校正 Cn]`）：规模实测
    **202 `\subsection` / 162 `\subsubsection` / 112 页 / 24 内容文件**（原文「36 节 / 250+ 节点」）；
    代号只到**两级**，`A1.1.3`/`A1.1.4` → `A1.3`/`A1.4`；R1 由「高」降为「低」（PyPI 实测有 wheel）；
    `\providecommand{\vizslot}` 要加父/子**两个** `main.tex`；R3 的三类结构分布核实
    （A1 段内 3 张 `tabular` + 1 处 `\boxed`，`\xrightarrow` 不在 A1，在 `C`/`de4`/`P`/`Q`/`li*`）。
  - **新发现（写进 `AGENTS.md` §6 校正 C5，是环境事实而非笔误）**：本仓与 `learning-viz` 都在父仓目录树内，
    `uv run` 会**向上发现父仓项目**并使用 `~/develop/learning/.venv` —— 所以两仓都**不建自己的
    `pyproject.toml`**，`uv add` 会把依赖写进**父仓** `pyproject.toml` + `uv.lock`（提交也要在父仓）。
    已实测确认：在 `learning-web/` 下 `uv run python -c ...` 用的是 `~/develop/learning/.venv/bin/python3`。
  - **公式渲染机制已定**（决策 D8）：`pandoc --katex` 是**客户端**渲染且默认指向 CDN，不符要求；
    改为**构建期用 Node 跑 `katex.renderToString` 预渲染**成静态 HTML（自带
    `<annotation encoding="application/x-tex">`，离线、无闪烁、可搜 TeX 源码）。
  - 仍未装依赖、未改 `notebook/`；`build/` 与 `graph/` 目前是空目录。

- **2026-10-20：新增 `src/linalg/so3.hpp` + ctest `so3`（51 项 FAIL 0；ctest 15/15 → 16/16）**。
  起因是用户问「叉积到底要不要、`[u]` 要不要实现、写在哪里」。定案三层否决：
  - **不进 `core/`** —— core 是「标量代数 + 哑容器」，叉积是几何运算；且 §3.7 刚立过
    「不为一次性用途扩 core 接口」（与 `Rotate` 就地写 trace 同理）。
  - **不进 `quat/`** —— D16 把四元数封死在本目录且不依赖 linalg；放这儿会让 so(3)
    反向依赖 quat，层次就乱。
  - **不塞 `operator.hpp`** —— 那文件自己写了「Dot/Abs2 是暂放，S6 派支线再定去留」，
    别往待搬迁的抽屉里再塞。
  ⇒ 新建 `src/linalg/so3.hpp`，只 include `core/matrix.hpp`、不碰任何分解。三个入口：
  `Cross(a,b)`（手写三分量）/ `Hat(u)`（R³→so(3)）/ `Vee(K)`（so(3)→R³，读 (2,1)(0,2)(1,0)，
  与 `_rot_core.py:71/:79` 同口径可对拍）。**`quat/` 不 include 它** —— `Rotate` 里那两个
  叉积继续手写三行，为 3 行算术引入跨层依赖不划算；这是**故意的重复**，已在头文件注释
  里钉死，防以后有人来"去重"把 D16 拆了。
  - **口径**：入参出参一律 3×1 列向量（D19）；形状错 `throw invalid_argument`；
    **`Hat` 不归一化**（`hat(λu)=λ·hat(u)` 是线性的，Rodrigues 需要 θ 单独乘在外面）。
  - **判据**：全部走「两条独立路径对拍」，[A] Hat 反对称结构 / [B] Hat·Vee 双向互逆 /
    [C] Cross 基元乘法表+反交换+正交 / [D] Cross vs `Hat(a)·b` / [E] `Hat³ = −|u|²·Hat`
    （Rodrigues 闭式成立的唯一关键，report.md §1）/ [F] 纯四元数积的向量部分 = 叉积 /
    [G] 形状契约 throw。分域 double + float 各一遍（float 阈值 1e-5）。
  - **注入对照（行为型，已回滚 + shasum 验过）**：把 `Hat` 整体翻号（等价左手系）⇒
    `[A3][B1][B2][D1][D2]` 共 10 项变红。**但 `[A1][E1][E2][E3]` 仍全绿** ——
    `Hat³ = −|u|²·Hat` 两边都是奇次，整体翻号互相抵消。
    ⇒ **判据盲区，已写进 `demo/linalg/so3.cpp` 注释：符号由 [A3]（六位置逐位）和
    [D1]（Cross 对拍）钉，不由 [E] 钉。** 与 `quaternion.hpp:110` 那条
    「叉积符号写反，往返测试抓不到」是同一类盲区。
- **2026-10-23：新增 `graphics/demo/pyapps/tensor-convention.md`（228 行，9 节）**。
  为「将来如果要做物理 demo」提前钉死的张量入场规约 —— **只有约定、没有代码**，
  **不进 C++ 库**（红线 #4 + `matrix.md` §6），**只活在 `demo/pyapps/`**。
  格式与数学口径跟随 `pyapps/sphere/_core.py` 头部（ASCII：`R^i_{jkl}`、`A' = R·A·R^T`）。
  - **升降按场景开关**（不再二选一）：场景 E（g=δ）**不写**升降函数（写了是空转自欺）；
    场景 G（球面非平凡 `g_ij`）/ 场景 R（Minkowski）**必须写**，R 还要先钉号差
    `(-,+,+,+)` vs `(+,-,-,-)`。每个脚本头部第一行必须声明场景（三选一）。
  - **上下标必须进变量名**（`V_up` / `V_down` / `A_up_down`）— numpy 数组只有轴没有上下标，
    这是全文唯一**靠命名而非靠测试**保证的规则（判据只在出错之后才报警）。
  - **反对称完整存储、只作判据不作存储约定**（否决"半存 + `get_F`"）：省不下东西，
    且符号一致地错时测试全绿 —— 与 `so3.hpp` 那次 `Hat` 整体翻号同盲区。
  - einsum 字符串即文档；**禁写 `contract(T,V,axes_T,axes_V)` 通用封装**
    （判据：封装是否携带语义）。手写 for 循环**鼓励写一遍**与 einsum 对拍。
  - 验收不用 `assert`（`python -O` 会关掉、失败无实测值），走 pyapps 现成口径
    「末尾打印数字 + PASS/FAIL」。`swapaxes` 必指明哪两轴 —— `.T` 对 ≥3 阶是全轴反转。
