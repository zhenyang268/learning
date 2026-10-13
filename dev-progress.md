# dev-progress.md — 进度与排期

> **读者**：人 + agent ｜ **最后更新**：2026-10-14 ｜ **本文是「进度」的唯一真源**
>
> **真源声明**
> - 本文是**总纲式进度真源**：回答「整体现在在哪、下一步做什么」。
> - `graphics/plan.md` 只负责任务定义与验收阈值（验收项记作 `V{阶段}.{序号}`），**不记录状态**。
> - **数学线的关卡明细，真源是 `learning-todo/ns-progress.md`**（逐次会话日志 + 关卡状态）。
>   本文只保留一行指针，**不重复记录**，以免两处漂移。
> - 冲突处理：工程线以本文为准；数学线以 `ns-progress.md` 为准。
>
> 分工：AI 排期 / 定门槛 / 验收，用户写代码。优先学习，工程增量推进。

---

## 0. 当前指针（每次会话先看这里）

> **数学线在 L1（测度与积分）中段；工程线在 `graphics` 阶段 1 中段，实际代码比 2026-09-29 日志领先（见第 3 节）。**

| 线 | 当前位置 | 下一个动作 |
|---|---|---|
| **数学（优先）** | L1 第 3 关进行中 —— **明细见 `learning-todo/ns-progress.md` 的「📍 当前位置」**（本表不重复记录） | Vitali 不可测集 → Carathéodory 判据 → 积分四步构造 |
| **工程** | `graphics` 阶段 1 中段：C++20 ✔、scalar/complex/matrix 重写 ✔、dual 修 ✔；**2026-10-13/14：`compute_graph` 基线成立（forward-mode）、接口层 bug 修完（标量契约 `fromReal` / Dual·Complex 路径 / 反向模式标记）、ctest 2/2 全绿，并新增 `graphics/run_tests.sh` 测试入口**。但 **V1.1–V1.7 一条都没通过**：linalg 528 行**从未被任何编译单元 include 过**（探针实测 32 个编译错误），`svd.hpp`/`expm.hpp` 仍是 0 字节，`demo/common.hpp` 不存在 | ① `demo/common.hpp` ② 修 `solve.hpp` ③ 修 `inverse.hpp` ④ 重写 `eigen.hpp` 为真 Jacobi 迭代 ⑤ `svd.hpp` ⑥ `expm` ⑦ 建 5 个验收 demo（清单见 §4） |
| **小项目** | 未解锁 | 等 Gate A（见第 2 节） |
| **笔记 web 化** | `learning-web` 已落**规约 + 设计 + 目录骨架**（`AGENTS.md` / `DESIGN.md` / `build/` / `graph/`），实现未开工 | 按 `learning-web/AGENTS.md` §12 第 1 步：`uv add --dev pypandoc-binary` |

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

| 门槛 | 数学 | 工程 | 解锁的小项目 |
|---|---|---|---|
| **Gate A** | L0 + ODE 数值解直觉（Euler/RK4/Verlet） | matrix 改造 ✔ + 最小 linalg（LU/solve）+ 一个积分器 | 动力系统可视化：单摆 / 行星轨道 / 双摆 / N 体 |
| **Gate B** | **L1 完成**（测度 + 积分 + 收敛定理 + 四步构造） | Gate A + FDM 基础 | PDE 可视化：热方程（图像模糊）、波动（鼓膜） |
| **Gate C** | L2 泛函 + Sobolev | 阶段 4 完整 | 弱解类（通往 NS 真目标） |

> **第一个小项目定在 Gate A**；L1 继续按 NS 路线推进，不互相阻塞。

---

## 3. 现状明细（截至 2026-10-14）

- **数学线**：**只以 `learning-todo/ns-progress.md` 为准**，本文不重复关卡明细。
- **工程线**（2026-10-13 实地盘点，代码在分支 `graphics_0`）：
  - 已落地（编译层验证过）：CMake C++20 ✔（`set(CMAKE_CXX_STANDARD 20)`）；
    `core/scalar.hpp`（Scalar concept + traits）✔；`core/complex.hpp`（130 行）✔；
    `core/matrix.hpp` 重写 ✔（静态工厂 Zero/Ones/Identity/Diagonal、扁平行主序、无缓存层）；
    `autodiff/dual.hpp` 修 ✔（multLine 零因子 bug 已按 matrix.md §7.1 修为前后缀积、
    operator/ 除零分支已删、`using std::` ADL 两步法已用）。
    **更正（2026-10-14 探针实测）**：上面记的「`linalg/solve.hpp` LU+Solve ✔（98 行）」
    不成立 —— 原 ✔ 的含义是"文件存在"，不是"能编译"。实例化实测报错，详见下面 linalg 小节。
    `linalg/` 下真正能通过实例化的只有 `cholesky.hpp`。
  - **原卡点已解（2026-10-13 会话）**：`demo/autodiff/compute_graph.cpp` 编译不过，
    **真实根因不是「`LnMultAll`/`Function` 未声明」** —— 这两个名字是 `GradOp` 的合法枚举值
    （`common.hpp`）。实际是三类问题：① `GradOp` 改成 `enum class` 后 `grad_graph.hpp` 的
    `case X:` 标签没加作用域前缀（35 处，连带 `SumAll`/`MultAll` 撞上 `jacobian.hpp` 同名函数，
    报成 `<overloaded function type> is not convertible to GradOp`）；② 类内 `Matrix(1, 1)`
    的 CTAD 失败（2 处，须写 `Matrix<T>(1, 1)`）；③ `Output`（透传）与 `Mult`（逐元素乘）
    两个分支缺失 —— 这类**能编译但静默出零**的问题比编译错误更危险。
    已全部修复，`ctest` 1/1 通过。
  - **基线边界（仍存在，见第二批修复说明）**：未实现算子（`CWisePow`/`LnMultAll`/
    `CrossEntropy`/`Sigmoid`/`Relu`/`Tanh`/`Softmax`）在 `Forward` 里统一 **throw「未实现」**
    （原来是静默 `break`，会返回零值）；反向模式整体未实现。
  - **linalg 实际状态（2026-10-14 探针实测，取代原先的"空壳文件"描述）**：
    自 10-13 起该目录有变动 —— `dot.hpp` 已删（`grep` 确认无引用残留），
    新增 `cholesky.hpp`(67) / `qr.hpp`(133) / `householder.hpp`(0) / `expm.hpp`(0)，
    `eigen.hpp` 38→89 行、`inverse.hpp` 111 行、`operator.hpp` 30 行。
    **但 `demo/` 只有 `scalar_traits` 与 `compute_graph` 两个 ctest，没有任何 TU include 过
    `linalg/`** —— 528 行代码一次都没进编译器。用 `-fsyntax-only` 实例化各入口实测，
    两轮共 32 个编译错误：
    - `solve.hpp`：`SwapRow` 定义 2 参而调用传 3 参；函数体用未声明的 `max_row`；
      `Tol` 未定义；`piv[j]` 应为 `result.piv[j]`；主元比较写裸 `>`（违反 D7，
      Complex/Dual 无 `operator>`）；**消元写成 `U(i,k) -= L(i,j)*U(i,k)`（应为 `U(j,k)`）**
      且整段被 `if (max_row != i)` 包住 —— 这两条即使编过数值也全错。
    - `inverse.hpp` ✅ **已修（2026-10-14，S2）**：`scalarTraits` 拼写；`Identity<T>(n)` 应为
      `Matrix<T>::Identity(n)`；`diagonal_inverse` 的 `tol` 未声明；`Cholesky`/`Solve` 未 include；
      `Inverse` 只有声明 → 已补分发。另有两条**编过会算错**的探测 bug：`is_symmetric` 比较对象
      写成 `AH(j,i)`（实方阵恒判 Symmetric）、`is_orthogonal` 只查对角线。详见 §6 同日条。
    - `eigen.hpp`：`is_symmetric(A)` 缺 `tol`；`Identity(n,n)` 只收 1 参；`ScalarTraits()::one()`
      语法错；`jacobi` 定义在使用点之后；**实现是 2x2 单次旋转，不是 Jacobi 迭代**。
    - `qr.hpp`：`real_t` 未定义、`aplha`/`Cows` 拼写、`QRComposition` 大小写、
      `ScalarTraits<double>::Zero` 不存在、`Dot` 调不到。
    - `operator.hpp`：`Dot` 是成员函数残骸（自由函数里用 `this->rows/cols/data`）。
    - `cholesky.hpp` **是唯一能通过实例化的**。
    - 空文件：`svd.hpp`（V1.4 必需）、`expm.hpp`（V1.6 必需）、`householder.hpp`。
    - 本阶段要求的 `demo/common.hpp`（`IsApprox` / 随机矩阵生成）**不存在**。
    完整缺口清单与收尾最小路径见 `graphics/plan.md`「待办 / 已知问题」2026-10-14 条。
  - **第二批修复（2026-10-13 同会话，基线的"接口层 bug"）**：把第一批列出的"已知边界"里
    属于**真 bug** 的三条改掉 ——
    ① `scalar.hpp` 的标量契约漂移：代码写的是设计文档里不存在的 `const_n(内层 T)`，
      而 §7.1 第 2 条要的是 `const_n → fromReal(real_t)`；后果是 `Scalar<Complex<double>>`
      与 `Scalar<Dual<float>>` **双双不成立**（`static_assert` 实测），任何
      `requires CommutativeScalar<T>` 的 linalg 函数都拒绝它们 —— **这是 V1.1 与 V1.5 的拦路石**。
      已按设计改名 + 补进 concept（顺带堵住"泛型代码裸写字面量"这个被 §7.1 点名的隐蔽假设）；
      连带修 `ScalarTraits<Dual<T>>::conj/abs2` 的成员名（`x.real`/`x.dual` → `x.value`/`x.derive`）；
    ② `jacobian.hpp` 的 `Derive(Matrix<T>&)` 接不了 `Apply(...)` 的临时量 → 改按值收参返回；
    ③ `dual.hpp` 手写的 `operator=` 触发 `-Wdeprecated-copy`（**触发条件是"有 user-provided
      拷贝赋值 + 无 user-declared 拷贝构造"，把参数改成 `const&` 并不解除** —— 已实测）
      → 删掉冗余的重载，交给编译器合成（拷贝/移动四个特化成员全在，`d = 0.5f` 语义不变，实测）；
    ④ `grad_graph.hpp` 补两处"明确失败"：二元算子 arity 不匹配、`Function` 节点 `matrixF`
      未设置（原来分别是 nullptr 解引用与 `matrixF[0][0]` 越界读）；
    ⑤ demo 里 `SumAll` 节点的占位形状 `1x2` → `1x1`（建图期形状与运行期不符的隐患，
      `y = sum(a)*k` 会抛尺寸不匹配，实测）。
  - **反向模式已明确标记未实现**：`ComputationGraph::Backward` 改为抛 `logic_error`
    （实测：节点 adjoint 从未按形状初始化，第一节点就抛 `Matrix dimension mismatch`，
    报错会把人误导到 Matrix 上）；`GradNode::Backward` 保留为草稿并注明未验证。
    按 `matrix.md` §7.5 这属于**待立项功能**而非 bug。
  - demo 注册：`scalar_traits`（新增，`demo/core/`，44 项检查）+ `compute_graph`（36 项检查），
    **ctest 2/2 通过**；`a_triangle` 无 CMakeLists 未注册。
  - **V1.1–V1.7 验收一个没跑。demo 全绿 ≠ V1 验收过（demo 是 forward-mode 计算图，
    与 V1.1–V1.4 的 linalg 验收是两条线）。**

---

## 4. 下一步

**数学线（优先）**

> 关卡明细与逐条状态**只以 `learning-todo/ns-progress.md` 为准**，本文不重复记录。
> （2026-09-23 移除本节原有的 `- [ ]` 勾选清单：它长期挂着「外测度未学」，
> 与 §0 指针及 `ns-progress.md` 直接矛盾 —— 属「两处记录」漂移，按真源约定删除。）

**工程线**（用户自定节奏；任务定义与阈值见 `graphics/plan.md` 阶段 1）

- [x] 修 `grad_graph.hpp` / `compute_graph.cpp` 编译错误 → **2026-10-13 完成**（真因见 §3 第 2 条：
      `enum class` 作用域前缀 35 处 + `Matrix(1,1)` CTAD 2 处 + 缺 `Output`/`Mult` 分支）；
      `ctest --test-dir graphics/build` 2/2 通过
- [x] 修「接口层」bug（同一会话第二批，见 §3）：标量契约 `const_n → fromReal` + concept 补条目、
      `ScalarTraits<Dual<T>>::conj/abs2` 成员名、`Derive` 收临时量、`Dual::operator=` 弃用告警、
      arity/`matrixF` 越界保护、`SumAll` 占位形状；反向模式明确拒绝执行。
      新增 ctest `scalar_traits`（44 项）+ `compute_graph`（36 项）
- [x] 落「基线」提交 → **2026-10-13/14 完成**：父仓 `graphics_0` 三次提交 ——
      `375f618`（基线修复）、`c82c7ed`（新增 `graphics/run_tests.sh` + AGENTS.md 登记）、
      `1f33b7a`（脚本补 `--build-only` / `--test-only`）；子仓 `learning-todo` `e733363`（`matrix.md`）。
      **均未 push**（推送由用户自行决定）。
- [ ] **V1.1 前置已清**（标量契约成立，`Scalar<Complex<double>>` / `Scalar<Dual<float>>` 均成立），
      但 V1.1 本身仍未通过 —— 要等 linalg 落地后，按「double 与 Complex 共用同一份源码」实跑。
- [ ] **V1 收尾总清单（2026-10-14 探针体检后重排；缺口明细见 `graphics/plan.md` 同日条）**：
      **进度：① ② ③ 已完成**（`demo/common.hpp` + `solve.hpp` 数值层 + `inverse.hpp` 分发；
      ctest 4/4、`solve` 34/34、`inverse` 48/48；见 §6 同日条）；余下 ④–⑧ 未动。原文如下：
      ① `demo/common.hpp`（`IsApprox` + 实/复/Dual 随机矩阵生成，plan 阶段 1 明确要求）；
      ② 修 `solve.hpp`（`SwapRow` 签名、`max_row`、`Tol`、`result.piv`、主元比较改走 `abs2`、
      **消元改为 `U(i,k) -= L(i,j)*U(j,k)`**、去掉 `if (max_row != i)` 包裹）→ 与 naive 高斯消元对拍；
      ③ ✅ 修 `inverse.hpp`（拼写 / `Matrix<T>::Identity` / `tol` 传参 / 补 include）+ 实现 `Inverse` 分发；
      ④ 重写 `eigen.hpp` 为**真 Jacobi 迭代**（收敛循环 + 非对角元判据 + 特征值排序）；
      ⑤ 实现 `svd.hpp`（经 A^H·A 特征分解）；⑥ 补 `expm`（V1.6 需要）；
      ⑦ 建 5 个验收 demo（V1.1+V1.2 / V1.3 / V1.4 / V1.5 / V1.6+V1.7）并注册 ctest；
      ⑧ 加一个「linalg 编译烟测」用例（只 include + 实例化各入口），防止"文件在、编译器从没跑过"复发。
      **进度更新（2026-10-14）**：`cholesky.hpp` **分解侧已完成**（S2a，见 §6 同日条，ctest 5/5、
      `cholesky` 98/98、V2.2 通过）。④ 之前先插一步 **S2a-2**：打开 `Classify` 的 SPD 档 +
      接上 `cholesky_inverse` / `is_positive_definite`（`inverse.hpp` 三处当前仍显式 throw）。
- [x] ~~检查 `Solve` 落盘版细节（消元循环里 `if (max_row != i)` 可疑）~~
      → **2026-10-14 探针证实"可疑"为真**：该 `if` 会整段跳过消元，且消元公式本身用错了行
      （`U(i,k)` 应为 `U(j,k)`）。详见 §3 的 linalg 小节。
- [x] ~~写 V1.2 验收 demo~~ / `inverse.hpp` / `eigen.hpp` / `svd.hpp` 四条 → 已并入上面的 V1 收尾总清单。
- [ ] `qr.hpp` / `operator.hpp` 是 10-13 后新建的半成品（`qr` 不在阶段 1 清单）：
      要么修到能编译再进阶段 2，要么先移出构建视野 —— 但现在它们连编译都不过（见 §3）。

---

## 5. 里程碑

- [ ] Gate A 达成 → 启动第一个小项目（动力系统可视化）
- [x] ~~第一个小项目完成 → 拉 `viz` / `web` 多仓基建~~
      → **已提前至 2026-09-22**：不等 Gate A，直接以「笔记 web 化」启动两仓（理由见第 6 节日志）
- [ ] 笔记 web 化 A1 垂直切片走完七步（`learning-web/AGENTS.md` §12 / `DESIGN.md` §5）

---

## 6. 进度日志（append-only：只追加，不修改历史条目）

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
  - **提交**：S2a 由用户放行后单独一刀（2026-10-14，未 push）。

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
