# dev-progress 现状与下一步快照（截止 2026-10-15）—— 已过期，仅供考古

> 从根 `dev-progress.md` 的 §3「现状明细」+ §4「下一步」整节迁出（2026-10-26）。
> **这不是真源**：现状看 `ctest` + `dev-progress.md` §0；现在做什么看根 `NOW.md`。
> 保留原因：里面有 V1.1–V1.7 验收落点的实测数字、S1–S6b 提交序、autodiff 域分裂的详细成因，
> 这些在同月日志里是散落的，此处是集中版。**不要再往里补内容。**

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

