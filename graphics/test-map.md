# test-map.md — 接口 ↔ 测试的可追溯路径（`graphics/`）

> **定位**：只回答「这个接口被什么测了、测到什么程度、用什么阈值」。
> **不写**怎么做（`acceptance.md`）｜**不写**数值契约（`numeric-contract.md`）｜
> **不写**进度与排期（根 `NOW.md`）。
>
> **为什么要有这份文件**：`demo/` 下的 17 个 ctest **不是悬垂的测试脚本**，
> 每一个都挂在 `src/` 的**具体接口**上，且头部注释里写着「多少算对」。
> 这条路径目前只存在于**每个文件自己的注释**里 —— 全局看不可见，
> 于是「哪个接口被漏测了」没人能一眼回答。本文就是那张全局图。

---

## 0. 三层结构：先认清追溯链

```
第 1 层  接口清单（真源：demo/linalg/smoke.cpp 的 EntryTable，40 条）
            ↓  取函数指针 ⇒ odr-use ⇒ 强制实例化
第 2 层  专项用例（16 个，每个覆盖一组接口，判据 18~266 项）
            ↓  头部注释即判据规格（「判据」一节，编号 ①②③…）
第 3 层  判据实现（文件内的 [A]/[B]/[C] 段与 [P] 占位段）
```

**任何一层的断点都是真问题**：
- 第 1 层漏项 → 接口从未被编译器见过（`acceptance.md` §2 闸门的意义）
- 第 2 层缺失 → 接口只有「跑得通」级验证，没有权威阈值
- 第 3 层是 `[P]` 占位 → **如实记录「它现在不可用」**，而不是假装可用

---

## 1. 全仓现状（2026-10-28 实测）

**17 个 ctest，合计 1175 项判据，全绿。**

| # | ctest | 判据数 | 覆盖模块 | 权威性 |
|---|---|---:|---|---|
| 1 | `scalar_traits` | 44 | `core/scalar.hpp` + `core/complex.hpp` + `autodiff/dual.hpp` | **权威**（concept 契约 + `a*b==b*a` 单测） |
| 2 | `tensor` | 24 | `core/tensor.hpp` | **占位**（`[P1]`–`[P24]` 全为「断言仍 throw」） |
| 3 | `compute_graph` | 36 | `autodiff/grad_graph.hpp` | **权威**（含三方对拍） |
| 4 | `apply` | 18 | `autodiff/jacobian.hpp` 的 `Apply` | 权威 |
| 5 | `backward` | 31 | `autodiff/grad_graph.hpp` 反向 | **权威**（最小闭环三方对拍） |
| 6 | `jacobian` | 18 | `autodiff/jacobian.hpp` | 权威 |
| 7 | `solve` | 34 | `linalg/solve.hpp` | **权威** |
| 8 | `inverse` | 85 | `linalg/inverse.hpp` | **权威** |
| 9 | `cholesky` | 98 | `linalg/cholesky.hpp` | **权威** |
| 10 | `eigen` | 119 | `linalg/eigen.hpp` | **权威** |
| 11 | `svd` | 266 | `linalg/svd.hpp` | **权威**（判据最多） |
| 12 | `qr` | 189 | `linalg/qr.hpp` | **权威**（三个分解各一遍） |
| 13 | `dual_complex` | 22 | `linalg` × `autodiff` 交叉 | 权威 |
| 14 | `expm` | 25 | `linalg/expm.hpp` | 权威（Taylor + Padé 互为对照） |
| 15 | `so3` | 51 | `linalg/so3.hpp` | **权威**（两条代码路径对拍） |
| 16 | `smoke` | 56 | **40 个公开入口** | **轻量**（阈值 `1e-10`，只判跑得通） |
| 17 | `quaternion` | 59 | `quat/quaternion.hpp` | **混合**（`[A]/[B]` 权威，`[P]` 占位） |

**读法**：`smoke` 的 56 项**不是**给某个算法定阈值 ——
它的 `[B]` 段每个入口只跑一次、阈值刻意放松（源码原话：「专项判据才是权威」）。

---

## 2. 接口 → 用例 映射表（40 个公开入口逐条）

**覆盖强度分级**：
- **★** 有专项用例，权威阈值
- **○** 只有 `smoke` 的轻量判据（跑得通 + 方向对）
- **▣** 占位（`[P]` 段，断言仍然 throw）

### 2.1 `linalg/solve.hpp`（4）

| 入口 | 覆盖 | 权威出处 | 关键判据 |
|---|---|---|---|
| `LU(A)` | ★ | `solve` ① ② | `P·A = L·U`，**`P` 由 `piv` 独立重建，不复用 LU 内部逻辑** |
| `Solve(A,B)` | ★ | `solve` ③ ④ ⑤ | `‖AX−B‖∞/(1+‖B‖∞) < 1e-10`；与「朴素列主元高斯消元」**逐元素对拍** |
| `SwapRow` | ★ | `solve`（结构判据内） | 行交换后 `L` 的乘子顺序必须同步换 |
| `SwapMultipliers` | ○ | `smoke [B]` only | **无专项判据** —— 见 §5 弱覆盖清单 |

### 2.2 `linalg/inverse.hpp`（11）

| 入口 | 覆盖 | 权威出处 | 关键判据 |
|---|---|---|---|
| `Inverse(A)` | ★ | `inverse` ① ② ③ | `‖A·A⁻¹−I‖∞ < 1e-8` **且** `‖A⁻¹·A−I‖∞ < 1e-8`（**两个方向都测**） |
| `Inverse(A, hint)` | ★ | `inverse` ⑤ | 与 `Inverse(A)` 结果一致 —— **`hint` 只跳过探测，不改结果** |
| `Classify(A, tol)` | ★ | `inverse` ④ | 对角/旋转/随机/对称分别落到对应 `kind` |
| `is_diagonal` | ★ | `inverse` ③ | 结构分支：对角 → 逐元倒数 |
| `is_symmetric` | ★ | `inverse` ④ | |
| `is_orthogonal` | ★ | `inverse` ③ | 正交（旋转）→ `A^T` |
| `is_positive_definite` | ★ | `cholesky` ⑥ | **判据本体 = Cholesky 走通** |
| `lu_inverse` | ★ | `inverse` ① | |
| `diagonal_inverse` | ○ | `smoke [B]` only | **无专项判据** —— 但它有独立奇异判据（`numeric-contract.md` §3.1） |
| `cholesky_inverse` | ★ | `inverse` ⑦ | SPD 档，与 LU 路径对拍 |
| `qr_inverse` | ★ | `inverse` ⑧ | `A⁻¹ = R⁻¹·Q^H`；**必查 `Q^H` 不是 `Q^T`**（实域等价，复域错） |

### 2.3 `linalg/cholesky.hpp`（2）

| 入口 | 覆盖 | 权威出处 | 关键判据 |
|---|---|---|---|
| `cholesky(A)` | ★ | `cholesky` ① ③ ④ | `‖L·L^H−A‖∞ < 1e-8`；严格上三角全 0；对角元为正实数 |
| `cholesky_solve` | ★ | `cholesky` ② | 与 LU 解之差 `< 1e-8`，自身残差 `< 1e-10` |

### 2.4 `linalg/eigen.hpp`（3）

| 入口 | 覆盖 | 权威出处 | 关键判据 |
|---|---|---|---|
| `eig(A)` | ★ | `eigen`（非对称分支） | 非对称走通用路径 |
| `eig_symmetric(A)` | ★ | `eigen` ① ② ③ ④ | 重构 `max|A−QΛQ^H| < 1e-8`；`‖Q^HQ−I‖∞ < 1e-12`；`Σλ = tr(A)` |
| `jacobi_rotation` | ○ | `smoke [B]` only | **无专项判据**（它是 `eig_symmetric` 的内部件，被间接覆盖） |

### 2.5 `linalg/svd.hpp`（1）

| 入口 | 覆盖 | 权威出处 | 关键判据 |
|---|---|---|---|
| `svd(A)` | ★ | `svd` ①–⑥ | 重构 `1e-8`；正交 `1e-12`；**独立交叉验证 `Σσ² = ‖A‖_F²`**（长方阵也成立）；`σ` 非负降序；秩亏 → `U` 走正交补 |

**注**：SVD 走 `A^HA`，**可信线 `σ/σ_max ≳ 1e-6`**，`< 1e-8` 当 0
（`numeric-contract.md` §4.3 ③）。

### 2.6 `linalg/qr.hpp`（3）

| 入口 | 覆盖 | 权威出处 | 关键判据 |
|---|---|---|---|
| `QR_Householder` | ★ | `qr` ①–⑤ | 重构 `1e-8`；`‖Q^HQ−I‖∞ < 1e-12`；`R` 严格下三角恒 0；**Q 是 m×n、R 是 n×n（thin）** |
| `QR_givens` | ★ | 同上 | 同上 |
| `QR_Shimit` (MGS) | ★ | 同上 | 同上 |

**为什么三条都测、且不直接比 `Q`**（`qr.cpp` 原话）：

> 同一个 `A`，三个实现给出的 `Q` **差一个对角相位**（`R` 对角符号约定各不相同），
> 所以 `Q` 本身不可比；但 `Q·R = A` 与 `Q^H·Q = I` **与约定无关**。

**这是「判据必须独立于实现路径」的教科书案例** —— 判据选的量要对**所有合法实现**都成立。

### 2.7 `linalg/expm.hpp`（**2 个入口，闸门只登了 1 个**）

| 入口 | 覆盖 | 权威出处 | 关键判据 |
|---|---|---|---|
| `expm(A)` | ★ | `expm` 25 项 `[A]`–`[D]` | 幂零恒等式（**判精确**）、`expm(0)`、`det(exp X) = exp(tr X)`、`RᵀR=I`、`tr R = 1+2cosθ` |
| `expm_pade(A)` | ★ | `expm` `[E]` 段（7 处） | 与 Taylor 版**互为对照基线** |
| **⚠️ `expm_pade` 不在 `smoke.cpp` EntryTable** | **漏** | — | **见 §4.1 已知漏项** |

### 2.8 `linalg/operator.hpp`（2）

| 入口 | 覆盖 | 权威出处 | 关键判据 |
|---|---|---|---|
| `Dot(a,b)` | ★ | `smoke [B12]` 是它的**第一个调用方** | `<v,v> = Σ|v_i|²`，**复数域必须共轭第一个参数** |
| `Abs2(x)` | ★ | `smoke [B12]` 同处 | 与 `Dot` 钉在一起 |

### 2.9 `linalg/so3.hpp`（3）

| 入口 | 覆盖 | 权威出处 | 关键判据 |
|---|---|---|---|
| `Cross(a,b)` | ★ | `so3 [C] [D] [F]` | 基元乘法表 + 反交换 + 正交；**手写三分量 vs 走 `Matrix` 乘法两条路径对拍**；**与四元数纯虚积对拍**（唯一同时穿过三条线的恒等式） |
| `Hat(u)` | ★ | `so3 [A] [B] [E]` | `K + K^T = 0`、对角**精确** 0；`Vee∘Hat = id` **且** `Hat∘Vee = id`；`Hat³ = −|u|²Hat` |
| `Vee(K)` | ★ | `so3 [B]` | 读 `(2,1)(0,2)(1,0)`，与 `_rot_core.py:71/:79` 同口径 |

**已知判据盲区**（`so3.cpp` 与 `acceptance.md` §2 都记了）：
`Hat³ = −|u|²Hat` 与 `Vee∘Hat = id` **抓不到 `Hat` 的整体符号翻转**
（两边都是奇次，翻号互相抵消）⇒ **符号必须靠六个非对角位置逐位比对 + `Cross` 对拍钉。**

### 2.10 `core/` 与 `autodiff/`

| 入口 | 覆盖 | 权威出处 |
|---|---|---|
| `Scalar` / `CommutativeScalar` concept | ★ | `scalar_traits`（**编译期静态断言**） |
| `ScalarTraits<T>::{zero,one,fromReal,conj,abs2,real}` | ★ | `scalar_traits` |
| `Complex<T>::operator/` | ★ | `scalar_traits`（Smith 1962 算法，**除零不判**） |
| `Matrix<T>` 13 个成员 | ○ | `smoke [A]`（取地址强制实例化）+ 各用例实际使用 |
| `Tensor<T>` 24 个入口 | ▣ | `tensor`（`[P1]`–`[P24]` 断言**仍然 throw**） |
| `Apply` / `Derive` | ★ | `apply` |
| 反向图 `Backward` | ★ | `backward`（含最小闭环三方对拍） |
| `Jacobian(f,x)` | ★ | `jacobian`（分子布局，一次 sweep 填第 j 列） |

### 2.11 `quat/`（1 类）

| 入口 | 覆盖 | 说明 |
|---|---|---|
| `+ - 一元- 标量乘`、`Conjugate`、`SquaredNorm`、`Norm`、`operator*` | ★ | `quaternion [A] [B]` 真判据 |
| `Inverse`、`Normalize`、`FromAxisAngle`、`ToMatrix`、`FromMatrix`、`Rotate`、`Slerp`、`Nlerp`、`AngleOf` | ▣ | `[P]` 段**断言仍然 throw**（10 处） |

**`[P]` 的结构是一套机制，不是偷懒**（`quaternion.cpp` 原话）：

> 你实现掉一条，对应的 `[P]` 就会**变红** —— 红就是提醒，提醒你把占位判据换成真判据。
> 留空的话，实现了也没人知道，没实现也没人知道。

---

## 3. 判据的六种类型（写新用例时照这个分类）

从 17 个用例归纳出来的，**不是每个都用满六种，但每种都要考虑**。

| 类型 | 问什么 | 例子 | 不受什么影响 |
|---|---|---|---|
| **① 结构性** | 形状/三角/对角是否**精确**成立 | `L` 严格上三角全 0；`K+K^T=0`；`R` 严格下三角恒 0 | 舍入（这些位置**必须精确 0**） |
| **② 重构** | `f` 的正逆是否闭合 | `‖Q·R−A‖∞`、`‖L·L^H−A‖∞`、`‖UΣV^H−A‖∞` | 实现路径 |
| **③ 独立不变量** | 有没有一条**与实现无关**的守恒量 | `Σσ²=‖A‖_F²`、`Σλ=tr(A)`、`det(exp X)=exp(tr X)` | **实现路径**（最强判据） |
| **④ 解析对拍** | 手算的闭式解对不对 | `2×2` 闭式特征值；`2×2` Hermitian 的 `L`；`|R11|,|R22|` | 数值方法（有闭式就用闭式） |
| **⑤ 跨实现对拍** | 两条**独立代码路径**是否一致 | SVD vs `eig(A^HA)`；`Cross` 手写三分量 vs `Matrix` 乘法；三个 QR 分解 | 单一实现的 bug |
| **⑥ 异常路径** | 该抛的是否抛、抛对类型 | 奇异/非方阵/`m<n`/空矩阵全 throw；**类型符合 `numeric-contract.md` §5** | — |

**判据强度排序**（`acceptance.md` §1 第 2 条的精神）：

> **③ 独立不变量 > ⑤ 跨实现对拍 > ④ 解析对拍 > ② 重构 > ① 结构性 > ⑥ 异常**

前四种**都不复用被测实现的逻辑**，是最可信的。而「拿实现 A 的输出当实现 A 的输入验」是最弱的。

**反面教材（本仓真实发生）**：`Hat³ = −|u|²Hat` 属于 ③，看起来很强 ——
但它对**奇次符号翻转不敏感**，所以 ③ 也有盲区，必须配 ① 逐位比对。

---

## 4. 已知断点（必须知道，改之前先看这里）

### 4.1 闸门漏项：`expm_pade`

`src/linalg/expm.hpp:72` 的 `expm_pade` 是**公开入口**，
`demo/linalg/expm.cpp` 里有 **7 处**使用（`[E]` 段 4 项判据），
但 **`smoke.cpp` 的 `EntryTable` 没有它**（`grep -c expm_pade demo/linalg/smoke.cpp` → **0**）。

**后果**：如果将来有人给 `expm_pade` 加模板分支、或它失去全部调用方后被改坏，
**`smoke` 闸门不会报警**。这违反 `smoke.cpp` 自己定的规矩：

> 新增公开入口时**必须**补一行（`EntryTable` + `GateCount` 两处，加 `kEntryCount` 计数），
> 否则闸门静默失效。

**修法**（三处，缺一处无效）：
1. `EntryTable` 加 `static constexpr auto expm_pade_all = &expm_pade<T>;`
2. `GateCount()` 加 `touch(E::expm_pade_all);`
3. `kEntryCount` 从 `40` 改成 `41`，注释「27 自由函数 + 13 成员」改成「28 + 13」

### 4.2 只有 `smoke` 级覆盖的入口（无专项判据）

| 入口 | 现状 | 风险 |
|---|---|---|
| `SwapMultipliers` | 仅 `smoke [B]` | 它是 `LU` 换行时同步乘子的内部件，**错了会静默给出错的 `L`** |
| `jacobi_rotation` | 仅 `smoke [B]` | `eig_symmetric` 的内部件，被间接覆盖（可接受） |
| `diagonal_inverse` | 仅 `smoke [B]` | ⚠️ **它有独立奇异判据**（`|d| <= tol` 抛），但**没有专项用例验证它真的抛** |

**这三条不一定要马上补**（`smoke [B]` 至少保证跑得通且方向对），
但**不能以为它们是"权威覆盖"** —— 它们的阈值是 `1e-10` 的轻量判据，不是算法精度判据。

### 4.3 占位覆盖（`[P]` 段）

| 用例 | 占位项数 | 对应 `src` 的 throw 数 |
|---|---:|---|
| `tensor` | 24 | 51 处 `throw` |
| `quaternion` | 剩余 10 个入口 | 10 处 `throw std::logic_error` |

**占位的语义是「如实钉住不可用」**，不是「测过了」。
**统计覆盖率时这两块必须单列**，不能算进「已实现」。

### 4.4 完全无代码的模块（无 ctest，因为无实现）

`src/math/`（`ode.hpp` / `pde.hpp` / `fourier.hpp`）**目录不存在**。
`src/linalg/operator.hpp` 的 `LinearOperator` 只有 40 行骨架。
`classify.hpp` 不存在（`Classify` 暂住在 `inverse.hpp`）。

---

## 5. 新增一个接口 / 一个用例的检查表

### 5.1 新增公开入口（`src/` 侧）

- [ ] **`numeric-contract.md` §8 九条全过**（形状检查/数值检查/容差参数/`abs2`/实复实例化/`Dual<Complex>` 实例化…）
- [ ] `smoke.cpp` **三处**：`EntryTable` 一行 + `GateCount()` 一行 + `kEntryCount` 计数
- [ ] 在本文 §2 的表里**补一行**
- [ ] 若有专项用例 → 标 **★**；若只有 `smoke` → 标 **○** 并进 §4.2 弱覆盖清单

### 5.2 新增一个 ctest 用例

- [ ] **头部注释先写判据**（`acceptance.md` §1「判据先行」）：
      编号 `①②③…` + **每条的具体阈值** + 「为什么是这个量级」（`numeric-contract.md` §4.2）
- [ ] 覆盖 §3 的六种类型中**能用的全部**，尤其是 **③ 独立不变量**
- [ ] **实/复各实例化一次**（`acceptance.md` §1 第 1 条，此坑已 5 次）
- [ ] **病态输入必列**：`n=0/1`、近奇异、零阵、尺度悬殊
- [ ] **判据函数与批量循环包 `try`**（异常会 terminate，退出码 134）
- [ ] `CMakeLists.txt` 注册 `add_test`
- [ ] 在本文 §1 表里补一行（填入**实测判据数**，不是估算）

### 5.3 判据数怎么填

**不要估算。** 跑一遍数 `[OK]` 行：

```bash
cd graphics && for t in core/scalar_traits linalg/svd ...; do
  echo "$t = $(./build/demo/$t 2>&1 | grep -c '\[OK\]')"
done
```

理由：`Checker` 有两种输出格式（`[OK]` 与 `got=/want=`）。
**用 `grep '检查'` 会漏掉 `scalar_traits` 和 `compute_graph`** —— 本轮实测踩过。

---

## 6. 怎么读失败输出

```bash
./graphics/run_tests.sh          # 全量：只看到哪个 ctest 红
./graphics/run_tests.sh -d       # 明细：看到逐项 [OK]/[FAIL]
./graphics/build/demo/linalg/svd # 单跑某个用例
```

**定位顺序**：

1. 先看 `[FAIL]` 那一行的**名字**（`Checker::check(bool, name)` 的名字是手写的，通常已含数值）
2. 回到该用例头部的「判据」编号，找到**它属于 ①②③ 哪一条**
3. 按那一类的**阈值依据**判断是「实现错」还是「判据病态」
   - ① 结构性（要求精确 0）红 → **几乎一定是实现错**
   - ③ 独立不变量红 → 先查判据的量级（`numeric-contract.md` §4.2）
4. 仍判断不了 → 跑 `run_tests.sh -w` 看 `-Wall -Wextra` 的模板告警

**退出码**：`0` 全过 / `1` 有 FAIL。`run_tests.sh` 用 ctest 聚合，
**单项用例的退出码 134 = SIGABRT = 判据函数没包 `try`**。

---

## 7. 与其它文档的关系

| 问题 | 去哪 |
|---|---|
| 「多少算对」「判据五条」「闸门与注入」 | `acceptance.md` |
| 「除零/容差/错误类型/布局」这些**硬约束** | `numeric-contract.md` |
| 「这个接口被什么测了、判据多少项」 | **本文** |
| 「接下来做什么」 | 根 `NOW.md` |
| 「数学路线上这个件对应哪一段」 | `docs/route-math.md` |
| 「这个件的工程水平算几档」 | `docs/baseline-levels.md` + `docs/eng-competency.md` |

**一句话**：`acceptance.md` 说**怎么验**，`numeric-contract.md` 说**必须满足什么**，
本文说**已经验了什么、还差什么**。
