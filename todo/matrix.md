# matrix.md — Matrix 库设计思路总结

> 本文汇总自 design 讨论，作为 graphics/src 矩阵库的设计备忘。
> 数学符号一律用 ASCII 写法（A^T 转置、A^H 共轭转置、[A,B]=AB-BA 对易子），便于终端阅读。

---

## 0. 一句话宪法

**容器哑、语法甜（[][] 派生不存储）、运算出（自由函数分文件）、
差异下沉 traits（算法层零分支）、算子走表示不走进元素、高维靠展平不造张量。**

---

## 1. 标量层：Scalar 概念 + ScalarTraits（地基之一）

### 1.1 设计动机

矩阵要同时支持 float / double / Complex<T> / Dual<T>，甚至组合类型
Dual<Complex<double>>（复系数 + 自动求导，服务波包演化）。
如果矩阵算法里到处 `if (是复数)`，代码会烂掉。解法：**类型要求 + 语义特化**。

### 1.2 两层机制

```cpp
// 主模板不定义 = 没特化的类型天然不满足约束
template <typename T> struct ScalarTraits;

template <typename T>
concept Scalar = requires (T a, T b) {
    { a + b } -> std::same_as<T>;        // -, *, / 同理
    { ScalarTraits<T>::zero() } -> std::same_as<T>;
    { ScalarTraits<T>::one()  } -> std::same_as<T>;
    { ScalarTraits<T>::conj(a) } -> std::same_as<T>;          // 共轭: 实数=恒等映射
    { ScalarTraits<T>::abs2(a) } -> std::same_as<typename ScalarTraits<T>::real_t>;
    // real_t 永远是 float/double: 一切"比较"只在实数域发生
};
```

三组特化（以后加新标量 = 只加一个特化，不动任何算法）：

| 特化 | conj | abs2 | 备注 |
|---|---|---|---|
| float / double | x | x*x | trivial |
| Complex<T> | 共轭 | 模方(不开方, 数值友好) | real_t = T |
| Dual<T> | 转发 T 的 conj | value 的 abs2 | 比较/选主元只看 value |

### 1.3 铁律

- 算法层（linalg/、math/）里**只允许**出现 `ScalarTraits<T>::abs2/conj/real_t`，
  禁止 `x*x`、禁止 `<` 直接比元素、禁止 `if constexpr (是复数)`。
- `if constexpr` 是逃生门，仅限数学上真正不同的场合（如复数 Householder 取相位），
  需注释说明理由。
- **Scalar 必须再加一个 CommutativeScalar 要求**（乘法交换，单测验证 a*b == b*a）。
  原因见 §5：LU/QR/Cholesky 的消元步 `m = a/p; x -= m*y` 隐式假设左右乘无差别。
  float/Complex/Dual 都交换；**算子不交换 → 从类型层面禁止算子做矩阵元素**。

---

## 2. Matrix 容器本体：该有什么、不该有什么

### 2.1 判定标准

"换一个存储布局（稀疏、GPU）后这个操作还成立吗？"
成立 → 不是容器的事。逆矩阵换算法还是逆矩阵，但换存储就不是 `T* operator[]` 了。

### 2.2 白名单（容器只保留这些）

- 存储：`std::vector<T> data`（行主序连续）+ rows/cols
- 访问：`at(i,j)`（带边界检查）+ `operator()` + `data()` + shape
- 构造：`(m,n)`、initializer_list（带校验：空表/行长一致/values.size()）
- 静态工厂：`Zero / Identity / Ones / Diagonal`
  （取代有歧义的 `Matrix(int m, T initVal)`——它与 `Matrix(int,int)` 在 float 下混淆）
- 逐元素四则：`+ - 标量乘`，逐元素乘改名 `cwiseProduct`
  （`operator*` 保留给矩阵乘，对齐业界惯例——注释与实现必须统一口径）
- 结构操作：`Transpose()`（纯转置）与 `Adjoint()`（共轭转置 = 复数世界的转置，
  Hermitian/酉矩阵、SVD、QR 全依赖它；对实数 traits 自动退化为 Transpose）

### 2.3 黑名单（明确移出容器的东西）

- `std::vector<T*>` 行指针层 —— **最大风险源，必须删**：
  1. 制造"matrix[i] == &data[i*cols]"这条只有人记得住的不变量；
     任何 resize/重新赋值 data 忘了重建就悬空（拷贝/移动构造都要手工重建，
     review 里已发现 moved-from 悬空指针 bug）。
  2. `operator[]` 返回的 T* 交出去后生命周期失控，`A[0][999]` 是无声 UB。
  3. 性能也是负资产：多一次指针间接（额外缓存行），不如 `i*cols+j` 一次地址计算，
     且扁平存储向量化友好。
  - **`A[i][j]` 语法不损失**：operator[] 改为即时派生
    `T* operator[](int i) { return data.data() + i*cols; }`
    ——指针每次从当前状态计算，永不过期。
- mutable 缓存（cachedTranspose/cachedInverse/decomposition）—— 失效逻辑守不住，
  整个删除；真需要缓存是上层事。
- 一切矩阵乘法以外的算法：multiply 移到 linalg；Det、求逆、分解全不进类。
- `SumLine / MultLine / MultDiagnal`：它们根本不是矩阵语义——"对一行累加/累乘"
  是 Dual 前向模式的传播规则，和"矩阵的行"只是存储巧合；
  且正是它们造成 friend 模板声明不匹配的编译级 P0 bug。
  → 迁到 `autodiff/jacobian.hpp`。

### 2.4 其它实现要点（原 review P0/P1 收口）

- `#pragma once`；析构在类外定义（unique_ptr 成员依赖完全类型，别挪回类内）。
- 移动构造初始化顺序与声明顺序一致（-Wreorder）；moved-from 清空。
- const 正确性全量清理；比较函数加 const。
- 容差用 `std::numeric_limits<real_t>::epsilon()` 的倍数，不用全局 EPSILON。
- 尺寸策略：本体保持动态 Matrix<T>；图形学热路径的 Vec<T,2/3/4> 独立小类；
  Matrix<T,M,N> 固定尺寸特化推迟到渲染出现真实性能瓶颈（学习目标是数学，不是性能）。
- 行主序自定，将来若上 OpenGL 记得上传时处理列主序。

---

## 3. 运算层：自由函数按能力分文件

### 3.1 为什么不用"抽象接口类"

"每次新运算加一个接口类"的正确落法是**头文件 + 自由函数模板**，不是 Java 式
`interface Invertible { ... }`：
- 代码全是 template<T>，虚函数与模板不能共存，接口类会退化成每个标量类型各继承一份；
- 运算输入输出都是值（Matrix 进、Matrix/结构体出），没有对象状态，基类是多余仪式；
- 插件式扩展同样满足：新算法 = 新头文件 + 新函数，不动老代码。

```cpp
// linalg/solve.hpp
template <typename T> struct LUResult { Matrix<T> L, U; std::vector<int> piv; };
template <CommutativeScalar T> LUResult<T> LU(const Matrix<T>& A);
template <CommutativeScalar T> Vec<T>     Solve(const Matrix<T>& A, const Vec<T>& b);

// linalg/inverse.hpp  —— 阶段 1 只此一种: LU + 逐列回代
template <CommutativeScalar T> Matrix<T> Inverse(const Matrix<T>& A);

// linalg/eigen.hpp    —— Jacobi 旋转法, 仅对称/Hermitian
// linalg/svd.hpp      —— 经 A^H·A 的 Jacobi 特征分解组装
```

### 3.2 求逆的最终形态（阶段 2+）

矩阵头文件里手写过的那张求逆方法表，可执行化为**枚举分发**，不需要多态类：

```cpp
enum class MatrixKind { Generic, Diagonal, Orthogonal/*酉*/, RigidTransform,
                        SymmetricPositiveDefinite/*Hermitian正定*/, Projection };
MatrixKind Classify(const Matrix<T>& A, real_t tol);   // 结构探测
Matrix<T>  Inverse(const Matrix<T>& A);                // Classify -> switch
Matrix<T>  Inverse(const Matrix<T>& A, MatrixKind hint);  // 已知结构跳过探测
```

学习价值：每个 Classify 判断本身就是一次数学练习——
正交检测 `A·A^H ≈ I` 用 traits 写后，对复数矩阵自动变成**酉矩阵检测**，免费；
SPD 检测 = Cholesky 不选主元能否走通。

### 3.3 阶段 1 方法选型原则

每种运算只实现一种"通用且简单"的方法：
- 求逆：LU 部分选主元（一切通用法的地基）
- 特征值：Jacobi（对称/Hermitian 情形简单稳定；非对称的 Hessenberg+QR 难一个量级，缓做）
- SVD：先经 A^H A 特征分解（够用，且顺便复习对称理论）

知道边界就好：对称/Hermitian 恰好覆盖 PCA、惯性张量、SVD、主方向、
以及量子力学的可观测量——价值密度最高的子集。

---

## 4. 硬验收标准：零分支双类型编译

**检验 traits 设计是否正确的唯一硬标准：**

> 同一份 linalg 源码，`Matrix<double>` 与 `Matrix<Complex<double>>`
> 零修改实例化、编译通过、结果正确（残差 ‖A·A⁻¹−I‖<ε、‖A−QΛQ^H‖<ε、‖UΣV^H−A‖<ε）。

编译不过 = 哪里漏了 traits（通常是某处直接比较元素或漏了 conj）。
这条同时验收了 Dual<Complex<double>> 的可行性（A5 验收项），
为薛定谔方程波包演化（复系数矩阵 + 自动求导）铺路。

---

## 5. 函数矩阵：两种对象，一个都不能含糊

### 5.1 类型 1：元素是函数（逐点语义，交换）

A(x) = [[sin x, e^x], [x^2, 1]]，乘法逐点、可交换——理论上就是
"交换环上的矩阵"，概念上兼容 Scalar 要求。但注意两个坑：
- 环 ≠ 域：函数有零点 → `x` 逐点不可除，LU 除法步骤产生极点；
  符号路线要进有理函数域，数值路线直接**采样**成网格值（每格点一个小矩阵）。
- 数值路线即 Field/Grid（plan.md 阶段 4），不是 Matrix<function>。

### 5.2 类型 2：元素是算子（乘法 = 复合，非交换）—— 微分算子在这里

D = d/dx 的"矩阵"是它在某组基下的**表示**，不是把 d/dx 塞进某个格子。
三个不可绕过的数学事实：

1. **复合不可交换**：[D, X] = DX − XD = I（乘积法则的直接改写）。
   所有消元类算法假设交换律 → 被 CommutativeScalar 约束挡在 Matrix 元素之外。
2. **无穷维 + 无界**：多项式基下 D 是次对角移位阵（x^n → n·x^{n−1}）；
   截断到 N 有截断误差/混叠；-d²/dx² 特征值 ~N²、条件数 ~N⁴，
   **容差必须随离散尺度缩放**，否则会把良定矩阵误判为奇异。
3. **"算子的矩阵"离开基和边界条件没有意义**：
   -d²/dx² 在 Dirichlet 正弦基下对角 diag(n²)、周期复指数基下对角、
   多项式基下稠密——同一算子三副面孔，**基是类型的一部分，不能靠口头约定**。

### 5.3 正确设计：LinearOperator（表示而非元素）

```cpp
// linalg/operator.hpp
template <typename T>
struct LinearOperator {
    Matrix<T> repr;        // 截断后的普通数矩阵 —— 有限维表示
    BasisKind basis;       // Monomial / SineTrig / ComplexExp / Chebyshev...
    int truncation;        // 基维度; 一切精度/收敛声明都挂它
    Vec<T,N> Apply(v)      const { return repr * v; }
    // 复合 = 有序矩阵乘: Compose(A,B) 与 Compose(B,A) 是两个不同结果
    // 伴随 = Adjoint(repr), 且这是"内积意义的伴随", 基与权函数必须记录
};
// 工厂: DifferentiationMatrix(basis, N) / LaplacianMatrix(bc, basis, N)
//       / MultiplicationMatrix(func, basis, N)
```

### 5.4 三分工（谁也不冒充谁）

| 需求 | 归属 | 说明 |
|---|---|---|
| 函数逐点运算 | Field/Grid（阶段 4） | 采样表示 + 每点切片 |
| 微分算子的计算 | LinearOperator | 基+BC+截断矩阵 |
| 符号表达式微分 | autodiff/表达式树（lexer.hpp 的位置） | 与 linalg 完全无关 |

配套 demo（与 learning.md 第二阶段 2.2、验收 B7 一一对应）：
- D 在多项式基下的矩阵 = 移位阵；expm(t·D) 作用 = Taylor 移位（与解析对拍）
- Dirichlet 正弦基下 Laplacian = diag(−n²)
- 对易子 Compose(D,X) − Compose(X,D) → I 的残差随 N 收敛到 0（记录收敛率）
  ——离散只能近似满足交换关系，这个 demo 本身就是"截断误差"的可视化

---

## 6. 矩阵只有二维吗？（高维数据的口径）

数学上"矩阵"严格指二维（线性映射的表出）。但应用中的"高维"几乎全能
**展平/切片归约为二维线性代数**——这不是妥协，这正是线性代数存在的原因
（有限维多重线性代数 ≅ 矩阵代数）：

| 场景 | 真实对象 | 处理 |
|---|---|---|
| f: R^n → R^m 的雅可比 | 二维 m×n | Matrix 直接胜任 |
| 向量场的 Hessian | 三维 | 按输出分量切成 m 个 n×n |
| 弹性张量 σ=C:ε | 四维 | Voigt 记法 → 6×6 矩阵 |
| 网格场 u(x,y)、图像 RGB | 数组多维 | 向量化：N 点 → 长度 N 向量；离散算子 → 大矩阵 |
| 分块（刚体求逆 [[R,t],[0,1]]） | 二维的视图 | Block 视图，不需要第三维 |

设计口径：Matrix（算子）/ Field、Grid（数据场，shape + flat buffer +
展平成列向量的视图，复用全部 linalg）/ Tensor<T,rank>（**不承诺，推迟**，
仅当出现展平处理不了的场景：GR 指标运算、einsum、连续介质细节）。

---

## 7. 与 Dual/autodiff 的接口备忘

- Dual<T> 自身也要满足 Scalar 语义后 Matrix<Dual> 才可用：
  `d.value==0` 的判零改走 `isZero`（经 traits）；除零分支返回 (0,0) 的策略需重审
  （前向模式传播 NaN 还是显式抛，取决于用途）。
- 求导函数签名 `std::function<T(T&)>` 改为 `T(const T&)`。
- sumLine/multLine 这类"一行 Dual 的传播规则"住 `autodiff/jacobian.hpp`，
  与矩阵的"行"再无关系。
- 矩阵函数 f(A)（阶段 2 expm 起）与"元素是 Dual 的矩阵"是两回事，
  前者是 Schur/特征分解的函数演算，后者只是 traits 的一次实例化——别混。

---

## 8. 决策一览（本主题全部拍板项）

| # | 决策 | 状态 |
|---|---|---|
| 1 | C++20 + concepts 定义标量约束 | 定 |
| 2 | 自实现 Complex<T>，不用 std::complex | 定 |
| 3 | Matrix 哑容器：扁平存储 + 派生行指针 + 逐元素 + Transpose/Adjoint | 定 |
| 4 | 运算 = 自由函数按能力分文件；不用抽象接口类 | 定 |
| 5 | 求逆/特征值/SVD 阶段 1 各只一种通用法（LU / Jacobi / A^H A 路线） | 定 |
| 6 | Classify 枚举分发做快速求逆 → 阶段 2 | 定 |
| 7 | 实复差异只在 traits 特化，算法层零分支；硬验收 = 双类型零修改编译 | 定 |
| 8 | Scalar 附加交换律要求；算子不做元素，走 LinearOperator | 定（本次新增） |
| 9 | 动态尺寸为主 + 独立 Vec<T,N>；Matrix<T,M,N> 与 Tensor 推迟 | 定 |
| 10 | 高维数据用展平/Field-Grid，不造张量类 | 定 |

---

## 9. 阶段 1 落地顺序（可直接派给 agent）

1. CMake 升 C++20；建 `core/scalar.hpp`（概念 + float/double traits + 交换律单测）
2. `complex/complex.hpp`：四则/模/共轭/exp（欧拉公式），特化 traits
3. 重写 `core/matrix.hpp`：§2 白名单 + P0 清理（对照原 review 逐条销账）
4. 修 `dual.hpp`：isZero 走 traits；`SumLine/MultLine/MultDiagnal` 迁 `jacobian.hpp`
5. `linalg/solve.hpp`(LU+Solve) → `inverse.hpp` → `eigen.hpp`(Jacobi) → `svd.hpp`
6. 验收 A1–A7（见 plan.md 阶段 1）；随后 LinearOperator 随阶段 2 的 B7 落地
