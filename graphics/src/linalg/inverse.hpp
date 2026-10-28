#pragma once

// ===========================================================================
// linalg/inverse.hpp —— 求逆 (按结构分发) + 结构探测
// ===========================================================================
//
// 【状态】已实现。
//
// ---------------------------------------------------------------------------
// 一、API 契约
// ---------------------------------------------------------------------------
// 结构探测 (全部要求 n×n; 非方阵一律返回 false, 不抛):
//   bool is_diagonal(A, tol)           所有非对角元 |a_ij| <= tol
//   bool is_symmetric(A, tol)          A == Aᵀ  (实域; Hermitian 走 is_hermitian)
//   bool is_orthogonal(A, tol)         Aᵀ·A == I
//   bool is_positive_definite(A, tol)  **Cholesky 能否走通** (见末尾)
//
// 求逆的两级入口:
//   Matrix<T> Inverse(A, MatrixKind kind)   按 kind 分发到对应快路径
//   Matrix<T> Inverse(A)                    先 Classify 探测, 再分发
//
// 五条并列的快路径 (**互不复用实现**):
//   InverseDiagonal / InverseSymmetric / InverseOrthogonal /
//   InverseSymmetricPositiveDefinite / InverseGeneral (LU)
//   外加 qr_inverse (走 A^{-1} = R^{-1}·Q^H) —— 它**不接** Inverse 分发。
//
// 异常:
//   · 非方阵 / 空矩阵            -> throw std::invalid_argument
//   · 奇异 (|pivot| <= tol)      -> throw std::runtime_error
//   · Projection 类 (A²=A≠I)     -> throw std::logic_error (不可逆)
//
// 关键口径: 全库**列向量** (D19); 容差由调用方显式给 (两层除零, D12)。
//
// ---------------------------------------------------------------------------
// 二、调用顺序 / 内部机制
// ---------------------------------------------------------------------------
//   `Inverse(A)` = Classify(A, DefaultTol) -> Inverse(A, kind)。
//   Classify 的**探测顺序有讲究**: 越特殊的档越先判。见末尾「为什么顺序错会
//   让一整个分支永远不可达」。
//   依赖方向: core -> solve / cholesky / qr -> inverse -> eigen (单向, 无循环)。
//
// ---------------------------------------------------------------------------
// 三、数学推导 / 踩坑
// ---------------------------------------------------------------------------
// 【SPD 判据 = "Cholesky 不选主元能否走通"】
//   存在性定理是**双向**的: A 正定 ⟺ Cholesky 走通。所以这不是"近似判断",
//   而是等价判据 —— 这也是"用算法本身当判据"少数正当的场合。
//   好处: 不必算特征值 (O(n³) 且更贵), 复用同一个分解。
//
// 【五条快路径为什么不复用实现】
//   每条路径的**数值性质不同**: 对角阵直接取倒数 (零误差), 正交阵用 A⁻¹ = Aᵀ
//   (零误差), SPD 走 Cholesky (半 flops), 通用走 LU (要选主元)。
//   强行抽公共层会把"零误差"降级成"微小误差", 而结构判据恰恰要靠这种精确性。
//   ⇒ 宁可有五份短代码, 不要一份"通用"的长代码。
//
// 【Classify 的探测顺序 —— 修过的 bug】
//   曾经把 is_symmetric 排在 is_diagonal 前面, 于是对角阵 (也是对称的) 全部
//   被对称分支吃掉, **对角分支永远不可达**。这类错不会报错, 只是某条快路径
//   成了死代码 —— 而它在测试里看起来"也对", 因为对称路径同样能算对。
//   ⇒ 判据: 每加一个 MatrixKind, 必须有**专门**用例证明它被单独命中,
//     不能只看"结果对"。见 demo/linalg/inverse.cpp 的 mid 分支覆盖。
//
// 【为什么 qr_inverse 不接 Inverse 分发】
//   QR 是**算法**, 不是**结构**。MatrixKind 描述的是"矩阵长什么样",
//   QR 描述的是"用什么方法算"。混进分发会逼着枚举里出现"我想用 QR"这种
//   非结构项, 语义就乱了。⇒ 求逆的默认路径是 LU; 需要 QR 就显式调
//   qr_inverse —— 这是调用方的选择, 不是容器的属性。
//
// 【A^H 与 A^T —— 复数域的陷阱】
//   本文件所有"正交/对称"判据在复数域必须用 Adjoint 而不是 Transpose。
//   用错时 double 域完全正常, Complex 域静默给出错误结论。
// ===========================================================================

#include <stdexcept>

#include "core/common.hpp"      // MatrixKind 已定义于此, 直接复用不重复定义
#include "core/matrix.hpp"
#include "core/scalar.hpp"

#include "linalg/solve.hpp"
//   LUResult{L,U,piv} / LU(A) / Solve(A,b) —— 列向量约定 (D19)
#include "linalg/cholesky.hpp"
//   SPD 判据 = "Cholesky 不选主元能否走通" —— 本文件只做分发, 不重复实现分解
#include "linalg/qr.hpp"
//   qr_inverse 走 A^{-1} = R^{-1}·Q^H; qr.hpp 不反向依赖本文件, 故无循环

// ---- 结构探测: 每个判断是一个数学练习 ----
// 比较一律走 abs2, 共轭走 conj, 算法层零实复分支; 容差由调用方给。
template <CommutativeScalar T>
bool is_diagonal(const Matrix<T>& A, typename ScalarTraits<T>::real_t tol)
{
    if (A.Rows() != A.Cols()) return false;
    for (int i = 0; i < A.Rows(); i++) {
        for (int j = 0; j < A.Cols(); j++) {
            if (i == j) continue;
            if (ScalarTraits<T>::abs2(A(i, j)) > tol * tol) return false;
        }
    }
    return true;
}

template <CommutativeScalar T>
bool is_symmetric(const Matrix<T>& A, typename ScalarTraits<T>::real_t tol)
{
    if (A.Rows() != A.Cols()) return false;
    Matrix<T> AH = A.Adjoint();          // AH(i,j) = conj(A(j,i))
    for (int i = 0; i < A.Rows(); i++) {
        for (int j = 0; j < i; j++) {
            // 对称判据 A ≈ A^H 即 A(i,j) ≈ AH(i,j) (= conj(A(j,i)))。
            if (ScalarTraits<T>::abs2(A(i, j) - AH(i, j)) > tol * tol) return false;
        }
    }
    return true;
}
// 酉对称, Hermitian 对称, 实对称是免费自带的, 不用if
// A ≈ A^H: 实数域=对称, 复数域=Hermitian (traits 写法免费升级)

template <CommutativeScalar T>
bool is_orthogonal(const Matrix<T>& A, typename ScalarTraits<T>::real_t tol)
{
    if (A.Rows() != A.Cols()) return false;
    Matrix<T> G = A.Adjoint() * A;       // 正交/酉 <=> A^H·A = I
    for (int i = 0; i < A.Rows(); i++) {
        for (int j = 0; j < A.Cols(); j++) {
            // 只查对角线会漏检非对角元: [[1,x],[x,1]] 也会被判成正交
            const T expect = (i == j) ? ScalarTraits<T>::one() : ScalarTraits<T>::zero();
            if (ScalarTraits<T>::abs2(G(i, j) - expect) > tol * tol) return false;
        }
    }
    return true;
}
// A·A^H ≈ I: 实数域=正交, 复数域=酉 (traits 写法免费升级)

// SPD 判据: "Cholesky 不选主元能否走通" —— 存在性定理 A 正定 <=> 分解走通 (理由见文件头「三」)。
// 只把 cholesky 抛出的 runtime_error (非正定 / 近奇异) 解读为"否"; 别的异常照常外传, 不吞。
template <CommutativeScalar T>
bool is_positive_definite(const Matrix<T>& A)
{
    if (A.Rows() != A.Cols() || A.Rows() == 0)
        return false;   // SPD 蕴含方阵; 非方阵拿去问 cholesky 只会换回一个入参错, 直接答"否"
    try {
        (void)cholesky(A);
        return true;
    } catch (const std::runtime_error&) {
        return false;
    }
}
// 口径: 近奇异(主元实部 <= 8n·eps)也算"否" —— 与 cholesky 同一把尺子, 不放近奇异矩阵过去。

template <CommutativeScalar T>
MatrixKind Classify(const Matrix<T>& A, typename ScalarTraits<T>::real_t tol)
{
    if (is_diagonal(A, tol)) return MatrixKind::Diagonal;
    if (is_orthogonal(A, tol)) return MatrixKind::Orthogonal;
    if (is_symmetric(A, tol)) {
        // SPD ⊂ Symmetric, 必须在这条**对称分支内部**筛子类。平铺成同级 if 的话,
        // 上一行先 return Symmetric 就把它吃掉了, 本档永远不可达。
        if (is_positive_definite(A)) return MatrixKind::SymmetricPositiveDefinite;
        return MatrixKind::Symmetric;
    }
    return MatrixKind::Generic;
}
// 探测顺序: Diagonal -> Orthogonal -> Symmetric{ -> SPD } -> Generic
// Diagonal 排在 SPD 之前是有意的: 对角矩阵即使正定也走 Diagonal(逐元倒数更省)。

// ---- 各分支专用求逆 (五个结构各一条路径; 近奇异/非正定在各自内部 throw, 不留静默 no-op) ----
template <CommutativeScalar T>
Matrix<T> lu_inverse(const Matrix<T>& A)
{
    Matrix<T> I = Matrix<T>::Identity(A.Rows());
    return Solve(A, I);
}
// 奇异判定: 选主元 |pivot|^2 <= tol 则抛


template <CommutativeScalar T>
Matrix<T> diagonal_inverse(const Matrix<T>& A, typename ScalarTraits<T>::real_t tol)
{
    if (A.Rows() != A.Cols()) throw std::invalid_argument("diagonal_inverse: A 必须是方阵");
    Matrix<T> B(A.Rows(), A.Cols());
    for (int i = 0; i < A.Rows(); i++) {
        if (ScalarTraits<T>::abs2(A(i, i)) <= tol * tol)   // 判据 |d| <= tol, 比较走 abs2
            throw std::runtime_error("diagonal_inverse: 对角元近奇异 (判据 |d| <= tol)");
        B(i, i) = ScalarTraits<T>::one() / A(i, i);
    }
    return B;
}
// 对角元倒数, 零/近零抛


// A = L·L^H, A^{-1} = (L^{-1})^H·L^{-1}   
// cholesky_solve(A, I) 展开正是 "L·Y = I 前代 + L^H·X = Y 回代", 与上式同一条路径,
// 于是直接复用那份压过 98 项验收的代码, 不另写一遍三角求逆 (少一份要维护的正确性)。
template <CommutativeScalar T>
Matrix<T> cholesky_inverse(const Matrix<T>& A)
{
    if (A.Rows() != A.Cols() || A.Rows() == 0)
        throw std::invalid_argument("cholesky_inverse: A 必须是非空方阵");
    return cholesky_solve(A, Matrix<T>::Identity(A.Rows()));   // A 非正定在此 throw
}

// A = Q·R (thin; 方阵时 Q、R 都是 n×n), A^{-1} = R^{-1}·Q^H。
// 只做两件**已有**的事: ① Householder QR (qr.hpp) ② 对上三角 R 的一次回代 (R·X = I ⇒ X = R^{-1})。
// 不写成 Solve(R, I): R 本身已经上三角, 再走一遍 LU 等于把 O(n³) 又花一遍;
// 也不显式求 R^{-1} 再去乘 —— 回代逐列解出 X 就是同一件事的紧凑写法 (与 qr.hpp 自带注释同口径)。
// 奇异判据沿用 lu_inverse / diagonal_inverse 的口径: |R_ii| <= tol ⇒ 抛。
// 注意它**不接矩阵分发**: MatrixKind 里没有 QR 档 (也不该有 —— QR 是算法不是结构),
// 所以本函数与 lu_inverse / cholesky_inverse 并列, 是一条"显式点名"的独立路径。
template <CommutativeScalar T>
Matrix<T> qr_inverse(const Matrix<T>& A)
{
    if (A.Rows() != A.Cols() || A.Rows() == 0)
        throw std::invalid_argument("qr_inverse: A 必须是非空方阵");

    const int n = A.Rows();
    const QRcomposition<T> qr = QR_Householder(A);      // m >= n 由方阵保证
    const typename ScalarTraits<T>::real_t tol = DefaultTol<T>(n);

    Matrix<T> X = Matrix<T>::Identity(n);               // 右端项 B = I, 原地被解成 R^{-1}
    for (int i = n - 1; i >= 0; i--) {
        if (ScalarTraits<T>::abs2(qr.R(i, i)) <= tol * tol)
            throw std::runtime_error("qr_inverse: R 对角元近奇异 (秩亏或近奇异, 判据 |R_ii| <= tol)");
        const T inv_diag = ScalarTraits<T>::one() / qr.R(i, i);
        for (int j = 0; j < n; j++) {
            T acc = X(i, j);                            // 此刻第 i 行还是单位阵的 δ_ij, 当右端项用
            for (int k = i + 1; k < n; k++)
                acc = acc - qr.R(i, k) * X(k, j);       // k > i 的行已回代完, 可直接引用
            X(i, j) = acc * inv_diag;
        }
    }
    return X * qr.Q.Adjoint();                          // R^{-1}·Q^H
}

// ---- 对外入口 (按 MatrixKind 枚举分发) ----

// 已知结构: 跳过探测直接走该结构对应的实现 (Orthogonal 直接 A.Adjoint(), A^{-1} = A^H)。
template <CommutativeScalar T>
Matrix<T> Inverse(const Matrix<T>& A, MatrixKind hint)
{
    if (A.Rows() != A.Cols() || A.Rows() == 0)
        throw std::invalid_argument("Inverse: A 必须是非空方阵");

    switch (hint) {
        case MatrixKind::Diagonal:                  return diagonal_inverse(A, DefaultTol<T>(A.Rows()));
        case MatrixKind::Orthogonal:                return A.Adjoint();
        case MatrixKind::Symmetric:
        case MatrixKind::Generic:                   return lu_inverse(A);
        case MatrixKind::SymmetricPositiveDefinite: return cholesky_inverse(A);
        case MatrixKind::RigidTransform:            throw std::logic_error("Inverse: RigidTransform 分支未实现");
        case MatrixKind::Projection:                throw std::logic_error("Inverse: 投影矩阵不可逆 (A^2 = A != I)");
    }
    throw std::logic_error("Inverse: 未覆盖的 MatrixKind");
}

// 通用入口: 先 Classify 探测结构, 再分发。
template <CommutativeScalar T>
Matrix<T> Inverse(const Matrix<T>& A)
{
    return Inverse(A, Classify(A, DefaultTol<T>(A.Rows())));
}
// ===========================================================================
// 四、数学公式与符号约定
// ===========================================================================
//
// **本文件所有公式集中在这里, 实现处只留"为什么这么写"的工程理由。**
//
// ---------------------------------------------------------------------------
// 符号表 (代码里出现的每个符号都在这里)
// ---------------------------------------------------------------------------
//   A              n×n 方阵
//   A⁻¹           A 的逆
//   Aᴴ             A 的共轭转置; 实域退化为 Aᵀ
//   I              n 阶单位阵
//   L, U           LU 分解的两个因子
//   Q, R           QR 分解的两个因子
//   tol            调用方给的容差
//   n              阶数
//   i, j           下标 (代码里是同名循环变量)
//   d_i            A 的第 i 个对角元
//   MatrixKind     结构枚举 (Diagonal / Symmetric / Orthogonal / SPD / …)
//
// ---------------------------------------------------------------------------
// (I1) 逆的定义与基本恒等式
// ---------------------------------------------------------------------------
//        A · A⁻¹ = A⁻¹ · A = I
//
//        两式**都要成立** —— 矩阵不交换, 只测一边漏一半代数。
//        (四元数那边同一个道理, 见 quat/quaternion.hpp。)
//
// ---------------------------------------------------------------------------
// (I2) 各结构下的逆 (五条快路径, 互不复用)
// ---------------------------------------------------------------------------
//        对角:    (A⁻¹)_ii = 1 / d_i,  其余为 0
//        正交:    A⁻¹ = Aᴴ = Aᵀ          (实域; 零误差, 不用任何除法)
//        对称:    A⁻¹ 也是对称的 → 可以只算一半 (实现上仍走 LU)
//        SPD:     A = L·Lᴴ  ⟹  A⁻¹ = (L⁻¹)ᴴ · L⁻¹   (式 C1/C5)
//        通用:    走 LU:  P·A = L·U  ⟹  A⁻¹ = U⁻¹·L⁻¹·P
//
//      **为什么不抽公共层**: 每条路径的数值性质不同 ——
//        对角取倒数零误差; 正交用 Aᵀ 零误差; SPD 走 Cholesky 半 flops。
//        强行统一会把"零误差"降级成"微小误差", 而结构判据恰恰靠这种精确性。
//
// ---------------------------------------------------------------------------
// (I3) 结构探测的判据 (每个都是一次数学练习)
// ---------------------------------------------------------------------------
//        对角:    |a_ij| <= tol  对所有 i ≠ j
//        对称:    |a_ij − conj(a_ji)| <= tol  对所有 i,j
//        正交:    ‖Aᴴ·A − I‖ <= tol
//        SPD:     Cholesky 不选主元能走通 (等价判据, 见式 C3)
//
//      ⚠ 复数域下"正交"必须用 **Aᴴ** 不是 Aᵀ —— 写错在 double 下完全正常,
//        只有 Complex 域静默给出错误结论。
//
// ---------------------------------------------------------------------------
// (I4) 探测顺序 (修过的 bug)
// ---------------------------------------------------------------------------
//        越特殊的档越先判, 且子类要判在父类**之前**:
//              Diagonal  ⊂  Symmetric  ⊂  General
//              SPD       ⊂  Symmetric
//
//        曾经把 is_symmetric 排在 is_diagonal 前面 ⇒ 对角阵全被对称分支吃掉,
//        **对角分支永远不可达**。这类错不报错, 只是某条快路径成了死代码 ——
//        而它在测试里看起来"也对", 因为对称路径同样能算对。
//
//      ⇒ 判据: 每加一个 MatrixKind, 必须有**专门**用例证明它被单独命中,
//        不能只看"结果对"。
//
// ---------------------------------------------------------------------------
// (I5) 奇异判据 (各条路径统一口径)
// ---------------------------------------------------------------------------
//        |pivot| <= tol  ⇒  throw
//
//        对角路径: |d_i| <= tol
//        LU 路径:  |u_kk| <= tol
//        QR 路径:  |R_ii| <= tol   (且要求 R 是对角占优才有保证)
//        SPD 路径: cholesky 内部自己判 (主元实部 <= tol)
//
//        比较走 abs2 (D7): abs2(pivot) <= tol²。
//
// ---------------------------------------------------------------------------
// (I6) 为什么 qr_inverse 不接 Inverse 分发
// ---------------------------------------------------------------------------
//        **QR 是算法, 不是结构**。MatrixKind 描述"矩阵长什么样",
//        QR 描述"用什么方法算"。混进分发会逼着枚举里出现"我想用 QR"
//        这种非结构项, 语义就乱了。
//
//      ⇒ 求逆的默认路径是 LU; 需要 QR 就显式调 qr_inverse ——
//        这是**调用方的选择**, 不是容器的属性。
//
// ---------------------------------------------------------------------------
// 实现口径 (跟公式直接相关)
// ---------------------------------------------------------------------------
//   · **列向量约定 (D19)**: 全库列向量, 逆作用于列向量右侧。
//   · **容差由调用方显式给** (两个入口: 带 tol 的与不带 tol 的)。
//     不写在 Matrix 里 —— 容差是算法的参数, 不是容器的属性。
//   · 投影矩阵 (A² = A ≠ I) **不可逆**, 明确 throw 而不是返回垃圾。
// ===========================================================================
