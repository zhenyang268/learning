// string_modes.cpp -- 弦的模态: 离散 Laplacian 的特征值 -> (k*pi)^2
//
// demo/apps 的第一个成员, 同时是 fan-plan 扇面 A1 切片的 C++ 版
// (Python 版按 fan-plan 另走 viz_modules; 这里是用已验收的 eig_symmetric 直接看现象)。
//
// 连续问题:   -u'' = lambda*u,  u(0) = u(1) = 0   ->   lambda_k = (k*pi)^2,  u_k = sin(k*pi*x)
// 离散化:     x_j = j*h,  h = 1/(N+1),  A = tridiag(-1, 2, -1)/h^2
// 矩阵特征值: lambda_k^h = 4*sin^2(k*pi/(2(N+1)))/h^2   ->   (k*pi)^2   (N -> inf)
//
// 现象 (不是判据): 比值 lambda_k^h/(k*pi)^2 随 N 增大逐列趋 1。
// 这就是 "矩阵 = 微分算子" 的第一次实锤 —— 矩阵的特征值序列, 认出了微分算子的谱。
//
// 附加区程序: 不注册 ctest (没有阈值断言, 只有现象与自检打印)。
// 运行: ./graphics/run_apps.sh string_modes

#include "../common.hpp"
#include "linalg/eigen.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using D = double;

// A = tridiag(-1, 2, -1) / h^2 —— Dirichlet 边界下一维 Laplacian 的中心差分
static Matrix<D> DiscreteLaplacian(int N, D h)
{
    Matrix<D> A = Matrix<D>::Zero(N, N);
    const D k = D(1) / (h * h);
    for (int i = 0; i < N; i++) {
        A(i, i) = D(2) * k;
        if (i > 0)     A(i, i - 1) = -k;
        if (i + 1 < N) A(i, i + 1) = -k;
    }
    return A;
}

int main()
{
    const D PI = std::acos(D(-1));
    const int Ns[] = {8, 16, 32, 64, 128};

    std::printf("离散 Laplacian A = tridiag(-1, 2, -1)/h^2,  h = 1/(N+1)\n");
    std::printf("比值 lambda_k^h / (k*pi)^2  (连续极限应为 1)\n\n");
    std::printf("%5s  %8s", "N", "h");
    for (int k = 1; k <= 5; k++) std::printf("  k=%d      ", k);
    std::printf("\n");
    std::printf("%s\n", std::string(66, '-').c_str());

    for (int N : Ns) {
        const D h = D(1) / D(N + 1);
        const Matrix<D> A = DiscreteLaplacian(N, h);

        const Eigendecomposition<D> e = eig_symmetric(A);

        // Jacobi 的输出顺序不保证升序 (按对角收敛顺序), 这里显式排序
        std::vector<D> lam;
        lam.reserve(N);
        for (int i = 0; i < N; i++) lam.push_back(e.Lamda(i, i));
        std::sort(lam.begin(), lam.end());

        std::printf("%5d  %8.5f", N, h);
        for (int k = 1; k <= 5; k++) {
            const D exact = D(k) * PI * D(k) * PI;
            std::printf("  %.6f", lam[k - 1] / exact);
        }
        std::printf("\n");
    }

    // ---- 一个模态的形状: 第一特征向量应该认出 sin(pi*x) ----
    const int N = 20;
    const D h = D(1) / D(N + 1);
    const Matrix<D> A = DiscreteLaplacian(N, h);
    const Eigendecomposition<D> e = eig_symmetric(A);

    // 第 1 模态 = 最小特征值所在的列
    int argmin = 0;
    for (int i = 1; i < N; i++)
        if (e.Lamda(i, i) < e.Lamda(argmin, argmin)) argmin = i;

    // 归一化到峰值为 1 (特征向量只在差一个常数倍的意义下确定)
    D peak = 0;
    for (int i = 0; i < N; i++) peak = std::max(peak, std::abs(e.Q(i, argmin)));
    std::vector<D> v(N);
    for (int i = 0; i < N; i++) v[i] = e.Q(i, argmin) / peak;

    std::printf("\n第 1 模态 (N=%d, 特征值 %.6f  vs  pi^2 = %.6f)\n", N,
                e.Lamda(argmin, argmin), PI * PI);
    std::printf("  x        数值 u_1        解析 sin(pi*x)     残差\n");
    D max_res = 0;
    for (int i = 0; i < N; i++) {
        const D x = D(i + 1) * h;
        const D exact = std::sin(PI * x);
        const D res = std::abs(v[i] - exact);
        max_res = std::max(max_res, res);
        if ((i + 1) % 2 == 0) {
            std::printf("  %.4f   %+.8f   %+.8f     %.2e\n", x, v[i], exact, res);
        }
    }
    std::printf("  最大残差 = %.3e  (中心差分在网格上的截断误差, 不是求解误差)\n", max_res);

    // ---- 自检: 结果不是抄出来的, 用两条独立判据现场验 ----
    const D recon = double(NormInf(e.Q * e.Lamda * e.Q.Adjoint() - A));
    D pair_max = 0;
    for (int j = 0; j < N; j++) {
        std::vector<D> col(N);
        for (int i = 0; i < N; i++) col[i] = e.Q(i, j);
        const Matrix<D> vj(N, 1, col);
        const D lam = e.Lamda(j, j);
        for (int i = 0; i < N; i++) col[i] *= lam;
        pair_max = std::max(pair_max, double(NormInf(A * vj - Matrix<D>(N, 1, col))));
    }
    std::printf("\n自检 (不是判据, 只是结果可信的证据):\n");
    std::printf("  ‖A - Q*L*Q^T‖inf = %.3e\n", recon);
    std::printf("  max_j ‖A*v_j - l_j*v_j‖inf = %.3e\n", pair_max);

    return 0;
}
