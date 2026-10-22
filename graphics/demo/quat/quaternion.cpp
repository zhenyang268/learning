// quaternion.cpp -- B1 四元数验收骨架 (plan.md 阶段 2 第 5 条; V2.5 / V2.6)
//
// ---------------------------------------------------------------------------
// 本文件现在是**骨架**, 不是验收
// ---------------------------------------------------------------------------
// 已实现的部分 ([A] / [B]) 是真判据, 定的是"多少算对";
// 待实现的部分 ([P]) 是「**断言它仍然 throw**」的占位判据 —— 这正是 smoke.cpp
// 对 qr_inverse 用的那条规矩: **不假装它可用**, 而是如实钉住"它现在不可用"。
//
// 为什么这样写而不是留空: 你实现掉一条, 对应的 [P] 就会**变红** —— 红就是提醒,
// 提醒你把占位判据换成真判据 (阈值已写在 plan.md V2.5 与 pyapps/rotation/report.md)。
// 留空的话, 实现了也没人知道, 没实现也没人知道。
//
// ---------------------------------------------------------------------------
// 四步 (与 dev-progress.md「下一版基线」B1 对齐)
// ---------------------------------------------------------------------------
//   ① 基本代数     + - 一元- 标量乘 ✔  |  Hamilton 积 operator* ⏳
//   ② 共轭/范数    Conjugate / SquaredNorm / Norm ✔  |  Inverse / Normalize ⏳
//   ③ 互转 + 旋转  FromAxisAngle / ToMatrix / FromMatrix / Rotate ⏳
//   ④ 插值         Slerp / Nlerp / AngleOf ⏳
//
// ---------------------------------------------------------------------------
// 分域纪律
// ---------------------------------------------------------------------------
// [A]/[B] 对 T = double 与 T = float **各跑一遍** (graphics 3.3 第 1 条:
// 模板 bug 只在实例化时暴露; 单域探针连续误判过 5 次, 正反两个方向都栽过)。
//
// 运行:
//   ./graphics/run_tests.sh                 # 全量 (本用例注册名: quaternion)
//   ./graphics/build/demo/quat/quaternion   # 单跑
//   ./graphics/run_tests.sh -w              # 顺带开 -Wall -Wextra

#include "../common.hpp"
#include "quat/quaternion.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

// ---- 小工具 ---------------------------------------------------------------

template <typename T>
std::string FmtQ(const Quaternion<T>& q)
{
    char buf[160];
    std::snprintf(buf, sizeof(buf), "(%.9g, %.9g, %.9g, %.9g)",
                  static_cast<double>(q.w), static_cast<double>(q.x),
                  static_cast<double>(q.y), static_cast<double>(q.z));
    return std::string(buf);
}

// 逐分量比较。**故意不写 operator==**: 四元数的"相等"有两个意思 ——
// 逐分量相等 (代数) vs |q0.q1| = 1 (同一旋转, 因为 M(q) = M(-q))。
// 把 == 定义成其中任何一个都会误导另一半。比较留在用例里, 写清是哪一种。
template <typename T>
bool SameQ(const Quaternion<T>& a, const Quaternion<T>& b, T tol)
{
    return std::fabs(a.w - b.w) <= tol && std::fabs(a.x - b.x) <= tol
        && std::fabs(a.y - b.y) <= tol && std::fabs(a.z - b.z) <= tol;
}

template <typename F>
bool ThrowsLogicError(F f)
{
    try {
        f();
    } catch (const std::logic_error&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

// ---- [A] ① 已实现: 构造 / 加减 / 一元负 / 标量乘 --------------------------

template <typename T>
void TestAlgebra(Checker& ck, const char* tag)
{
    using Q = Quaternion<T>;
    const T tol = T(1e-6);

    // 默认构造 = 单位元 —— 与 Dual 的默认 (0,0) 故意不同:
    // (0,0,0,0) 没有逆、不是旋转, 默认建造出它是个陷阱。
    const Q q0;
    ck.check(q0.w == T(1) && q0.x == T(0) && q0.y == T(0) && q0.z == T(0),
             std::string("[A1] ") + tag + " 默认构造 = 单位元 (1,0,0,0)  实测 "
                 + FmtQ(q0));

    ck.check(SameQ(q0, Q::Identity(), T(0)),
             std::string("[A2] ") + tag + " Identity() 与默认构造一致");

    const Q a(T(1), T(2), T(3), T(4));
    ck.check(a.w == T(1) && a.x == T(2) && a.y == T(3) && a.z == T(4),
             std::string("[A3] ") + tag + " 四参构造逐位就位  实测 " + FmtQ(a));

    const Q b(T(0.5), T(-1.5), T(2), T(-0.25));
    ck.check(SameQ(a + b, Q(T(1.5), T(0.5), T(5), T(3.75)), tol),
             std::string("[A4] ") + tag + " 加法逐位  a+b = " + FmtQ(a + b));

    ck.check(SameQ(a - b, Q(T(0.5), T(3.5), T(1), T(4.25)), tol),
             std::string("[A5] ") + tag + " 减法逐位  a-b = " + FmtQ(a - b));

    ck.check(SameQ(a + (-a), Q(T(0), T(0), T(0), T(0)), T(0)),
             std::string("[A6] ") + tag + " 一元负: a + (-a) = 0 (逐位精确)");

    // 标量乘两侧都要: q*2 与 2*q 必须一致
    ck.check(SameQ(a * T(2), T(2) * a, T(0)) && (a * T(2)).x == T(4),
             std::string("[A7] ") + tag + " 标量乘 q*2 == 2*q, 且逐位翻倍  "
                 + FmtQ(a * T(2)));

    // 分配律 —— 不是凑数: Hamilton 积实现之后这条仍须成立 (标量乘与乘法相容)
    ck.check(SameQ((a + b) * T(2), a * T(2) + b * T(2), tol),
             std::string("[A8] ") + tag + " 标量乘对加法分配");
}

// ---- [B] ② 已实现: 共轭 / 范数 -------------------------------------------

template <typename T>
void TestConjugateNorm(Checker& ck, const char* tag)
{
    using Q = Quaternion<T>;
    const Q a(T(1), T(2), T(3), T(4));

    const Q c = a.Conjugate();
    ck.check(c.w == T(1) && c.x == T(-2) && c.y == T(-3) && c.z == T(-4),
             std::string("[B1] ") + tag + " 共轭只翻向量部分  q* = " + FmtQ(c));

    // 1 + 4 + 9 + 16 = 30, 整数运算在 float/double 下都精确
    ck.check(a.SquaredNorm() == T(30),
             std::string("[B2] ") + tag + " |q|^2 = 30 (逐位精确)");

    ck.check(std::fabs(a.Norm() - T(std::sqrt(30.0))) <= T(1e-6),
             std::string("[B3] ") + tag + " |q| = sqrt(30)");
}

// ---- [C] ① 的魂: Hamilton 积 (2026-10-19 已实现, 原 [P1] 占位判据升级到这里) ----
//
// 为什么这一节的判据要这么密: 叉积**符号写反**, 往返测试 (q * q^-1) 和姿态往返
// 都**抓不到** (正反抵消), 只有
//   ① 基元乘法表 (i*j=+k 等 9 条)  —— 钉死 ij=k 还是 ji=k
//   ② 共轭反自同构 (ab)* = b* a*   —— 钉死乘法**顺序**
//   ③ 非交换 [i,k] ≠ 0             —— 钉死"确实不是复数"
// 三者合起来才能唯一确定 Hamilton 约定。

template <typename T>
void TestHamilton(Checker& ck, const char* tag)
{
    using Q = Quaternion<T>;
    const T tol = T(1e-6);

    const Q i(T(0), T(1), T(0), T(0));
    const Q j(T(0), T(0), T(1), T(0));
    const Q k(T(0), T(0), T(0), T(1));
    const Q one(T(1), T(0), T(0), T(0));

    auto chkMul = [&](const Q& a, const Q& b, const Q& want, const std::string& name) {
        ck.check(SameQ(a * b, want, tol), name + "  实测 " + FmtQ(a * b));
    };

    // 乘法表 (Hamilton 约定 ij = k, ji = -k)
    chkMul(i, j,  k,  std::string("[C1] ") + tag + " i*j = +k");
    chkMul(j, k,  i,  std::string("[C2] ") + tag + " j*k = +i");
    chkMul(k, i,  j,  std::string("[C3] ") + tag + " k*i = +j");
    chkMul(j, i,  -k, std::string("[C4] ") + tag + " j*i = -k  (反向)");
    chkMul(i, k,  -j, std::string("[C5] ") + tag + " i*k = -j  (反向)");

    // 非交换: ik = -j 而 ki = +j, 差 2j —— 量级 O(1), 不是舍入残渣
    ck.check(SameQ(i * k + k * i, Q(T(0), T(0), T(0), T(0)), tol),
             std::string("[C6] ") + tag + " 不同基元反交换: ik + ki = 0");
    ck.check((i * k - k * i).SquaredNorm() == T(4),
             std::string("[C7] ") + tag + " 非交换量 |ik - ki|^2 = 4 (O(1), 不是舍入)");

    // 单位元
    const Q a(T(1), T(2), T(3), T(4));
    ck.check(SameQ(a * one, a, T(0)) && SameQ(one * a, a, T(0)),
             std::string("[C8] ") + tag + " 单位元: a*1 = 1*a = a (逐位精确)");

    // 结合律 / 分配律
    const Q b(T(0.5), T(-1.5), T(2), T(-0.25));
    const Q c(T(-2), T(0.75), T(1), T(3));
    ck.check(SameQ((a * b) * c, a * (b * c), tol * T(10)),
             std::string("[C9] ") + tag + " 结合律 (ab)c = a(bc)");
    ck.check(SameQ(a * (b + c), a * b + a * c, tol * T(10)),
             std::string("[C10] ") + tag + " 分配律 a(b+c) = ab + ac");

    // 共轭是**反**自同构: (ab)* = b* a* —— 顺序反了就不成立, 专抓乘法写反
    ck.check(SameQ((a * b).Conjugate(), b.Conjugate() * a.Conjugate(), tol * T(10)),
             std::string("[C11] ") + tag + " 共轭反自同构 (ab)* = b* a*");

    // 范数乘性 (Lagrange 恒等式 |ab| = |a||b|): 判据独立于四个分量各自的值
    ck.check(std::fabs((a * b).SquaredNorm() - a.SquaredNorm() * b.SquaredNorm())
                 <= tol * a.SquaredNorm() * b.SquaredNorm() * T(10),
             std::string("[C12] ") + tag + " |ab|^2 = |a|^2 |b|^2");

    // 同轴旋转相加: q(θ1)·q(θ2) = q(θ1+θ2)  —— 半角逻辑的正确性
    auto qy = [](T theta) {
        const T h = theta / T(2);
        return Q(std::cos(h), T(0), std::sin(h), T(0));
    };
    ck.check(SameQ(qy(T(0.8)) * qy(T(0.6)), qy(T(1.4)), tol),
             std::string("[C13] ") + tag + " q(0.8)*q(0.6) = q(1.4) (同轴相加)");

    // 双覆盖: 转 2π 得 -q (转 4π 才回到 q)
    const T pi = T(3.14159265358979323846);
    ck.check(SameQ(qy(T(0.5) + T(2) * pi), -qy(T(0.5)), tol),
             std::string("[C14] ") + tag + " q(θ+2π) = -q(θ)  (SU(2) 双覆盖)");
}

// ---- [P] 待实现: 如实断言它们仍然 throw ------------------------------------

void TestPending(Checker& ck)
{
    using Q = Quaternion<double>;
    const Q a(1.0, 0.0, 0.0, 0.0);
    const Q b(0.7071067811865476, 0.0, 0.7071067811865476, 0.0);

    // 原 [P1] (operator* 仍是占位) 已于 2026-10-19 兑现 —— 判据已升级到 [C1]–[C14],
    // 这里不再断言它 throw。其余编号保持不变 (与 dev-progress 日志里的记录对得上)。

    ck.check(ThrowsLogicError([&] { (void)a.Inverse(); }),
             "[P2] Inverse 仍是占位 —— 待第 ② 步");
    ck.check(ThrowsLogicError([&] { (void)a.Normalize(); }),
             "[P3] Normalize 仍是占位 —— 待第 ② 步");

    ck.check(ThrowsLogicError([&] { (void)Q::FromAxisAngle(0.0, 1.0, 0.0, 1.2); }),
             "[P4] FromAxisAngle 仍是占位 —— 待第 ③ 步");
    ck.check(ThrowsLogicError([&] { (void)a.ToMatrix(); }),
             "[P5] ToMatrix 仍是占位 —— 待第 ③ 步");
    ck.check(ThrowsLogicError([&] { (void)Q::FromMatrix(Matrix<double>::Identity(3)); }),
             "[P6] FromMatrix 仍是占位 —— 待第 ③ 步");
    ck.check(ThrowsLogicError([&] { (void)a.Rotate(Matrix<double>(3, 1)); }),
             "[P7] Rotate 仍是占位 —— 待第 ③ 步");

    ck.check(ThrowsLogicError([&] { (void)a.AngleOf(); }),
             "[P8] AngleOf 仍是占位 —— 待第 ④ 步");
    ck.check(ThrowsLogicError([&] { (void)Slerp(a, b, 0.5); }),
             "[P9] Slerp 仍是占位 —— 待第 ④ 步");
    ck.check(ThrowsLogicError([&] { (void)Nlerp(a, b, 0.5); }),
             "[P10] Nlerp 仍是占位 —— 待第 ④ 步");
}

int main()
{
    std::printf("== B1 四元数 (plan.md V2.5 / V2.6) —— 骨架期 ==\n");
    std::printf("   已实现: 构造/单位元/加减/一元负/标量乘/共轭/平方范数/范数\n");
    std::printf("   待实现: operator* / Inverse / Normalize / ③ 互转与旋转 / ④ 插值\n");

    Checker ck;

    std::printf("\n-- [A] ① 基本代数 (已实现; double 与 float 各一遍) --\n");
    TestAlgebra<double>(ck, "double");
    TestAlgebra<float>(ck, "float ");

    std::printf("\n-- [B] ② 共轭与范数 (已实现; double 与 float 各一遍) --\n");
    TestConjugateNorm<double>(ck, "double");
    TestConjugateNorm<float>(ck, "float ");

    std::printf("\n-- [C] ① Hamilton 积 (已实现; double 与 float 各一遍) --\n");
    TestHamilton<double>(ck, "double");
    TestHamilton<float>(ck, "float ");

    std::printf("\n-- [P] 待实现部分 (占位判据: 断言仍 throw) --\n");
    TestPending(ck);

    return CheckReport(ck.ok, ck.fail, "quaternion");
}
