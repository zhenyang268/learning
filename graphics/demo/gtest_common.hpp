#pragma once
// ---------------------------------------------------------------------------
// demo/gtest_common.hpp -- gtest 用例侧的公共工具
//
// 与 demo/common.hpp 的区别:
//   common.hpp       量度 (MaxAbs2 / NormInf / IsApprox / RandMatrix) —— **两套框架共用**
//   gtest_common.hpp 断言适配层 —— 只给 gtest 用例用
//
// 为什么不把 gtest 头塞进 common.hpp: 还有未迁移的老用例不链接 gtest,
// 一旦 common.hpp 依赖 <gtest/gtest.h>, 那些 TU 会直接编不过 (缺 include 路径)。
// 迁移期两套并存, 所以公共量度与断言层必须分开。
//
// ---------------------------------------------------------------------------
// 为什么需要 EXPECT_EXACT_THROW 而不是直接用 EXPECT_THROW
// ---------------------------------------------------------------------------
// std::invalid_argument 与 std::out_of_range **都派生自 std::logic_error**。
// 所以 EXPECT_THROW(stmt, std::logic_error) 对这两者都成立 —— 契约段会**假绿**:
// 实现把 invalid_argument 误写成 out_of_range, 用例照样全绿。
// 契约判据的全部价值就在"抛哪一种", 所以必须判**精确类型**。
// 全仓只准通过本文件的 EXPECT_EXACT_THROW 来判异常类型, 不要直接写 EXPECT_THROW。
// ---------------------------------------------------------------------------

#include <gtest/gtest.h>

#include "common.hpp"

#include <exception>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <type_traits>
#include <utility>

#if defined(__has_include)
#  if __has_include(<cxxabi.h>)
#    include <cxxabi.h>
#    include <cstdlib>
#    define GRAPHICS_HAVE_CXA_DEMANGLE 1
#  endif
#endif

namespace graphics_test
{

// typeid(...).name() 是可读性极差的修饰名, 能 demangle 就 demangle
inline std::string Demangle(const char* mangled)
{
#ifdef GRAPHICS_HAVE_CXA_DEMANGLE
    int status = 0;
    char* out = abi::__cxa_demangle(mangled, nullptr, nullptr, &status);
    if (status == 0 && out != nullptr) {
        std::string s(out);
        std::free(out);
        return s;
    }
    return mangled;
#else
    return mangled;
#endif
}

struct ExactThrowOutcome
{
    bool threw = false;         // 是否抛了
    bool exact = false;         // 抛了且类型**恰好**是期望的那个
    std::string actual_type;    // 实际抛出的类型名 (未抛时为空)
};

// 捕获并判定"是否恰好抛出 E"
//
// 判 E 而不是判基类: 见文件头。这里先按 E 抓, 再比 typeid ——
// 只 catch (const E&) 不够, 因为派生类也会被这条 catch 命中。
template <typename E, typename F>
ExactThrowOutcome CheckExactThrow(F&& f)
{
    static_assert(std::is_polymorphic<E>::value,
                  "E 必须是多态异常类型, 否则 typeid 取不到动态类型");
    ExactThrowOutcome r;
    try {
        f();
    } catch (...) {
        r.threw = true;
        std::exception_ptr ep = std::current_exception();
        try {
            std::rethrow_exception(ep);
        } catch (const E& e) {
            r.exact = (typeid(e) == typeid(E));
            r.actual_type = Demangle(typeid(e).name());
        } catch (const std::exception& e) {
            r.actual_type = Demangle(typeid(e).name());
        } catch (...) {
            r.actual_type = "(非 std::exception)";
        }
    }
    return r;
}

inline std::string ToString(const ExactThrowOutcome& o)
{
    if (!o.threw) return "没有抛出";
    return "抛出了 " + o.actual_type;
}

}  // namespace graphics_test

// 断言: statement 恰好抛出 ExceptionType (不允许派生类冒充)
#define EXPECT_EXACT_THROW(statement, ExceptionType)                                  \
    do {                                                                              \
        const auto _gfx_outcome =                                                     \
            ::graphics_test::CheckExactThrow<ExceptionType>([&] { statement; });      \
        EXPECT_TRUE(_gfx_outcome.exact)                                               \
            << "期望恰好抛出 " #ExceptionType ", 实际: "                              \
            << ::graphics_test::ToString(_gfx_outcome);                               \
    } while (0)

// 断言: statement 不抛
#define EXPECT_NO_THROW_STMT(statement)                                               \
    do {                                                                              \
        const auto _gfx_outcome =                                                     \
            ::graphics_test::CheckExactThrow<std::exception>([&] { statement; });     \
        EXPECT_FALSE(_gfx_outcome.threw)                                              \
            << "不该抛, 实际: " << ::graphics_test::ToString(_gfx_outcome);           \
    } while (0)

// 断言: 两个 Matrix 在 ∞ 范数下相近 (D7: 比较走 abs2, 阈值取 tol^2, 这里只是开方打印)
#define EXPECT_MATRIX_APPROX(A, B, tol)                                               \
    do {                                                                              \
        const auto& _gfx_a = (A);                                                     \
        const auto& _gfx_b = (B);                                                     \
        EXPECT_TRUE(IsApprox(_gfx_a, _gfx_b, (tol)))                                  \
            << "‖A-B‖∞ = " << NormInf(_gfx_a - _gfx_b) << " > " << (tol);             \
    } while (0)
