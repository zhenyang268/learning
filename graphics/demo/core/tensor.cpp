// tensor.cpp -- Tensor<T> / TensorView<T> 验收骨架
//               (matrix.md §6 第三次订正 / 决策表 #21, 2026-10-23)
//
// ---------------------------------------------------------------------------
// 本文件现在是**骨架**, 不是验收
// ---------------------------------------------------------------------------
// 接口契约已经全部写在 src/core/tensor.hpp 的注释里, 但实现体都是 throw。
// 所以这里对每条占位都写「**断言它仍然 throw**」的用例 —— 这正是 smoke.cpp
// 对 qr_inverse、quaternion.cpp 对 9 个入口用的那条规矩:
//   **不假装它可用**, 而是如实钉住"它现在不可用"。
//
// 为什么这样写而不是留空: 你实现掉一条, 对应的 [P] 就会**变红** ——
//   红就是提醒, 提醒你把占位判据换成真判据。
//   留空的话, 实现了也没人知道, 没实现也没人知道。
//
// 真判据的完整规格见文件末尾「升级路径」: 实现到哪一步, 就把对应 [P] 换成
// 那一节。**别一次全换** —— 一条一条兑, 每一步都跑一遍。
//
// ---------------------------------------------------------------------------
// 实现顺序 (与 tensor.hpp 头部的 ①–⑥ 一一对应)
// ---------------------------------------------------------------------------
//   ① 几何工具  ComputeStrides / ShapeSize / TensorOffset / CheckBoundsAt  [P1]–[P4]
//   ② Tensor    构造 / IsContiguous / operator() / at                      [P5]–[P9]
//   ③ 视图      View / Slice / Permute / Transpose / Reshape / Flatten     [P10]–[P17]
//   ④ 深拷贝    Tensor(const TensorView<U>&)                               [P18]
//   ⑤ 逐元素    + - 一元- 标量乘 +=                                        [P19]–[P24]
//   ⑥ 契约      [G] 段代码已写好 (在文件里但不注册), 实现 ② 后打开
//
// ---------------------------------------------------------------------------
// 分域纪律 (3.3 第 1 条)
// ---------------------------------------------------------------------------
// 模板 bug 只在实例化时暴露: double 过不代表 float 过 (本项目连续误判过 5 次)。
// [G] 契约段用 double 与 float **各跑一遍**。
//
// 运行:
//   ./graphics/run_tests.sh              # 全量 (本用例注册名: tensor)
//   ./graphics/build/demo/core/tensor    # 单跑
//   ./graphics/run_tests.sh -w           # 顺带开 -Wall -Wextra

#include "../common.hpp"
#include "core/tensor.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

// ---- 小工具 ---------------------------------------------------------------
//
// ⚠ 这两条不能混用: std::invalid_argument 与 std::out_of_range **都派生自
//   std::logic_error**。所以"断言它 throw logic_error"在契约段会**假绿**。
//   想精确区分类型必须用下面第二个。这个坑写在头注释里, 免得以后有人为了
//   "让用例变绿"把 [G] 改成抓 logic_error, 那等于什么都没测。

template <typename F>
bool ThrowsLogicError(F f)
{
    try { f(); } catch (const std::logic_error&) { return true; } catch (...) { return false; }
    return false;
}

template <typename F>
bool ThrowsInvalidArgument(F f)
{
    try { f(); } catch (const std::invalid_argument&) { return true; } catch (...) { return false; }
    return false;
}

// ---- [P] 待实现: 如实断言它们仍然 throw ------------------------------------

void TestPending(Checker& ck)
{
    using D = double;

    // ---- ① 几何工具: 四个自由函数 ----
    ck.check(ThrowsLogicError([&] { (void)ComputeStrides({2, 3, 4}); }),
             "[P1]  ComputeStrides 仍是占位 —— 待 ①");
    ck.check(ThrowsLogicError([&] { (void)ShapeSize({2, 3, 4}); }),
             "[P2]  ShapeSize 仍是占位 —— 待 ①");
    ck.check(ThrowsLogicError([&] {
                 const std::size_t i[3] = {1, 2, 3};
                 const std::size_t s[3] = {12, 4, 1};
                 (void)TensorOffset(i, s, 3);
             }),
             "[P3]  TensorOffset 仍是占位 —— 待 ①");
    ck.check(ThrowsLogicError([&] {
                 const std::size_t i[1] = {5};
                 const std::size_t s[1] = {3};
                 CheckBoundsAt(i, s, 1, "probe");
             }),
             "[P4]  CheckBoundsAt 仍是占位 —— 待 ①");

    // ---- ② Tensor 构造与形状 ----
    ck.check(ThrowsLogicError([&] { (void)Tensor<D>(); }),
             "[P5]  Tensor() 默认构造仍是占位 —— 待 ②");
    ck.check(ThrowsLogicError([&] { (void)Tensor<D>({2, 3}); }),
             "[P6]  Tensor({2,3}) 仍是占位 —— 待 ②");
    ck.check(ThrowsLogicError([&] {
                 (void)Tensor<D>(std::vector<std::size_t>{2, 3},
                                 std::vector<D>{1, 2, 3, 4, 5, 6});
             }),
             "[P7]  Tensor(shape, values) 仍是占位 —— 待 ②");

    // IsContiguous 与 operator() 要等构造可用才能测; 现在连构造都 throw,
    // 只能钉"整条链路不通"。实现 ② 之后把这两条换成 [A3][B1]。
    ck.check(ThrowsLogicError([&] { (void)Tensor<D>({2, 3}).IsContiguous(); }),
             "[P8]  构造不通 => IsContiguous 取不到 (待 ② 后改写为 [A3])");
    ck.check(ThrowsLogicError([&] {
                 Tensor<D> t({2, 3});
                 (void)t(0, 0);
             }),
             "[P9]  构造不通 => operator() 取不到 (待 ② 后改写为 [B1])");

    // ---- ③ 视图 ----
    ck.check(ThrowsLogicError([&] { (void)TensorView<D>(); }),
             "[P10] TensorView() 默认构造仍是占位 —— 待 ③");
    ck.check(ThrowsLogicError([&] {
                 D buf[6] = {0, 1, 2, 3, 4, 5};
                 (void)TensorView<D>(buf, {2, 3}, {3, 1});
             }),
             "[P11] TensorView(ptr,shape,strides) 仍是占位 —— 待 ③");
    ck.check(ThrowsLogicError([&] {
                 Tensor<D> t({2, 3});
                 (void)t.View();
             }),
             "[P12] Tensor::View 仍是占位 —— 待 ③ (且 ② 要先通)");
    ck.check(ThrowsLogicError([&] {
                 D buf[6] = {0, 1, 2, 3, 4, 5};
                 (void)TensorView<D>(buf, {2, 3}, {3, 1}).Slice(0, 0, 1);
             }),
             "[P13] Slice 仍是占位 —— 待 ③");
    ck.check(ThrowsLogicError([&] {
                 D buf[6] = {0, 1, 2, 3, 4, 5};
                 (void)TensorView<D>(buf, {2, 3}, {3, 1}).Permute({1, 0});
             }),
             "[P14] Permute 仍是占位 —— 待 ③");
    ck.check(ThrowsLogicError([&] {
                 D buf[6] = {0, 1, 2, 3, 4, 5};
                 (void)TensorView<D>(buf, {2, 3}, {3, 1}).Transpose();
             }),
             "[P15] Transpose 仍是占位 —— 待 ③");
    ck.check(ThrowsLogicError([&] {
                 D buf[6] = {0, 1, 2, 3, 4, 5};
                 (void)TensorView<D>(buf, {2, 3}, {3, 1}).Reshape({6});
             }),
             "[P16] Reshape 仍是占位 —— 待 ③");
    ck.check(ThrowsLogicError([&] {
                 D buf[6] = {0, 1, 2, 3, 4, 5};
                 (void)TensorView<D>(buf, {2, 3}, {3, 1}).Flatten();
             }),
             "[P17] Flatten 仍是占位 —— 待 ③");

    // ---- ④ 深拷贝 ----
    ck.check(ThrowsLogicError([&] {
                 D buf[6] = {0, 1, 2, 3, 4, 5};
                 const TensorView<const D> v(buf, {2, 3}, {3, 1});
                 (void)Tensor<D>(v);
             }),
             "[P18] Tensor(view) 深拷贝仍是占位 —— 待 ④");

    // ---- ⑤ 逐元素 (连构造都不通, 只能钉"链路不通") ----
    ck.check(ThrowsLogicError([&] {
                 Tensor<D> a({2, 2});
                 Tensor<D> b({2, 2});
                 (void)(a + b);
             }),
             "[P19] operator+ 取不到 —— 待 ⑤ (且 ② 要先通)");
    ck.check(ThrowsLogicError([&] {
                 Tensor<D> a({2, 2});
                 Tensor<D> b({2, 2});
                 (void)(a - b);
             }),
             "[P20] operator- 取不到 —— 待 ⑤ (且 ② 要先通)");
    ck.check(ThrowsLogicError([&] {
                 Tensor<D> a({2, 2});
                 (void)(-a);
             }),
             "[P21] 一元 - 取不到 —— 待 ⑤ (且 ② 要先通)");
    ck.check(ThrowsLogicError([&] {
                 Tensor<D> a({2, 2});
                 (void)(a * D(2));
             }),
             "[P22] 标量乘取不到 —— 待 ⑤ (且 ② 要先通)");
    ck.check(ThrowsLogicError([&] {
                 Tensor<D> a({2, 2});
                 (void)(D(2) * a);
             }),
             "[P23] 左标量乘取不到 —— 待 ⑤ (且 ② 要先通)");
    ck.check(ThrowsLogicError([&] {
                 Tensor<D> a({2, 2});
                 Tensor<D> b({2, 2});
                 a += b;
             }),
             "[P24] += 取不到 —— 待 ⑤ (且 ② 要先通)");
}

// ---- [G] 契约: 定义在此, 但骨架期**不注册** (见 main) ---------------------
//
// 为什么现在不跑: 契约段断言的是"抛 invalid_argument / out_of_range",
// 而骨架期抛的是 logic_error("待实现"), 于是它们会**全部变红**。
// red-by-design 的用例一多, "全绿"这个信号就失去意义了。
// 所以代码留着, 等 ② 实现后把 main 里那三行注释放开即可。
//
// 注意 [G2] 用的是 ThrowsLogicError 而不是 ThrowsOutOfRange: 因为
// out_of_range 派生自 logic_error —— 这样写在骨架期是绿的 (占位抛的也是
// logic_error), 但实现后如果**没抛**就会红, 仍然有守住"越界要报错"的作用。
// 只是它分不清 out_of_range 与别的 logic_error, 精度差一档, 记在这里。

template <typename T>
void TestContract(Checker& ck, const char* tag)
{
    ck.check(ThrowsInvalidArgument([&] {
                 Tensor<T> t({2, 3, 4});
                 (void)t(0, 1);
             }),
             std::string("[G1] ") + tag + " 下标个数 2 != 秩 3 => invalid_argument");

    ck.check(ThrowsLogicError([&] {
                 Tensor<T> t({2, 3});
                 (void)t.at(0, 5);
             }),
             std::string("[G2] ") + tag + " at 越界 => out_of_range (派生自 logic_error)");

    ck.check(ThrowsInvalidArgument([&] {
                 (void)Tensor<T>(std::vector<std::size_t>{2, 2},
                                 std::vector<T>{T(1), T(2), T(3)});
             }),
             std::string("[G3] ") + tag + " values 个数 != prod(shape) => invalid_argument");

    ck.check(ThrowsInvalidArgument([&] {
                 Tensor<T> t({2, 3});
                 (void)t.View().Permute({0, 0});
             }),
             std::string("[G4] ") + tag + " Permute 重复轴 => invalid_argument");

    ck.check(ThrowsInvalidArgument([&] {
                 Tensor<T> t({2, 3});
                 (void)t.View().Slice(0, 1, 5);
             }),
             std::string("[G5] ") + tag + " Slice 越过该维长度 => invalid_argument");

    ck.check(ThrowsInvalidArgument([&] {
                 Tensor<T> t({2, 3});
                 (void)t.View().Reshape({5});
             }),
             std::string("[G6] ") + tag + " Reshape 元素个数不符 => invalid_argument");

    ck.check(ThrowsLogicError([&] {
                 Tensor<T> t({2, 3});
                 (void)t.View().Permute({1, 0}).Reshape({6});
             }),
             std::string("[G7] ") + tag + " 非连续 view 上 Reshape => throw (不偷偷拷贝)");

    ck.check(ThrowsInvalidArgument([&] {
                 Tensor<T> t({2, 3, 4});
                 (void)t.View().Transpose();
             }),
             std::string("[G8] ") + tag + " 秩 3 上 Transpose => invalid_argument");
}

// ---- main -----------------------------------------------------------------

int main()
{
    std::printf("== Tensor / TensorView —— 骨架期 ==\n");
    std::printf("   契约已在 src/core/tensor.hpp 写全; 实现体全 throw\n");
    std::printf("   实现顺序: ① 几何工具 ② Tensor ③ 视图 ④ 深拷贝 ⑤ 逐元素\n");
    std::printf("   实现掉一条, 对应 [P] 就变红 —— 那就把 [P] 换成真判据\n");

    Checker ck;

    std::printf("\n-- [P] 待实现部分 (占位判据: 断言仍 throw) --\n");
    TestPending(ck);

    // ---- [G] 骨架期故意不跑 ----
    // 打开条件: ② Tensor 构造 + operator()/at 实现之后。
    //
    // std::printf("\n-- [G] 契约 (实现 ② 后打开) --\n");
    // TestContract<double>(ck, "double");
    // TestContract<float>(ck, "float ");

    return CheckReport(ck.ok, ck.fail, "tensor");
}

// ===========================================================================
// 升级路径 —— 实现到位后, 把对应 [P] 换成下面这些真判据
// ===========================================================================
//
// 【A】形状与 strides (换掉 [P1][P2][P6][P8]) —— 结构性, 逐位精确
//   [A1] Tensor<double>({2,3,4}): Rank=3, Size=24
//   [A2] Strides() == {12,4,1}          ← 行主序, 单位是元素个数
//   [A3] IsContiguous() == true
//   [A4] 秩 0: Rank=0, Size=1, Strides().empty()
//   [A5] 含 0 维: Tensor({2,0,3}) -> Size=0, Rank=3 (合法空张量, 不 throw)
//   [A6] Tensor({1}) -> Strides() == {1}
//
// 【B】索引与行主序 (换掉 [P9]) —— **两条独立路径**
//   用 0,1,2,...,n-1 顺序填充, 让"扁平下标 == 元素值", 索引错会肉眼可见。
//   [B1] t(i,j,k) == i*12 + j*4 + k     (手算 vs 实现)
//   [B2] 三层循环顺次遍历 == Data() 线性扫描 (内存布局 vs 索引算术)
//   [B3] operator() 返回**引用**: 写 t(1,2,3)=99 后 Data()[1*12+2*4+3]==99
//   [B4] 秩 0: t() = 标量本身
//
// 【C】视图别名 (换掉 [P12]–[P15]) —— **本文件最核心, 只有它能抓住"偷偷拷贝"**
//   一个"看起来像视图"的拷贝实现能通过 [A][B][F] 全部。
//   [C1] 改宿主 => view 看见
//   [C2] 改 view => 宿主看见
//   [C3] View().Data() == 宿主 Data()      (零拷贝, 指针相等)
//   [C4] Slice(1,1,3) => shape (2,2,4)
//   [C5] Slice 指针偏移 == begin*strides[axis]
//   [C6] sl(i,j,k) == host(i,j+1,k)
//   [C7] Permute({2,0,1}) => shape (4,2,3)
//   [C8] Permute 后 strides == {1,12,4}   (shape 与 strides 一起换)
//   [C9] pm(k,i,j) == host(i,j,k)
//   [C10] 秩 2 Transpose: mt(j,i) == m(i,j), 零拷贝
//
// 【D】非连续性 (换掉 [P16][P17])
//   [D1] 原始 view 连续
//   [D2] Permute 后**不连续**
//   [D3] 非连续 view 上 Reshape 必须 throw
//   [D4] 连续 view 上 Reshape 正常
//   [D5] 非连续 view 上 Flatten 必须 throw
//   [D6] Slice(0,...) 仍连续          ← 只削外层
//   [D7] Slice(1,...) **不连续**      ← 内层被跳过
//
// 【E】显式深拷贝 (换掉 [P18])
//   [E1] Tensor(view) 保留 view 的 shape
//   [E2] 拷贝结果连续 (strides 重算)
//   [E3] 拷贝后指针不同
//   [E4] dense(k,i,j) == host(i,j,k)   ← odometer 写对了的证明
//   [E5] 改宿主后 dense 不变           ← 真深拷贝
//   [E6] 拷贝结果可 Reshape
//
// 【F】逐元素 (换掉 [P19]–[P24])
//   [F1]–[F6] + - 一元- 标量乘(两侧) +=
//   [F7] 形状不一致的 + 必须 throw (不静默广播)
//
// 【G】契约: 见 TestContract, 实现 ② 之后把 main 里那三行注释放开。
//
// ---------------------------------------------------------------------------
// 注入对照 (3.3 第 6 条) —— 实现完至少做这两条, 否则不知道闸门有没有牙
// ---------------------------------------------------------------------------
//   ① 把 TensorOffset 的和式写成 idxs[d]*strides[n-1-d] => [B1][C6][C9] 变红
//   ② 把 View() 改成"拷贝一份再给 view 指过去" =>
//      **只有 [C1][C2][C3] 变红, 其余全绿** —— 这正是 [C] 段存在的理由。
//      与 so3 那次「Hat 整体翻号, 立方恒等式仍全绿」是同一类盲区。
// ===========================================================================