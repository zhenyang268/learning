#pragma once

// ===========================================================================
// core/tensor.hpp —— Tensor<T> (占有) + TensorView<T> (非占有)  【框架】
//                    matrix.md §6 第三次订正 / 决策表 #21 (2026-10-23)
// ===========================================================================
//
// ⚠ 本文件现在是**框架**: 契约写全了, 实现体全是 throw。
//   你实现掉一条, demo/core/tensor.cpp 里对应的 [P] 用例就会变红 ——
//   红就是提醒你把占位判据换成真判据 (与 quaternion.hpp 的规矩完全一样)。
//
// ---------------------------------------------------------------------------
// 红线口径 (改之前先读这条, 它决定这个文件里"不该有什么")
// ---------------------------------------------------------------------------
//   · 容器**只准哑**: 存 (data) + 索引 (operator()) + 形状查询 (Shape/Strides/Size)。
//   · **一切指标收缩 / 乘法 / 求逆 / einsum 一律是自由函数**, 绝不进这两个类。
//     (§3.7: 能力用自由函数表达, 容器保持哑的; 与 Multiply(A,B,hint) 同性质)
//   · 不进 Matrix: 这是**独立类型**, 不是 Matrix 的 N 维版, 也不与它隐式互转。
//
// 为什么现在放行 (2026-10-23 用户定案, 理由见 matrix.md §6 第三次订正):
//   上一版结论是"存储层免费 —— 张量的价值全在运算层"。用户据此指出:
//   **这个理由恰恰支持写容器** —— 存储层既然只有 ~100 行且无数学内容, 那它
//   不值得写"完整容器类"(容器+算子+结构标签+表达式模板), 但值得写"哑容器"。
//   放行的是容器那一半; 运算那一半的红线一字未动。
//
// ---------------------------------------------------------------------------
// 单线程, 最简单版本 —— 这是**设计目标**, 不是待补的欠债
// ---------------------------------------------------------------------------
// 明确不做 (写下来防范围蠕变):
//   · 并发读写 / 原子引用计数 / 写时复制 (COW)
//   · 表达式模板 —— 先要正确, 再谈少一次临时对象
//   · 运行期结构标签 (MatrixKind 那种) —— §3.7 已判过: 缓存标签必然过期
//   · 稀疏存储 / SIMD 对齐 / 非行主序布局
//   · 编译期 rank (Tensor<T,rank>) —— 运行期 shape_ 就够; 把 rank 写成模板参数
//     会把"不同 rank 的 view 互转"这条最常用的路径堵死 (Slice 掉一维就变 rank-1)
//
// ---------------------------------------------------------------------------
// 布局: 行主序 (C order), strides_ 以"元素个数"计
// ---------------------------------------------------------------------------
//   shape   {2, 3, 4}          -> Size() = 24
//   strides {12, 4, 1}         -> strides[d] = prod(shape[d+1 ..])
//   (i,j,k) 的扁平下标 = i*12 + j*4 + k
//
// 为什么存 strides 而不是每次乘出来: 这是 TensorView **存在的前提** ——
//   Slice / Permute 只改 shape 与 strides, **一个字节都不搬**。
//   若只存 shape 而每次现算 strides, 就只能表达"稠密行主序"一种布局,
//   view 立刻退化成"拷贝一份", 那就白造了。
//
// ---------------------------------------------------------------------------
// 两个类的关系
// ---------------------------------------------------------------------------
//   Tensor<T>      占有 (vector<T>), 拷贝 = 深拷贝
//   TensorView<T>  非占有 (T* + shape + strides), 拷贝 = 浅拷贝, 不可 resize
//
//   Tensor::View()         整块 view (不复制)
//   View.Slice(axis,b,e)   沿一轴取 [b,e), 指针偏移 b*strides[axis]
//   View.Permute(axes)     换轴 (只重排 shape/strides)
//   View.Reshape(...)      **仅当 IsContiguous()**, 否则 throw
//
// **生命周期契约**: view 不得长于其宿主。与 std::string_view 同款风险 ——
//   文档写清, **不做运行期检查** (检查要引用计数, 与"单线程最简单"冲突)。
//   特别注意: 不要把临时 Tensor 的 View() 存下来用:
//       auto v = MakeTensor().View();   // ← 悬空!
//
// **Reshape 在非连续 view 上 throw, 而不是偷偷拷贝**:
//   numpy 的 reshape 需要时静默复制, 于是 `a.T.reshape(...)` 与 `a.reshape(...)`
//   一个是视图一个是拷贝, 性能与别名语义全不同却不报错。本项目纪律
//   「宁可明确报错, 也不要静默算错」⇒ 要复制请显式写 `Tensor<T>(view)`。
//
// ---------------------------------------------------------------------------
// 实现顺序 (依赖序, 单线程最简单版; 与 demo/core/tensor.cpp 编号对齐)
// ---------------------------------------------------------------------------
//   ① 几何工具  ComputeStrides / ShapeSize / TensorOffset / CheckBoundsAt   ↔ 全员依赖
//   ② Tensor    构造 / Rank / Size / Shape / Strides / IsContiguous /
//               operator() / at                                          ↔ [A][B]
//   ③ 视图与索引 View / operator() / Slice / Permute / Transpose /
//               Reshape / Flatten                                        ↔ [C][D]
//   ④ 深拷贝    Tensor(const TensorView<U>&)                             ↔ [E]
//   ⑤ 逐元素    + - 一元- 标量乘 += / CheckSameShape                     ↔ [F]
//   ⑥ 契约      throw 的位置与类型 (§G 用例已写好, 实现完自然通过)
//
// 已实现 (不用你改): CheckIndexCount / 各一行 getter / Data() /
//                    friend operator*(T, Tensor)
// ===========================================================================

#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

template <typename T> class TensorView;

// ===========================================================================
// ① 几何工具 —— 两个类共用, 保证索引算术只有**一处**实现
// ===========================================================================
//
// 为什么抽成自由函数而不是两个类各写一遍: 索引算术只准有一处实现, 否则
// Tensor 与 TensorView 迟早算得不一样 (这正是 3.3 第 2 条「判据独立于实现路径」
// 要防的东西 —— 两处实现会互相掩护对方的错)。

// 行主序 strides: strides[rank-1] = 1, strides[d] = strides[d+1] * shape[d+1]
//
// 【待实现】例: {2,3,4} -> {12,4,1}。注意从**后往前**算。
inline std::vector<std::size_t> ComputeStrides(const std::vector<std::size_t>& shape)
{
    (void)shape;
    throw std::logic_error("ComputeStrides: 待实现 (tensor 框架 ①)");
}

// 元素个数 = prod(shape)。任一维为 0 => 0 (合法的空张量)
// 注意: 秩 0 的 prod 是 **1** (空乘积约定), 不是 0 —— 那是"标量有一个元素"。
//
// 【待实现】
inline std::size_t ShapeSize(const std::vector<std::size_t>& shape)
{
    (void)shape;
    throw std::logic_error("ShapeSize: 待实现 (tensor 框架 ①)");
}

// 扁平下标 = Σ idxs[d] * strides[d]
//
// 【待实现】
inline std::size_t TensorOffset(const std::size_t* idxs,
                                const std::size_t* strides,
                                std::size_t n)
{
    (void)idxs; (void)strides; (void)n;
    throw std::logic_error("TensorOffset: 待实现 (tensor 框架 ①)");
}

// 形状校验: 下标个数必须与秩一致。
//
// ⚠ 这条**不是可选优化**: 个数少给时 TensorOffset 会读超 strides_ 尾部 (越界读),
//   个数多给时多余下标被静默丢掉 (静默算错)。两者都是本项目明令禁止的。
//   而个数检查是 O(1)、不进内层循环, 不影响性能 —— 所以 operator() 也要查。
//
// 【已实现, 不用你改】
inline void CheckIndexCount(std::size_t got, std::size_t want, const char* who)
{
    if (got != want) {
        throw std::invalid_argument(
            std::string(who) + ": 下标个数 (" + std::to_string(got) +
            ") 与秩 (" + std::to_string(want) +
            ") 不一致 —— 多余下标会被静默忽略、缺少下标会越界读, 所以直接报错");
    }
}

// 逐维越界检查 (只给 at() 用; operator() 故意不查 —— 与 Matrix 同口径)
//
// 【待实现】抛 std::out_of_range, 消息里带上"第几维、下标值、该维长度"。
inline void CheckBoundsAt(const std::size_t* i, const std::size_t* shape,
                          std::size_t n, const char* who)
{
    (void)i; (void)shape; (void)n; (void)who;
    throw std::logic_error("CheckBoundsAt: 待实现 (tensor 框架 ①)");
}

// ===========================================================================
// ② Tensor<T> —— 占有内存的哑容器
// ===========================================================================
template <typename T>
class Tensor
{
public:
    // 默认构造 = **秩 0 的标量张量**: Rank()==0, Size()==1, strides 空。
    // 不是空张量 —— 秩 0 是合法张量 (标量); 要真空张量让某一维为 0:
    //   Tensor<double>({0})  -> Size()==0
    //
    // 【待实现】提示: data_ 要有 1 个元素 (值 T{}), shape_/strides_ 都空。
    Tensor() { throw std::logic_error("Tensor::Tensor(): 待实现 (② 构造)"); }

    // 按形状构造, 元素零初始化 (用 T{} —— 对 Dual / Complex 也成立)
    //
    // 【待实现】注意初始化列表顺序要与**成员声明顺序**一致, 否则 -Wreorder 告警。
    explicit Tensor(std::vector<std::size_t> shape)
    {
        (void)shape;
        throw std::logic_error("Tensor::Tensor(shape): 待实现 (② 构造)");
    }

    // 【待实现】委托给上面那个: : Tensor(std::vector<std::size_t>(shape)) {}
    Tensor(std::initializer_list<std::size_t> shape)
    {
        (void)shape;
        throw std::logic_error("Tensor::Tensor(init_list): 待实现 (② 构造)");
    }

    // 按形状 + 现成数据构造; **个数必须等于 prod(shape)**, 否则 throw。
    //
    // 【待实现】
    Tensor(std::vector<std::size_t> shape, std::vector<T> values)
    {
        (void)shape; (void)values;
        throw std::logic_error("Tensor::Tensor(shape, values): 待实现 (② 构造)");
    }

    // 【待实现】委托
    Tensor(std::initializer_list<std::size_t> shape, std::vector<T> values)
    {
        (void)shape; (void)values;
        throw std::logic_error("Tensor::Tensor(init_list, values): 待实现 (② 构造)");
    }

    // 从 view 深拷贝 —— 唯一合法的"非连续 -> 稠密"通道, 且必须**显式**调用。
    // 用 U 而不是 const T: 非 const 宿主的 View() 给的是 TensorView<T>,
    // 两者都要能吃 (否则 Tensor<double>(t.View()) 编不过)。
    //
    // 【待实现】做法: 按行主序**里程表 (odometer)** 遍历自己的 idx,
    //   用 v.Strides() 算出源偏移后取值。伪码:
    //     idx 全 0;
    //     for n in [0, Size()):
    //         data_[n] = v.Data()[TensorOffset(idx, v.Strides(), rank)];
    //         从最后一维开始 ++idx, 满了进位并归零;
    //   秩 0 单独处理: data_[0] = v();
    template <typename U>
    explicit Tensor(const TensorView<U>& v)
    {
        (void)v;
        throw std::logic_error("Tensor::Tensor(view): 待实现 (④ 深拷贝)");
    }

    // ---- 形状查询 --------------------------------------------------------
    // 这五个是**一行 getter**, 已经实现 (它们没有可错的逻辑)

    std::size_t Rank() const { return shape_.size(); }
    std::size_t Size() const { return data_.size(); }
    const std::vector<std::size_t>& Shape() const { return shape_; }
    const std::vector<std::size_t>& Strides() const { return strides_; }

    // 稠密行主序 <=> strides_ 等于"按 shape 现算"的结果
    //
    // 【待实现】一行: return strides_ == ComputeStrides(shape_);
    //   (注意: 连续性是**布局**性质, 与指针无关。Slice 出来的 view 可能连续
    //    也可能不连续, 见 demo [D6][D7] —— 别去检查内存地址。)
    bool IsContiguous() const
    {
        throw std::logic_error("Tensor::IsContiguous: 待实现 (② 形状)");
    }

    // ---- 元素访问 --------------------------------------------------------
    //
    // 口径 (与 Matrix 一致, 但 operator() **更严一档**):
    //   operator()  只查**下标个数** (O(1)), 不查每维越界
    //   at()        查个数 + 每次越界
    //
    // 【待实现, 两处都用变参展开】写法:
    //     CheckIndexCount(sizeof...(Idx), shape_.size(), "Tensor::operator()");
    //     const std::size_t i[] = { static_cast<std::size_t>(idx)... };
    //     return data_[TensorOffset(i, strides_.data(), sizeof...(Idx))];
    //   ⚠ 秩 0 时 sizeof...(Idx)==0, `const std::size_t i[] = {}` 是零长数组
    //   —— **GNU 扩展**, -Wpedantic 会告警。用 if constexpr 把 0 参单独分出来。

    template <typename... Idx>
    T& operator()(Idx... idx)
    {
        (void)std::initializer_list<int>{ (static_cast<void>(idx), 0)... };
        throw std::logic_error("Tensor::operator(): 待实现 (② 索引)");
    }

    template <typename... Idx>
    const T& operator()(Idx... idx) const
    {
        (void)std::initializer_list<int>{ (static_cast<void>(idx), 0)... };
        throw std::logic_error("Tensor::operator() const: 待实现 (② 索引)");
    }

    // 【待实现】同 operator(), 但多一次 CheckBoundsAt
    template <typename... Idx>
    T& at(Idx... idx)
    {
        (void)std::initializer_list<int>{ (static_cast<void>(idx), 0)... };
        throw std::logic_error("Tensor::at: 待实现 (② 索引)");
    }

    template <typename... Idx>
    const T& at(Idx... idx) const
    {
        (void)std::initializer_list<int>{ (static_cast<void>(idx), 0)... };
        throw std::logic_error("Tensor::at const: 待实现 (② 索引)");
    }

    // 裸数据出口 (行主序一维连续) —— 已实现
    T*       Data()       { return data_.data(); }
    const T* Data() const { return data_.data(); }

    // ---- 视图 (定义在文件末尾, 需要 TensorView 完整) ----------------------
    TensorView<T>       View();
    TensorView<const T> View() const;

    // ---- 逐元素运算 ------------------------------------------------------
    //
    // **不给 operator***: 张量乘法的语义太多 (缩并哪两个指标? 外积? Hadamard?),
    // 一个 operator* 装不下。全部走自由函数 (einsum 那条线)。
    //
    // 【待实现】形状不一致必须 throw, **不静默广播** (广播是另一层语义,
    //   要的话做成显式函数, 别藏在 + 里)。

    Tensor operator+(const Tensor& o) const
    {
        (void)o;
        throw std::logic_error("Tensor::operator+: 待实现 (⑤ 逐元素)");
    }

    Tensor operator-(const Tensor& o) const
    {
        (void)o;
        throw std::logic_error("Tensor::operator-: 待实现 (⑤ 逐元素)");
    }

    // 【待实现】注意用 T(0) - x 而不是 -x: 一元负不在 Scalar concept 里保证
    //   (Dual 有, 但不是所有 T 都有) —— 用减法更稳。
    Tensor operator-() const
    {
        throw std::logic_error("Tensor::operator- (一元): 待实现 (⑤ 逐元素)");
    }

    Tensor operator*(const T& s) const
    {
        (void)s;
        throw std::logic_error("Tensor::operator*(scalar): 待实现 (⑤ 逐元素)");
    }

    // 【已实现】左标量乘走 friend, 转发给成员版; 不用你改
    friend Tensor operator*(const T& s, const Tensor& t) { return t * s; }

    Tensor& operator+=(const Tensor& o)
    {
        (void)o;
        throw std::logic_error("Tensor::operator+=: 待实现 (⑤ 逐元素)");
    }

private:
    // 形状一致性检查 (给 + - += 用)
    //
    // 【待实现】throw std::invalid_argument(std::string(who) + ": 形状不一致");
    void CheckSameShape(const Tensor& o, const char* who) const
    {
        (void)o; (void)who;
        throw std::logic_error("Tensor::CheckSameShape: 待实现 (⑤)");
    }

    std::vector<T>           data_;     // 一维连续存储, 行主序
    std::vector<std::size_t> shape_;    // 各维长度, 如 {2,3,4}
    std::vector<std::size_t> strides_;  // 各维步长, 单位是"元素个数"
};

// ===========================================================================
// ③ TensorView<T> —— 非占有视图
// ===========================================================================
//
// 拷贝 view **不复制数据**, 只是再指一次同一块内存 (与 std::string_view 同)。
// 宿主销毁后 view 悬空 —— 见头注释的生命周期契约。
//
// T 可以是 const T (来自 const Tensor 的 View()): 此时只能读不能写,
// 由**编译器**保证, 不额外做运行期检查。
template <typename T>
class TensorView
{
public:
    // 【待实现】默认构造: ptr_ = nullptr, shape_/strides_ 空。
    //   ⚠ 不能 = default —— ptr_ 是裸指针, = default 不初始化它 (读它就是 UB)。
    TensorView() { throw std::logic_error("TensorView::TensorView(): 待实现 (③)"); }

    // 【待实现】存下 ptr/shape/strides; **shape 与 strides 长度不等必须 throw**。
    //   (长度不等 = 布局自相矛盾, 之后所有索引都会错)
    TensorView(T* ptr, std::vector<std::size_t> shape, std::vector<std::size_t> strides)
    {
        (void)ptr; (void)shape; (void)strides;
        throw std::logic_error("TensorView::TensorView(ptr,shape,strides): 待实现 (③)");
    }

    // ---- 形状查询: 与 Tensor 完全同口径 --------------------------------
    // 一行 getter 已实现

    std::size_t Rank() const { return shape_.size(); }
    std::size_t Size() const { return ShapeSize(shape_); }
    const std::vector<std::size_t>& Shape() const { return shape_; }
    const std::vector<std::size_t>& Strides() const { return strides_; }

    // 【待实现】一行
    bool IsContiguous() const
    {
        throw std::logic_error("TensorView::IsContiguous: 待实现 (③)");
    }

    // ---- 元素访问: 与 Tensor 同口径 (operator() 查个数, at() 再查越界) ----

    template <typename... Idx>
    T& operator()(Idx... idx)
    {
        (void)std::initializer_list<int>{ (static_cast<void>(idx), 0)... };
        throw std::logic_error("TensorView::operator(): 待实现 (③)");
    }

    template <typename... Idx>
    const T& operator()(Idx... idx) const
    {
        (void)std::initializer_list<int>{ (static_cast<void>(idx), 0)... };
        throw std::logic_error("TensorView::operator() const: 待实现 (③)");
    }

    template <typename... Idx>
    T& at(Idx... idx)
    {
        (void)std::initializer_list<int>{ (static_cast<void>(idx), 0)... };
        throw std::logic_error("TensorView::at: 待实现 (③)");
    }

    template <typename... Idx>
    const T& at(Idx... idx) const
    {
        (void)std::initializer_list<int>{ (static_cast<void>(idx), 0)... };
        throw std::logic_error("TensorView::at const: 待实现 (③)");
    }

    T*       Data()       { return ptr_; }
    const T* Data() const { return ptr_; }

    // ---- 切片: 沿 axis 取 [begin, end), **只改指针偏移与 shape, 不搬数据** ---
    //
    // 契约: 0 <= begin <= end <= shape[axis], 否则 throw。
    // 【待实现】新指针 = ptr_ + begin * strides_[axis]; 新 shape 只改第 axis 维。
    //   ⚠ strides **不变** —— 切一维不会改变剩下各维的步长。
    TensorView Slice(std::size_t axis, std::size_t begin, std::size_t end) const
    {
        (void)axis; (void)begin; (void)end;
        throw std::logic_error("TensorView::Slice: 待实现 (③)");
    }

    // ---- 换轴: 只重排 shape/strides 的顺序, **不搬数据** -------------------
    //
    // axes[d] = "新第 d 轴取原第几轴"。
    // 【待实现】必须是 0..rank-1 的**排列**, 否则 throw —— 重复或缺失会让
    //   Size() 与真实内存对不上, 是静默越界的源头。
    //   提示: 用 vector<bool> seen 查重; 新 shape[d] = 旧 shape[axes[d]]。
    TensorView Permute(const std::vector<std::size_t>& axes) const
    {
        (void)axes;
        throw std::logic_error("TensorView::Permute: 待实现 (③)");
    }

    // 便捷: 转置 (仅秩 2)
    //
    // 为什么只做秩 2: numpy 的 `.T` 对 >=3 阶是**全轴反转**而不是转置前两轴,
    // 那个歧义本项目不复制。一般换轴请用 Permute。
    //
    // 【待实现】秩 != 2 就 throw, 否则 return Permute({1, 0});
    TensorView Transpose() const
    {
        throw std::logic_error("TensorView::Transpose: 待实现 (③)");
    }

    // ---- 变形: **仅当连续**, 否则 throw (不偷偷拷贝) -----------------------
    //
    // 【待实现】顺序很关键:
    //   1. !IsContiguous() -> throw std::logic_error (不是 invalid_argument:
    //      入参本身没错, 是"这个 view 做不到这个操作")
    //   2. ShapeSize(newShape) != Size() -> throw std::invalid_argument
    //   3. 重算 strides 后再构造。⚠ **先把 strides 算成局部变量再 move newShape**:
    //        TensorView(ptr_, std::move(newShape), ComputeStrides(newShape))
    //      是错的 —— 函数实参求值顺序**未指定**, std::move 若先绑到形参就把
    //      newShape 搬空了, ComputeStrides 于是算一个空 vector。
    //      (2026-10-23 实测就是这么 SIGABRT 的。)
    TensorView Reshape(std::vector<std::size_t> newShape) const
    {
        (void)newShape;
        throw std::logic_error("TensorView::Reshape: 待实现 (③)");
    }

    // 【待实现】一行: return Reshape({Size()});
    TensorView Flatten() const
    {
        throw std::logic_error("TensorView::Flatten: 待实现 (③)");
    }

private:
    // 下标 -> 偏移; 口径与 Tensor 同 (查个数, 不查越界)
    // 【待实现】与 Tensor::operator() 同一套写法 (含秩 0 的 if constexpr 分支)
    template <typename... Idx>
    std::size_t OffsetOf(Idx... idx) const
    {
        (void)std::initializer_list<int>{ (static_cast<void>(idx), 0)... };
        throw std::logic_error("TensorView::OffsetOf: 待实现 (③)");
    }

    T*                       ptr_;
    std::vector<std::size_t> shape_;
    std::vector<std::size_t> strides_;
};

// ===========================================================================
// 依赖 TensorView 的成员, 定义放最后 (此时 TensorView 才完整)
// ===========================================================================

// 【待实现】一行: return TensorView<T>(data_.data(), shape_, strides_);
template <typename T>
TensorView<T> Tensor<T>::View()
{
    throw std::logic_error("Tensor::View: 待实现 (③)");
}

template <typename T>
TensorView<const T> Tensor<T>::View() const
{
    throw std::logic_error("Tensor::View const: 待实现 (③)");
}