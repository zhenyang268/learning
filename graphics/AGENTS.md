# AGENTS.md — graphics 项目宪法

本仓库是"用计算机图形学/物理可视化消化数学自学内容"的学习项目。
AI 助手在此项目中的工作方式被本文件严格约束：**只做增量填格子，禁止生成新项目、禁止无关重构。**

## 1. 目录职责（固定，不得随意新增顶层目录）

```
src/
├── core/        第0层 | 无内部依赖
│                scalar.hpp : Scalar/BuiltinScalar/DiffScalar 概念 + ScalarTraits<T>
│                matrix.hpp : 纯容器。存储/访问/构造/逐元素四则/Transpose/Adjoint
│                vec.hpp    : Vec<T,2/3/4> 固定小向量（图形学热路径, 零堆分配）
├── complex/     第0.5层 | 依赖 core/scalar
│                自实现 Complex<T>（禁用 std::complex）, 特化 ScalarTraits
├── autodiff/    第0.5层 | 依赖 core/scalar
│                dual.hpp   : Dual<T>, 前向模式
│                grad_graph.hpp / jacobian.hpp : 计算图与雅可比传播
├── linalg/      第1层 | 依赖 core + (traits 语义)
│                按能力分文件: solve.hpp / inverse.hpp / eigen.hpp / svd.hpp
│                qr.hpp / cholesky.hpp / expm.hpp / classify.hpp
│                全部自由函数 + 结果结构体, 不做成员方法
├── quat/        第2层 | 依赖 core, linalg
│                Quaternion<T> + 旋转矩阵互转 + Slerp
├── math/        第2层 | 依赖 core, linalg
│                ode.hpp (Euler/RK4/Verlet/辛) / pde.hpp (FDM/迭代解法)
│                mc.hpp (RNG/低差异序列/蒙特卡洛) / fourier.hpp (DFT)
├── gfx/         第3层 | 依赖以上全部
│                窗口(SDL2)/帧缓冲/光栅化/相机/raytracer/scene/ImGui 面板
demo/
├── <模块名>/    每个模块一个子目录, 内含若干 demo 程序
│                demo 的 main(): 计算 + assert 校验 + return 0 == 测试通过
third_party/     只允许 SDL2 / imgui 这类 IO 库, 禁止任何数学库(Eigen/glm/std::complex)
build/           CMake 产物, 不提交内容
```

依赖方向只允许从上到下（core → complex/autodiff → linalg → quat/math → gfx）。
横向依赖（如 linalg 内部文件互相依赖）允许，逆向依赖禁止。

## 2. 硬规则

1. **Matrix 容器保持哑**：新的数学运算一律以自由函数加入 linalg/math，
   禁止向 Matrix 类添加成员方法（求逆/分解/行列式/apply 等一概不收）。
2. **算法层零类型分支**：实数/复数/Dual 的差异只允许存在于 `ScalarTraits<T>` 特化里。
   `linalg/`、`math/` 代码中禁止 `if constexpr (is_complex)` 之类的分支；
   比较、容差、选主元一律用 `ScalarTraits<T>::real_t` 与 `abs2()`。
   逃生门 `if constexpr` 仅限数学上真正不同的场合，需在 PR/注释中说明理由。
3. **不引入现成计算实现**：矩阵/复数/微分方程/采样/傅里叶全部自写。
   Python 侧 matplotlib 仅用于画验证曲线，禁止 numpy/scipy 参与正式实现。
4. **demo 即测试**：每个功能必须配 demo 并注册进 `demo/<模块>/CMakeLists.txt`
   与 ctest。新增功能未配测试不算完成。
5. **每次会话只做当前阶段 plan.md 中的任务**，不顺手"优化"其他模块；
   发现其他模块 bug 记录到 plan.md 的「待办/已知问题」区，不擅自修。
6. **禁止删除或改写本文件与 plan.md 的验收标准**，只允许追加「待办/已知问题」。

## 3. 工程约定

- C++20（concepts 定义 Scalar）；头文件 `#pragma once`；模板定义放头文件。
- 命名：类型 PascalCase，函数/变量 snake_case，常量 kXxx；文件级自由函数不带类包装。
- 构造 Matrix 失败（尺寸不符/奇异）抛 `std::invalid_argument`；数值容差用
  `std::numeric_limits<real_t>::epsilon()` 的倍数，不用全局魔数。
- 每个头文件顶部一行注释写明：数学出处（notebook 对应章节）+ 使用的算法。
- 行主序存储；`operator[]` 返回**即时派生**的 `T*`（不存储行指针数组）。

## 4. 常用命令

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## 5. 给 agent 的任务模板

提交给 agent 的任务应形如：
「实现 `linalg/xxx.hpp` 中的 `Fxx()`，签名如下…，依据 AGENTS.md 第 2 条走 traits，
新增 `demo/xxx/` 测试并注册 ctest，通过标准见 plan.md 阶段 N 验收项 M。」
