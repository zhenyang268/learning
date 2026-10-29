# 可复用 prompt：把剩余 ctest 用例迁移到 googletest

> 用法：把本文件**全文**交给负责迁移的 agent（一次只让它做 1–2 个文件）。
> 里面的路径全部是绝对路径，可跨会话、跨工具直接使用。
> 读完本文件前**不要动任何文件**。

---

## 你的任务

把 `graphics/demo/` 下**还没迁移**的 ctest 用例从手写 `Checker` 框架改成 googletest，
**判据一条都不许删、阈值一个都不许改**。

已经迁完的样板（**先读这两份，照抄它们的风格**）：

- `/Users/apple/develop/learning/graphics/demo/quat/quaternion.cpp` —— 双域 `TYPED_TEST` 的样板（663 行 → 14 个用例，72 条判据编号全保留）
- `/Users/apple/develop/learning/graphics/demo/core/scalar_traits.cpp` —— 单域 `TEST` 的样板（11 条 `static_assert` 留在文件作用域）
- `/Users/apple/develop/learning/graphics/demo/gtest_common.hpp` —— 断言适配层，**只准通过它判异常类型**

## 待迁移清单（12 个文件；`solve` / `expm` 已迁完）

数字是**迁移前** `./graphics/run_tests.sh -d` 报的断言数，用来核对「一条没丢」。迁移后 `-d` 会改报「用例 N 个」，**不能拿它直接比**。

> **开工前先跑一次** `grep -rn graphics_add_gtest graphics/demo/*/CMakeLists.txt` 看当前真实进度 ——
> 这个仓可能有别的 agent 同时在迁，本表可能已经过期。

| 文件 | 行数 | 断言数 | 备注 |
|---|---:|---:|---|
| `demo/linalg/inverse.cpp` | 282 | 85 | 含 `Classify` 的 SPD 档分发，契约段要精确到子类 |
| `demo/linalg/cholesky.cpp` | 283 | 98 | |
| `demo/linalg/eigen.cpp` | 247 | 119 | |
| `demo/linalg/svd.cpp` | 299 | 266 | 断言最多的一份 |
| `demo/linalg/qr.cpp` | 251 | 189 | |
| ~~`demo/linalg/solve.cpp`~~ | 293 | 34 | ✅ 已迁（另一 agent，2026-10-30） |
| ~~`demo/linalg/expm.cpp`~~ | 299 | 25 | ✅ 已迁（另一 agent，2026-10-30） |
| `demo/linalg/so3.cpp` | 367 | 51 | **有一处已知盲区**：`Hat³` 恒等式抓不到整体翻号，见文件头 |
| `demo/linalg/dual_complex.cpp` | 319 | 22 | 双域（实/复 × Dual） |
| `demo/linalg/smoke.cpp` | 353 | 56 | **特殊**：`[A]` 段用「取函数指针、刻意不调用」强制实例化，这个手法必须原样保留，不许改成 `TYPED_TEST` |
| `demo/autodiff/apply.cpp` | 394 | 18 | |
| `demo/autodiff/backward.cpp` | 413 | 31 | |
| `demo/autodiff/jacobian.cpp` | 361 | 18 | |
| `demo/autodiff/compute_graph.cpp` | 321 | 36 | |
| ~~`demo/core/tensor.cpp`~~ | 364 | 24 | **禁止迁移**，理由见下方「硬规则 7」 |

## 开工前必读

按顺序读，不要跳：

1. `/Users/apple/develop/learning/graphics/demo/README.md` —— **§5.7 两套框架并存**（迁移规矩）+ **§4.3.1 判据边界表**（三列：阈值 / 病态输入 / 专抓什么错）
2. `/Users/apple/develop/learning/graphics/demo/gtest_common.hpp` —— 文件头解释了「为什么判异常必须用 `EXPECT_EXACT_THROW`」
3. `/Users/apple/develop/learning/graphics/demo/quat/quaternion.cpp`（样板）
4. 你要迁的那个文件**的头部注释** —— 判据边界表、已知盲区、分域纪律都写在那里
5. `/Users/apple/develop/learning/docs/pitfalls.md` §九（注入对照的两个隐藏陷阱）+ §十二（工具链与脚本的坑）

## 硬规则（违反任何一条即视为未完成）

1. **一个文件只准用一种框架。** 迁完必须改该子目录 `CMakeLists.txt` 的注册行。
   两种混在同一文件里 = 两套输出 + `-d` 统计失真。

   **注册行现在有三种形态，都要认**（2026-10-30 实测：`linalg/` 里三种并存）：

   | 现状 | 你要做的 |
   |---|---|
   | `graphics_add_legacy_test(<名>)` | 换成 `graphics_add_gtest(<名>)` |
   | **原始四行块**（`add_executable` + `target_include_directories` + `add_test`） | **整块删掉**，换成一行 `graphics_add_gtest(<名>)` |
   | 已经是 `graphics_add_gtest(<名>)` | 不用动（可能已被别的 agent 迁过，见下） |

   **⚠ 编辑 `CMakeLists.txt` 前必须重新读一遍当前内容。** 这个仓可能同时有多个 agent 在迁移
   不同文件，而 `linalg/CMakeLists.txt` 是所有人共用的一个文件 ——
   凭记忆改会覆盖掉别人刚加的行。改完 `git diff <该文件>` 确认只动了你自己那一块。
2. **判据一条都不许删、不许合并、不许放宽阈值。** 原文件写 `1e-6` 就还是 `1e-6`。
   如果你认为某个阈值不对，**停下来写进报告**，不要自己改。
3. **判据编号必须保留**（`[A1]`、`[G3]`…）。用 `SCOPED_TRACE("[A1] …")` 携带，
   且**标签与断言必须在同一个块内**：
   ```cpp
   { STEP("[C1] i*j = +k"); EXPECT_TRUE(SameQ(i * j, k, tol)); }   // 对
   STEP("[C1] …"); EXPECT_TRUE(...);                              // 错：标签会累积
   ```
   理由：`SCOPED_TRACE` 的生命周期是整个块，不是下一条语句（`docs/pitfalls.md` §十二）。
4. **同一个模板要跑多个域就用 `TYPED_TEST`**；只有单一类型才用 `TEST`。
   老写法里 `TestXxx<double>(ck, "double")` + `TestXxx<float>(ck, "float ")` 这种手写两遍的，
   一律换成 `TYPED_TEST_SUITE`。
5. **异常判据只准用 `EXPECT_EXACT_THROW`**，禁止 `EXPECT_THROW(stmt, std::logic_error)`。
   `invalid_argument` / `out_of_range` 都派生自 `logic_error`，判基类会假绿 ——
   契约判据的全部价值就在「抛哪一种」。
6. **不许 `FetchContent`、不许 `gtest_discover_tests`。**
   依赖由 `graphics/scripts/fetch_deps.sh` 预下载（已 pin sha256）；
   一个可执行文件必须仍是 1 个 ctest 用例（`gtest_discover_tests` 会把 17 项拆成几百项，
   `graphics/test-map.md` 的映射表立刻失效）。
7. **`tensor.cpp` 禁止单独迁移。** 它的 24 条 `[P]` 是「断言实现仍 throw」的占位判据，
   先迁到 gtest 再实现等于把占位写两遍 —— 它的迁移必须和 `tensor.hpp` 的实现合并做。
8. **不许改 `graphics/src/`。** 迁移是纯测试层工作。如果你在迁移中发现
   实现有 bug 或契约与实现不符，**停下来**：写出「哪条判据、期望抛什么、实际抛什么」，
   交给用户决定，不要顺手改实现。

## 机械映射表

| 老写法 | 新写法 |
|---|---|
| `ck.check(cond, "名字")` | `EXPECT_TRUE(cond) << "名字";` 或 `CHECK_TRUE("名字", cond)` |
| `check(name, got, want, eps)`（手写 printf） | `EXPECT_NEAR(got, want, eps) << name;` |
| `check_true(name, cond)` | `EXPECT_TRUE(cond) << name;` |
| `ThrowsLogicError(f)` | （仅占位期用）`EXPECT_THROW(f(), std::logic_error)` |
| `ThrowsInvalidArgument(f)` | `EXPECT_EXACT_THROW(f(), std::invalid_argument)` |
| `ThrowsDomainError(f)` | `EXPECT_EXACT_THROW(f(), std::domain_error)` |
| `MatDiff(A,B) <= tol` | `EXPECT_MATRIX_APPROX(A, B, tol)`（在 `gtest_common.hpp`）或 `EXPECT_LE(MatDiff(A,B), tol)` |
| `std::printf("\n-- [C] … --\n")` 分节 | 拆成独立的 `TEST`/`TYPED_TEST`，节名进测试名 |
| `main()` + `CheckReport` | **删掉**，`GTest::gtest_main` 提供 main |
| `#include "../common.hpp"` | `#include "gtest_common.hpp"`（它已经 include 了 `common.hpp`，量度函数照旧可用） |

参考实现见 `graphics/demo/gtest_common.hpp` 里的 `CHECK_TRUE` / `EXPECT_EXACT_THROW` /
`EXPECT_MATRIX_APPROX` 宏定义。

## 每个文件的流程（严格按顺序，不许跳步）

1. **记录基线**：`./graphics/run_tests.sh -d`，抄下该文件的断言数（用上面的表核对）。
2. **通读**目标文件头部注释，把「判据边界表」和「已知盲区」记下来 —— 第 8 步要用。
3. **改写**，按上面的映射表与硬规则。
4. **改 CMakeLists**（先重读一遍当前内容，别凭记忆）：把该用例的注册改/换成
   `graphics_add_gtest(<名>)`，形态见硬规则 1 的三种情况表。改完 `git diff` 确认
   只动了你这一块，没有碰别人刚加的行。
5. **构建 + 跑**：`./graphics/run_tests.sh`，必须 **17/17**。
6. **开告警复查**：`./graphics/run_tests.sh -w -c`，必须 **0 告警**。
7. **核对判据没丢**：数新文件里的 `EXPECT_*` / `CHECK_*` 条数，**必须 >= 基线断言数**
   （允许变多，不允许变少；变少就说明合并了判据）。
8. **注入对照（本步不许省）**：按该文件头注释里记的已知盲区，往 `src/` 里**故意注入**那个 bug，
   跑该用例，确认**仍然是红的**；然后 `git checkout -- graphics/src/` 还原，再跑一次确认 17/17。
   - 注入期间改了 `src/` 是**唯一允许**的例外，还原必须干净（用 `git status` 确认）。
   - **新框架抓不住旧 bug ⇒ 就地回退这次迁移**，不要继续下一个文件。
9. **覆盖率看一眼**：`./graphics/run_tests.sh --coverage`，找到你迁的这个用例相关的未覆盖行，
   确认它们能归入三类之一（未实例化 / 契约性 throw / 死代码）。归不进去的，写进报告。

## 完成判据（三条同时满足才算完）

1. `./graphics/run_tests.sh` → **17/17**；`-w -c` → **0 告警**
2. 断言数**不比迁移前少**，判据编号**全部保留**
3. **注入对照重跑过**，且报告里写明「注入了什么 / 哪几条判据变红 / 是否还原」

## 交付格式

每完成一个文件，在 `docs/logs/2026-10.md`（或当月文件）**最上面**追加一条，写清：

- 迁了哪个文件、迁移前后断言数
- 做了哪个注入对照、结果（哪几条判据变红）
- 有没有发现问题（阈值可疑 / 契约与实现不符 / 未解释行）
- 提交号

提交信息用中文短句说明**做了什么**，例如：
`迁移 linalg/cholesky 到 gtest（98 条断言保留，注入 SPD 档漏判已验）`

## git 纪律（2026-10-30 起沙箱不再拦截危险命令，必须自律）

> 详细禁令与事故档案 → `/Users/apple/develop/learning/docs/env-git.md` §三。

1. **一律 `git add <具体路径>`，禁止 `git add -A` / `git add .`。**
   同一个仓里可能有别的 agent 在改别的文件，一次 `-A` 就会把它半途的改动一起提交。
2. **提交前 `git status --short` 复核。看到不认识的改动就停下报告**，不要擅自提交它。
3. **每个文件一个提交**，提交信息用中文短句说明做了什么。
4. **`git clean` 永远不加 `-x`。** 父仓 `.gitignore` 里登记了 6 个独立子仓
   （`learning-todo/`、`learning-viz/`、`learning-web/`、`learning-benchmark/`、
   `11408-notebook/`、`learning-408/`）+ `graphics/third_party/` + `.venv/`，
   `-x` 会把它们连同各自的历史一起删掉，父仓无法恢复。
5. **还原注入只准按文件**：`git checkout -- graphics/src/<具体文件>`。
   **禁止** `git checkout .` / `git restore .` / `git reset --hard` ——
   那会把你自己的未提交改动一起冲掉，而 **`reflog` 救不了未提交的改动**。
6. **不要 `git push`。**

## 明确不要做的事

- 不要顺手「优化」判据：不合并、不减少迭代次数、不把 `1e-6` 调成 `1e-5`。
- 不要给 `src/` 补注释（那是另一条线的工作，会撞车）。
- 不要改 `graphics/test-map.md` 以外的文档来「解释」你的改动 —— 该文件 §4.3 有一列「框架」，
  迁移完成后把对应行改成 `gtest` 即可。
- 不要一次迁移超过 2 个文件。注入对照要一个个做，批量做必然漏。
