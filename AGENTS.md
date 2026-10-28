# AGENTS.md — 项目协作规约（每个 session 自动加载）

> **读者**：编码 agent ｜ **最后更新**：2026-10-26
> **本文是入口文件**。展开内容按下方「文档地图」**按需读取，不要一次性全读**（`learning.md` 约 34 KB）。

---

## 一、文档地图

| 文件 | 内容 | 什么时候读 |
|---|---|---|
| `NOW.md` | **近期要做什么**（全仓唯一滚动文件，1–3 周一版，过期整篇重写） | **每次会话开始时** |
| `dev-progress.md` | **进度总览**（§0 现状表 + §2 Gate + §5 里程碑）。**不含日志正文** | 每次会话开始时 |
| `docs/logs/<YYYY-MM>.md` | **进度日志**（按月一个文件，append-only，最新在上） | 想知道"哪天做过什么"时，按月份查 |
| `docs/logs/dev-progress-status-2026-10-15.md` | 2026-10-15 的现状/下一步快照（**已过期，仅供考古**） | 只在追 V1.1–V1.7 实测数字、S1–S6b 提交序时 |
| `README.md` | 人读总纲：仓库职责、目录关系、样板 | 需要理解项目意图，或改动跨仓时 |
| `docs/route-math.md` | 数学**地图**（十段、依赖顺序、每段"看懂线"、不做的边界、Lee 读法） | 讲数学、判断先后顺序时 |
| `learning.md` | 数学**手推题库**（12 阶段经典问题）。**只是题，不含路线** | 讲数学题时，**只读相关阶段** |
| `graphics/acceptance.md` | **验收判据规范**（怎么算通过；含已退役 `V*` 编号对照） | 写/改任何 ctest 用例前 |
| `graphics/numeric-contract.md` | **数值契约**（`D7`/`D12`/`D16`/`D19` 完整内容：除零两层、容差体系、错误类型三分、布局形状）。接管已删除 `matrix.md` 的契约部分 | 写任何 `src/` 数值代码前；看到源码里裸 `D12`/`§10.8` 引用时 |
| `graphics/test-map.md` | **接口 ↔ 测试可追溯图**（40 个公开入口 → 权威用例 → 判据编号；含已知断点） | 想知道"某接口被什么测了"、新增入口/用例、或怀疑覆盖有洞时 |
| `docs/pitfalls.md` | **口径与踩坑**（判据病态、共轴掩盖、holonomy 口径、泛型域、实例化才叫落地） | 判据出怪数、或怀疑"实现错"时 |
| `docs/eng-competency.md` | **工程能力 60/80 基线**（四赛道检查点、赛道重叠矩阵、本仓现状对照）。**只管"能做出来"** | 评估"这个件做到什么程度算过关"，或对岗位要求时 |
| `docs/baseline-levels.md` | **四线 × 三档（60/80/90）× 两把尺子**（学习 vs 工程）。覆盖 **数学 / 数学物理 / CG / CS**；含**反例判据（§0.6 一票否决）**与**验证机制（§0.7）** | 自我评估"我到哪一档了"、区分"读懂"与"做出来"、找**降档信号**时 |
| `docs/interview-map.md` | **三方向 × 学科线的档位要求**（C++ 后端 / 算法测开 / CG）。回答"为了这个方向的面试，每条线要到哪一档"；含 Raft 等具体资产的现状与缺口 | 准备某个方向的面试、判断"这块要不要补"时 |
| `graphics/demo/pyapps/tensor-convention.md` | pyapps 层张量记号约定（轴序 = 指标序、上下标进变量名、einsum 即文档） | 写 pyapps 张量脚本时 |
| `docs/env-git.md` | git 写操作被沙箱拒绝的成因与处理（锁文件清单、worktree 补法） | git 报 `Operation not permitted` / `File exists` 时 |
| `docs/doc-checkup.md` | 文档体检清单（低频，每月一次） | 大改文档后 |
| `learning-web/AGENTS.md` | 笔记 web 化的**入口规约**（真源约定、四问判据、硬规则、红线）｜ 在子仓内，只能按绝对路径访问 | 动 `learning-web` 任何东西之前 |
| `learning-web/DESIGN.md` | 笔记 web 化的架构与契约（决策 D1–D8、数据流、契约样例） | 动契约 / 渲染器 / `build/` 之前 |
| `11408-notebook/LEARNING.md` | **11408 考研学习引擎**（冲刺日历、刷题自主权、判卷/遗忘/归档工作流） | 会话涉及考研学习时，学哪个模块读对应 `modules/` |
| `11408-notebook/PROGRESS.md` | 11408 学习进度仪表盘（唯一进度真源） | 11408 学习会话开场 |
| `11408-notebook/AGENTS.md` | 11408 归档规约（对话/照片 → ipynb，只追加） | 归档或判卷照片时 |
| `learning-408/README.md` | 408 实践代码导航（**主题 → 代码映射真源**）｜ 独立子仓，只能按绝对路径访问 | 复习 408 要引实践代码时，**先读它再引代码，不得凭记忆编造路径** |

> **已删除（2026-10-26）** —— 下面这些路径**已不存在**，内容都已迁走，别再按旧路径找：
>
> | 已删除 | 内容现在在哪 |
> |---|---|
> | `learning-todo/matrix.md`、`docs/matrix.md` | **已彻底删除**。三处接管：红线 → `acceptance.md` §0；**数值契约（除零/容差/错误类型/布局）→ `graphics/numeric-contract.md`**；张量契约 → `tensor.hpp` 头注释。**源码里的 `matrix.md §N` 是历史出处，不是可跳转链接**，解码表见 `numeric-contract.md` §9 |
> | `learning-todo/route-math.md` | `docs/route-math.md` |
> | `learning-todo/interview-depth.md` | `docs/interview-depth.md` |
> | `learning-todo/ns-progress.md` | `docs/math-log.md` |
> | `learning-todo/{math-plan,fan-plan,lie-plan}.md` | `docs/route-math.md`（口径警告在 `docs/pitfalls.md`） |
> | `learning-todo/stage1-wrapup.md` | `NOW.md` + `docs/interview-depth.md` |
> | `graphics/plan.md` | 红线 → `graphics/acceptance.md` §0；`V*` 编号对照 → §4 |
>
> **`docs/logs/` 与 `docs/math-log.md` 的历史条目仍按这些旧路径引用** —— 那些引用现在指向
> 不存在的文件，**属于历史记录，不要修**（改历史等于造假）。

**真源约定（避免重复维护）**：同一条事实只在一个文件里定义，别处只引用不重述。
各主题的真源已在上表「内容」列标出，此处不再逐条复述。只额外强调三条：
**引用 408 实践代码前必读 `learning-408/README.md`**（不许凭记忆编造路径）；
**git 环境坑只写在 `docs/env-git.md`**；**判据与踩坑只写在 `graphics/acceptance.md` + `docs/pitfalls.md`**。

### 文档分三类（2026-10-26 定 —— 这是防腐化条款）

用户原话：「大阶段的 plan 没有意义，我也不会按着他来。每段时间都要整理一次文档，没有意义。」
据此把文档按**是否会腐化**分三类，越线即违规：

| 类 | 文件 | 允许写什么 |
|---|---|---|
| **滚动** | `NOW.md`（**唯一**） | 可以写"接下来做什么"。1–3 周一版，过期**整篇重写**，不追加、不归档、不留历史 |
| **契约/参考** | `graphics/acceptance.md`、`docs/pitfalls.md`、`docs/route-math.md`、`docs/interview-depth.md`、`tensor-convention.md`、`docs/*.md`、`README.md` | **只写"是什么"。禁止写：下一步、阶段编号、月份、时间表、进度、验收关卡** |
| **历史** | `docs/logs/<YYYY-MM>.md`（按月的进度日志）、`docs/math-log.md`（数学线日志）、`math-notebook/*.tex`、`learning-todo/imgs/` | append-only，可写日期与既成事实 |

> **产物不进父仓**：二进制/大文件（PNG、会话 ipynb）→ `learning-todo/`（独立仓）。
> 判断标准：**能进 git 的文本 → 父仓；二进制/大产物 → learning-todo。**

**禁止事项（硬规则）**：

1. **不得新建任何大阶段计划** —— 「阶段 1–6」「八刀路线」「扇面切片清单」「六阶段读书表」
   这类东西一律不再产生。要做的事只进 `NOW.md`，1–3 周为限。
2. **不得周期性"整理文档"** —— 整理只在两种情况下发生：触犯 `docs/doc-checkup.md` 的某一条，
   或用户明确要求。**"定期梳理一遍"本身不是理由。**
3. **旧路径不保留跳转页**（2026-10-26 用户改口径）。删/迁后 `docs/logs/` 里的历史引用会指向不存在的
   文件 —— **那是历史记录，不许修，也不必补兼容层**。
   **删前唯一必做的检查**：这个文件是否兼任**契约真源**（红线 / 判据 / 阈值 / 编号）？
   是则先把那部分原样搬走（`graphics/acceptance.md` §0）再删。
   实例：`graphics/plan.md` 看着是"阶段计划"，实际 §0 是**设计红线 #1–#5 的唯一真源**，
   已抢救进 `graphics/acceptance.md` §0。⇒ **"看起来是计划"≠"只有计划"。**

---

## 二、背景（为什么这么要求）

用户是**非数学专业出身**（有高等数学基础）在自学数学系，目标是拿到
**PDE 与计算机图形学所需的工程数学基础**。

**目标的口径（2026-10-14 用户定）**：NS 方程只是**打开整个学习路线的钥匙**，
用来选题与排序，**不代表这些基础不重要** —— 基础本身就是目标之一。
因此 agent 在排期建议里不得把「暂时用不到」当成长期降低某块基础的理由；
推迟只意味着「排到它的消费点再学」，不意味着「可以不学」。

因此默认取向是：**严谨优先于速度，直觉优先于形式化。**
特别复杂或纯技巧性的习题可以降级，但基础概念不允许含糊。

### 工程落脚点与口径

- **多仓结构与进度真源**：见 §一 文档地图 与 §七（子仓各自独立、父仓不跟踪、须绝对路径访问）。
  用户原话：这些子仓「都是外部的子仓，只是放到一起方便编程」。
- **不再提起的决定**：`graphics/` 定位偏窄这件事用户已明确决定「保持现状」——
  **不要再提**，也不要在排期建议里隐含否定它。
- **两条通用工程事故的本仓实例**（通用规则在用户级 `MEMORY.md`，实例见 `.workbuddy/memory/2026-10-17.md`）：
  ① 「读—改—写」脚本不幂等 ⇒ `graphics/demo/linalg/expm.cpp` 的 `[E]` 段被插两遍、`redefinition` 构建失败；
  ② 整文件 Write 覆盖未提交 WIP 后，靠 VSCode 本地历史恢复 `graphics/src/linalg/solve.hpp`（2026-10-14）。

---

## 三、命令

```bash
uv sync                                     # 安装/同步 Python 依赖
uv run pytest                               # Python 测试
cmake -S graphics -B graphics/build -DCMAKE_BUILD_TYPE=Debug
cmake --build graphics/build -j             # 构建 C++ 库与 demo
ctest --test-dir graphics/build --output-on-failure   # 跑 C++ 测试
./graphics/run_tests.sh                     # C++ 测试一键入口(配置+构建+ctest)
                                            #   -c 从零构建 | -d 打各用例内部明细 | -w 开 -Wall -Wextra | -h 帮助
latexmk -pdfxe -interaction=nonstopmode main.tex   # 编译笔记（在 math-notebook/ 下执行；引擎坑见 §七）
bl image generate --model qwen-image-3.0-pro --watermark false --size "2048*2048"
```

---

## 四、代号速查

项目文档里大量使用缩写引用，释义如下：

- `B1` / `G3` / `G4` / `H5` / `E1` … = `math-notebook/<X>_*.tex` 内的 `\subsection{X#.}` 编号
  （`A_analysis`→A、`B_algebra`→B、`C_partial`→C、`D_classic`→D、`E_Field`→E、
  `F_quantum`→F、`G_numeric`→G、`H_cg`→H、`I_fluid`→I、`J_rtx`→J、
  `K_Feyman`→K、`L_Feyman`→L、`M_addition`→M、`P_computation`→P、
  `Q_topology`→Q、`R_complex`→R）
- `de1`–`de4` = `math-notebook/diffrential/de*.tex`（拉格朗日/哈密顿篇的「第 N 层」）
  > ⚠️ 目录名 **`diffrential` 就是这么拼的**（少一个 `e`），是真实存在的目录，
  > **不要「修正」成 `differential`** —— 改了会直接找不到文件。
- `lie1`–`lie3` = `math-notebook/li/lie*.tex`（李群篇）
  > 历史改名（`B_algbra`→`B_algebra`、`M_addtion`→`M_addition`、计算必修 `G1`–`G9`→`P1`–`P9`）
  > 见 `docs/doc-checkup.md` §历史改名记录；见到旧名按新名理解。
- `§NN.N` = **编译后**的全局章节号（`main.pdf` 里的编号），见下方警告
- `L0`–`L5` = **已废弃**的 NS 学习路线关卡（2026-10-23 被十段表取代，十段表见 `route-math.md`）。
  旧文里见到按层理解即可，**不要再使用**。
- `V1.1`–`V6.4` = **已退役**的验收项编号（2026-10-26）。对照表见 `graphics/acceptance.md` §4，
  **只为读懂旧日志存在，不表示任何计划**。
- `Gate A/B/C` = `dev-progress.md` 的小项目启动门槛（数学与工程的**双条件**）

> **两个必须区分的东西**
> 1. **验收项写 `V*`，笔记代号写字母+数字**。历史上验收项用过 `A1`、`B1`、`E1` 这类裸字母编号，
>    与笔记代号（`B1`、`E1`…）撞名，是 agent 误读的主要来源，**已废弃**。
> 2. `§NN.N` 是全局编译编号，**插入新的 `.tex` 就会整体位移**。引用时优先写文件内编号（`B1`、`G1.1`），
>    只有必须跨篇交叉引用时才用 `§`，并注明可能漂移。

---

## 五、符号与书写作用域

| 场合 | 写法 |
|---|---|
| 聊天回复 | LaTeX 渲染：行内 `$…$`、独立 `$$…$$`（2026-10-26 用户定：行内禁 Unicode/ASCII 混排公式，如 `J^i_j`）；不输出裸 LaTeX 源码 |
| 数学内容文档：`.ipynb`、`learning.md` | KaTeX：行内 `$…$`、行间 `$$…$$` |
| 工程文档：`graphics/acceptance.md`、`docs/pitfalls.md`、`dev-progress.md`、`NOW.md`、`AGENTS.md`、`README.md`、`docs/*.md` | ASCII：`A^T`（转置）、`A^H`（共轭转置）、`[A,B] = AB − BA`（对易子） |
| C++ 代码与注释 | 同上一行（ASCII） |

**工程文档的统一细则（2026-09-22 起生效）**

- **上标一律用 `^`**，不用 Unicode 上标字符：写 `A^T`、`A^H`、`A^-1`、`n^2`、`h^4`、`O(n^2)`、`|ψ|^2`、`χ^2`，
  不写 `Aᵀ`、`Aᴴ`、`A⁻¹`、`n²`、`h⁴`、`O(n²)`、`|ψ|²`、`χ²`；`½` 写 `(1/2)`。
- 其他符号（`∫ Σ ∇ Λ Φ θ ω μ π δ ε ∞ ≈ ≤ ≥ ⟨⟩ ‖·‖ · −`）保持现状 —— 它们不是"两套写法"，是标准数学排版字符。
- **历史条目不动**：`ns-progress.md` 的会话日志是 append-only，旧条目里的 Unicode 写法保留原样；
  **新写入的内容**按上面的细则来。

> 判定规则很简单：**给 agent / 终端读的用 ASCII，给 Jupyter 渲染的用 KaTeX。**
> 同一份文档内不要混用两套写法。

---

## 六、硬规则

### 数学

- **正确性第一。** 一个符号或函数写错就会浪费大量时间。**不确定就明确标注「待核」，不要含糊过去。**
- 先给直觉，再上定理；英文定理一律用括号附中文；希腊/拉丁字母附读法；集合类符号附解释。
- 复杂表达式必须拆解注解，不能只丢一个式子。
- 明确标出「谁是谁的基础」「谁解决什么问题」「历史上是被什么逼出来的」。
- 纯技巧性或过于严苛的习题可降级：**每节最多保留 1 道硬题**。

### 输出

- terminal 输出当课堂讲义用；`.ipynb` 当黑板与笔记用。
- **会话 ipynb 不再归档**（`learning-todo/chats/` 已于 2026-10-26 删除）。需要留档时，
  由用户当次指定路径，不要自建目录。
- **生成产物**一律进 `learning-todo/imgs/`（按主题分子目录，最多两级）。
- 只有明确要求写代码时才输出代码，不要附带多余的验证性代码。

### 图像生成（阿里云百炼 CLI）

- 目标是**自学导读 / 知识图谱**，不是海报：知识密度要高、字体不要大、不要醒目标题。
- 一律使用 `bl image generate --model qwen-image-3.0-pro`（走 Token Plan 订阅）。
  `qwen-image-3.0` 非 pro 版只扣免费额度且余额少，**不要默认使用**。
- 产物保存到 `learning-todo/` 及子目录（图片默认 `learning-todo/imgs/`，可按主题分类，**最多两级子目录**）。
- 尺寸参数含 `*` 时必须加引号（zsh 通配符）：`--size "2048*2048"`。
- 生成后**必须读图校验**数学内容（公式、坐标轴刻度、箭头指向）；有错就重出，或明确告知用户瑕疵。
- AI 生成图的坐标轴刻度数字经常乱码 → **讲义级精确配图改用 matplotlib 代码绘制**。

### 编程

- 介绍语言特性或设计写法时，**先讲它背后解决了什么问题**，再讲语法。
- C++ 使用 C++20；矩阵与线代的既定设计以**各头文件顶部注释** + `graphics/acceptance.md` §0 为准，不要擅自改动（红线见第九节）。

---

## 七、环境与已知坑

- **Python**：项目用 uv 管理（`.venv` + `pyproject.toml` + `uv.lock`，Python 3.14，见 `.python-version`）。
  跑命令用 `uv run xxx`，装包用 `uv add xxx`。
  **禁止**使用系统 `python3`（`/usr/bin` 的 3.9.6）或裸 `pip install`。
  全局 `python3` 是 uv 托管的 3.14（`~/.local/bin`）；Homebrew python 已卸载，**不要再装回**。
- **四子仓**：`learning-todo/`、`learning-viz/`、`learning-web/`、`learning-408/` 是**独立 git 仓**，各有自己的 remote，
  父仓不跟踪（无 `.gitmodules`、无 gitlink）。
  → 它们在 `git worktree` 与新的 clone 中**都不出现**，必须用绝对路径访问
  （如 `~/develop/learning/learning-todo/`）；**改动须在各自仓内提交**。
- **LaTeX 必须走 xelatex**：`latexmk -pdfxe`。若用 `-pdf`（走 pdfLaTeX），ctex 会报
  `fontset 'mac' unavailable` 并**卡在交互输入**等键盘，命令永远不返回（2026-09-24 实测）。
- **git 写操作会被沙箱拒绝**：`commit` / `update-index` / `branch -D` / `worktree remove` 都可能报
  `unable to unlink .../.git/index.lock: Operation not permitted`，并级联出
  `Another git process seems to be running`（看着像随机失败，其实是前一次失败的回声）。
  **处理**：用 shell `rm -f` 删掉**具体那一个**锁文件 → **立刻重跑**，循环 3–8 次基本必过。
  细节（成因、锁文件清单、worktree 补法）→ `docs/env-git.md`；skill `gitsandbox-write-workaround`。
- 全仓递归 `grep` / `ls -R` 容易被沙箱拒绝 → 改用结构化的文件检索工具。

---

## 八、目录用途

> 各仓**职责与状态**见 `README.md` §仓库职责（真源）；本节只记 agent 操作时需要的内部结构。

- `math-notebook/`：LaTeX 笔记（`.tex`，由 `main.tex` 统一 `\input`，已编译 `main.pdf`）
- `graphics/`：C++ 图形/线代引擎（`src/` 头文件库；`demo/<模块>/` 为 ctest 验收用例，
  一键入口 `run_tests.sh`；`demo/apps/` 为**附加区**——现象演示，不注册 ctest，用 `run_apps.sh` 跑；
  `test/` 是早期 C 实验与 notebook，与 ctest 无关）
- `learning-todo/`：生成产物与参考文档（`chats/`、`imgs/`、`ns-progress.md`（数学日志）、
  `route-math.md`（数学地图）、`interview-depth.md`、`*_plan.md` 退役跳转页）
- `learning-viz/`：职责**收窄为计算/生成层**（Flask `/api/viz/*` + ImGui/matplotlib 本地渲染）；
  **不再承担笔记转换**（`notes_pipeline` 已移出）
- `learning-web/`：笔记 web 化的全部落地：`build/`（转换，Python）+ `graph/`（语义边真源）
  + `src/`（Vite 外壳）+ `dist/`（产物，gitignore）
- `11408-notebook/`：11408（政治/英一/数一/408）考研学习体系，五件套分工
  （`AGENTS.md` 归档规约 + `LEARNING.md` 学习引擎 + `TOOLKIT.md` 图技能 + `modules/` 模块路线，
  真源约定见其内部文件地图）。**只追加 ipynb，不改已有单元**
- `learning-408/`：408 实践代码汇总仓。本体只保留一份**主题 → 代码导航**（`README.md`），
  三个实践子仓 `jos/`（操作系统）、`cs144/`（网络）、`csapp/`（组成）被其 `.gitignore` 忽略、各自独立管理
- `docs/`：外移文档（`env-git.md` 环境坑、`doc-checkup.md` 低频体检清单）
- **四个子仓与父仓的关系、路径坑**见第七节；依赖约定：**子仓都不建自己的 `pyproject.toml`**，
  `uv add` 会写进**父仓** `pyproject.toml` + `uv.lock`

---

## 九、不要修改（Do Not Modify）

- **四个子仓不属于父仓**：`learning-todo/`、`learning-viz/`、`learning-web/`、`learning-408/` 的改动只能在各自仓内提交。
- **`Matrix<T>` 的既定设计**（哑容器、不加缓存层、不把算法塞进类、实复差异只走 traits）：
  这些已在**各头文件顶部注释**里定稿。需要改设计，**先改注释或 `acceptance.md` §0 并写清理由**，再动代码。
- **已写好的 `.ipynb` 单元**：只追加新单元，不重写旧单元。
- **`math-notebook/` 下的笔记**：非明确指示，不要修改。

---

## 十、提交

- **未经明确要求不要 `git push`**；在工作区内可以正常 `git add` / `git commit`。
- 提交信息用中文短句，说明**做了什么**，而不是罗列改了哪些文件。

---

## 十一、安全

- 不在聊天或文件中回显完整 API Key，也不要输出 `.env` 的内容。
- `bl` 已配置好鉴权，**无需**在命令里传 `--api-key`。

---

## 十二、文档体检

**低频**（每次大改文档后，或每月一次），不必每个 session 做。7 类历史漂移、对应改法、
以及历史改名记录（`B_algbra`、`G1`–`G9`、`notebook/` 等）→ `docs/doc-checkup.md`。
