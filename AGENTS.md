# AGENTS.md — 项目协作规约（每个 session 自动加载）

> **读者**：编码 agent ｜ **最后更新**：2026-09-22
> **本文是入口文件**。展开内容按下方「文档地图」**按需读取，不要一次性全读**（`learning.md` 有 43 KB）。

---

## 一、文档地图

| 文件 | 内容 | 什么时候读 |
|---|---|---|
| `README.md` | 人读总纲：仓库职责、目录关系、样板 | 需要理解项目意图，或改动跨仓时 |
| `dev-progress.md` | **进度的唯一真源** + Gate A/B/C 门槛 | 每次会话开始时 |
| `graphics/plan.md` | C++ 工程任务与验收阈值（阶段 1–6） | 做图形/线代工程任务时 |
| `learning.md` | 数学手推题库（12 阶段经典问题 + 读书路线） | 讲数学题时，**只读相关阶段** |
| `learning-todo/ns-progress.md` | NS 学习路线进度（L0–L5） | 推进 NS 路线时 |
| `learning-todo/matrix.md` | Matrix 库设计的**理由**（决策与权衡） | 动 `graphics/src` 矩阵/线代代码之前 |
| `learning-web/AGENTS.md` | 笔记 web 化的**入口规约**（真源约定、四问判据、硬规则、红线）｜ 在子仓内，只能按绝对路径访问 | 动 `learning-web` 任何东西之前 |
| `learning-web/DESIGN.md` | 笔记 web 化的架构与契约（决策 D1–D8、数据流、契约样例、A1 切片验收） | 动契约 / 渲染器 / `build/` 之前 |

**真源约定（避免重复维护）**：同一条事实只在一个文件里定义，别处只引用不重述。
进度 → `dev-progress.md`；工程验收 → `graphics/plan.md`；库设计理由 → `matrix.md`；
笔记 web 化 → `learning-web` 子仓（规约 `AGENTS.md`、设计 `DESIGN.md`）。

---

## 二、背景（为什么这么要求）

用户是**非数学专业出身**（有高等数学基础）在自学数学系，目标是拿到
**PDE 与计算机图形学所需的工程数学基础**。

因此默认取向是：**严谨优先于速度，直觉优先于形式化。**
特别复杂或纯技巧性的习题可以降级，但基础概念不允许含糊。

---

## 三、命令

```bash
uv sync                                     # 安装/同步 Python 依赖
uv run pytest                               # Python 测试
cmake -S graphics -B graphics/build -DCMAKE_BUILD_TYPE=Debug
cmake --build graphics/build -j             # 构建 C++ 库与 demo
ctest --test-dir graphics/build --output-on-failure   # 跑 C++ 测试
latexmk -pdf main.tex                       # 编译笔记（在 notebook/ 下执行）
bl image generate --model qwen-image-3.0-pro --watermark false --size "2048*2048"
```

---

## 四、代号速查

项目文档里大量使用缩写引用，释义如下：

- `B1` / `G3` / `G4` / `H5` / `E1` … = `notebook/<X>_*.tex` 内的 `\subsection{X#.}` 编号
  （`A_analysis`→A、`B_algebra`→B、`C_partial`→C、`D_classic`→D、`E_Field`→E、
  `F_quantum`→F、`G_numeric`→G、`H_cg`→H、`I_fluid`→I、`J_rtx`→J、
  `K_Feyman`→K、`L_Feyman`→L、`M_addition`→M、`P_computation`→P、
  `Q_topology`→Q、`R_complex`→R）
- `de1`–`de4` = `notebook/diffrential/de*.tex`（拉格朗日/哈密顿篇的「第 N 层」）
- `lie1`–`lie3` = `notebook/li/lie*.tex`（李群篇）
  > 注：`B_algbra.tex`、`M_addtion.tex` 原为拼写错误，**已于 2026-09-22 改名为
  > `B_algebra.tex`、`M_addition.tex`**（`main.tex` 的 `\input` 已同步）。见到旧名按新名理解。
  > 同类改名还有一处：**「计算必修清单」的九组原记 `G1`–`G9`，因与 `G_numeric.tex`
  > 的 `G1`–`G5` 撞名，已统一改记为 `P1`–`P9`**（黑板里的旧编号 `G$n$` = 笔记里的 `P$n$`）。
- `§NN.N` = **编译后**的全局章节号（`main.pdf` 里的编号），见下方警告
- `L0`–`L5` = NS 学习路线关卡（`learning-todo/ns-progress.md`）
- `V1.1`–`V6.4` = `graphics/plan.md` 的验收项，记作 `V{阶段}.{序号}`
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
| 聊天回复 | Unicode（∑ ∫ √ λ ᵀ）——**不输出裸 LaTeX** |
| 数学内容文档：`.ipynb`、`learning.md` | KaTeX：行内 `$…$`、行间 `$$…$$` |
| 工程文档：`matrix.md`、`graphics/plan.md`、`dev-progress.md`、`AGENTS.md`、`README.md` | ASCII：`A^T`（转置）、`A^H`（共轭转置）、`[A,B] = AB − BA`（对易子） |
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
- 正式输出写入 `learning-todo/chats/<session-id>.ipynb`：**每次回复追加新的 markdown 单元，
  不修改已经写好的单元**（禁止大规模重写历史内容）。
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
- C++ 使用 C++20；矩阵与线代的既定设计以 `matrix.md` 为准，不要擅自改动（红线见第九节）。

---

## 七、环境与已知坑

- **Python**：项目用 uv 管理（`.venv` + `pyproject.toml` + `uv.lock`，Python 3.14，见 `.python-version`）。
  跑命令用 `uv run xxx`，装包用 `uv add xxx`。
  **禁止**使用系统 `python3`（`/usr/bin` 的 3.9.6）或裸 `pip install`。
  全局 `python3` 是 uv 托管的 3.14（`~/.local/bin`）；Homebrew python 已卸载，**不要再装回**。
- **三子仓**：`learning-todo/`、`learning-viz/`、`learning-web/` 是**独立 git 仓**，各有自己的 remote，
  父仓不跟踪（无 `.gitmodules`、无 gitlink）。
  → 它们在 `git worktree` 与新的 clone 中**都不出现**，必须用绝对路径访问
  （如 `~/develop/learning/learning-todo/`）；**改动须在各自仓内提交**。
- **git 的写操作会被执行层拒绝 unlink**，典型输出是
  `warning: unable to unlink '.../.git/index.lock': Operation not permitted` 加 `fatal: Unable to write new index file`。
  - **不是仓库问题**：`.git` 没有 `uchg/schg` 标志、没有 ACL、属主可写，
    且锁文件的 mode/owner 与 shell 新建的文件完全一致；同一路径下 shell 的 `rm` / `mv`
    （含覆盖已存在文件）都成功 —— **只有 `git` 进程被拒**（可复现）。
  - **也与"是否免沙箱"无关**：在标记为 `Sandbox bypassed` 的运行里同样出现该报错。
    环境里存在 `CODEBUDDY_SANDBOX_BROKER_TRACE_ID`，但**具体是哪一层拦的、规则是什么，
    从进程内看不到** —— 不要臆断成因，直接按下面的办法处理。
  - **真正的坑是级联**：失败会留下 `*.lock`，之后所有 git 命令报
    `Unable to create ...: File exists` / `Another git process seems to be running`——
    看起来像随机失败，其实是前一次失败的回声。
  - **实测有效的处理**：用 shell `rm -f` 删掉**具体那几个**锁文件
    （`.git/index.lock`、`.git/packed-refs.lock`、`.git/refs/heads/**/*.lock`、`.git/ORIG_HEAD.lock`），
    然后**立刻重跑**该 git 命令；循环 3–8 次基本必过（通常第 1 次就过）。
    **不要用通配符删** —— zsh 在无匹配时会中断整条命令。
  - `GIT_OPTIONAL_LOCKS=0` 只对纯读命令（`git status`、`git diff`）有意义，
    对 `commit` / `update-index` / `branch -D` 无效。
  - **删 worktree 会撞同一个坑**：`git worktree remove` 可能已经删掉了工作目录，却卡在
    `.git/worktrees/<name>` 上。补法：手动 `rm -rf .git/worktrees/<name>` → `git worktree prune`
    → 再 `git branch -D <branch>`。
- 全仓递归 `grep` / `ls -R` 容易被沙箱拒绝 → 改用结构化的文件检索工具。

---

## 八、目录用途

- `notebook/`：LaTeX 笔记（`.tex`，由 `main.tex` 统一 `\input`，已编译 `main.pdf`）
- `graphics/`：C++ 图形/线代引擎（`src/`、`demo/` 为 ctest 测试程序、`test/`）
- `linear_algebra/`：线代笔记（`.ipynb`）
- `learning-todo/`：**独立子仓**，生成产物与待办输出（`chats/`、`imgs/`、`ns-progress.md`、`matrix.md`）
- `learning-viz/`：**独立子仓**，2026-09-22 起有提交。职责**收窄为计算/生成层**（Flask `/api/viz/*` +
  ImGui/matplotlib 本地渲染）；**不再承担笔记转换**（`notes_pipeline` 已移出）
- `learning-web/`：**独立子仓**，2026-09-22 起有提交。笔记 web 化的全部落地：
  `build/`（转换，Python）+ `graph/`（语义边真源）+ `src/`（Vite 外壳）+ `dist/`（产物，gitignore）
- 三者的关系与路径坑见第七节；依赖约定：两个子仓都不建自己的 `pyproject.toml`，
  `uv add` 会写进**父仓** `pyproject.toml` + `uv.lock`

---

## 九、不要修改（Do Not Modify）

- **三个子仓不属于父仓**：`learning-todo/`、`learning-viz/`、`learning-web/` 的改动只能在各自仓内提交。
- **`Matrix<T>` 的既定设计**（哑容器、不加缓存层、不把算法塞进类、实复差异只走 traits）：
  这些已在 `matrix.md` 里定稿并给出理由。需要改设计，**先改 `matrix.md` 并写清理由**，再动代码。
- **已写好的 `.ipynb` 单元**：只追加新单元，不重写旧单元。

---

## 十、输出与提交

- 生成产物一律进 `learning-todo/`：图片进 `imgs/`，会话记录进 `chats/`。
- **未经明确要求不要 `git push`**；在工作区内可以正常 `git add` / `git commit`。
- 提交信息用中文短句，说明**做了什么**，而不是罗列改了哪些文件。

---

## 十一、安全

- 不在聊天或文件中回显完整 API Key，也不要输出 `.env` 的内容。
- `bl` 已配置好鉴权，**无需**在命令里传 `--api-key`。
