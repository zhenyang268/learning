# dev-progress.md — 进度与排期

> **读者**：人 + agent ｜ **最后更新**：2026-09-22 ｜ **本文是「进度」的唯一真源**
>
> **真源声明**
> - 本文是**总纲式进度真源**：回答「整体现在在哪、下一步做什么」。
> - `graphics/plan.md` 只负责任务定义与验收阈值（验收项记作 `V{阶段}.{序号}`），**不记录状态**。
> - **数学线的关卡明细，真源是 `learning-todo/ns-progress.md`**（逐次会话日志 + 关卡状态）。
>   本文只保留一行指针，**不重复记录**，以免两处漂移。
> - 冲突处理：工程线以本文为准；数学线以 `ns-progress.md` 为准。
>
> 分工：AI 排期 / 定门槛 / 验收，用户写代码。优先学习，工程增量推进。

---

## 0. 当前指针（每次会话先看这里）

> **数学线在 L1（测度与积分）中段；工程线停在 `graphics` 阶段 1 的起点，尚未开工。**

| 线 | 当前位置 | 下一个动作 |
|---|---|---|
| **数学（优先）** | L1 第 3 关进行中 —— **明细见 `learning-todo/ns-progress.md` 的「📍 当前位置」**（本表不重复记录） | Vitali 不可测集 → Carathéodory 判据 → 积分四步构造 |
| **工程** | `graphics` 阶段 1 未动：matrix 未做实数/复数 scalar trait 改造，无 `linalg`，CMake 仍是 C++17 | CMake 升 C++20；`core/scalar.hpp` + `core/matrix.hpp` 重写 |
| **小项目** | 未解锁 | 等 Gate A（见第 2 节） |
| **笔记 web 化** | 两仓骨架已建（`learning-viz` / `learning-web` 各一个 init 提交），A1 垂直切片未开工 | 按 `learning-viz/python/notes_pipeline/DESIGN.md` §12 第 1 步：`uv add --dev pypandoc-binary` |

> **2026-09-22 修正**：本表原先写「L1 第 3 关 —— 外测度已学、勒贝格积分刚学」，
> 与 `ns-progress.md` 的「外测度 / Carathéodory **未接触**」**互相矛盾**。
> 已改为**只留指针不做重复记录** —— 数学线的关卡状态一律以 `ns-progress.md` 为准。
> 同日用户确认：**外测度确已学**，滞后的是 `ns-progress.md` 一侧，该文件已一并更新。

**更新方式**：每次会话用户报「学到哪 / 写到哪」→ 改上表 + 在第 6 节日志追加一条。
上表只保留**当前**状态（覆盖式），历史一律进日志（追加式，不删改）。

---

## 1. 协作模式

- 数学线（慢，优先）与工程线（快，够用即做）并行。
- 每次会话：用户报「学到哪 / 写到哪」 → AI 更新本文件并给下一步。
- **matrix 改造**与**自动微分**由用户自行抽时间完成，AI 只给验收点。

---

## 2. 小项目启动门槛（Gate）

| 门槛 | 数学 | 工程 | 解锁的小项目 |
|---|---|---|---|
| **Gate A** | L0 + ODE 数值解直觉（Euler/RK4/Verlet） | matrix 改造 ✔ + 最小 linalg（LU/solve）+ 一个积分器 | 动力系统可视化：单摆 / 行星轨道 / 双摆 / N 体 |
| **Gate B** | **L1 完成**（测度 + 积分 + 收敛定理 + 四步构造） | Gate A + FDM 基础 | PDE 可视化：热方程（图像模糊）、波动（鼓膜） |
| **Gate C** | L2 泛函 + Sobolev | 阶段 4 完整 | 弱解类（通往 NS 真目标） |

> **第一个小项目定在 Gate A**；L1 继续按 NS 路线推进，不互相阻塞。

---

## 3. 现状明细（截至 2026-09-21）

- **数学线**：见 `learning-todo/ns-progress.md`（据 2026-09-22 用户确认：L1 第 3 关进行中，
  已过 σ-代数 / 测度性质 / 三大收敛定理 / **外测度**；差 Vitali → Carathéodory → 积分四步构造）。
- **工程线**：`graphics` 阶段 1 未动 —— matrix 未做实数/复数 scalar trait 改造；无 `linalg`；CMake 仍 C++17。

---

## 4. 下一步

**数学线（优先）**

- [ ] 外测度（方盒覆盖、次可加）
- [ ] Vitali 不可测集
- [ ] Carathéodory 判据 → 定理
- [ ] Lebesgue 积分四步构造
- [ ]（顺带）Vitali 收敛定理

**工程线**（用户自定节奏；任务定义与阈值见 `graphics/plan.md` 阶段 1）

- [ ] CMake 升 C++20 ｜ 对应 `plan.md` V1.x 前置
- [ ] matrix 实数/复数 scalar trait 改造（用户自行）
- [ ] linalg：LU + Solve
- [ ] 积分器：Euler / RK4 / Verlet

---

## 5. 里程碑

- [ ] Gate A 达成 → 启动第一个小项目（动力系统可视化）
- [x] ~~第一个小项目完成 → 拉 `viz` / `web` 多仓基建~~
      → **已提前至 2026-09-22**：不等 Gate A，直接以「笔记 web 化」启动两仓（理由见第 6 节日志）
- [ ] 笔记 web 化 A1 垂直切片走完五步（`learning-viz/python/notes_pipeline/DESIGN.md` §7）

---

## 6. 进度日志（append-only：只追加，不修改历史条目）

### 2026-09-22

- 文档治理：统一四层结构（入口 `AGENTS.md` / 真源 `README.md`+本文 / 展开 `plan.md`+`learning.md`+`matrix.md` / 产物 `learning-todo/`）。
- 本次改动不动学习进度，仅治理文档；`plan.md` 的验收项编号由 `A1`–`F4` 改为 `V{阶段}.{序号}`，
  以消除与笔记代号（`B1`、`E1` 等）的撞名。
- 删除仓库根目录的 `requirements`（与 uv 工作流重复，且内容是 pip 源说明，会误导 agent 走 pip 路径）。
- 清理工作区：删除多余的 worktree `graphics-0-2fe0063d` 及其分支 `workbuddy/graphics-0-2fe0063d`
  （确认无独有提交，分支 `12dfde4` 已合并进 `graphics_0`）。现在只剩主仓与 `graphics-0-0c774982`。
  被删 worktree 的 `.workbuddy/memory/` 已备份到
  `graphics-0-0c774982/.workbuddy/memory/archive/2fe0063d-2026-09-22/`（含 md5 校验）。
- `notebook/` 笔记体检与修补（commit `b899d33`，18 文件）：
  - 修 bug：`diffrential/de4.tex` 重复的空标题「知识网络」（原 135/136 行），后续小节重编号为 25.2–25.5。
  - 内容缺口：`A_analysis.tex` 补 σ-代数「可数并而非任意并」的理由；「可数可加性」由三行定义
    扩为完整论证（有限可加太弱的三反例表 + 失效定理清单 + 不可数可加太强 + 集合列连续性 + Vitali 三难）。
  - 消冗余：新增 11 个 `\label` 与 6 处指路句，打通 F→K、G→de4、H→I、B→lie1/lie2、A→de3 的重复链路。
  - 新增 `partN.tex` §15.5「章节与教材对照索引」（教材定位表 + 精读/略读分级）。
  - `main.tex` 给两个专题加 `\part` 分界（专题一 p.41 / 专题二 p.61），与 A.–O. 字母编号拉开层次。
  - 顺手修两处既有缺陷：全篇 38 处 ASCII 直引号规范为中文引号；9 处被静默丢字的字符
    （7×`↔`、2×`✓` → `\leftrightarrow` / `\checkmark`，Latin Modern 无此字形）。
  - 校验：`latexmk -xelatex` 通过，76 页，0 缺字 / 0 Overfull / 0 未定义引用 / 0 重复标签；
    另用 pdfpages+sips 导出关键页目视复核（p.11–12 新论证、p.39 索引表、p.41 `\part`、p.48/53 符号修复）。
  - `.gitignore` 补 `*.xdv`（XeLaTeX 中间产物）。
- `notebook/` 内容补缺（第二次会话）：用户要求把 `learning-todo/chats/` 里的「计算必修清单」落成正式章节，
  并盘查其余缺口。新增/扩写五块，`main.pdf` 由 76 页增至 **112 页**（A–R 字母连续，两个专题顺延到 §19 起）：
  - **P 计算必修**（`P_computation.tex`，P1–P10）：由黑板 `ses_f5e6…ipynb` md-24～md-34 转写为 LaTeX。
    组号 **G1–G9 → P1–P9**（`G1`–`G5` 已被 `G_numeric.tex` 占用，撞名已消）；自检题按约定留在聊天/黑板，不入笔记。
  - **C6–C10 概率深化**（扩写 `C_partial.tex`，不重编号 C1–C5，避免 A_analysis 里「→ C4」「C2 整节失效」的引用失效）：
    收敛模式与层级 + 反例博物馆 + 均匀可积性 + Egorov/Scheffé/Vitali 收敛/单调类定理 + Borel–Cantelli + 两种大数定律、
    特征函数与 CLT 证明、马尔可夫链与 MCMC/MLT、蒙特卡洛与方差缩减、Itô 微积分与 Fokker–Planck。
  - **A4 Sobolev 空间与弱解**（扩写 `A_analysis.tex`）：弱导数、Sobolev 嵌入与紧嵌入、Poincaré、Lax–Milgram、正则性提升。
  - **Q 拓扑学**（新建 `Q_topology.tex`）：点集拓扑、紧性（含与弱紧性/Arzelà–Ascoli/Rellich–Kondrachov 的对接）、
    商空间、同伦与 $\pi_1$、单纯同调与 Euler 特征、de Rham 上同调；与 `lie2` 的「拓扑结构」明确分工。
  - **R 复分析**（新建 `R_complex.tex`）：Cauchy 理论、留数（补上 `de3` 缺的 Laplace 反演一环）、共形映射与共形参数化、
    解析延拓（接 K4 的 Wick 转动）、复势与二维流体。
  - 顺带：表格定宽列统一为 `L{}`（`\raggedright`，定义在两个 `main.tex`），消除 CJK 窄列拉伸；
    `partN.tex` 教材对照索引补 P/Q/A4/C6–C10 五行；`AGENTS.md` 代号速查补 P/Q/R。
  - 校验：0 缺字 / 0 Overfull / 0 Underfull / 0 未定义引用 / 0 重复标签；子文档 `diffrential/main.tex` 也能独立编译（0 报错）。
- 仍缺、**尚未动手**的内容（下次可从这些里挑）：抽象代数（群环域、模、Galois）、图论与组合、
  凸优化与变分法系统化、数值 PDE 理论（Lax 等价定理、von Neumann 稳定性分析）、有限元理论、连续介质力学/弹性、
  复几何与辛几何、测度论进阶（Radon–Nikodym 的完整证明、乘积测度与 Fubini 的构造）、偏微分方程弱解的存在唯一（Evans Ch.6 级）。
- **多仓基建提前启动**（本日第三次会话）：`learning-viz` 与 `learning-web` 各落一个 init 提交，
  不再等「第一个小项目完成」。起因是拿到一份「笔记 web 化」设计（另一 session 产出，已存为
  `learning-viz/python/notes_pipeline/DESIGN.md` 并作为该模块的设计真源）：
  - `learning-viz`（commit `f782f4c`）：`README.md`（定位/红线/跨仓路径坑）、`.gitignore`（`dist/` 不入库）、
    `python/` 骨架（`notes_pipeline/vendored`、`app`、`viz_modules`、`contracts`、`tests`）、
    `python/notes_pipeline/DESIGN.md`。
  - `learning-web`（commit `3e44024`）：`README.md`、`.gitignore`、`src/` 骨架（`shell`、`renderers`、`viz`、`modules`）+ `public/`。
  - **不做** CDN、不做前端计算；代号（`A1`、`G1.2`）作路由 ID，不用 `§NN.N`。
  - 设计经评审后**就地校正 10 处**（文中标 `[校正 Cn]`），要点：`pandoc --katex` 实为**客户端**渲染须改机制
    （改走构建期 Node 预渲染 KaTeX）；filter 正则 `[A-O]` → `[A-R]`；插槽代号 `A1.1.3`/`A1.1.4` →
    `A1.3`/`A1.4`（不存在三级代号）；`\providecommand{\vizslot}` 需加父/子**两个** `main.tex`；
    R1「Python 3.14 缺 wheel」经 PyPI 核实**已解除**（`pypandoc-binary 1.17` 有 `py3-none-macosx_11_0_arm64.whl`、
    `flask 3.1.3` 为 `py3-none-any`）；笔记规模由「21 文件 / 72 页 / 33 节」更正为
    **24 内容文件 / 112 页 / 202 个 `\subsection`**。
  - 本次只做骨架与文档治理，未装任何依赖、未改 `notebook/`；`\vizslot` 兜底与步骤 1 尚未开始。
  - 父仓文档同步：`README.md` 两仓状态「空仓」→「开发中」；`AGENTS.md` 文档地图新增 `DESIGN.md` 一行、
    第八节目录用途改写；本节即第 5 节里程碑提前的说明。
