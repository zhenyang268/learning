# learning — 项目总纲

> **读者**：人（项目全景）｜ **最后更新**：2026-10-26
> 给 agent 的协作规约、目录分层与硬规则全部在 **`AGENTS.md`**（每 session 自动加载），本文不重复。
>
> **真源约定**：**现在做什么 → `NOW.md`**（唯一滚动文件）｜ 进度与历史 → `dev-progress.md` ｜
> 验收判据 → `graphics/acceptance.md` ｜ 数学地图 → `docs/route-math.md` ｜
> 红线 → `graphics/acceptance.md` §0，设计理由 → 各头文件顶部注释 ｜ 踩坑与口径 → `docs/pitfalls.md`。
>
> **不再维护大阶段计划**（阶段 1–6 / 八刀 / 扇面切片 / 六阶段读书表已于 2026-10-26 全部退役）。
> 要做的事只进 `NOW.md`，1–3 周为限，过期整篇重写。理由见 `AGENTS.md` §一「文档分三类」。

> 自学数学系（→ PDE / 图形学工程数学）的多仓学习工程：**后端算、前端看。**
> 目标口径：NS / PDE 是**打开路线的钥匙**（用于选题与排序），基础本身也是目标 —— 不是可跳过的脚手架。

---

## 快速开始

```bash
# 1. Python 环境（uv 管理，Python 3.14）
uv sync

# 2. 构建 C++ 引擎与 demo，跑测试
cmake -S graphics -B graphics/build -DCMAKE_BUILD_TYPE=Debug
cmake --build graphics/build -j
ctest --test-dir graphics/build --output-on-failure
# 或一键：./graphics/run_tests.sh    （-d 明细 / -c 从零 / -w 开 -Wall -Wextra）

# 3. 编译数学笔记（产出 math-notebook/main.pdf）
cd math-notebook && latexmk -pdfxe -interaction=nonstopmode main.tex && cd ..
```

> ⚠️ **笔记必须用 `-pdfxe`（xelatex）。** 用 `-pdf`（pdfLaTeX）时 ctex 会报
> `fontset 'mac' unavailable` 并**卡在交互输入等键盘，命令永远不返回**（2026-09-24 实测）。
> 目录是 `math-notebook/`，**不是** `notebook/`。

> 注意：**四个子仓（`learning-todo/`、`learning-viz/`、`learning-web/`、`learning-408/`）
> 不在本仓跟踪范围内**（无 `.gitmodules`、无 gitlink），新 clone 或 `git worktree` 里都不会出现。
> 绝对路径通常是 `~/develop/learning/<仓名>/`，**必须按绝对路径访问，改动在各自仓内提交**。

---

## 仓库职责

> `viz` 与 `web` 已建仓并各有提交（骨架 + 文档，2026-09-22）；
> 首个真实内容是「笔记 web 化」，**已完成**（`learning-web` 独立仓 318/318 单元入库，2026-10-04 收工）。
> 规约与设计真源在 `learning-web` 子仓（`AGENTS.md` + `DESIGN.md`，只能按绝对路径访问）。
>
> **职责切分（2026-09-22 起）**：笔记的**转换链路整体归 `learning-web/build/`**；
> `learning-viz` 收窄为**计算模块的后端**（`/api/viz/*` + 本地渲染），不再有 `/api/notes/*`。

| 仓 | 地址 | 职责 | 状态 |
|---|---|---|---|
| `learning` | https://github.com/zhenyang268/learning.git （镜像 `git@gitee.com:we_we_we/learning.git`） | 主仓/后端。C++ 单一实现，单元测试保接口与数据正确。含数学、图形学、笔记。 | 开发中 |
| `learning-todo` | https://github.com/zhenyang268/learning-todo.git | **产物存档仓**：只装不适合进父仓的二进制/大文件 —— `imgs/`（59 MB / 47 张图）+ `chats/`（会话 ipynb）。**文档不放这里**（2026-10-26 已全部搬回父仓 `docs/`），详见该仓 `README.md`。 | 使用中 |
| `learning-viz` | https://github.com/zhenyang268/learning-viz.git | 计算/生成层。调 `learning` 接口，产完整结果（轨道、DAG）并本地渲染（ImGui / matplotlib）。 | 开发中 |
| `learning-web` | https://github.com/zhenyang268/learning-web.git | 展示平台 + 笔记转换。统一外壳 + 模块注册；`build/` 把笔记转成知识图谱与阅读面。 | **已完成**（318/318） |
| `learning-408` | 独立子仓（被 `.gitignore` 忽略，各自管理） | 408 实践代码汇总。本体只保留一份**主题 → 代码导航**（`README.md`），下辖 `jos/`（操作系统）、`cs144/`（网络）、`csapp/`（组成）三个实践子仓。 | 使用中 |

仓库内目录（**非子仓**）：`11408-notebook/` —— 11408 考研学习体系（五件套：归档规约
`AGENTS.md` / 学习引擎 `LEARNING.md` / 图技能 `TOOLKIT.md` / 模块路线 `modules/`、进度仪表盘 `PROGRESS.md`），
学习对话按规约归档为 `notes/` 下的 ipynb。**它只服务考研，与 408 实践代码仓是两回事**；
agent 默认不碰（用户 2026-10-17 明确「不用管考研线」）。

---

## 项目设计

- **数据流**：`learning`(纯函数) ← `viz`(生成+渲染) ← `web`(外壳+模块)。
- **例外**：**笔记 web 化**不走这条链 —— 笔记是内容不是计算，其转换链路
  （`.tex` → JSON / standalone HTML）整体在 `learning-web/build/`，只有**计算模块的演示**
  才经 `viz` 出契约。见子仓 `learning-web/AGENTS.md` 决策 D1 / D6。
- **原则**：
  - 单一实现 + 绑定（Python 复用 C++，仅测试层允许朴素对拍）；
  - 不重复实现；
  - `scene` 归 `viz`、`learning` 无状态；
  - **时间轴 = 渲染层、物理参数 = 计算层**；
  - 每功能独立数据契约 + version。

---

## 目录关系

```
learning/  ──构建──▶  .so + python绑定  ──被调用──▶  viz/  ──Flask+契约──▶  web/
```

```
learning/                        viz/                        web/
  graphics/    C++ 引擎            cpp/    ImGui/OpenGL         build/      笔记转换(Python)
  math-notebook/ LaTeX 笔记        python/ app(Flask) / 模块      graph/      语义边真源(yaml)
  learning-todo/  (子仓) 存档       viz_modules/ 计算模块          src/        shell / renderers / viz / modules
  learning-408/   (子仓) 408 代码   contracts/ 契约 schema        dist/       产物(gitignore)
  NOW.md       近期要做什么
  learning.md  数学手推题库
```

---

## 工程示例（各部件承担什么）

| 示例 | `learning`（后端：算） | `viz`（生成 + 本地渲染） | `web`（渲染 + 交互） |
|---|---|---|---|
| **太阳系** | 轨道要素 + Kepler 方程求解 → 行星位置序列 | 生成整条轨道与太阳方向，出服务端点 | 3D 渲染器：画星球/轨道、相机、时间轴播放 |
| **自动微分** | `grad_graph`/`dual` → 计算 DAG + 每节点前值/微分 | 生成 DAG 结构（含数值） | 图渲染器：节点-边展示，可点节点看值/导数 |
| **TCP 传输** | （由 `learning-408/cs144` 承担） | （同左） | 时间轴渲染器：重放报文/状态动画 |

> 注：TCP 的后端是 `learning-408` 子仓（承担 `learning`/`viz` 的角色），`web` 复用同一时间轴渲染器
> ——这正是 web「多后端复用」的验证。

### 地球（唯一展开的样板）

**前后端责任划分**

| 目标 | 后端（`learning` + `viz`） | 前端（`web`） |
|---|---|---|
| 数据准确度 | 轨道要素、Kepler 求解、地轴倾角、自转角、太阳方向、时间系统 | 不参与，只把后端给的方向/角度喂给光照 |
| 视觉逼真度 | 不参与 | 贴图、PBR 材质、昼夜混合、大气、云层、后处理、性能优化 |

**接口**：后端随时间/参数吐出

```
{ time, earth: {position, rotation_angle, axis_tilt}, sun: {direction} }
```

**准确度分级**：T1 理想开普勒（低） / **T2 JPL 近似轨道要素表（推荐起步）** / T3 DE440 星历 + 岁差章动（高）。

**资源清单**

- 贴图（进 `web`）：Blue Marble 地表、地形法线、夜灯、云层、星空 cubemap —— 来源 NASA Visible Earth / Earth at Night / solarsystemscope。
- 天文数据（进 `learning`/`viz`）：JPL 近似轨道要素表（T2）、DE440 星历（T3）、IAU 常数。

**仓库约定**：大贴图 → `web` 仓走 Git LFS/CDN；星历大文件 → 不提交、脚本下载 + 缓存；每个资源记来源 + 授权。

---

## 文档导航

| 文件 | 一句话 |
|---|---|
| `AGENTS.md` | **协作规约与入口**：文档地图、文档三分类（防腐化条款）、命令、代号速查、硬规则、环境坑 |
| `NOW.md` | **近期要做什么**（唯一滚动文件，1–3 周一版，过期整篇重写） |
| `dev-progress.md` | **进度总览**：§0 现状表 + §2 Gate + §5 里程碑（**不含日志正文**） |
| `docs/logs/<YYYY-MM>.md` | **进度日志**（按月一个文件，append-only，最新在上） |
| `docs/route-math.md` | 数学**地图**：十段、双线依赖、每段"看懂线"、不做的边界、Lee 的正确读法 |
| `learning.md` | 数学**手推题库**（12 阶段经典问题；只是题，不含路线） |
| `graphics/acceptance.md` | **验收判据规范**：判据五条、闸门与注入、历史 `V*` 编号对照 |
| `docs/pitfalls.md` | **口径与踩坑**：判据病态、共轴掩盖 bug、holonomy 口径、泛型域、实例化才叫落地 |
| `docs/eng-competency.md` | **工程能力 60/80 基线**：四赛道检查点、赛道重叠矩阵、本仓现状对照（**不含顺序与期限**） |
| `docs/baseline-60.md` | **四线 × 两把 60 分尺子**（学习 vs 工程）：数学 / 数学物理 / CG / CS；补工程侧未覆盖的两条线 |
| `docs/math-log.md` | 数学线**会话日志**（append-only） |
| `docs/interview-depth.md` | 面试深度清单 + 本仓实况校正 + 明确不做的边界 |
| `docs/env-git.md` / `docs/doc-checkup.md` | git 沙箱坑 / 文档体检清单（低频） |
| `graphics/demo/pyapps/tensor-convention.md` | pyapps 层张量记号约定 |

> **已删除（2026-10-26 / 10-28）**：`graphics/plan.md`、`learning-todo/{math-plan,fan-plan,lie-plan,stage1-wrapup}.md`、
> **`docs/matrix.md`（设计文档，1197 行）**。
> 内容已迁走（红线 → `graphics/acceptance.md` §0；张量契约 → `tensor.hpp` 头注释；
> 数学路线 → `docs/route-math.md`），**不要再按旧路径找内容**；
> `docs/logs/` 里的历史条目仍按旧路径引用，属历史记录。
