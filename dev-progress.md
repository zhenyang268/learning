# dev-progress.md — 进度与排期

> **读者**：人 + agent ｜ **最后更新**：2026-09-22 ｜ **本文是「进度」的唯一真源**
>
> **真源声明**
> - 全项目「现在在哪、下一步做什么」**只以本文为准**。
> - `graphics/plan.md` 只负责任务定义与验收阈值（验收项记作 `V{阶段}.{序号}`），**不记录状态**。
> - `learning-todo/ns-progress.md` 负责 NS 路线的逐次会话日志与关卡状态。
> - 冲突时以本文为准；发现冲突请直接改本文，并同步删掉别处的重复表述。
>
> 分工：AI 排期 / 定门槛 / 验收，用户写代码。优先学习，工程增量推进。

---

## 0. 当前指针（每次会话先看这里）

> **数学线在 L1（测度与积分）中段；工程线停在 `graphics` 阶段 1 的起点，尚未开工。**

| 线 | 当前位置 | 下一个动作 |
|---|---|---|
| **数学（优先）** | L1 第 3 关 —— 外测度已学、勒贝格积分刚学 | Carathéodory 判据 → 积分四步构造 |
| **工程** | `graphics` 阶段 1 未动：matrix 未做实数/复数 scalar trait 改造，无 `linalg`，CMake 仍是 C++17 | CMake 升 C++20；`core/scalar.hpp` + `core/matrix.hpp` 重写 |
| **小项目** | 未解锁 | 等 Gate A（见第 2 节） |

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

- **数学线**：L1 第 3 关 —— 外测度已学、勒贝格积分刚学；差 Vitali → Carathéodory → 积分四步构造。
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
- [ ] 第一个小项目完成 → 拉 `viz` / `web` 多仓基建

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
