# dev-progress.md — 开发进度与排期

> 分工：AI 排期 / 定门槛 / 验收，用户写代码。优先学习，工程增量推进。
> 关联：README.md（架构总纲）· learning.md（数学手推）· graphics/plan.md（工程里程碑）· learning-todo/ns-progress.md（NS 进度）

---

## 一、协作模式

- 数学线（慢，优先）与工程线（快，够用即做）并行。
- 每次会话：用户报"学到哪 / 写到哪" → AI 更新本文件并给下一步。
- **matrix 改造**与**自动微分**由用户自行抽时间完成，AI 只给验收点。

---

## 二、当前现状（2026-09-21）

- **数学线**：L1 第 3 关 —— 外测度已学、勒贝格积分刚学；差 Vitali → Carathéodory → 积分四步构造。
- **工程线**：graphics 阶段 1 未动 —— matrix 未做实数/复数 scalar trait 改造；无 linalg；CMake 仍 C++17。

---

## 三、小项目启动门槛（Gate）

| 门槛 | 数学 | 工程 | 解锁的小项目 |
|---|---|---|---|
| **Gate A** | L0 + ODE 数值解直觉（Euler/RK4/Verlet） | matrix 改造 ✔ + 最小 linalg（LU/solve）+ 一个积分器 | 动力系统可视化：单摆 / 行星轨道 / 双摆 / N 体 |
| **Gate B** | **L1 完成**（测度 + 积分 + 收敛定理 + 四步构造） | Gate A + FDM 基础 | PDE 可视化：热方程（图像模糊）、波动（鼓膜） |
| **Gate C** | L2 泛函 + Sobolev | 阶段 4 完整 | 弱解类（通往 NS 真目标） |

> **第一个小项目定在 Gate A**；L1 继续按 NS 路线推进，不互相阻塞。

---

## 四、下一步

**数学线（优先）**
- [ ] 外测度（方盒覆盖、次可加）
- [ ] Vitali 不可测集
- [ ] Carathéodory 判据 → 定理
- [ ] Lebesgue 积分四步构造
- [ ]（顺带）Vitali 收敛定理

**工程线**（用户自定节奏；任务细节见 `graphics/plan.md`）
- [ ] CMake 升 C++20
- [ ] matrix 实数/复数 scalar trait 改造（用户自行）
- [ ] linalg：LU + Solve
- [ ] 积分器：Euler / RK4 / Verlet

---

## 五、里程碑

- [ ] Gate A 达成 → 启动第一个小项目（动力系统可视化）
- [ ] 第一个小项目完成 → 拉 viz / web 多仓基建