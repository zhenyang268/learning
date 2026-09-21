# learning — 项目总纲

> 自学数学系（→ PDE / 图形学工程数学）的多仓学习工程：**后端算、前端看。**

---

## 仓库职责

| 仓 | 地址 | 职责 |
|---|---|---|
| `learning` | https://github.com/zhenyang268/learning.git （镜像 `git@gitee.com:we_we_we/learning.git`） | 主仓/后端。C++ 单一实现，出 `.so` + Python 绑定；单元测试保证接口与数据正确。含数学、图形学、笔记。 |
| `todo` | https://github.com/zhenyang268/learning-todo/.git | 子仓。session 存档备份（进度管理 + 图片输出），整体不拆。 |
| `viz` | https://github.com/zhenyang268/learning-viz.git （待创建） | 接口/生成层。调 `learning` 接口，产完整结果（轨道、DAG）并本地渲染（ImGui / matplotlib）。 |
| `web` | https://github.com/zhenyang268/learning-web.git （待创建） | 展示平台。统一外壳 + 模块注册，渲染知识图谱与时间轴演示；可被多个后端仓复用。 |
| *(future) 408* | 待创建 | 另一后端内容仓，端点自备，`web` 一并托管。 |

---

## 项目设计

- **数据流**：`learning`(纯函数) ← `viz`(生成+渲染) ← `web`(外壳+模块)。
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
learning/                    viz/                        web/
  graphics/    C++ 引擎        cpp/    ImGui/OpenGL         shell/     外壳(导航/路由/主题)
  notebook/    LaTeX 笔记      python/ Flask / matplotlib    renderers/ graph / timeline / 专用
  learning-todo/  (子仓) session 存档   generate/ 产出完整结果          modules/   功能注册(端点+契约+渲染器)
  learning.md  数学手推清单      (契约 schema)                 assets/    贴图等(LFS/CDN)
```

---

## 文档导航

- `AGENTS.md`（协作规约：每 session 自动加载的分层与硬规则）
- `learning.md`（数学手推清单：12 阶段经典问题 + 6 阶段几何物理读书路线）
- `graphics/plan.md`（C++ 工程实现里程碑：阶段 1–6 与验收节点）
- `learning-todo/ns-progress.md`（NS 路线进度：L0–L5 与当前位置）
- `learning-todo/matrix.md`（Matrix 库设计备忘：容器/算法/算子分层）

---

## 工程示例（各部件承担什么）

| 示例 | `learning`（后端：算） | `viz`（生成 + 本地渲染） | `web`（渲染 + 交互） |
|---|---|---|---|
| **太阳系** | 轨道要素 + Kepler 方程求解 → 行星位置序列 | 生成整条轨道与太阳方向，出服务端点 | 3D 渲染器：画星球/轨道、相机、时间轴播放 |
| **自动微分** | `grad_graph`/`dual` → 计算 DAG + 每节点前值/微分 | 生成 DAG 结构（含数值） | 图渲染器：节点-边展示，可点节点看值/导数 |
| **TCP 传输** | （由未来 `408` 仓承担） | （由未来 `408` 仓承担） | 时间轴渲染器：重放报文/状态动画 |

> 注：TCP 的后端是 `408` 仓（承担 `learning`/`viz` 的角色），`web` 复用同一时间轴渲染器——这正是 web "多后端复用"的验证。

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