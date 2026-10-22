# 三维旋转 Python 数值验证报告

> 位置：`graphics/demo/pyapps/rotation/`（纯 Python 现象演示区，不进 ctest、不进 CMake）
> 运行方式：`uv run python graphics/demo/pyapps/rotation/<名字>.py`，每个文件独立可跑
> 图片输出：`learning-todo/imgs/rot/`
> 日期：2026-10-17 ｜ 依赖：numpy + matplotlib（未用 scipy）
> 对应切片：`fan-plan.md` **A5「四元数 + Slerp」** 的 pyapps 预演（原标 C++、缺 `quat/`）
> 编号：**A5.1–A5.6**（与球面目录 G1–G6 并列）

---

## 0. 主线：五种表示、一个对象

**轴角（罗德里格斯）→ 四元数 → 欧拉角 → 矩阵：同一个旋转的五张皮，各有各的坏点。**

| 步 | 文件 | 一句话 | 验收数字 |
|---|---|---|---|
| A5.1 | `rodrigues.py` | 无穷级数收成三行闭式，靠一条恒等式 | 闭式 vs 级数 6.7e-16 |
| A5.2 | `quaternion.py` | 半角来自"平方"，双覆盖来自"对径点" | M(q)=M(−q) 精确 0 |
| A5.3 | `euler_gimbal.py` | 万向锁是参数化缺陷，不是数学奇点 | 塌缩等式 ≤3.3e-16 |
| A5.4 | `slerp.py` | 在 S³ 上走大圆才匀速 | 中点精确平分 47.5940° |
| A5.5 | `representations.py` | 四种互转的坏点各不相同 | 往返全部 <1e-12 |
| A5.6 | `compose.py` | 左乘/右乘、一阶复合、顺序 bug | 变轴顺序差 O(1) |

**建议阅读顺序**：A5.1（从旋转移到矩阵）→ A5.2（四元数，最常用）→ A5.4（插值，立刻能用）
→ A5.3（欧拉角与万向锁，理解"为什么不用它"）→ A5.5（互转与坏点）→ A5.6（实战陷阱）。

---

## 1. A5.1 罗德里格斯公式（`rodrigues.py`）

叉积写成矩阵（hat 映射），**反对称**即李代数 so(3)：

$$[u]_\times = \begin{pmatrix}0&-u_z&u_y\\ u_z&0&-u_x\\ -u_y&u_x&0\end{pmatrix},\qquad a\times b = [a]_\times b$$

指数映射的闭式（罗德里格斯）：

$$\exp(\theta[u]_\times) = I + \sin\theta\,[u]_\times + (1-\cos\theta)\,[u]_\times^2$$

**成立的唯一关键是** $[u]_\times^3 = -[u]_\times$（实测 2.2e-16）——它把幂次降回 1/2 次，
于是级数只剩三类项，系数实测反解为 $A=1,\ B=\sin\theta,\ C=1-\cos\theta$。

**两个等价形式**（实测互差 0）：

$$I + \sin\theta\,[u]_\times + (1-\cos\theta)[u]_\times^2
\;\equiv\; \cos\theta\,I + \sin\theta\,[u]_\times + (1-\cos\theta)\,uu^T$$

互换靠恒等式 $[u]_\times^2 = uu^T - I$。**注意**：写成 $\cos\theta I + \cdots + (1-\cos\theta)[u]_\times^2$
是**错的**（把 $I$ 系数改掉时必须同时换 $[u]_\times^2 \to uu^T$）。

验收：闭式 vs Taylor 级数（36 组轴×角）最大偏差 **6.7e-16**；随机 2000 组全在 SO(3)
（正交 1.6e-15、|det−1| 1.8e-15）；小角极限 $(R-I)/\theta \to [u]_\times$ 误差 ∝ θ。

## 2. A5.2 四元数（`quaternion.py`）

$$q = \cos\tfrac\theta2 + \sin\tfrac\theta2\,u,\qquad v' = q\,v\,q^*$$

**半角的来历**：四元数乘法 = 旋转复合（实测 $q(\theta_1)\otimes q(\theta_2) = q(\theta_1+\theta_2)$，
$q_{90}\cdot q_{90} = q_{180}$ 差 0）。所以"转 θ"要**两次**半角。

**双覆盖**：$M(q) = M(-q)$（1000 组精确 0）。几何上 q 与 −q 是 S³ 上的**对径点**。
转 360° 时 $q\to -q$，转 720° 才回到 $q$ —— 旋量（spinor）的数学原型。

**非交换**：$(w_1+\mathbf v_1)(w_2+\mathbf v_2) = (w_1w_2-\mathbf v_1\cdot\mathbf v_2)
+ (w_1\mathbf v_2+w_2\mathbf v_1+\mathbf v_1\times\mathbf v_2)$，
最后那一项 $\mathbf v_1\times\mathbf v_2$ 是唯一来源。实测 $q_xq_z \ne q_zq_x$（差 1.0），
与矩阵 $R_xR_z\ne R_zR_x$ 完全一致 ⇒ **复数表示不了三维旋转**。

**漂移（实测修正了一个常见误解）**：

| 复合次数 | 四元数 \|q\|−1 | 矩阵 \|RᵀR−I\| |
|---|---|---|
| 1 000 | 1.8e-15 | 4.1e-15 |
| 10 000 | 6.9e-15 | 2.6e-15 |
| 1 000 000 | 4.5e-14 | 6.7e-14 |

两者都是 **√n·eps 级、量级相当** —— 不存在"矩阵漂得快、四元数不漂"。
真差别在**修复成本**：四元数 q←q/|q| 是 O(1) 且仍在单位球面上；矩阵必须重正交化（投影，改变原矩阵）。

## 3. A5.3 欧拉角与万向锁（`euler_gimbal.py`）

12 种约定（6 轴序列 × 内外旋）。三条实测恒等式（≤2.2e-16）：

① 内旋(seq,a) = 外旋(seq 倒序, a 倒序)
② 内旋(seq,a) = 外旋(seq, −a)ᵀ
③ 内旋(seq,a)ᵀ = 内旋(seq 倒序, −(a 倒序))

**反面**：内旋(seq,a) ≠ 外旋(seq,a)ᵀ（实测差 1.4）。这条最容易写错。

**万向锁（ZYX 内旋）**：$R = R_z(y)R_y(p)R_x(r)$，分解式 $p=\arcsin(-R_{20})$、
$r=\text{atan2}(R_{21},R_{22})$、$y=\text{atan2}(R_{10},R_{00})$。

|pitch| = 90° 时 $R_{21}=R_{22}=0$ ⇒ atan2(0,0) 无定义。实测塌缩等式：

$$R_z(y)R_y({+90°})R_x(r) = R_z(y-r)R_y({+90°})R_x(0)$$
$$R_z(y)R_y({-90°})R_x(r) = R_z(y+r)R_y({-90°})R_x(0)$$

⇒ **yaw +θ ≡ roll −θ**（实测 5.6e-17）。几何本体：pitch 转 90° 后，**yaw 轴（世界 z）
与 roll 轴（物体 x）共线（夹角 180°）**。

**准确表述**（纠正常见误解）：不是某个旋钮"失效"——实测 yaw 与 roll 各自转 30°
都确实造成 30° 姿态变化；**丢失的是可分辨性**（两轴共线 ⇒ 只剩 1 个独立自由度）。
同一姿态有**无穷多组**欧拉角（$y-r$ 不变即可，四组实测差 1.7e-16）。

## 4. A5.4 Slerp（`slerp.py`）

$$\text{slerp}(q_0,q_1,t) = \frac{\sin((1-t)\Omega)}{\sin\Omega}q_0
+ \frac{\sin(t\Omega)}{\sin\Omega}q_1,\qquad \cos\Omega = q_0\cdot q_1$$

权重是 sin 而非线性，因为 $|q(t)|$ 必须恒为 1（线性权重会掉进球内，实测 $|q|$ 最低 0.857）。

| 方式 | 端点 | 角速度 |
|---|---|---|
| slerp | 精确 | **匀速**（相邻增量 ptp 2.0e-15，良态口径）|
| nlerp | 精确 | 不匀（ptp 2.0e-2，两端慢中间快）|
| 欧拉角线性插值 | 精确 | 不匀且中点不平分（三轴联动时差 0.259）|

**短弧**：$q_0\cdot q_1<0$ 时必须翻 $q_1$ 符号，否则走 350° 长弧（实测：绕 z 转 350°
姿态上只离 10°，不翻号会走长弧）。这是**必需**步骤，不是优化。

**踩坑（良态判据）**：匀速性若用 $\arccos((\text{tr}-1)/2)$ 量角度，误差会假性涨到 3e-8；
换 $\text{atan2}(\|\text{vee}(R-R^T)\|/2,(\text{tr}-1)/2)$ 后回到 **2.0e-15**。
**是判据病态（θ→0 时 arccos 导数发散），不是 slerp 错。**

## 5. A5.5 四种表示互转（`representations.py`）

| 转换路径 | 往返最大误差 |
|---|---|
| 矩阵 ↔ 四元数（Shepperd） | 机器精度 |
| 矩阵 ↔ 欧拉(ZYX) | 机器精度 |
| 轴角 ↔ 矩阵（trace / vee） | 机器精度 |
| 四元数 ↔ 矩阵 | 机器精度 |

**坏点才是难点**：

- **四元数**：矩阵→四元数**非单值**（1000 组中 170 组取回 $-q$）。用法上必须比较 $|q_0\cdot q_1|$
  或相邻帧强制同号，否则插值会跳。
- **轴角**：θ=0 时轴无定义（任意轴都行）；θ=π 时轴定向二义。
  **踩坑**：`rotation_angle_of` 原用 `acos((tr−1)/2)`，在 θ=π（tr=−1）处导数发散，
  给出 179.999998°，**令"反解轴+角重建矩阵"差 2e-8**。换 atan2 良态口径后
  θ=π 重建回到 **6.7e-16**。与 A5.4 同族：**判据病态 ≠ 实现错**。
- **欧拉角**：pitch 89.999° → 90.001°（姿态差 0.002°）时，反解 yaw 从 0° **跳到 180°**。

选型表：矩阵 9 数冗余 6 约束 / 四元数 4 数有 ± / 轴角 4 数不能直接复合 /
欧拉角 3 数有 12 种约定 + 万向锁。

## 6. A5.6 复合旋转（`compose.py`）

**左乘/右乘**：增量在左 = 绕**世界**轴；在右 = 绕**物体**轴。两者不同（实测差 1.0，
因绕世界 x 与绕物体 x 不是一回事）。统一口诀：

$$R_0 R_\delta \;\;\equiv\;\; \tilde R_\delta R_0,\qquad \tilde R_\delta = \exp(\theta[R_0 u]_\times)$$

（把增量轴先经 $R_0$ 搬运，再做世界系旋转。）四元数乘法顺序与矩阵**一一对应**。

**关键洞察（本演示最值钱的一条）**：

> 如果增量矩阵**都相同**（常轴 ω），那么 $R\leftarrow R_pR$ 迭代 n 次与
> $R\leftarrow RR_p$ 迭代 n 次**都是 $R_p^n$** —— **顺序无关**！
> 且 $\exp(\theta_1K)\exp(\theta_2K)=\exp((\theta_1+\theta_2)K)$ 精确成立 ⇒ 连误差都没有。

| 情形 | 左乘 vs 右乘 |
|---|---|
| 常轴（增量相同） | 舍入级（4.8e-15）⇒ **测不出顺序 bug** |
| 变轴（增量随 t 变） | **O(1)**（姿态跑偏 25.53°），细步长救不回来 |

**"一阶复合"的确切位置**：不是"小步复合"本身，而是**近似增量** $I+\omega\,dt$ 代替 $\exp(\omega\,dt)$：

| 步数 n | exp 增量 | $I+\omega\,dt$ 增量 |
|---|---|---|
| 200 | 5.0e-15 | 5.4e-02 |
| 800 | 1.0e-14 | 1.3e-02 |

后者误差 ∝ dt（一阶），且矩阵偏离 SO(3)。**纪律**：验证复合顺序必须喂
**非共轴 + 时变**序列（与 graphics 3.3「判据要独立于实现路径」同族）。

**欧拉角速度陷阱**：$\omega = W(\alpha,\beta,\gamma)\cdot(\dot\alpha,\dot\beta,\dot\gamma)$，
$W$ 依赖姿态、在万向锁附近奇异。实测：物体系纯绕一根轴转 1 rad/s，
欧拉角**三个分量都在变**（如绕 y 时 $(-11.9, 56.2, -3.5)$ °/s）⇒ 姿态动力学不该用欧拉角做状态量。

---

## 7. 方法学备忘（本次踩过的真坑）

1. **判据病态 ≠ 实现错**（连续踩两次）：`rotation_angle_of` 与 slerp 的匀速性判据
   都因 $\arccos$ 在 θ→0/π 处导数发散而给出假误差（2e-8 / 3e-8）。
   换 $\text{atan2}(\|\text{vee}\|/2,(\text{tr}-1)/2)$ 后双双回到机器精度。
   → 与 graphics 3.3「浮点残渣必然非零 ⇒ 判量级」同族，但**更细一层**：
   **不仅要判量级，还要检查判据函数本身的条件数**。
2. **"顺序很重要"必须先验证前提**：常轴时左乘≡右乘（都是 $R_p^n$），
   我第一版用常轴实验"证明"顺序重要，实测差 0 —— 造实验前先想清**什么条件下才该有差异**。
3. **Python 字面量会静默吃掉 LaTeX 控制符**：普通三引号里 `\t \v \b \f`
   被转义（`\theta`→TAB+`heta`）。**但 `.py` 里中文双引号写成 ASCII `"` 会直接语法错误**，
   本次 6 个脚本共修 30+ 处（统一改「」）。写 pyapps 时**中文引号一律用「」**。
4. **"内旋 = 外旋转置"是错的**（差 1.4）。正确：内旋(seq,a) = 外旋(倒序, 倒序 a)。
   我连续三次推断错，最后靠**枚举实测**定案 —— 猜公式不如跑一遍。
5. `rotation_axis_of` 在 θ=π 用 $R+I$ 的列取方向是**确定的一支**（不随机），
   所以重建可靠；但 θ=0 时确实无信息（返回兜底值），调用方要先判 θ。

## 8. 待做（下一轮）

- A5.7 刚体动力学：$\dot R = R[\omega_{body}]_\times$ 与欧拉方程（惯性张量）的数值积分
- A5.8 SO(3) 上的优化：姿态残差最小化（接 A3 点云配准，SVD 与 Gauss-Newton 对拍）
- 与 C++ 侧挂接时机：一旦做 `quat/`（A5 的 C++ 正课），本目录可作为**判据来源**

## 9. 运行一览

```bash
cd /Users/apple/develop/learning
uv run python graphics/demo/pyapps/rotation/rodrigues.py        # A5.1
uv run python graphics/demo/pyapps/rotation/quaternion.py       # A5.2
uv run python graphics/demo/pyapps/rotation/euler_gimbal.py     # A5.3
uv run python graphics/demo/pyapps/rotation/slerp.py            # A5.4
uv run python graphics/demo/pyapps/rotation/representations.py  # A5.5
uv run python graphics/demo/pyapps/rotation/compose.py          # A5.6
```

公共数学层 `_rot_core.py`：hat/vee、罗德里格斯（闭式 + 级数双实现）、SO(3) 判定、
四元数全套（乘法/共轭/矩阵互转/旋转向量/slerp/nlerp）、欧拉角（12 约定 + ZYX 分解含锁死分支）、
旋转的轴角反解（**atan2 良态口径**）、3D 画图辅佐。

**自检**：`_rot_core` 有 15 组断言（见 `.workbuddy/memory/2026-10-17.md`），
覆盖罗德里格斯双实现、$[u]_\times^3=-[u]_\times$、四元数同态、双覆盖、内外旋三恒等式、
12 序列链式一致性、ZYX 全域往返（含锁死点）、短路 slerp、轴角反解。