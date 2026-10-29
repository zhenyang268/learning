# env-git.md — git 环境：remote 口径 + 危险操作自我审查 + 沙箱拒绝的备用处理

> 从根 `AGENTS.md` §七 外移（2026-10-17），入口只留摘要与指针。
> **什么时候读**：① 要推/拉、或配置/检查 remote 时（§一）；
> ② **要执行任何破坏性 git 命令前**（`reset --hard` / `clean` / `checkout .` / `push --force` /
> `branch -D` / `stash clear`）—— **必读 §三**；
> ③ `git commit` / `update-index` / `branch -D` / `worktree remove` 报
> `Operation not permitted`，或报 `File exists` / `Another git process seems to be running` 时（§二，已降级为备用）。

---

# 一、remote 统一口径（2026-10-28 定）

**全部 8 个仓一律两端**，命名固定，**不许再各仓各叫法**：

| remote | 指向 | 用途 |
|---|---|---|
| **`origin`** | **GitHub**（`git@github.com:zhenyang268/<repo>.git`） | 留档 / 对外 |
| **`gitee`** | **Gitee**（`git@gitee.com:we_we_we/<repo>.git`） | **日常推送** |

**已统一的仓**（2026-10-28 实测）：

```
learning                     origin=zhenyang268  gitee=we_we_we
learning-todo                origin=zhenyang268  gitee=we_we_we
learning-viz                 origin=zhenyang268  gitee=we_we_we
learning-web                 origin=zhenyang268  gitee=we_we_we
learning-408                 origin=zhenyang268  gitee=we_we_we
11408-notebook               origin=zhenyang268  gitee=we_we_we
learning-benchmark           origin=zhenyang268  gitee=we_we_we
raft（~/develop/raft）        origin=zhenyang268  gitee=we_we_we
```

**本次修的三处历史不一致**（别再改回去）：

1. 父仓 `learning` 曾同时有 `github` 和 `origin` **指向同一地址** ⇒ 已删 `github`；
   分支只跟踪 `origin/*`，删除不影响。
2. `raft` 的 `origin` 原本指向 **Gitee**（与其他仓相反）⇒
   `origin` 改为 GitHub，另加 `gitee` 别名。
3. `11408-notebook` 只有 `gitee` ⇒ 补 `origin`（GitHub）。

**推送惯例**：日常 `git push gitee <branch>`；GitHub 手动推，不设自动同步
（详见 `learning-benchmark/docs/bench-workflow.md` §3 的同口径表述）。

**批量检查命令**（在新会话里先跑这个，别凭记忆）：

```bash
cd ~/develop
for d in learning learning/learning-todo learning/learning-viz learning/learning-web \
         learning/11408-notebook learning/learning-408 learning/learning-benchmark raft; do
  printf "%-22s" "$(basename $d)"
  echo "origin=$(git -C $d remote get-url origin 2>/dev/null || echo MISSING)  \
gitee=$(git -C $d remote get-url gitee 2>/dev/null || echo MISSING)"
done
```

---

# 二、git 写操作被沙箱拒绝：成因与处理（**2026-10-30 起降级为备用**）

> 用户已于 2026-10-30 **无条件允许**危险 git 操作，沙箱不再拦截写操作。
> 本节保留是因为：① 万一环境回退可查；② 级联锁文件的成因仍成立。
> **但真正的风险位置换了** —— 从「命令被拒」换成「命令成功且不可逆」，见 §三。

## 典型症状

```
warning: unable to unlink '.../.git/index.lock': Operation not permitted
fatal: Unable to write new index file
```

## 不是仓库的问题

- `.git` 没有 `uchg` / `schg` 标志、没有 ACL、属主可写；
- 锁文件的 mode / owner 与 shell 新建的文件完全一致；
- 同一路径下 shell 的 `rm` / `mv`（含覆盖已存在文件）都成功 —— **只有 `git` 进程被拒**（可复现）。

## 也不是「免沙箱」就能解决

在标记为 `Sandbox bypassed` 的运行里同样出现该报错。环境里存在
`CODEBUDDY_SANDBOX_BROKER_TRACE_ID`，但**具体是哪一层拦的、规则是什么，从进程内看不到**
—— **不要臆断成因**，直接按下面的办法处理。

## 真正的坑是级联

一次失败会留下 `*.lock`，之后所有 git 命令报
`Unable to create ...: File exists` / `Another git process seems to be running`。
看起来像随机失败，**其实是前一次失败的回声**。

## 实测有效的处理

1. 用 shell `rm -f` 删掉**具体那几个**锁文件：
   - `.git/index.lock`
   - `.git/packed-refs.lock`
   - `.git/refs/heads/**/*.lock`
   - `.git/ORIG_HEAD.lock`
2. **立刻重跑**该 git 命令。循环 3–8 次基本必过（通常第 1 次就过）。

**不要用通配符删** —— zsh 在无匹配时会中断整条命令。

## `GIT_OPTIONAL_LOCKS` 的适用范围

`GIT_OPTIONAL_LOCKS=0` 只对**纯读命令**（`git status`、`git diff`）有意义；
对 `commit` / `update-index` / `branch -D` **无效**。

## 删 worktree 会撞同一个坑

`git worktree remove` 可能已经删掉了工作目录，却卡在 `.git/worktrees/<name>` 上。补法：

```bash
rm -rf .git/worktrees/<name>
git worktree prune
git branch -D <branch>
```

## 相关

- 同内容封装成的 skill：`gitsandbox-write-workaround`（现在只在环境回退时用得上）
- 入口摘要：根 `AGENTS.md` §七

---

# 三、危险操作的自我审查（2026-10-30 起**每条都必须走**）

## 3.0 前提变了，先说清

2026-10-30 用户**无条件允许**了危险 git 操作。这句话的实际含义不是「随便跑」，
而是**唯一的机械安全网被拆掉了**：

| | 之前 | 现在 |
|---|---|---|
| 跑 `git reset --hard` | 沙箱拒绝 → 命令失败 → 工作区完好 | **成功执行** → 未提交改动没了 |
| 代价 | 烦，要重试 | **不可逆** |
| 谁在拦 | 环境 | **只有这份文档** |

所以本节不是「建议」，是**替代沙箱的那道闸**。闸在自己脑子里，写下来才算存在。

> **一句话规则：机械上允许 ≠ 被授权。**
> 允许（approved）解决的是「能不能跑」，不是「该不该跑」。
> 判「该不该」的依据只有一条：**用户在这个会话里有没有要我做这件事。**
> 没要 = 不做；拿不准 = 先问一句，问的成本远低于恢复的成本。

## 3.1 绝对禁止（除非用户**点名**要求，并复述后果后确认）

| 命令 | 为什么 |
|---|---|
| `git push --force`（尤其 `main` / `master`） | 覆盖远端历史，别人的 clone 全部失配。**对主分支一律不做**，用户要求时也要先警告 |
| `git clean -fdx` | **本仓最危险的命令**，见 §3.2 |
| `git reset --hard`（工作区有未提交改动时） | 未提交的 WIP 直接蒸发，`reflog` 救不回来 |
| `git checkout .` / `git restore .` | 同上，且它**不在 reflog 里** |
| `git stash clear` / `git stash drop` | stash 里的东西同样不进 reflog |
| `git branch -D` | 未合并的提交会失去唯一引用（reflog 有时限） |
| `git rebase` / `git commit --amend`（已推送的提交上） | 改写公共历史 |
| `git filter-branch` / `git gc --prune=now` | 不可逆地重写/清理对象库 |

## 3.2 本仓特有的高危：`git clean` 会删掉 6 个独立子仓

父仓 `.gitignore` 里登记了 **6 个独立 git 仓 + 第三方源码 + 虚拟环境 + 产物仓**：

```
learning-todo/  learning-viz/  learning-web/  learning-benchmark/
11408-notebook/ learning-408/  graphics/third_party/  .venv/  build/
```

`git clean -fdx` 的 `-x` 表示**连 ignored 一起删** ⇒ 它会：

- 删掉 **6 个独立仓库的整个工作树**（每个都有自己独立的提交历史，父仓无法恢复）
- 删掉 `learning-todo/imgs/`（59 MB 产物图，独立仓）
- 删掉 `third_party/googletest` 与 `third_party/imgui` 源码（要重新下载，需网络）
- 删掉 `.venv/`（要重跑 `uv sync`）

**`-x` 这个字母，本仓永远不加。** 需要清理未跟踪文件就写具体的 `-e` 排除，或一个目录一个目录地删。

## 3.3 执行前的四步仪式（破坏性操作一律走）

```bash
# ① 记锚点：出事了至少知道回到哪
git rev-parse --short HEAD

# ② 看现场：有未提交改动就是高危，先判断这些改动值不值得丢
git status --short

# ③ 零成本快照（不切分支、不改变任何东西，可事后单文件取回）
git branch backup/pre-op-$(date +%H%M%S)

# ④ 把「我要执行 X，会丢失 Y」在会话里**写出来**，再执行
```

第 ④ 步不是形式。**默默跑掉一条 `reset --hard` 和明确说出「这会丢掉 solve.hpp 的改动」再跑，
是两种完全不同的行为** —— 后者几乎总会让人发现「其实不用 reset」。

第 ③ 步的价值：`backup/*` 分支不占工作区、不影响任何命令，
事后一个 `git checkout backup/pre-op-<时间> -- <文件>` 就能捞回单个文件。
这是替代沙箱的那张网，成本接近于零。

## 3.4 `reflog` 能救什么、不能救什么

| 丢的东西 | `reflog` 能救？ |
|---|---|
| 被 `reset --hard` / `rebase` 移走的**已提交**内容 | **能**（`git reflog` 找旧 HEAD） |
| **未提交**的工作区改动（`checkout .` / `reset --hard` 之后） | **不能** |
| stash 里被 `drop` / `clear` 掉的 | **不能**（stash 是独立引用，不在 reflog 里） |
| `git clean` 删掉的 untracked 文件 | **不能**（从未进过对象库） |

⇒ 「有 reflog 兜底」**只对已提交内容成立**。本仓两次真事故都是「未提交的 WIP 被覆盖」，
reflog 一次都没用上（`AGENTS.md` §二 记着这两起）。

## 3.5 多 agent 并行时的额外约束

现在有多个 agent 在同一个仓里干活（例如按 `docs/prompt-gtest-migration.md` 迁移测试用例）。
git 的索引是**共享状态**，两个 agent 同时 `git add` 会互相截胡：

- **同一个仓，同一时刻只允许一个写操作序列**（`add` → `commit` 之间不许别人插进来）。
- 写操作前先 `git status --short`：**看到不认识的改动就停下**，
  那可能是另一个 agent 正在改的文件，`git add -A` 会把它一起提交。
- **禁止 `git add -A` / `git add .`**：一律按路径显式 `git add <file>`。
  这既是防截胡，也是防把 `third_party/` 之类意外内容带进索引（`.gitignore` 改过之后尤其要复核）。
- 冲突处理：不认识的改动 → **报告，不动**。

## 3.6 事故档案（本仓真实发生过的，都是「命令成功了」）

1. **「读—改—写」脚本不幂等** ⇒ `graphics/demo/linalg/expm.cpp` 的 `[E]` 段被插两遍，
   构建 `redefinition` 失败。**脚本成功返回了**，错的是它执行了两次。
2. **整文件 Write 覆盖未提交 WIP** ⇒ `graphics/src/linalg/solve.hpp` 的工作区改动被冲掉，
   靠 VSCode 本地历史才恢复。`git` 层面的任何兜底在这里都**没生效**。

两条的共性：**不是命令被拒，是命令按预期执行了，而预期本身是错的。**
这正是 §3.3 第 ④ 步（先说出来）存在的理由。
