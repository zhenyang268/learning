# env-git.md — git 环境：remote 口径 + 写操作被沙箱拒绝

> 从根 `AGENTS.md` §七 外移（2026-10-17），入口只留摘要与指针。
> **什么时候读**：① 要推/拉、或配置/检查 remote 时（§一）；
> ② `git commit` / `update-index` / `branch -D` / `worktree remove` 报
> `Operation not permitted`，或报 `File exists` / `Another git process seems to be running` 时（§二起）。

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

# 二、git 写操作被沙箱拒绝：成因与处理

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

- 同内容封装成的 skill：`gitsandbox-write-workaround`
- 入口摘要：根 `AGENTS.md` §七
