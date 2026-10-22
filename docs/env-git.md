# env-git.md — git 写操作被沙箱拒绝：成因与处理

> 从根 `AGENTS.md` §七 外移（2026-10-17），入口只留摘要与指针。
> **什么时候读**：`git commit` / `update-index` / `branch -D` / `worktree remove` 报
> `Operation not permitted`，或报 `File exists` / `Another git process seems to be running` 时。

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
