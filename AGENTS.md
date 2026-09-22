# 项目约定（每个 session 自动加载）

## 学习
- 我是非数学系有高数基础, 在自学数学系, 目的是能够有pde以及计算机图形学等相关的工程学的数学基础, 所以帮我补完基础, 可以特别严谨或者特别复杂的习题的优先级可以降低
- 一定一定要保证正确性, 一点符号书写函数错误都会浪费我很多时间
- 要保证一定美观和可读性, 如果碰到复杂的表达, 要进行解释和注解
- 介绍数学知识时要先符合直觉, 再去解释相关的定理, 英文定理可以都加上括号写上中文, 希腊/拉丁符号可以附上读法, 复杂的数学符号比如集合类的符号也可以加上
- 可以展示哪部分是哪部分的基础, 谁是用来解决什么问题的, 比如历史上碰到什么问题这些
- 介绍编程时先介绍这种语言特性或者设计写法背后解决的问题
- **NS 学习路线**：进度统一记录在 `learning-todo/ns-progress.md`，每完成一个话题追加会话日志；每天结束说"收尾"时汇总更新路线图状态；用户跳转时先检查停车场和路线偏差

## 图像生成（阿里云百炼 CLI）
- 图片优先保证正确性, 知识密度要高, 是自学导读, 字体不要太大, 不需要醒目标题之类的东西, 不是海报, 是知识图谱
- 一律使用 `bl image generate --model qwen-image-3.0-pro`（走 Token Plan 订阅；qwen-image-3.0 非 pro 只扣免费额度，剩余额度少，不要默认用它）
- 生成到imgs目下的时候可以考虑分类, 最多允许创建两级子目录
- 所有生成产物保存到 `learning-todo/` 及其子目录（图片默认 `learning-todo/imgs/`），`--watermark false`
- 尺寸参数含 `*` 时必须加引号（zsh 通配符问题），如 `--size "2048*2048"`
- 生成后必须读图校验数学内容正确性（公式、坐标、箭头指向），发现错误重出或明确告知用户瑕疵
- AI 生成图的坐标轴刻度数字经常乱码：讲义级精确配图改用 matplotlib 代码绘制

## 数学表达
- terminal输出指导当老师用, ipynb当黑板和笔记用
- 聊天回复中公式用 Unicode（∑ ∫ √ λ ᵀ 等），不输出裸 LaTeX
- 将正式输出指向到learning-todo/chats目录下的ipynb文件中, 可以按session id确定不同文件的归属session, 每次回复都追加新的markdown格式框, 不要大规模修改已经写好的文本框
- 需要正式笔记/渲染公式时写 `.ipynb`（Jupyter 用 KaTeX 渲染 `$$...$$`）
- 不需要额外输出验证的代码, 只有明确输出代码时才需要

## Python 环境
- 项目用 uv 管理：`.venv` + `pyproject.toml` + `uv.lock`，Python 3.14（`.python-version`）
- 跑命令用 `uv run xxx`；装包用 `uv add xxx`；禁止用系统 python3（/usr/bin 的 3.9.6）或 pip 直接装
- 全局 `python3` 是 uv 托管的 3.14（~/.local/bin），Homebrew python 已卸载，不要再装回

## 目录用途
- `notebook/`：LaTeX 笔记（.tex，已编译 main.pdf）
- `learning-todo/`：生成产物与待办输出
- `linear_algebra/`、`graphics/`：学习笔记目录
- 注：`learning-todo/`、`learning-viz/`、`learning-web/` 是同级独立 git 仓（各有自己的 remote），父仓不跟踪、只是物理放在一起，新 clone/worktree 中不出现，改动需在各自仓内提交

## 安全
- 不在聊天/文件中回显完整 API Key；bl 已配置鉴权，无需在命令里传 --api-key
