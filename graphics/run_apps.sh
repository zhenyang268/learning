#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# graphics/run_apps.sh -- 附加区 (demo/apps) 的构建 + 运行入口
#
#   ./run_apps.sh                列出所有 app
#   ./run_apps.sh <name>         构建并运行
#   ./run_apps.sh <name> -c      先清构建目录, 从零配置再构建
#   ./run_apps.sh <name> -- 3 7  把 -- 之后的参数原样透传给 app
#   ./run_apps.sh -w             打开 -Wall -Wextra
#   ./run_apps.sh -j 4           并行度 (默认: 本机核数)
#   ./run_apps.sh -b /tmp/gb     指定构建目录 (默认: graphics/build)
#   ./run_apps.sh -h             帮助
#
# 与 run_tests.sh 的分工: 那个跑 ctest 验收 (绿灯 = 实现对不对);
# 本脚本跑 apps 演示 (看现象, 不判对错, 退出码 = app 自己的退出码)。
#
# 新增 app: 在 demo/apps/ 放 <name>.cpp, 并在 demo/apps/CMakeLists.txt 里加
#   add_executable(<name> EXCLUDE_FROM_ALL <name>.cpp) + target_include_directories。
# 约定: 文件名 = CMake 目标名 (本脚本靠扫目录发现, 不用 ctest 自省)。
#
# 退出码: 0 成功 | 1 配置/构建失败 | 64 参数错误
# ---------------------------------------------------------------------------

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APPS_DIR="$SCRIPT_DIR/demo/apps"
BUILD_DIR="$SCRIPT_DIR/build"
JOBS="$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)"
CLEAN=0; WARN=0; APP=""; PASS_ARGS=()

if [ -t 1 ]; then
    C_STEP=$'\033[1;34m'; C_OK=$'\033[32m'; C_ERR=$'\033[31m'; C_OFF=$'\033[0m'
else
    C_STEP=''; C_OK=''; C_ERR=''; C_OFF=''
fi
step() { printf '\n%s==> %s%s\n' "$C_STEP" "$*" "$C_OFF"; }
ok()   { printf '%s%s%s\n' "$C_OK" "$*" "$C_OFF"; }
err()  { printf '%s%s%s\n' "$C_ERR" "$*" "$C_OFF" >&2; }

usage() {
    # 打印文件头两条分隔线之间的说明块 (避免行号硬编码)
    awk '
        /^# -{10,}/ { n++; next }
        n == 1 && /^#/ { sub(/^# ?/, ""); print; next }
        n >= 2 { exit }
    ' "${BASH_SOURCE[0]}"
}

list_apps() {
    [ -d "$APPS_DIR" ] || return 0
    for f in "$APPS_DIR"/*.cpp; do
        [ -e "$f" ] || continue          # 无匹配时 glob 保持字面量
        basename "$f" .cpp
    done
}

while [ $# -gt 0 ]; do
    case "$1" in
        -c|--clean)     CLEAN=1 ;;
        -w|--warn)      WARN=1 ;;
        -j)             [ $# -ge 2 ] || { err "-j 需要参数 (并行度)"; exit 64; }; JOBS="$2"; shift ;;
        -j*)            JOBS="${1#-j}" ;;
        -b|--build-dir) [ $# -ge 2 ] || { err "-b 需要参数 (构建目录)"; exit 64; }; BUILD_DIR="$2"; shift ;;
        -b*)            BUILD_DIR="${1#-b}" ;;
        -h|--help)      usage; exit 0 ;;
        --)             shift; PASS_ARGS=("$@"); break ;;
        -*)             err "未知参数: $1"; echo; usage; exit 64 ;;
        *)              if [ -n "$APP" ]; then err "只接受一个 app 名 (多余: $1)"; exit 64; fi
                        APP="$1" ;;
    esac
    shift
done

# 构建目录规范化为绝对路径 (相对路径按当前 cwd 解释)
case "$BUILD_DIR" in
    /*) ;;
    *)  BUILD_DIR="$PWD/$BUILD_DIR" ;;
esac

# --- 0) 不带名字: 只列清单 --------------------------------------------------
if [ -z "$APP" ]; then
    step "demo/apps 下可运行的 app"
    n=0
    while IFS= read -r a; do
        printf '  %s\n' "$a"
        n=$((n + 1))
    done < <(list_apps)
    if [ "$n" -eq 0 ]; then
        printf '  (空) 在 demo/apps/ 放 <name>.cpp 并在 demo/apps/CMakeLists.txt 登记\n'
    fi
    printf '\n用法: ./run_apps.sh <name>    (共 %d 个)\n' "$n"
    exit 0
fi

# --- 1) 校验 app 存在 ------------------------------------------------------
if [ ! -f "$APPS_DIR/$APP.cpp" ]; then
    err "没有这个 app: $APP"
    printf '可用:' >&2
    while IFS= read -r a; do printf ' %s' "$a" >&2; done < <(list_apps)
    printf '\n' >&2
    exit 64
fi

LOG="$(mktemp -t graphics-apps.XXXXXX)"
trap 'rm -f "$LOG"' EXIT
tail_log() { tail -30 "$LOG"; }

# --- 2) 清理 (可选) --------------------------------------------------------
if [ "$CLEAN" = 1 ]; then
    case "$BUILD_DIR" in
        ""|"/"|".") err "拒绝清理构建目录: '$BUILD_DIR'"; exit 64 ;;
    esac
    step "清理构建目录 $BUILD_DIR"
    rm -rf "$BUILD_DIR"
fi

# --- 3) 配置 ---------------------------------------------------------------
CMAKE_ARGS=(-S "$SCRIPT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Debug)
if [ "$WARN" = 1 ]; then
    CMAKE_ARGS+=(-DCMAKE_CXX_FLAGS="-Wall -Wextra")
fi

step "配置 (cmake ${CMAKE_ARGS[*]})"
if ! cmake "${CMAKE_ARGS[@]}" >"$LOG" 2>&1; then
    err "cmake 配置失败:"; tail_log; exit 1
fi
ok "配置 OK"

# --- 4) 构建 (只构建这一个 target; EXCLUDE_FROM_ALL 要求显式指定) ----------
step "构建目标 $APP (cmake --build --target $APP -j $JOBS)"
if ! cmake --build "$BUILD_DIR" --target "$APP" -j "$JOBS" >"$LOG" 2>&1; then
    err "构建失败:"; tail_log; exit 1
fi
ok "构建 OK"

BIN="$BUILD_DIR/demo/apps/$APP"
if [ ! -x "$BIN" ]; then
    err "构建产物不存在: $BIN"; exit 1
fi

# --- 5) 运行 ---------------------------------------------------------------
step "运行 $APP${PASS_ARGS[0]+ (带参)}"
# bash 3.2 + set -u 下空数组展开会报 unbound, 用 ${arr[@]+...} 惯用法
"$BIN" ${PASS_ARGS[@]+"${PASS_ARGS[@]}"}
rc=$?
if [ "$rc" -eq 0 ]; then
    ok "退出码 0"
else
    err "退出码 $rc"
fi
exit "$rc"
