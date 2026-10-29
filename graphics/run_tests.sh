#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# graphics/run_tests.sh -- C++ 测试一键入口 (配置 + 构建 + ctest)
#
#   ./run_tests.sh              增量构建并跑全部用例
#   ./run_tests.sh -c           删掉构建目录从零构建 (怀疑缓存时用)
#   ./run_tests.sh -d           跑完后再逐用例打印内部检查明细 (末行计数)
#   ./run_tests.sh -w           打开 -Wall -Wextra 编译告警
#   ./run_tests.sh -j 4         并行度 (默认: 本机核数)
#   ./run_tests.sh -b /tmp/gb   指定构建目录 (默认: graphics/build)
#   ./run_tests.sh -l           只列出已注册用例, 不构建不跑
#   ./run_tests.sh --build-only 只配置 + 构建, 不跑测试
#   ./run_tests.sh --test-only  只跑测试, 不重新配置/构建 (需先构建过)
#   ./run_tests.sh --coverage   带覆盖率构建并跑测试, 出 gcovr 报告
#                               (默认构建目录换成 build-cov, 不污染日常构建;
#                                分母只算 src/, 见 demo/README.md §5)
#   ./run_tests.sh -h           帮助
#
# 退出码: 0 全绿 | 1 配置/构建失败 | 2 有用例失败 | 64 参数错误
#
# 新增用例的做法: 在 demo/<子目录>/ 放 xxx.cpp + CMakeLists.txt
# (照抄 demo/core/CMakeLists.txt 的 add_executable + add_test 四行),
# 再在 demo/CMakeLists.txt 里 add_subdirectory(子目录)。本脚本自动带出。
# 只做演示、不需要进 ctest 的程序放 demo/apps/ (附加区), 用 ./run_apps.sh 跑。
# ---------------------------------------------------------------------------

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

BUILD_DIR="$SCRIPT_DIR/build"
JOBS="$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)"
CLEAN=0; DETAIL=0; WARN=0; LIST_ONLY=0; BUILD_ONLY=0; TEST_ONLY=0
COVERAGE=0; BUILD_DIR_GIVEN=0

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

while [ $# -gt 0 ]; do
    case "$1" in
        -c|--clean)     CLEAN=1 ;;
        -d|--detail)    DETAIL=1 ;;
        -w|--warn)      WARN=1 ;;
        -l|--list)      LIST_ONLY=1 ;;
        --build-only)   BUILD_ONLY=1 ;;
        --test-only)    TEST_ONLY=1 ;;
        --coverage)     COVERAGE=1 ;;
        -j)             [ $# -ge 2 ] || { err "-j 需要参数 (并行度)"; exit 64; }; JOBS="$2"; shift ;;
        -j*)            JOBS="${1#-j}" ;;
        -b|--build-dir) [ $# -ge 2 ] || { err "-b 需要参数 (构建目录)"; exit 64; }; BUILD_DIR_GIVEN=1; BUILD_DIR="$2"; shift ;;
        -b*)            BUILD_DIR_GIVEN=1; BUILD_DIR="${1#-b}" ;;
        -h|--help)      usage; exit 0 ;;
        *)              err "未知参数: $1"; echo; usage; exit 64 ;;
    esac
    shift
done

# --coverage 默认用独立构建目录: 插桩 + -O0 会污染日常构建的 .o
if [ "$COVERAGE" = 1 ] && [ "$BUILD_DIR_GIVEN" = 0 ]; then
    BUILD_DIR="$SCRIPT_DIR/build-cov"
    # 并且**默认从零重建**。.gcda 会跨次累积; 一旦被测代码或框架变过,
    # 复用旧目录会给出失真数字 —— 2026-10-30 实测: 框架迁移后
    # src/quat/quaternion.hpp 的覆盖率从 96% 「掉」成 5%, 数字是假的。
    # 想增量跑覆盖率就显式给 -b <dir>, 那种情况由你自己负责口径。
    CLEAN=1
fi

# 构建目录规范化为绝对路径 (相对路径按当前 cwd 解释)
case "$BUILD_DIR" in
    /*) ;;
    *)  BUILD_DIR="$PWD/$BUILD_DIR" ;;
esac

LOG="$(mktemp -t graphics-test.XXXXXX)"
trap 'rm -f "$LOG"' EXIT
tail_log() { tail -30 "$LOG"; }

if [ "$BUILD_ONLY" = 1 ] && [ "$TEST_ONLY" = 1 ]; then
    err "--build-only 与 --test-only 互斥"
    exit 64
fi

if [ "$TEST_ONLY" = 1 ] && [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
    err "构建目录尚未配置: $BUILD_DIR (先跑一次 ./run_tests.sh)"
    exit 1
fi

# --- 0) 只列用例 -----------------------------------------------------------
if [ "$LIST_ONLY" = 1 ]; then
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        err "构建目录尚未配置: $BUILD_DIR (先跑一次 ./run_tests.sh)"
        exit 1
    fi
    ctest --test-dir "$BUILD_DIR" -N
    exit $?
fi

# --- 1) 清理 ---------------------------------------------------------------
if [ "$TEST_ONLY" = 0 ] && [ "$CLEAN" = 1 ]; then
    case "$BUILD_DIR" in
        ""|"/"|".") err "拒绝清理构建目录: '$BUILD_DIR'"; exit 64 ;;
    esac
    step "清理构建目录 $BUILD_DIR"
    rm -rf "$BUILD_DIR"
fi

if [ "$TEST_ONLY" = 0 ]; then

# --- 2) 配置 ---------------------------------------------------------------
EXTRA_CXX_FLAGS=""
if [ "$WARN" = 1 ]; then EXTRA_CXX_FLAGS="$EXTRA_CXX_FLAGS -Wall -Wextra"; fi
# -O0 -fno-inline: 让行号归属尽量贴近源码, 否则内联会把判据行算到别的行上
if [ "$COVERAGE" = 1 ]; then EXTRA_CXX_FLAGS="$EXTRA_CXX_FLAGS -O0 -g --coverage -fno-inline"; fi

CMAKE_ARGS=(-S "$SCRIPT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Debug)
if [ -n "$EXTRA_CXX_FLAGS" ]; then
    CMAKE_ARGS+=(-DCMAKE_CXX_FLAGS="${EXTRA_CXX_FLAGS# }")
fi
if [ "$COVERAGE" = 1 ]; then
    CMAKE_ARGS+=(-DCMAKE_EXE_LINKER_FLAGS="--coverage")
fi

step "配置 (cmake ${CMAKE_ARGS[*]})"
if ! cmake "${CMAKE_ARGS[@]}" >"$LOG" 2>&1; then
    err "cmake 配置失败:"; tail_log; exit 1
fi
ok "配置 OK"

# --- 3) 构建 ---------------------------------------------------------------
step "构建 (cmake --build -j $JOBS)"
if ! cmake --build "$BUILD_DIR" -j "$JOBS" >"$LOG" 2>&1; then
    err "构建失败:"; tail_log; exit 1
fi
n_warn="$(grep -ci 'warning:' "$LOG" || true)"
n_err="$(grep -ci 'error:' "$LOG" || true)"
ok "构建 OK (告警 $n_warn, 错误 $n_err)"
if [ "$WARN" = 1 ] && [ "$n_warn" -gt 0 ]; then
    grep -i 'warning:' "$LOG" | head -10
fi

fi   # TEST_ONLY = 0: 配置 + 构建至此结束

if [ "$BUILD_ONLY" = 1 ]; then
    ok "只构建完成 (未跑测试)"
    exit 0
fi

# --- 4) 测试 ---------------------------------------------------------------
rc_cov=0
step "测试 (ctest --output-on-failure)"
ctest --test-dir "$BUILD_DIR" --output-on-failure
rc_test=$?

# --- 4.5) 覆盖率 (可选) ----------------------------------------------------
# 口径: 分母只算 src/。demo/ 是判据本身, 算进去会把"测试写了多少行"当成覆盖率。
# **不要把这个百分比当质量指标**: 模板未实例化的分支天然是 0%,
# 正确用法是看 demo/README.md §5 的「未解释行」——每一行未覆盖必须能说出理由。
if [ "$COVERAGE" = 1 ]; then
    step "覆盖率 (gcovr; 分母只算 src/)"
    GCOVR_BIN="$(command -v gcovr || true)"
    if [ -z "$GCOVR_BIN" ] && [ -x "$HOME/.local/bin/gcovr" ]; then
        GCOVR_BIN="$HOME/.local/bin/gcovr"
    fi
    if [ -z "$GCOVR_BIN" ]; then
        err "找不到 gcovr。先跑 ./graphics/scripts/fetch_deps.sh"
        exit 1
    fi

    mkdir -p "$BUILD_DIR/coverage"
    if ( cd "$BUILD_DIR" && "$GCOVR_BIN" \
            --root "$SCRIPT_DIR" \
            --filter "$SCRIPT_DIR/src/" \
            --exclude-unreachable-branches \
            --txt --txt-metric line \
            --print-summary \
            --html-details "coverage/index.html" ); then
        ok "覆盖率报告: $BUILD_DIR/coverage/index.html"
    else
        err "gcovr 失败 —— 常见原因: AppleClang 的 .gcda 与 gcov 不兼容。"
        err "排查: cd $BUILD_DIR && gcovr --root $SCRIPT_DIR --filter $SCRIPT_DIR/src/ -v"
        rc_cov=1
    fi
fi

# --- 5) 明细 (可选) --------------------------------------------------------
if [ "$DETAIL" = 1 ]; then
    step "各用例内部检查明细"
    while IFS= read -r exe; do
        [ -n "$exe" ] && [ -x "$exe" ] || continue
        out="$("$exe" 2>&1)"; rc=$?
        # 两套框架的输出格式不同, 两个都要认:
        #   老用例(Checker):  "  [OK]   xxx" / "  [FAIL] xxx"
        #   gtest:            "[       OK ] Suite.Test (0 ms)" / "[  FAILED  ] Suite.Test"
        # 只写 '\[OK\]' 会让 gtest 用例恒为 0 项 —— 静默失效, 且看起来"跑过了"。
        n_ok="$(printf '%s\n' "$out" | grep -cE '\[[[:space:]]*OK[[:space:]]*\]')"
        n_bad="$(printf '%s\n' "$out" | grep -cE '\[[[:space:]]*(FAIL|FAILED)[[:space:]]*\]')"
        # 两套框架数的东西**不是一回事**, 不能混在一列:
        #   老 Checker  : 每行 [OK] 是一条断言  -> 数的是「断言数」
        #   gtest       : 每级 [       OK ] 是一个 TEST() -> 数的是「用例数」
        # 拿 gtest 的用例数去和老的断言数比较会得出「迁移后判据变少了」的错误结论。
        if printf '%s\n' "$out" | grep -qE '^\[==========\]'; then
            metric="用例 ${n_ok} 个 (FAIL ${n_bad})"
        else
            metric="断言 ${n_ok} 项 (FAIL ${n_bad})"
        fi
        printf '  %-16s rc=%d  %-26s %s\n' \
            "$(basename "$exe")" "$rc" "$metric" \
            "$(printf '%s\n' "$out" | tail -1)"
    done < <(ctest --test-dir "$BUILD_DIR" -N -V 2>/dev/null | sed -n 's/^.*Test command: //p')
fi

# --- 6) 汇总 ---------------------------------------------------------------
# 覆盖率失败也必须反映到退出码 —— 否则"要求了覆盖率却没生成"会被当成绿。
if [ "$rc_test" -ne 0 ]; then
    err "有用例失败 (ctest exit=$rc_test)"
    exit 2
fi
if [ "$rc_cov" -ne 0 ]; then
    err "用例通过, 但覆盖率报告未生成 (见上面的排查提示)"
    exit 1
fi
ok "全部用例通过"
exit 0
