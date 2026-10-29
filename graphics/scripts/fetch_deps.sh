#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# graphics/scripts/fetch_deps.sh —— graphics 全部外部依赖的一键准备
#
#   ./fetch_deps.sh             下载 C++ 依赖 + 装 Python 工具 + 同步 Python 库
#   ./fetch_deps.sh --verify    只校验（不下载、不安装），缺失/校验不符则非 0 退出
#   ./fetch_deps.sh --force     强制重新下载 C++ 依赖
#   ./fetch_deps.sh --no-python 跳过 Python 部分（只想跑 C++ 时）
#   ./fetch_deps.sh --clean     删掉本清单管理的 third_party 子目录（不碰 imgui）
#   ./fetch_deps.sh -h          帮助
#
# 退出码：0 全部就绪 | 1 下载/校验/安装失败 | 64 参数错误
#
# 设计约束（不是偏好，是事故换来的）：
#   1. **构建期不联网**。CMake 里不用 FetchContent，一律"先跑到位再构建"。
#      否则 `cmake` 的成败取决于网络，且同一份源码在别人机器上可能编出不同东西。
#   2. **pin sha256**，不只 pin tag。tag 可被上游移动，sha256 不能。
#   3. **幂等**：靠 third_party/<name>/.deps-stamp 判断，重复跑是空操作。
#   4. **不碰非本清单目录**：third_party/imgui/ 是手工放的，脚本对它零操作。
# ---------------------------------------------------------------------------

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GRAPHICS_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$GRAPHICS_DIR/.." && pwd)"
MANIFEST="$GRAPHICS_DIR/deps.manifest"
THIRD_PARTY="$GRAPHICS_DIR/third_party"
CACHE="$THIRD_PARTY/.cache"

FORCE=0; VERIFY_ONLY=0; NO_PYTHON=0; CLEAN=0

if [ -t 1 ]; then
    C_STEP=$'\033[1;34m'; C_OK=$'\033[32m'; C_ERR=$'\033[31m'; C_OFF=$'\033[0m'
else
    C_STEP=''; C_OK=''; C_ERR=''; C_OFF=''
fi
step() { printf '\n%s==> %s%s\n' "$C_STEP" "$*" "$C_OFF"; }
ok()   { printf '%s  %s%s\n' "$C_OK" "$*" "$C_OFF"; }
err()  { printf '%s  %s%s\n' "$C_ERR" "$*" "$C_OFF" >&2; }

usage() {
    awk '
        /^# -{10,}/ { n++; next }
        n == 1 && /^#/ { sub(/^# ?/, ""); print; next }
        n >= 2 { exit }
    ' "${BASH_SOURCE[0]}"
}

while [ $# -gt 0 ]; do
    case "$1" in
        --force)     FORCE=1 ;;
        --verify)    VERIFY_ONLY=1 ;;
        --no-python) NO_PYTHON=1 ;;
        --clean)     CLEAN=1 ;;
        -h|--help)   usage; exit 0 ;;
        *)           err "未知参数: $1"; echo; usage; exit 64 ;;
    esac
    shift
done

[ -f "$MANIFEST" ] || { err "依赖清单不存在: $MANIFEST"; exit 1; }

# sha256 工具：macOS 是 shasum，Linux 是 sha256sum
if command -v shasum >/dev/null 2>&1; then
    sha256_of() { shasum -a 256 "$1" | awk '{print $1}'; }
elif command -v sha256sum >/dev/null 2>&1; then
    sha256_of() { sha256sum "$1" | awk '{print $1}'; }
else
    err "找不到 shasum / sha256sum，无法校验"; exit 1
fi

mkdir -p "$THIRD_PARTY" "$CACHE"

# ---------------------------------------------------------------------------
# --clean：只删清单里点名的目录（imgui 不在清单里，因此安全）
# ---------------------------------------------------------------------------
if [ "$CLEAN" = 1 ]; then
    step "清理清单管理的第三方目录"
    while read -r name version hash url; do
        case "$name" in ''|\#*) continue ;; esac
        [ -n "${url:-}" ] || continue          # 2 字段行是 Python 工具，跳过
        if [ -d "$THIRD_PARTY/$name" ]; then
            rm -rf "${THIRD_PARTY:?}/$name"
            ok "已删除 third_party/$name"
        fi
    done < "$MANIFEST"
    ok "清理完成（third_party/imgui 未触碰）"
    exit 0
fi

rc=0
n_cpp=0; n_py=0

# ---------------------------------------------------------------------------
# 1) C++ 源码包
# ---------------------------------------------------------------------------
step "C++ 第三方依赖（${THIRD_PARTY}）"

while read -r name version hash url; do
    case "$name" in ''|\#*) continue ;; esac
    [ -n "${url:-}" ] || continue
    n_cpp=$((n_cpp + 1))

    dest="$THIRD_PARTY/$name"
    stamp="$dest/.deps-stamp"
    want="$version $hash"

    if [ ! -f "$stamp" ] || [ "$(cat "$stamp" 2>/dev/null)" != "$want" ]; then
        :   # 需要下载
    else
        ok "$name $version —— 已就绪（stamp 匹配）"
        continue
    fi

    if [ "$VERIFY_ONLY" = 1 ]; then
        err "$name $version —— 缺失或版本/校验值不匹配"
        rc=1
        continue
    fi

    tarball="$CACHE/$name-$version.tar.gz"
    if [ "$FORCE" = 1 ] || [ ! -f "$tarball" ]; then
        printf '  下载 %s %s ...\n' "$name" "$version"
        if ! curl -fsSL --max-time 300 -o "$tarball.part" "$url"; then
            err "$name —— 下载失败：$url"
            rm -f "$tarball.part"; rc=1; continue
        fi
        mv "$tarball.part" "$tarball"
    fi

    got="$(sha256_of "$tarball")"
    if [ "$got" != "$hash" ]; then
        err "$name —— sha256 不匹配"
        err "    期望 $hash"
        err "    实得 $got"
        err "    （上游 tag 可能被移动，或下载被替换。确认无误后更新 deps.manifest）"
        rm -f "$tarball"; rc=1; continue
    fi
    ok "$name $version —— sha256 校验通过"

    printf '  解包 ...\n'
    tmp="$(mktemp -d "${TMPDIR:-/tmp}/deps.XXXXXX")"
    if ! tar xzf "$tarball" -C "$tmp"; then
        err "$name —— 解包失败"; rm -rf "$tmp"; rc=1; continue
    fi
    # 压缩包顶层目录名通常是 <name>-<version> 或 <name>-<version-with-v>
    top="$(find "$tmp" -mindepth 1 -maxdepth 1 -type d | head -1)"
    if [ -z "$top" ]; then
        err "$name —— 压缩包结构异常（无顶层目录）"; rm -rf "$tmp"; rc=1; continue
    fi
    rm -rf "$dest"
    mv "$top" "$dest"
    rm -rf "$tmp"
    printf '%s' "$want" > "$stamp"
    ok "$name $version —— 安装到 third_party/$name"
done < "$MANIFEST"

if [ "$n_cpp" -eq 0 ]; then
    err "清单里没有 C++ 依赖，检查 $MANIFEST"
fi

# ---------------------------------------------------------------------------
# 2) Python 工具 + 3) Python 库
# ---------------------------------------------------------------------------
if [ "$NO_PYTHON" = 1 ]; then
    step "跳过 Python 部分（--no-python）"
elif [ "$VERIFY_ONLY" = 1 ]; then
    step "Python 工具（只校验）"
    have="$(uv tool list 2>/dev/null || true)"
    while read -r name version hash url; do
        case "$name" in ''|\#*) continue ;; esac
        [ -n "${version:-}" ] || continue
        [ -z "${hash:-}" ] || continue          # 4 字段行是 C++ 依赖，不是 Python 工具
        if printf '%s' "$have" | grep -q "^$name "; then
            ok "$name —— 已安装"
        else
            err "$name —— 未安装（跑 ./fetch_deps.sh）"; rc=1
        fi
    done < "$MANIFEST"
else
    step "Python 工具（uv tool install）"
    while read -r name version hash url; do
        case "$name" in ''|\#*) continue ;; esac
        [ -n "${version:-}" ] || continue
        [ -z "${hash:-}" ] || continue          # 4 字段行是 C++ 依赖，不是 Python 工具
        n_py=$((n_py + 1))
        have="$(uv tool list 2>/dev/null || true)"
        if printf '%s' "$have" | grep -q "^$name "; then
            ok "$name —— 已安装，跳过（升级请 --force）"
            continue
        fi
        if uv tool install "$name==$version" >/dev/null 2>&1; then
            ok "$name $version —— 安装完成"
        else
            err "$name $version —— 安装失败"; rc=1
        fi
    done < "$MANIFEST"

    step "Python 库（父仓 uv sync）"
    if [ -f "$REPO_ROOT/pyproject.toml" ]; then
        if ( cd "$REPO_ROOT" && uv sync >/dev/null 2>&1 ); then
            ok "uv sync 完成（${REPO_ROOT}）"
        else
            err "uv sync 失败，手动跑一次看报错：cd $REPO_ROOT && uv sync"
            rc=1
        fi
    else
        ok "父仓无 pyproject.toml，跳过"
    fi
fi

# ---------------------------------------------------------------------------
# 汇总
# ---------------------------------------------------------------------------
echo
if [ "$rc" -eq 0 ]; then
    ok "全部依赖就绪（C++ $n_cpp 项，Python 工具 $n_py 项）"
    ok "现在可以构建：./graphics/run_tests.sh"
    exit 0
else
    err "有依赖未就绪，构建会失败。按上面的提示处理后重跑本脚本。"
    exit 1
fi
