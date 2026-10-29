"""pytest 的路径准备。

pyapps 的脚本不是包（没有 __init__.py，本来就是"能直接运行看图"的形态），
所以这里用 sys.path 把要测的目录挂进去，而不是把整个 pyapps 改造成包 ——
改造成包会让 `python sphere/G1_xxx.py` 这种用法失效。

只挂**数值核**所在目录。绘图脚本不进 pytest。
"""

import sys
from pathlib import Path

PYAPPS = Path(__file__).resolve().parents[1]      # graphics/demo/pyapps
for sub in ("sphere", "rotation"):
    p = PYAPPS / sub
    if p.is_dir() and str(p) not in sys.path:
        sys.path.insert(0, str(p))
