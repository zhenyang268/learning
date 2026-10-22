"""
G1 · 大圆 vs 等角航线：最短路径不等于看起来直的线
================================================================
现象：北京 → 纽约，大圆航线（最短）在经纬度展开图上是一条**弯向北方的弧**，
      而"图上直"的等角航线（rhumb line）反而更长 —— 飞机绕北极飞不是绕远。

一句话在说什么：
    地图（墨卡托投影）为了保方向角，把纬线越往高纬拉得越开 ——
    高纬地区的长度在图上被**放大**了。所以球面上"真直"的大圆
    在图上显示为弯的，而图上直的恒向线在球面上是一条越走越偏北的螺旋。
    **投影会撒谎，度规不会** —— 这正是"度规定长度"的第一课。

数值内容：
    1. haversine 大圆距离（与几何法对拍）
    2. 等角航线长度（墨卡托空间的精确公式）
    3. 两者的 3D 轨迹与平面展开对比
    4. 墨卡托平面上恒向线是**直线**（斜率恒定）的实测证明 —— 与前一张图的
       "在 lat-lon 平面上不是直线"正好互为镜像

运行： uv run python graphics/demo/pyapps/sphere/great_circle.py
输出： learning-todo/imgs/geo/great_circle_vs_rhumb.png
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _core as C  # noqa: E402

plt = C.setup_style()

R_EARTH = 6371.0               # 地球平均半径 km
BEIJING = (39.9, 116.4)        # (纬度, 经度) 度
NEWYORK = (40.7, -74.0)


def latlon_to_cart(lat_deg, lon_deg, R=R_EARTH):
    """地理纬度（赤道=0，北正）+ 经度 -> 笛卡尔。注意与极角 θ 互余。"""
    return np.stack(C.sph_to_cart(np.pi / 2 - np.radians(lat_deg),
                                  np.radians(lon_deg), R), axis=-1)


def latlon_path_to_cart(lat_deg, lon_deg, R=R_EARTH):
    """点列版：整条路径一次转成笛卡尔 (N,3)，避免逐点 stack 出错。"""
    return np.stack(C.sph_to_cart(np.pi / 2 - np.radians(lat_deg),
                                  np.radians(lon_deg), R), axis=1)


def mercator_y(lat_deg):
    """墨卡托纵坐标 ψ = ln tan(π/4 + φ/2)：把纬度"拉高"的映射。"""
    lat = np.radians(lat_deg)
    return np.log(np.tan(np.pi / 4 + lat / 2))


def rhumb_line(lat1, lon1, lat2, lon2, R=R_EARTH, n=500):
    """等角航线（恒向线）：墨卡托图上的直线。返回 (点列 (n,3), 长度 km, 经度序列)。

    长度精确公式：
        Δψ = ψ2 − ψ1,  Δλ = λ2 − λ1（取短方向）
        q = Δφ/Δψ（Δψ≈0 时取 cosφ）
        d = R·√(Δφ² + q²Δλ²)

    长度用**球面度规**算（ds² = R²dφ² + R²cos²φ dλ²），不是把图上长度加起来 ——
    图上长度是被放大过的，这正是本脚本要教学的点。
    """
    dlon = np.radians(lon2 - lon1)
    if dlon > np.pi:
        dlon -= 2 * np.pi
    if dlon < -np.pi:
        dlon += 2 * np.pi
    dphi = np.radians(lat2 - lat1)
    dpsi = mercator_y(lat2) - mercator_y(lat1)
    q = dphi / dpsi if abs(dpsi) > 1e-12 else np.cos(np.radians(lat1))
    dist = R * np.hypot(dphi, q * dlon)

    # 轨迹：墨卡托空间直线插值 -> 逆映射回球面
    t = np.linspace(0, 1, n)
    psis = mercator_y(lat1) + t * dpsi
    lons = lon1 + np.degrees(t * dlon)
    lats = np.degrees(2 * np.arctan(np.exp(psis)) - np.pi / 2)   # ψ -> φ
    path = latlon_path_to_cart(lats, lons, R)
    return path, dist, lats, lons


def main():
    lat1, lon1 = BEIJING
    lat2, lon2 = NEWYORK
    p1 = latlon_to_cart(lat1, lon1)
    p2 = latlon_to_cart(lat2, lon2)

    print("=" * 70)
    print(f"【航线对比】北京({lat1}N,{lon1}E) → 纽约({lat2}N,{abs(lon2)}W)")
    print("=" * 70)

    # 大圆：短弧点列（单位球）+ 中心角
    gc_unit, gc_ang = C.great_circle_path(p1, p2)
    gc_len = gc_ang * R_EARTH
    hav = C.haversine(lat1, lon1, lat2, lon2)
    err_hav = abs(gc_len - hav)
    print(f"  大圆距离（几何法）   = {gc_len:10.1f} km")
    print(f"  大圆距离（haversine）= {hav:10.1f} km   两者差 = {err_hav:.2e} km")
    print("  → haversine 就是大圆中心角公式的小距离稳定版，两个算法对拍。")

    # 短弧判定：终点必须落在目标点上（单位球）
    p2n = p2 / np.linalg.norm(p2)
    err_end = np.linalg.norm(gc_unit[-1] - p2n)
    print(f"  短弧判定：路径终点与纽约偏差 = {err_end:.2e}  （原版此处拿未归一化的 p2")
    print("    相比，距离恒 ≈ 6371，判据形同虚设；现在拿归一化点比，才真的能择向）")

    rh_path, rh_len, rh_lat, rh_lon = rhumb_line(lat1, lon1, lat2, lon2)
    print(f"  等角航线长度         = {rh_len:10.1f} km")
    print(f"  等角航线多飞         = {rh_len - gc_len:10.1f} km "
          f"（约 {(rh_len / gc_len - 1) * 100:.1f}%）")
    ok = (err_hav < 1e-6) and (err_end < 1e-9) and (rh_len > gc_len)
    print(f"  [{'PASS' if ok else 'FAIL'}]  三条判据：haversine 对拍 < 1e-6 km、"
          f"短弧择向 < 1e-9、恒向线确实更长")

    # 大圆最高点
    i_max = int(np.argmax(gc_unit[:, 2]))
    th_m, ph_m = C.cart_to_sph(gc_unit[i_max])
    print(f"  大圆最高到达纬度     = {90 - np.degrees(th_m):6.1f}°N "
          f"（出发地只有 ~{lat1}°N —— 为了抄近路，飞机往北极方向绕）")

    # 墨卡托平面：恒向线是直线（斜率恒定）的实测
    psi_rh = mercator_y(rh_lat)
    lon_rh_unwrapped = C.unwrap_lon_deg(rh_lon)
    slope_mc = np.gradient(lon_rh_unwrapped, psi_rh)
    slope_ll = np.gradient(lon_rh_unwrapped, rh_lat)
    rel_mc = np.ptp(slope_mc) / abs(np.median(slope_mc))
    rel_ll = np.ptp(slope_ll) / abs(np.median(slope_ll))
    print("\n" + "=" * 70)
    print("【平面性检验】一条曲线是直线 ⟺ 斜率处处相同（相对波动 → 0）")
    print("=" * 70)
    print(f"  {'平面':<22}{'斜率中位':>16}{'相对波动':>14}")
    print(f"  {'经纬度平面 dλ/dφ':<22}{np.median(slope_ll):>16.4f}{rel_ll:>14.2e}")
    print(f"  {'墨卡托平面 dλ/dψ':<20}{np.median(slope_mc):>16.4f}{rel_mc:>14.2e}"
          f"  ← 舍入级，即恒为常数")
    print("  → 墨卡托把纬度做了非线性拉伸 z=ψ(φ)，正好把恒向线拉成直线 ——")
    print("    代价是**长度被拉伸**（高纬放大），所以 3D 里它比大圆长 30.9%。")

    # ---------- 画图 ----------
    fig = plt.figure(figsize=(16, 5.8))

    # (a) 3D 球面（全部归一到单位球，路径与底球同尺度）
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    C.draw_sphere_base(ax)
    ax.plot(*gc_unit.T, color="#c0392b", lw=2.5,
            label=f"大圆（最短）{gc_len:,.0f} km")
    ax.plot(*(rh_path / R_EARTH).T, color="#2980b9", lw=2,
            label=f"等角航线 {rh_len:,.0f} km")
    b1 = p1 / np.linalg.norm(p1)
    b2 = p2 / np.linalg.norm(p2)
    ax.scatter(*b1, color="k", s=40)
    ax.scatter(*b2, color="k", s=40)
    ax.text(*b1 * 1.2, "北京", fontsize=11)
    ax.text(*b2 * 1.2, "纽约", fontsize=11)
    ax.view_init(elev=32, azim=30)
    ax.set_box_aspect((1, 1, 1))
    ax.set_title("① 球面视角：大圆（红）比恒向线（蓝）短")
    ax.legend(fontsize=9, loc="lower left")

    # (b) 经纬度平面展开（经度必须解卷，否则跨 180° 拉假跳线）
    ax = fig.add_subplot(1, 3, 2)
    th_gc, ph_gc = C.cart_to_sph(gc_unit)
    lat_gc = 90 - np.degrees(th_gc)
    lon_gc = C.unwrap_lon_deg(np.degrees(ph_gc))
    # 把整段平移到与北京(116)同一支，避免和 rhumb 的经度区间错开
    lon_gc = lon_gc - 360 * round((lon_gc[0] - lon1) / 360)
    ax.plot(lon_gc, lat_gc, color="#c0392b", lw=2.5, label="大圆：地图上的弯路")
    ax.plot(lon_rh_unwrapped, rh_lat, color="#2980b9", lw=2,
            label="等角航线：地图上的直路")
    ax.scatter([lon1, lon2], [lat1, lat2], color="k", s=35, zorder=5)
    ax.text(lon1 + 4, lat1 - 3, "北京", fontsize=10)
    ax.text(lon2 - 14, lat2 + 2, "纽约", fontsize=10)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("经度（度）")
    ax.set_ylabel("纬度（度）")
    ax.set_title("② 经纬度平面：短的弯，直的贵\n（高纬在图上被拉宽）")
    ax.legend(fontsize=9, loc="lower center")
    ax.grid(alpha=0.3)
    ax.set_ylim(30, 90)

    # (c) 墨卡托平面：恒向线是直线
    ax = fig.add_subplot(1, 3, 3)
    psi_gc = mercator_y(lat_gc)
    ax.plot(lon_gc, psi_gc * 180 / np.pi, color="#c0392b", lw=2.5,
            label="大圆：在有拉伸的图上弯了")
    ax.plot(lon_rh_unwrapped, psi_rh * 180 / np.pi, color="#2980b9", lw=2,
            label="等角航线：**直线**")
    ax.scatter([lon1, lon2], [mercator_y(lat1) * 180 / np.pi,
                              mercator_y(lat2) * 180 / np.pi],
               color="k", s=35, zorder=5)
    ax.set_xlabel("经度（度）")
    ax.set_ylabel("墨卡托纵坐标 ψ（度）")
    ax.set_title("③ 墨卡托平面：蓝线是直线\n（斜率恒定 —— 航海图的历史理由）")
    ax.legend(fontsize=9, loc="lower center")
    ax.grid(alpha=0.3)

    fig.suptitle("G1 · 大圆 vs 等角航线：为什么飞机往北极绕", fontsize=13.5, y=0.985)
    fig.subplots_adjust(top=0.80, wspace=0.26)
    path = C.savefig(fig, "great_circle_vs_rhumb.png")
    plt.close(fig)

    print("\n图：", path)
    print("\n这一步在说什么：")
    print("  「最短」由度规说了算：ds² = R²dφ² + R²cos²φ dλ² 里 cosφ 告诉你")
    print("  高纬度的经度圈越来越「瘦」，所以往高处绕反而近。")
    print("  地图投影保了角度（墨卡托 = 航海图的历史理由），代价是撒了长度的谎。")


if __name__ == "__main__":
    main()