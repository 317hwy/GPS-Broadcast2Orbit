
# GPS_Coordination.py
# -------------------
# 本文件用于基于GPS广播星历文件进行卫星轨迹计算、精密星历误差对比及可视化。
# 主要功能包括：
#   - 轨迹计算与导出
#   - 轨迹三维可视化
#   - 广播星历与精密星历误差对比
#
# 变量/常量说明：
#   MU: 地球引力常数 (m^3/s^2)
#   OMEGA_EARTH: 地球自转角速度 (rad/s)
#   SECONDS_PER_WEEK: 一周的秒数
#   SECONDS_PER_DAY: 一天的秒数
#   HALF_WEEK: 半周的秒数
#
#   ephemeris_list: 星历字典列表，每个元素为一颗卫星的星历参数
#   prn: 卫星PRN号（如G01、G02）
#   t_sow: GPS周内秒 (seconds of week)
#   points: 轨迹点列表，每个元素为一时刻的卫星坐标
#   trajectories: 所有卫星的轨迹点字典
#   rows: 误差或轨迹等数据的行列表
#
from __future__ import annotations
import argparse
import csv
import math
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from GPS_IO import read_nav_body_file, save_nav_body_file


MU = 3.986005e14
OMEGA_EARTH = 7.2921151467e-5
SECONDS_PER_WEEK = 604800.0
SECONDS_PER_DAY = 86400.0
HALF_WEEK = 302400.0


def normalize_prn(prn: str) -> str:
    """
    归一化PRN号：去除前缀G，补零到两位。
    参数：prn - 原始PRN字符串
    返回：两位数字字符串，如'01', '02'
    """
    p = prn.strip().upper().replace("G", "")
    return p.zfill(2)


def normalize_tk(tk: float) -> float:
    """
    归一化tk（周内秒差），保证其在[-HALF_WEEK, HALF_WEEK]区间。
    参数：tk - 原始时间差
    返回：归一化后的时间差
    """
    if tk > HALF_WEEK:
        tk -= SECONDS_PER_WEEK
    elif tk < -HALF_WEEK:
        tk += SECONDS_PER_WEEK
    return tk


def gps_seconds_of_week( 
    year: int,
    month: int,
    day: int,
    hour: int,
    minute: int,
    second: int,
) -> float:
    '''
    计算指定UTC时间对应的GPS周内秒
    参数：年、月、日、时、分、秒
    返回：GPS周内秒（0~604800）
    '''
    gps_epoch = datetime(1980, 1, 6, 0, 0, 0)
    t = datetime(year, month, day, hour, minute, second)
    dt = (t - gps_epoch).total_seconds()
    return dt % SECONDS_PER_WEEK


def _filter_prn_records(ephemeris_list: Iterable[Dict[str, Any]], prn: str) -> List[Dict[str, Any]]:
    """
    从星历列表中过滤出指定PRN号的所有记录。
    参数：ephemeris_list - 星历字典列表
         prn - 目标卫星PRN号
    返回：该PRN的所有星历记录列表
    """
    target = normalize_prn(prn)
    records: List[Dict[str, Any]] = []
    for eph in ephemeris_list:
        raw_prn = eph.get("PRN")
        if raw_prn is None:
            continue
        if normalize_prn(str(raw_prn)) == target:
            records.append(eph)
    return records


def find_best_ephemeris_by_toe(
        
    

    ephemeris_list: Iterable[Dict[str, Any]],
    prn: str,
    t_sow: float,
) -> Dict[str, Any]:
    """
    查找最接近指定时刻（t_sow）的星历记录。
    参数：ephemeris_list - 星历字典列表
         prn - 目标卫星PRN号
         t_sow - 目标周内秒
    返回：包含tk/abs_tk的最佳星历字典
    """
    records = _filter_prn_records(ephemeris_list, prn)
    if not records:
        raise ValueError(f"未找到 PRN={prn} 的星历记录")

    best: Optional[Dict[str, Any]] = None
    best_abs_tk = float("inf")
    best_tk = 0.0

    for eph in records:
        toe = eph.get("toe")
        if toe is None:
            continue
        tk = normalize_tk(float(t_sow) - float(toe))
        abs_tk = abs(tk)
        if abs_tk < best_abs_tk:
            best_abs_tk = abs_tk
            best_tk = tk
            best = eph

    if best is None:
        raise ValueError(f"PRN={prn} 的记录中未读取到 toe")

    selected = dict(best)
    selected["tk"] = best_tk
    selected["abs_tk"] = best_abs_tk
    return selected


def get_toc_from_record(eph: Dict[str, Any]) -> datetime:
    """
    从星历记录中提取TOC（钟参考历元）为datetime对象。
    """
    return datetime(
        int(eph["YEAR"]),
        int(eph["MONTH"]),
        int(eph["DAY"]),
        int(eph["HOUR"]),
        int(eph["MINUTE"]),
        int(eph["SECOND"]),
    )


def _solve_kepler(Mk: float, e: float, max_iter: int = 30, tol: float = 1e-13) -> float:
    """
    解Kepler方程，返回偏近点角Ek。
    参数：Mk - 平近点角，e - 偏心率
    """
    Ek = Mk
    for _ in range(max_iter):
        f = Ek - e * math.sin(Ek) - Mk
        fp = 1.0 - e * math.cos(Ek)
        delta = -f / fp
        Ek += delta
        if abs(delta) < tol:
            break
    return Ek


def compute_satellite_ecef(eph: Dict[str, Any], t_sow: float) -> Tuple[float, float, float]:
    """
    根据星历参数和时刻计算卫星在地心地固坐标系(ECEF)下的三维坐标。
    参数：eph - 星历字典
         t_sow - 目标周内秒
    返回：(x, y, z)坐标，单位米
    """
    # tk: 卫星时刻与toe的时间差（归一化）
    # A: 轨道长半轴
    # n0: 平均角速度
    # n: 改正后角速度
    # Mk: 平近点角
    # e: 偏心率
    # Ek: 偏近点角
    # vk: 真近点角
    # phi: 轨道幅角
    # du, dr, di: 改正项
    # uk, rk, ik: 改正后参数
    # xk_orb, yk_orb: 轨道平面坐标
    # omega_k: 卫星升交点赤经
    tk = normalize_tk(t_sow - float(eph["toe"]))
    A = float(eph["sqrt_a"]) ** 2
    n0 = math.sqrt(MU / (A**3))
    n = n0 + float(eph["delta_n"])
    Mk = float(eph["m0"]) + n * tk

    e = float(eph["eccentricity"])
    Ek = _solve_kepler(Mk, e)

    sin_vk = math.sqrt(1.0 - e * e) * math.sin(Ek) / (1.0 - e * math.cos(Ek))
    cos_vk = (math.cos(Ek) - e) / (1.0 - e * math.cos(Ek))
    vk = math.atan2(sin_vk, cos_vk)

    phi = vk + float(eph["omega"])
    two_phi = 2.0 * phi

    du = float(eph["cus"]) * math.sin(two_phi) + float(eph["cuc"]) * math.cos(two_phi)
    dr = float(eph["crs"]) * math.sin(two_phi) + float(eph["crc"]) * math.cos(two_phi)
    di = float(eph["cis"]) * math.sin(two_phi) + float(eph["cic"]) * math.cos(two_phi)

    uk = phi + du
    rk = A * (1.0 - e * math.cos(Ek)) + dr
    ik = float(eph["i0"]) + float(eph["idot"]) * tk + di

    xk_orb = rk * math.cos(uk)
    yk_orb = rk * math.sin(uk)

    omega_k = (
        float(eph["omega0"])
        + (float(eph["omega_dot"]) - OMEGA_EARTH) * tk
        - OMEGA_EARTH * float(eph["toe"])
    )

    x = xk_orb * math.cos(omega_k) - yk_orb * math.cos(ik) * math.sin(omega_k)
    y = xk_orb * math.sin(omega_k) + yk_orb * math.cos(ik) * math.cos(omega_k)
    z = yk_orb * math.sin(ik)

    return x, y, z


def _day_start_sow_for_prn(ephemeris_list: Sequence[Dict[str, Any]], prn: str) -> float:
    """
    获取该PRN卫星当天0点对应的GPS周内秒。
    """
    records = _filter_prn_records(ephemeris_list, prn)
    if not records:
        raise ValueError(f"未找到 PRN={prn} 的星历记录")

    toc = get_toc_from_record(records[0])
    toc_sow = gps_seconds_of_week(
        toc.year,
        toc.month,
        toc.day,
        toc.hour,
        toc.minute,
        toc.second,
    )
    seconds_in_day = toc.hour * 3600 + toc.minute * 60 + toc.second
    return (toc_sow - seconds_in_day) % SECONDS_PER_WEEK


def generate_daily_trajectory(
    
    ephemeris_list: Sequence[Dict[str, Any]],
    prn: str,
    step_seconds: int = 60,
    t_start: int = 0,
    t_end: int = 86400,
) -> List[Dict[str, Any]]:
    """
    生成指定卫星一天内的轨迹点（ECEF坐标）。
    参数：ephemeris_list - 星历字典列表
         prn - 目标卫星PRN号
         step_seconds - 采样步长（秒）
         t_start, t_end - 起止秒数
    返回：轨迹点字典列表
    """
    if step_seconds <= 0:
        raise ValueError("step_seconds 必须为正整数")
    if t_start < 0 or t_end > 86400 or t_start > t_end:
        raise ValueError("要求 0 <= t_start <= t_end <= 86400")

    day_start_sow = _day_start_sow_for_prn(ephemeris_list, prn)
    points: List[Dict[str, Any]] = []

    for t_day in range(t_start, t_end + 1, step_seconds):
        t_sow = (day_start_sow + t_day) % SECONDS_PER_WEEK
        eph = find_best_ephemeris_by_toe(ephemeris_list, prn, t_sow)
        x, y, z = compute_satellite_ecef(eph, t_sow)
        points.append(
            {
                "PRN": f"G{normalize_prn(prn)}",
                "t_day_s": t_day,
                "t_sow_s": round(t_sow, 3),
                "X_m": round(x, 3),
                "Y_m": round(y, 3),
                "Z_m": round(z, 3),
            }
        )
    return points


def get_all_gps_prns(ephemeris_list: Sequence[Dict[str, Any]]) -> List[str]:
    """
    获取星历文件中所有GPS卫星的PRN号列表。
    """
    prns = set()
    for eph in ephemeris_list:
        raw = eph.get("PRN")
        if raw is None:
            continue
        prn_num = normalize_prn(str(raw))
        if prn_num.isdigit():
            prns.add(f"G{prn_num}")
    return sorted(prns)


def generate_all_gps_trajectories(
    ephemeris_list: Sequence[Dict[str, Any]],
    step_seconds: int = 60,
    t_start: int = 0,
    t_end: int = 86400,
) -> Dict[str, List[Dict[str, Any]]]:
    """
    计算所有GPS卫星的轨迹。
    返回：{prn: 轨迹点列表}
    """
    trajectories: Dict[str, List[Dict[str, Any]]] = {}
    for prn in get_all_gps_prns(ephemeris_list):
        try:
            trajectories[prn] = generate_daily_trajectory(
                ephemeris_list=ephemeris_list,
                prn=prn,
                step_seconds=step_seconds,
                t_start=t_start,
                t_end=t_end,
            )
        except ValueError:
            continue
    if not trajectories:
        raise ValueError("未生成任何GPS卫星轨迹，请检查广播星历文件")
    return trajectories


def save_trajectory_csv(points: Sequence[Dict[str, Any]], output_csv: Path) -> None:
    """
    将轨迹点写入CSV文件。
    """
    if not points:
        raise ValueError("没有可写入的轨迹点")

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    headers = ["PRN", "t_day_s", "t_sow_s", "X_m", "Y_m", "Z_m"]
    with output_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(points)


def plot_trajectory_3d(points: Sequence[Dict[str, Any]], prn: str, output_html: Path) -> None:
    """
    绘制单颗卫星三维轨迹并输出为HTML。
    """
    if not points:
        raise ValueError("没有可绘图的轨迹点")

    x = [p["X_m"] for p in points]
    y = [p["Y_m"] for p in points]
    z = [p["Z_m"] for p in points]

    fig = go.Figure()
    fig.add_trace(
        go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode="lines",
            name=f"G{normalize_prn(prn)}",
            line={"width": 5, "color": "#1f77b4"},
        )
    )
    fig.add_trace(
        go.Scatter3d(
            x=[x[0]],
            y=[y[0]],
            z=[z[0]],
            mode="markers",
            name="Start",
            marker={"size": 6, "color": "green"},
        )
    )
    fig.add_trace(
        go.Scatter3d(
            x=[x[-1]],
            y=[y[-1]],
            z=[z[-1]],
            mode="markers",
            name="End",
            marker={"size": 6, "color": "red"},
        )
    )

    fig.update_layout(
        title=f"GPS Satellite Trajectory in ECEF (G{normalize_prn(prn)})",
        scene={
            "xaxis_title": "X (m)",
            "yaxis_title": "Y (m)",
            "zaxis_title": "Z (m)",
            "aspectmode": "data",
        },
    )

    output_html.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(output_html), include_plotlyjs="cdn")


def _build_time_labels(points: Sequence[Dict[str, Any]], label_step: int) -> Tuple[List[float], List[float], List[float], List[str]]:
    """
    构建轨迹点的时间标签（每隔label_step个点）。
    返回：坐标与标签文本列表
    """
    if label_step <= 0:
        raise ValueError("label_step 必须为正整数")
    xs: List[float] = []
    ys: List[float] = []
    zs: List[float] = []
    texts: List[str] = []
    for i, p in enumerate(points):
        if i % label_step != 0:
            continue
        t_day = int(p["t_day_s"])
        hh = t_day // 3600
        mm = (t_day % 3600) // 60
        ss = t_day % 60
        xs.append(float(p["X_m"]))
        ys.append(float(p["Y_m"]))
        zs.append(float(p["Z_m"]))
        texts.append(f"{hh:02d}:{mm:02d}:{ss:02d}")
    return xs, ys, zs, texts


def plot_trajectory_3d_with_time_labels(
    points: Sequence[Dict[str, Any]],
    prn: str,
    output_html: Path,
    label_step: int,
) -> None:
    """
    绘制带时间标签的单颗卫星三维轨迹。
    """
    if not points:
        raise ValueError("没有可绘图的轨迹点")

    x = [p["X_m"] for p in points]
    y = [p["Y_m"] for p in points]
    z = [p["Z_m"] for p in points]
    lx, ly, lz, ltext = _build_time_labels(points, label_step)

    fig = go.Figure()
    fig.add_trace(
        go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode="lines",
            name=f"G{normalize_prn(prn)}",
            line={"width": 5, "color": "#1f77b4"},
        )
    )
    fig.add_trace(
        go.Scatter3d(
            x=lx,
            y=ly,
            z=lz,
            mode="markers+text",
            name="Time Labels",
            text=ltext,
            textposition="top center",
            marker={"size": 3, "color": "#333333"},
            textfont={"size": 10},
        )
    )

    fig.update_layout(
        title=f"GPS Satellite Trajectory in ECEF with Time Labels (G{normalize_prn(prn)})",
        scene={
            "xaxis_title": "X (m)",
            "yaxis_title": "Y (m)",
            "zaxis_title": "Z (m)",
            "aspectmode": "data",
        },
    )

    output_html.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(output_html), include_plotlyjs="cdn")


def plot_all_gps_trajectories_3d(
    trajectories: Dict[str, Sequence[Dict[str, Any]]],
    output_html: Path,
) -> None:
    """
    绘制所有GPS卫星的三维轨迹。
    """
    if not trajectories:
        raise ValueError("没有可绘制的卫星轨迹")

    fig = go.Figure()
    for prn in sorted(trajectories.keys()):
        points = trajectories[prn]
        if not points:
            continue
        x = [p["X_m"] for p in points]
        y = [p["Y_m"] for p in points]
        z = [p["Z_m"] for p in points]
        fig.add_trace(
            go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode="lines",
                name=prn,
                line={"width": 2},
            )
        )

    fig.update_layout(
        title="All GPS Satellite Trajectories in ECEF",
        scene={
            "xaxis_title": "X (m)",
            "yaxis_title": "Y (m)",
            "zaxis_title": "Z (m)",
            "aspectmode": "data",
        },
        legend={"itemsizing": "constant"},
    )
    output_html.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(output_html), include_plotlyjs="cdn")


def plot_all_gps_trajectories_3d_with_time_labels(
    trajectories: Dict[str, Sequence[Dict[str, Any]]],
    output_html: Path,
    label_step: int,
) -> None:
    """
    绘制所有GPS卫星带时间标签的三维轨迹。
    """
    if not trajectories:
        raise ValueError("没有可绘制的卫星轨迹")

    fig = go.Figure()
    for prn in sorted(trajectories.keys()):
        points = trajectories[prn]
        if not points:
            continue
        x = [p["X_m"] for p in points]
        y = [p["Y_m"] for p in points]
        z = [p["Z_m"] for p in points]
        lx, ly, lz, ltext = _build_time_labels(points, label_step)

        fig.add_trace(
            go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode="lines",
                name=prn,
                line={"width": 2},
            )
        )
        fig.add_trace(
            go.Scatter3d(
                x=lx,
                y=ly,
                z=lz,
                mode="markers+text",
                name=f"{prn}-time",
                text=ltext,
                textposition="top center",
                marker={"size": 2, "color": "#222222"},
                textfont={"size": 8},
                showlegend=False,
            )
        )

    fig.update_layout(
        title="All GPS Satellite Trajectories in ECEF with Time Labels",
        scene={
            "xaxis_title": "X (m)",
            "yaxis_title": "Y (m)",
            "zaxis_title": "Z (m)",
            "aspectmode": "data",
        },
        legend={"itemsizing": "constant"},
    )
    output_html.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(output_html), include_plotlyjs="cdn")


def flatten_trajectories(trajectories: Dict[str, Sequence[Dict[str, Any]]]) -> List[Dict[str, Any]]:
    """
    将所有卫星轨迹点展平成单一列表。
    """
    rows: List[Dict[str, Any]] = []
    for prn in sorted(trajectories.keys()):
        rows.extend(trajectories[prn])
    return rows


def compute_all_gps_positions_at_t(
    ephemeris_list: Sequence[Dict[str, Any]],
    t_day_s: int,
) -> List[Dict[str, Any]]:
    """
    计算所有GPS卫星在某一时刻的坐标。
    参数：t_day_s - 当天秒数
    返回：所有卫星坐标点列表
    """
    if t_day_s < 0 or t_day_s > 86400:
        raise ValueError("t_day_s 必须在 0~86400 秒")

    rows: List[Dict[str, Any]] = []
    for prn in get_all_gps_prns(ephemeris_list):
        try:
            day_start_sow = _day_start_sow_for_prn(ephemeris_list, prn)
            t_sow = (day_start_sow + t_day_s) % SECONDS_PER_WEEK
            eph = find_best_ephemeris_by_toe(ephemeris_list, prn, t_sow)
            x, y, z = compute_satellite_ecef(eph, t_sow)
            rows.append(
                {
                    "PRN": prn,
                    "t_day_s": t_day_s,
                    "t_sow_s": round(t_sow, 3),
                    "X_m": round(x, 3),
                    "Y_m": round(y, 3),
                    "Z_m": round(z, 3),
                }
            )
        except ValueError:
            continue

    if not rows:
        raise ValueError("未计算出任意卫星在该时刻的位置")
    return rows


def _gps_sow_from_datetime(dt: datetime) -> float:
    """
    将datetime对象转为GPS周内秒。
    """
    return gps_seconds_of_week(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)


def parse_sp3_gps_positions(sp3_file: Path, interval_seconds: int = 300) -> List[Dict[str, Any]]:
    """
    读取SP3精密星历文件，提取GPS卫星坐标。
    interval_seconds: 采样间隔（秒）
    返回：包含PRN、坐标等的点列表
    """
    """读取SP3精密星历中的GPS卫星坐标，单位转换为米。"""
    if interval_seconds <= 0:
        raise ValueError("interval_seconds 必须为正整数")

    points: List[Dict[str, Any]] = []
    current_epoch: Optional[datetime] = None
    first_day: Optional[datetime] = None

    with sp3_file.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("*"):
                parts = line.split()
                if len(parts) < 7:
                    continue
                year = int(parts[1])
                month = int(parts[2])
                day = int(parts[3])
                hour = int(parts[4])
                minute = int(parts[5])
                sec = int(float(parts[6]))
                current_epoch = datetime(year, month, day, hour, minute, sec)
                if first_day is None:
                    first_day = current_epoch.replace(hour=0, minute=0, second=0, microsecond=0)
                continue

            if current_epoch is None or first_day is None:
                continue

            if not line.startswith("PG"):
                continue

            t_day = int((current_epoch - first_day).total_seconds())
            if t_day < 0 or t_day > int(SECONDS_PER_DAY):
                continue
            if t_day % interval_seconds != 0:
                continue

            prn = line[1:4].strip()
            x_km = float(line[4:18])
            y_km = float(line[18:32])
            z_km = float(line[32:46])

            if abs(x_km) > 99999 or abs(y_km) > 99999 or abs(z_km) > 99999:
                continue

            points.append(
                {
                    "PRN": prn,
                    "epoch": current_epoch,
                    "t_day_s": t_day,
                    "t_sow_s": _gps_sow_from_datetime(current_epoch),
                    "X_precise_m": x_km * 1000.0,
                    "Y_precise_m": y_km * 1000.0,
                    "Z_precise_m": z_km * 1000.0,
                }
            )

    if not points:
        raise ValueError("SP3 文件中未解析到GPS卫星坐标")
    return points


def compare_broadcast_with_precise(
    ephemeris_list: Sequence[Dict[str, Any]],
    sp3_points: Sequence[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """
    对比广播星历与精密星历在同一历元下的坐标误差。
    返回：误差明细行列表
    """
    rows: List[Dict[str, Any]] = []
    for p in sp3_points:
        prn = str(p["PRN"])
        t_sow = float(p["t_sow_s"])
        try:
            eph = find_best_ephemeris_by_toe(ephemeris_list, prn, t_sow)
        except ValueError:
            continue

        xb, yb, zb = compute_satellite_ecef(eph, t_sow)
        dx = xb - float(p["X_precise_m"])
        dy = yb - float(p["Y_precise_m"])
        dz = zb - float(p["Z_precise_m"])
        err = math.sqrt(dx * dx + dy * dy + dz * dz)

        rows.append(
            {
                "PRN": prn,
                "epoch": p["epoch"].strftime("%Y-%m-%d %H:%M:%S"),
                "t_day_s": int(p["t_day_s"]),
                "t_sow_s": round(t_sow, 3),
                "X_brdc_m": round(xb, 3),
                "Y_brdc_m": round(yb, 3),
                "Z_brdc_m": round(zb, 3),
                "X_precise_m": round(float(p["X_precise_m"]), 3),
                "Y_precise_m": round(float(p["Y_precise_m"]), 3),
                "Z_precise_m": round(float(p["Z_precise_m"]), 3),
                "dX_m": round(dx, 3),
                "dY_m": round(dy, 3),
                "dZ_m": round(dz, 3),
                "err_3d_m": round(err, 3),
            }
        )

    if not rows:
        raise ValueError("未生成误差记录，请检查PRN覆盖或输入文件")
    return rows


def summarize_errors(rows: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    统计每颗卫星的误差均值、RMSE、最大值。
    返回：每颗卫星的误差统计字典列表
    """
    grouped: Dict[str, List[float]] = {}
    for r in rows:
        grouped.setdefault(str(r["PRN"]), []).append(float(r["err_3d_m"]))

    summary: List[Dict[str, Any]] = []
    for prn in sorted(grouped.keys()):
        vals = grouped[prn]
        n = len(vals)
        mean_v = sum(vals) / n
        rmse_v = math.sqrt(sum(v * v for v in vals) / n)
        max_v = max(vals)
        summary.append(
            {
                "PRN": prn,
                "count": n,
                "mean_err_m": round(mean_v, 3),
                "rmse_m": round(rmse_v, 3),
                "max_err_m": round(max_v, 3),
            }
        )
    return summary


def summarize_constellation_overall(rows: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    """
    统计全星座（所有卫星合并）3D误差总体指标：RMSE与95%分位。
    返回：总体统计字典
    """
    if not rows:
        raise ValueError("没有可统计的误差数据")

    vals = [float(r["err_3d_m"]) for r in rows]
    n = len(vals)
    rmse_v = math.sqrt(sum(v * v for v in vals) / n)
    vals_sorted = sorted(vals)
    # 线性插值计算95%分位，避免样本数较小时出现阶跃偏差
    pos = 0.95 * (n - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        p95_v = vals_sorted[lo]
    else:
        w = pos - lo
        p95_v = vals_sorted[lo] * (1.0 - w) + vals_sorted[hi] * w

    return {
        "PRN": "ALL_GPS",
        "count": n,
        "rmse_m": round(rmse_v, 3),
        "p95_err_m": round(p95_v, 3),
    }


def save_rows_csv(rows: Sequence[Dict[str, Any]], output_csv: Path, headers: Sequence[str]) -> None:
    """
    将任意数据行写入CSV文件。
    """
    if not rows:
        raise ValueError("没有可写入的数据")
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(headers))
        writer.writeheader()
        writer.writerows(rows)


def plot_error_timeseries(rows: Sequence[Dict[str, Any]], output_html: Path) -> None:
    """
    绘制所有卫星误差时序图。
    """
    if not rows:
        raise ValueError("没有可绘制的误差数据")

    grouped: Dict[str, List[Dict[str, Any]]] = {}
    for r in rows:
        grouped.setdefault(str(r["PRN"]), []).append(r)

    overall = summarize_constellation_overall(rows)

    fig = go.Figure()
    for prn in sorted(grouped.keys()):
        records = sorted(grouped[prn], key=lambda x: int(x["t_day_s"]))
        x_hour = [int(r["t_day_s"]) / 3600.0 for r in records]
        y_err = [float(r["err_3d_m"]) for r in records]
        fig.add_trace(
            go.Scatter(
                x=x_hour,
                y=y_err,
                mode="lines",
                name=prn,
                line={"width": 1.5},
            )
        )

    fig.update_layout(
        title=(
            "Broadcast vs Precise Orbit Error Time Series (5 min)<br>"
            f"ALL_GPS: RMSE={overall['rmse_m']} m, P95={overall['p95_err_m']} m"
        ),
        xaxis_title="Time of Day (hour)",
        yaxis_title="3D Position Error (m)",
        legend_title="PRN",
        hovermode="x unified",
    )
    output_html.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(output_html), include_plotlyjs="cdn")


def plot_error_bar_by_satellite(summary_rows: Sequence[Dict[str, Any]], output_html: Path) -> None:
    """
    绘制所有GPS卫星广播星历误差柱状图（按RMSE）。
    """
    if not summary_rows:
        raise ValueError("没有可绘制的误差统计数据")

    records = sorted(summary_rows, key=lambda x: str(x["PRN"]))
    prns = [str(r["PRN"]) for r in records]
    rmse_vals = [float(r["rmse_m"]) for r in records]
    mean_vals = [float(r["mean_err_m"]) for r in records]
    max_vals = [float(r["max_err_m"]) for r in records]

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=prns,
            y=rmse_vals,
            name="RMSE",
            marker={"color": "#1f77b4"},
            customdata=list(zip(mean_vals, max_vals)),
            hovertemplate=(
                "PRN=%{x}<br>RMSE=%{y:.3f} m"
                "<br>Mean=%{customdata[0]:.3f} m"
                "<br>Max=%{customdata[1]:.3f} m<extra></extra>"
            ),
        )
    )

    fig.update_layout(
        title="Broadcast vs Precise Orbit Error by Satellite (RMSE)",
        xaxis_title="PRN",
        yaxis_title="3D Position Error (m)",
        bargap=0.2,
    )
    output_html.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(output_html), include_plotlyjs="cdn")


def plot_error_xyz_timeseries_subplots(rows: Sequence[Dict[str, Any]], output_html: Path) -> None:
    """
    绘制每颗卫星三维坐标分量误差(dX/dY/dZ)时序子图。
    """
    if not rows:
        raise ValueError("没有可绘制的误差数据")

    grouped: Dict[str, List[Dict[str, Any]]] = {}
    for r in rows:
        grouped.setdefault(str(r["PRN"]), []).append(r)

    overall = summarize_constellation_overall(rows)
    prns = sorted(grouped.keys())
    cols = 3
    n = len(prns)
    subplot_rows = (n + cols - 1) // cols

    fig = make_subplots(
        rows=subplot_rows,
        cols=cols,
        subplot_titles=prns,
        horizontal_spacing=0.05,
        vertical_spacing=0.08,
    )

    for idx, prn in enumerate(prns):
        r = idx // cols + 1
        c = idx % cols + 1
        records = sorted(grouped[prn], key=lambda x: int(x["t_day_s"]))
        x_hour = [int(rec["t_day_s"]) / 3600.0 for rec in records]
        dx_vals = [float(rec["dX_m"]) for rec in records]
        dy_vals = [float(rec["dY_m"]) for rec in records]
        dz_vals = [float(rec["dZ_m"]) for rec in records]

        fig.add_trace(
            go.Scatter(x=x_hour, y=dx_vals, mode="lines", name="dX", line={"width": 1.0, "color": "#d62728"}, showlegend=(idx == 0)),
            row=r,
            col=c,
        )
        fig.add_trace(
            go.Scatter(x=x_hour, y=dy_vals, mode="lines", name="dY", line={"width": 1.0, "color": "#2ca02c"}, showlegend=(idx == 0)),
            row=r,
            col=c,
        )
        fig.add_trace(
            go.Scatter(x=x_hour, y=dz_vals, mode="lines", name="dZ", line={"width": 1.0, "color": "#1f77b4"}, showlegend=(idx == 0)),
            row=r,
            col=c,
        )

        fig.update_xaxes(title_text="Hour", row=r, col=c)
        fig.update_yaxes(title_text="Coord Error (m)", row=r, col=c)

    fig.update_layout(
        title=(
            "Broadcast vs Precise XYZ Coordinate Error by Satellite (5 min)<br>"
            f"ALL_GPS: RMSE={overall['rmse_m']} m, P95={overall['p95_err_m']} m"
        ),
        height=max(340 * subplot_rows, 720),
        width=1450,
        hovermode="x unified",
    )
    output_html.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(output_html), include_plotlyjs="cdn")


def main() -> None:
    """
    命令行主入口，解析参数并调度各功能。
    """
    parser = argparse.ArgumentParser(
        description="广播星历坐标解算与精密星历误差评估"
    )
    parser.add_argument("--mode", choices=["trajectory", "compare", "both"], default="both", help="运行模式")
    parser.add_argument("--nav-file", default="brdc0380.26n", help="RINEX广播星历文件")
    parser.add_argument("--body-file", default="nav_body.txt", help="仅含星历体部分的中间文件")
    parser.add_argument("--sp3-file", default="WUM0MGXFIN_20260380000_01D_05M_ORB.SP3", help="SP3精密星历文件")
    parser.add_argument("--prn", default="G01", help="目标卫星PRN，如 G01")
    parser.add_argument(
        "--all-gps",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="是否计算并绘制所有GPS卫星轨迹（默认开启，可用 --no-all-gps 切换为单星）",
    )
    parser.add_argument("--show-time-labels", action="store_true", help="在轨迹点上显示时间标注")
    parser.add_argument("--time-label-step", type=int, default=30, help="时间标注点间隔（按采样点计）")
    parser.add_argument("--step", type=int, default=60, help="时间步长（秒）")
    parser.add_argument("--t-start", type=int, default=0, help="起始时刻（秒）")
    parser.add_argument("--t-end", type=int, default=86400, help="结束时刻（秒）")
    parser.add_argument("--csv", default="trajectory.csv", help="输出CSV文件")
    parser.add_argument("--html", default="trajectory_3d.html", help="输出3D轨迹图HTML文件")
    parser.add_argument("--all-csv", default="all_gps_trajectories.csv", help="全部GPS轨迹CSV")
    parser.add_argument("--all-html", default="all_gps_trajectories_3d.html", help="全部GPS轨迹图HTML")
    parser.add_argument("--query-t", type=int, default=None, help="输出当天任意时刻t(0~86400s)所有GPS卫星坐标")
    parser.add_argument("--query-csv", default="all_gps_positions_at_t.csv", help="任意时刻坐标输出CSV")
    parser.add_argument("--compare-interval", type=int, default=300, help="误差对比时间间隔（秒）")
    parser.add_argument("--error-csv", default="orbit_error_5min.csv", help="误差明细CSV")
    parser.add_argument("--error-summary-csv", default="orbit_error_summary.csv", help="误差统计CSV")
    parser.add_argument("--error-html", default="orbit_error_timeseries.html", help="误差时序图HTML")
    parser.add_argument("--error-bar-html", default="orbit_error_bar.html", help="误差柱状图HTML")
    parser.add_argument(
        "--error-xyz-subplots-html",
        default="orbit_error_xyz_timeseries_subplots.html",
        help="各卫星三维坐标分量误差子图HTML",
    )
    args = parser.parse_args()

    nav_file = Path(args.nav_file)
    body_file = Path(args.body_file)
    sp3_file = Path(args.sp3_file)
    csv_file = Path(args.csv)
    html_file = Path(args.html)
    all_csv_file = Path(args.all_csv)
    all_html_file = Path(args.all_html)
    query_csv_file = Path(args.query_csv)
    error_csv_file = Path(args.error_csv)
    error_summary_csv_file = Path(args.error_summary_csv)
    error_html_file = Path(args.error_html)
    error_bar_html_file = Path(args.error_bar_html)
    error_xyz_subplots_html_file = Path(args.error_xyz_subplots_html)

    if not nav_file.exists():
        raise FileNotFoundError(f"未找到广播星历文件: {nav_file}")

    save_nav_body_file(str(nav_file), str(body_file))
    ephemeris_list = read_nav_body_file(str(body_file))

    if args.query_t is not None:
        query_rows = compute_all_gps_positions_at_t(ephemeris_list, args.query_t)
        save_trajectory_csv(query_rows, query_csv_file)
        print(f"已输出 t={args.query_t}s 时刻全部GPS卫星坐标（mm精度）")
        print(f"时刻坐标CSV输出: {query_csv_file.resolve()}")

    if args.mode in ("trajectory", "both"):
        if args.all_gps:
            all_trajectories = generate_all_gps_trajectories(
                ephemeris_list=ephemeris_list,
                step_seconds=args.step,
                t_start=args.t_start,
                t_end=args.t_end,
            )
            all_rows = flatten_trajectories(all_trajectories)
            save_trajectory_csv(all_rows, all_csv_file)
            if args.show_time_labels:
                plot_all_gps_trajectories_3d_with_time_labels(
                    all_trajectories,
                    all_html_file,
                    args.time_label_step,
                )
            else:
                plot_all_gps_trajectories_3d(all_trajectories, all_html_file)

            print("已完成全部GPS卫星轨迹解算")
            print(f"卫星数量: {len(all_trajectories)}")
            print(f"轨迹总点数: {len(all_rows)}")
            print(f"全部轨迹CSV输出: {all_csv_file.resolve()}")
            print(f"全部轨迹图输出: {all_html_file.resolve()}")
        else:
            points = generate_daily_trajectory(
                ephemeris_list=ephemeris_list,
                prn=args.prn,
                step_seconds=args.step,
                t_start=args.t_start,
                t_end=args.t_end,
            )
            save_trajectory_csv(points, csv_file)
            if args.show_time_labels:
                plot_trajectory_3d_with_time_labels(points, args.prn, html_file, args.time_label_step)
            else:
                plot_trajectory_3d(points, args.prn, html_file)

            print(f"已完成 PRN G{normalize_prn(args.prn)} 轨迹解算")
            print(f"轨迹点数量: {len(points)}")
            print(f"轨迹CSV输出: {csv_file.resolve()}")
            print(f"轨迹图输出: {html_file.resolve()}")
            print(f"首点坐标(m): ({points[0]['X_m']}, {points[0]['Y_m']}, {points[0]['Z_m']})")
            print(f"末点坐标(m): ({points[-1]['X_m']}, {points[-1]['Y_m']}, {points[-1]['Z_m']})")

    if args.mode in ("compare", "both"):
        if not sp3_file.exists():
            raise FileNotFoundError(f"未找到SP3精密星历文件: {sp3_file}")

        sp3_points = parse_sp3_gps_positions(sp3_file, interval_seconds=args.compare_interval)
        error_rows = compare_broadcast_with_precise(ephemeris_list, sp3_points)
        summary_rows = summarize_errors(error_rows)
        overall_row = summarize_constellation_overall(error_rows)

        detail_headers = [
            "PRN",
            "epoch",
            "t_day_s",
            "t_sow_s",
            "X_brdc_m",
            "Y_brdc_m",
            "Z_brdc_m",
            "X_precise_m",
            "Y_precise_m",
            "Z_precise_m",
            "dX_m",
            "dY_m",
            "dZ_m",
            "err_3d_m",
        ]
        summary_headers = ["PRN", "count", "mean_err_m", "rmse_m", "max_err_m"]
        save_rows_csv(error_rows, error_csv_file, detail_headers)
        save_rows_csv(summary_rows, error_summary_csv_file, summary_headers)
        plot_error_timeseries(error_rows, error_html_file)
        plot_error_bar_by_satellite(summary_rows, error_bar_html_file)
        plot_error_xyz_timeseries_subplots(error_rows, error_xyz_subplots_html_file)

        print(f"已完成广播星历与精密星历误差评估，间隔: {args.compare_interval}s")
        print(f"误差明细CSV: {error_csv_file.resolve()}")
        print(f"误差统计CSV: {error_summary_csv_file.resolve()}")
        print(f"误差时序图: {error_html_file.resolve()}")
        print(f"误差柱状图: {error_bar_html_file.resolve()}")
        print(f"三轴误差子图: {error_xyz_subplots_html_file.resolve()}")
        print(f"参与对比卫星数量: {len(summary_rows)}")
        print("全星座总体统计:")
        print(f"  ALL_GPS RMSE(m): {overall_row['rmse_m']}")
        print(f"  ALL_GPS P95(m): {overall_row['p95_err_m']}")
        print("部分卫星RMSE(m):")
        for row in summary_rows[:10]:
            print(f"  {row['PRN']}: {row['rmse_m']}")


if __name__ == "__main__":
    main()
