from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Any, Dict, List, Sequence

import plotly.graph_objects as go
from flask import Flask, render_template_string, request

from GPS_Coordination import (
    compare_broadcast_with_precise,
    flatten_trajectories,
    generate_all_gps_trajectories,
    parse_sp3_gps_positions,
    summarize_constellation_overall,
    summarize_errors,
)
from GPS_IO import read_nav_body_file, save_nav_body_file


app = Flask(__name__)


PAGE = """
<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>GPS 星历交互分析</title>
  <style>
    body { font-family: 'Segoe UI', sans-serif; margin: 0; background: #f4f8ff; color: #14213d; }
    .wrap { max-width: 1280px; margin: 0 auto; padding: 20px; }
    .hero { background: linear-gradient(135deg, #d9f0ff, #fff5d6); border-radius: 14px; padding: 20px; }
    .grid { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; margin-top: 16px; }
    .card { background: #fff; border-radius: 12px; padding: 16px; box-shadow: 0 8px 20px rgba(0,0,0,0.06); }
    h1, h2, h3 { margin: 0 0 12px 0; }
    label { display: block; margin-top: 10px; font-size: 14px; }
    input[type='number'], input[type='file'] { width: 100%; margin-top: 6px; }
    button { margin-top: 14px; padding: 8px 14px; border: none; background: #1565c0; color: #fff; border-radius: 8px; cursor: pointer; }
    button:hover { background: #0d47a1; }
    .metrics { display: flex; gap: 10px; flex-wrap: wrap; margin: 10px 0 16px 0; }
    .metric { background: #eef6ff; border-radius: 10px; padding: 8px 12px; }
    .full { margin-top: 16px; }
    table { border-collapse: collapse; width: 100%; background: #fff; }
    th, td { border: 1px solid #e3e3e3; padding: 6px; text-align: center; }
    .err { color: #b00020; font-weight: 600; }
    @media (max-width: 900px) { .grid { grid-template-columns: 1fr; } }
  </style>
</head>
<body>
  <div class="wrap">
    <div class="hero">
      <h1>GPS 广播星历与精密星历交互分析</h1>
      <div>上传文件后可在线计算轨迹与误差图。</div>
      {% if error %}<div class="err">{{ error }}</div>{% endif %}
    </div>

    <div class="grid">
      <div class="card">
        <h2>全天轨迹（所有 GPS 卫星）</h2>
        <form method="post" action="/run_trajectory" enctype="multipart/form-data">
          <label>广播星历（RINEX NAV）</label>
          <input name="nav_file" type="file" required>
          <label>采样步长（秒）</label>
          <input name="step" type="number" min="60" max="1800" value="300" required>
          <label>起始时刻（秒）</label>
          <input name="t_start" type="number" min="0" max="86400" value="0" required>
          <label>结束时刻（秒）</label>
          <input name="t_end" type="number" min="0" max="86400" value="86400" required>
          <button type="submit">计算并绘制轨迹</button>
        </form>
      </div>

      <div class="card">
        <h2>广播-精密误差分析</h2>
        <form method="post" action="/run_compare" enctype="multipart/form-data">
          <label>广播星历（RINEX NAV）</label>
          <input name="nav_file" type="file" required>
          <label>精密星历（SP3）</label>
          <input name="sp3_file" type="file" required>
          <label>对比间隔（秒）</label>
          <input name="compare_interval" type="number" min="60" max="1800" value="300" required>
          <button type="submit">计算并绘制误差</button>
        </form>
      </div>
    </div>

    {% if metrics %}
      <div class="metrics">
        {% for m in metrics %}
          <div class="metric"><strong>{{ m.label }}</strong>: {{ m.value }}</div>
        {% endfor %}
      </div>
    {% endif %}

    {% if plot1 %}<div class="card full">{{ plot1|safe }}</div>{% endif %}
    {% if plot2 %}<div class="card full">{{ plot2|safe }}</div>{% endif %}
    {% if table_html %}<div class="card full">{{ table_html|safe }}</div>{% endif %}
  </div>
</body>
</html>
"""


def _save_uploaded_file(file_obj: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    file_obj.save(path)


def _build_all_gps_trajectory_figure(trajectories: Dict[str, Sequence[Dict[str, Any]]]) -> go.Figure:
    fig = go.Figure()
    for prn in sorted(trajectories.keys()):
        points = trajectories[prn]
        if not points:
            continue
        fig.add_trace(
            go.Scatter3d(
                x=[p["X_m"] for p in points],
                y=[p["Y_m"] for p in points],
                z=[p["Z_m"] for p in points],
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
    return fig


def _build_error_timeseries_figure(rows: Sequence[Dict[str, Any]]) -> go.Figure:
    grouped: Dict[str, List[Dict[str, Any]]] = {}
    for r in rows:
        grouped.setdefault(str(r["PRN"]), []).append(r)

    overall = summarize_constellation_overall(rows)
    fig = go.Figure()
    for prn in sorted(grouped.keys()):
        records = sorted(grouped[prn], key=lambda x: int(x["t_day_s"]))
        fig.add_trace(
            go.Scatter(
                x=[int(r["t_day_s"]) / 3600.0 for r in records],
                y=[float(r["err_3d_m"]) for r in records],
                mode="lines",
                name=prn,
                line={"width": 1.4},
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
    return fig


def _build_error_bar_figure(summary_rows: Sequence[Dict[str, Any]]) -> go.Figure:
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
    return fig


def _summary_table_html(summary_rows: Sequence[Dict[str, Any]]) -> str:
    headers = ["PRN", "count", "mean_err_m", "rmse_m", "max_err_m"]
    lines = ["<h3>卫星误差统计</h3>", "<table><thead><tr>"]
    for h in headers:
        lines.append(f"<th>{h}</th>")
    lines.append("</tr></thead><tbody>")
    for row in summary_rows:
        lines.append("<tr>")
        for h in headers:
            lines.append(f"<td>{row[h]}</td>")
        lines.append("</tr>")
    lines.append("</tbody></table>")
    return "".join(lines)


def _render(
    error: str = "",
    metrics: Sequence[Dict[str, Any]] | None = None,
    plot1: str = "",
    plot2: str = "",
    table_html: str = "",
):
    return render_template_string(
        PAGE,
        error=error,
        metrics=metrics or [],
        plot1=plot1,
        plot2=plot2,
        table_html=table_html,
    )


@app.route("/", methods=["GET"])
def index():
    return _render()


@app.route("/run_trajectory", methods=["POST"])
def run_trajectory():
    nav_file = request.files.get("nav_file")
    if nav_file is None or nav_file.filename == "":
        return _render(error="请上传广播星历文件。")

    try:
        step = int(request.form.get("step", 300))
        t_start = int(request.form.get("t_start", 0))
        t_end = int(request.form.get("t_end", 86400))
        if t_start > t_end:
            raise ValueError("起始时刻不能大于结束时刻")

        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            nav_path = td_path / nav_file.filename
            body_path = td_path / "nav_body.txt"
            _save_uploaded_file(nav_file, nav_path)

            save_nav_body_file(str(nav_path), str(body_path))
            ephemeris_list = read_nav_body_file(str(body_path))
            trajectories = generate_all_gps_trajectories(
                ephemeris_list=ephemeris_list,
                step_seconds=step,
                t_start=t_start,
                t_end=t_end,
            )
            all_rows = flatten_trajectories(trajectories)

        fig = _build_all_gps_trajectory_figure(trajectories)
        metrics = [
            {"label": "卫星数量", "value": len(trajectories)},
            {"label": "轨迹总点数", "value": len(all_rows)},
        ]
        return _render(
            metrics=metrics,
            plot1=fig.to_html(full_html=False, include_plotlyjs="cdn"),
        )
    except Exception as e:
        return _render(error=f"轨迹计算失败: {e}")


@app.route("/run_compare", methods=["POST"])
def run_compare():
    nav_file = request.files.get("nav_file")
    sp3_file = request.files.get("sp3_file")
    if nav_file is None or nav_file.filename == "" or sp3_file is None or sp3_file.filename == "":
        return _render(error="请同时上传广播星历和精密星历文件。")

    try:
        compare_interval = int(request.form.get("compare_interval", 300))
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            nav_path = td_path / nav_file.filename
            sp3_path = td_path / sp3_file.filename
            body_path = td_path / "nav_body.txt"
            _save_uploaded_file(nav_file, nav_path)
            _save_uploaded_file(sp3_file, sp3_path)

            save_nav_body_file(str(nav_path), str(body_path))
            ephemeris_list = read_nav_body_file(str(body_path))
            sp3_points = parse_sp3_gps_positions(sp3_path, interval_seconds=compare_interval)
            error_rows = compare_broadcast_with_precise(ephemeris_list, sp3_points)
            summary_rows = summarize_errors(error_rows)
            overall = summarize_constellation_overall(error_rows)

        fig_ts = _build_error_timeseries_figure(error_rows)
        fig_bar = _build_error_bar_figure(summary_rows)
        metrics = [
            {"label": "参与对比卫星数", "value": len(summary_rows)},
            {"label": "ALL_GPS RMSE(m)", "value": overall["rmse_m"]},
            {"label": "ALL_GPS P95(m)", "value": overall["p95_err_m"]},
        ]
        return _render(
            metrics=metrics,
            plot1=fig_ts.to_html(full_html=False, include_plotlyjs="cdn"),
            plot2=fig_bar.to_html(full_html=False, include_plotlyjs=False),
            table_html=_summary_table_html(summary_rows),
        )
    except Exception as e:
        return _render(error=f"误差评估失败: {e}")


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=8501, debug=False)
