"""Self-contained interactive HTML report for completed Spacer outputs."""

from __future__ import annotations

import html
import json
from pathlib import Path


def write_html_report(
    path: Path,
    *,
    analysis_rows: list[dict[str, str]],
    sensitivity_rows: list[dict[str, str]],
    interference_rows: list[dict[str, str]],
    agreement_rows: list[dict[str, str]],
    discordant_rows: list[dict[str, str]],
    validation_rows: list[dict[str, str]],
    plot_rows: list[dict[str, object]],
) -> None:
    group = next(row for row in analysis_rows if row["summary_level"] == "single_group")
    matched = [row for row in agreement_rows if row["subset"] == "all_matched"]
    correlations = [float(row["pearson_r"]) for row in matched if row["pearson_r"]]
    median_difference = [float(row["median_signed_difference_pp"]) for row in matched if row["median_signed_difference_pp"]]
    pd_spectrum_mapped = [
        row for row in agreement_rows if row["subset"] == "pd_spectrum_mapped"
    ]
    pd_identified = [
        row for row in agreement_rows if row["subset"] == "identified"
    ]
    mapped_scorable = sum(int(row["scorable_ms2"]) for row in pd_spectrum_mapped)
    mapped_likely = sum(int(row["likely_chimeric"]) for row in pd_spectrum_mapped)
    identified_scorable = sum(int(row["scorable_ms2"]) for row in pd_identified)
    identified_likely = sum(int(row["likely_chimeric"]) for row in pd_identified)
    group_distribution = [
        row for row in interference_rows if row["summary_level"] == "single_group"
    ]
    median_interference = float(group_distribution[0]["median_interference_fraction"])
    interpretation = [
        "Chimericity is quantified continuously as MS1 isolation-window interference: 0% means all window signal is assigned to the transmitted target envelope; higher values mean a larger non-target share of the window signal.",
        f"The median continuous interference estimate was {median_interference * 100:.1f}% across scorable MS2; the Interference degree tab shows the full binned distribution.",
        f"Spacer classified {float(group['likely_chimeric_fraction']) * 100:.1f}% of scorable MS2 as likely chimeric at the configured threshold.",
        f"Technical-replicate likely-chimeric fractions ranged from {float(group['replicate_likely_fraction_min']) * 100:.1f}% to {float(group['replicate_likely_fraction_max']) * 100:.1f}%.",
        f"{float(group['indeterminate_fraction']) * 100:.1f}% of all MS2 were indeterminate and are not silently treated as low interference.",
        (
            f"Among scorable MS2 with an exact PD MS/MS-spectrum record, {mapped_likely}/{mapped_scorable} ({100 * mapped_likely / mapped_scorable:.1f}%) were likely chimeric."
            if mapped_scorable
            else "No scorable MS2 had an exact PD MS/MS-spectrum record."
        ),
        (
            f"Among the narrower q-value-passing PD-identified subset, {identified_likely}/{identified_scorable} ({100 * identified_likely / identified_scorable:.1f}%) were likely chimeric."
            if identified_scorable
            else "No scorable q-value-passing PD-identified MS2 were available."
        ),
        (
            f"Exact-scan PD agreement had mean Pearson r={sum(correlations) / len(correlations):.3f} and mean signed Spacer minus PD difference={sum(median_difference) / len(median_difference):.2f} percentage points."
            if correlations and median_difference
            else "No finite PD interference agreement values were available."
        ),
        "PD values are descriptive reference context only; they did not modify Spacer scores, classifications, or thresholds.",
    ]
    payload = json.dumps(
        {
            "analysis": analysis_rows,
            "sensitivity": sensitivity_rows,
            "interference": interference_rows,
            "agreement": agreement_rows,
            "discordant": discordant_rows,
            "validation": validation_rows,
            "plots": plot_rows,
        },
        ensure_ascii=False,
    ).replace("</", "<\\/")
    bullets = "".join(f"<li>{html.escape(item)}</li>" for item in interpretation)
    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Spacer Single-group report</title>
<style>
body {{ font-family:system-ui,sans-serif; margin:2rem auto; max-width:1200px; color:#17212b; }}
.note {{ color:#53616e; }} .cards {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(180px,1fr)); gap:1rem; margin:1.25rem 0; }}
.card {{ background:#f4f7f9; border-radius:8px; padding:1rem; }} .card b {{ font-size:1.35rem; display:block; }}
button {{ padding:.55rem .8rem; border:1px solid #b8c2cc; background:#fff; cursor:pointer; }} button.active {{ background:#245b8c; color:#fff; }}
.panel {{ display:none; margin-top:1rem; }} .panel.active {{ display:block; }} input {{ padding:.5rem; width:min(400px,100%); margin:.5rem 0; }}
table {{ border-collapse:collapse; width:100%; font-size:.85rem; }} th,td {{ border:1px solid #d8dee4; padding:.45rem; text-align:left; vertical-align:top; }}
th {{ background:#edf2f6; }} .scroll {{ overflow:auto; max-height:560px; }} .warning {{ background:#fff8e8; padding:.8rem; border-left:4px solid #dd9a24; }}
.plot-grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(420px,1fr)); gap:1rem; }} figure {{ margin:0; }} svg {{ width:100%; height:auto; border:1px solid #d8dee4; background:#fff; }} figcaption {{ font-weight:600; margin:.5rem 0; }}
</style></head><body>
<h1>Spacer Single-group report</h1><p class="note">Interactive, self-contained report. All score calculations use mzML evidence; PD is reference context only.</p>
<div class="cards"><div class="card"><span>Runs</span><b>{html.escape(group["run_count"])}</b></div><div class="card"><span>Scorable MS2</span><b>{html.escape(group["scorable_ms2"])}</b></div><div class="card"><span>Median interference</span><b>{median_interference * 100:.1f}%</b></div><div class="card"><span>Likely chimeric</span><b>{float(group["likely_chimeric_fraction"]) * 100:.1f}%</b></div><div class="card"><span>Indeterminate</span><b>{float(group["indeterminate_fraction"]) * 100:.1f}%</b></div></div>
<h2>Interpretation</h2><ul>{bullets}</ul><p class="warning">Agreement is not accuracy: PD exports are not ground truth and never recalibrate Spacer.</p>
<nav><button class="tab active" data-panel="run">Run summary</button><button class="tab" data-panel="interference">Interference degree</button><button class="tab" data-panel="sensitivity">Sensitivity</button><button class="tab" data-panel="agreement">PD agreement</button><button class="tab" data-panel="discordant">Discordant scans</button><button class="tab" data-panel="plots">PD-matched plots</button><button class="tab" data-panel="validation">Validation</button></nav>
<input id="filter" placeholder="Filter active table"><div id="run" class="panel active"></div><div id="interference" class="panel"></div><div id="sensitivity" class="panel"></div><div id="agreement" class="panel"></div><div id="discordant" class="panel"></div><div id="plots" class="panel"></div><div id="validation" class="panel"></div>
<script>
const data = {payload};
const maps = {{run:data.analysis,interference:data.interference,sensitivity:data.sensitivity,agreement:data.agreement,discordant:data.discordant,validation:data.validation}};
let active = "run";
function esc(value) {{ const node=document.createElement("span"); node.textContent=value ?? ""; return node.innerHTML; }}
function render(name) {{
  if (name === "plots") {{ renderPlots(); return; }}
  const rows=maps[name], filter=document.getElementById("filter").value.toLowerCase();
  const visible=rows.filter(function(row) {{ return Object.values(row).join(" ").toLowerCase().includes(filter); }});
  const columns=rows.length?Object.keys(rows[0]):[];
  let body="<div class=\\"scroll\\"><table><thead><tr>";
  columns.forEach(function(column) {{ body+="<th>"+esc(column)+"</th>"; }}); body+="</tr></thead><tbody>";
  visible.forEach(function(row) {{ body+="<tr>"; columns.forEach(function(column) {{ body+="<td>"+esc(row[column])+"</td>"; }}); body+="</tr>"; }});
  body+="</tbody></table></div><p>"+visible.length+" of "+rows.length+" rows shown.</p>";
  document.getElementById(name).innerHTML=body;
}}
function sticks(points, width, height, color, transform) {{
  if (!points.length) return "";
  const max=Math.max(...points.map(function(point) {{ return Number(point.intensity); }}),1);
  return points.map(function(point) {{
    const xy=transform(point, max);
    return "<line x1='"+xy.x+"' x2='"+xy.x+"' y1='"+height+"' y2='"+xy.y+"' stroke='"+(xy.color || color)+"' stroke-width='1'/>";
  }}).join("");
}}
function renderPlots() {{
  const target=document.getElementById("plots");
  if (!data.plots.length) {{ target.innerHTML="<p>No exact PD-matched discordant scans were available for plotting.</p>"; return; }}
  const index=Math.min(Number(target.dataset.index || 0),data.plots.length-1), plot=data.plots[index];
  const ms1=plot.ms1, ms2=plot.ms2, width=530, height=250, pad=36;
  const ms1Min=Math.min(...ms1.map(function(p){{return Number(p.mz);}})), ms1Max=Math.max(...ms1.map(function(p){{return Number(p.mz);}}));
  const ms2Min=Math.min(...ms2.map(function(p){{return Number(p.mz);}})), ms2Max=Math.max(...ms2.map(function(p){{return Number(p.mz);}}));
  function x(mz,min,max) {{ return pad+(Number(mz)-min)*(width-2*pad)/(max-min || 1); }}
  const ms1Sticks=sticks(ms1,width-pad,height-pad,"#303030",function(p,max) {{
    const targetPeak=p.is_target_envelope==="true", transmitted=p.is_transmitted_target_signal==="true";
    return {{x:x(p.mz,ms1Min,ms1Max),y:height-pad-(Number(p.intensity)/max)*(height-2*pad),color:transmitted ? "#00a86b" : targetPeak ? "#2a9d8f" : p.in_isolation_window==="true" ? "#303030" : "#b5b5b5"}};
  }});
  const ms2Sticks=sticks(ms2,width-pad,height-pad,"#305f9f",function(p,max) {{
    return {{x:x(p.mz,ms2Min,ms2Max),y:height-pad-(Number(p.intensity)/max)*(height-2*pad)}};
  }});
  const lower=x(plot.isolation_lower_mz,ms1Min,ms1Max), upper=x(plot.isolation_upper_mz,ms1Min,ms1Max);
  target.innerHTML="<p><button id='previous-plot'>Previous</button> <button id='next-plot'>Next</button> "+(index+1)+" of "+data.plots.length+
    " &mdash; <b>"+esc(plot.run_basename)+" scan "+esc(plot.scan_id)+"</b>; Spacer interference "+(100*Number(plot.interference_fraction)).toFixed(1)+"%; PD interference "+Number(plot.pd_isolation_interference_percent).toFixed(1)+"%; absolute difference "+Number(plot.absolute_difference_percent_points).toFixed(1)+" percentage points; "+esc(plot.classification)+
    ".</p><p class='note'>These are the exact PD-matched scans with the largest Spacer-versus-PD differences. The amber band is the recorded isolation window; green is transmitted target-envelope signal; teal is target-envelope signal outside that window.</p>"+
    "<div class='plot-grid'><figure><figcaption>MS1 isolation window</figcaption><svg viewBox='0 0 "+width+" "+height+"' role='img'><rect x='"+lower+"' y='"+pad+"' width='"+(upper-lower)+"' height='"+(height-2*pad)+"' fill='#f4a261' fill-opacity='.25'/>"+ms1Sticks+"<line x1='"+pad+"' x2='"+(width-pad)+"' y1='"+(height-pad)+"' y2='"+(height-pad)+"' stroke='#17212b'/><text x='"+pad+"' y='"+(height-8)+"'>m/z "+ms1Min.toFixed(2)+"–"+ms1Max.toFixed(2)+"</text></svg></figure>"+
    "<figure><figcaption>MS2 spectrum (relative intensity)</figcaption><svg viewBox='0 0 "+width+" "+height+"' role='img'>"+ms2Sticks+"<line x1='"+pad+"' x2='"+(width-pad)+"' y1='"+(height-pad)+"' y2='"+(height-pad)+"' stroke='#17212b'/><text x='"+pad+"' y='"+(height-8)+"'>m/z "+ms2Min.toFixed(2)+"–"+ms2Max.toFixed(2)+"</text></svg></figure></div>";
  document.getElementById("previous-plot").onclick=function() {{ target.dataset.index=(index-1+data.plots.length)%data.plots.length; renderPlots(); }};
  document.getElementById("next-plot").onclick=function() {{ target.dataset.index=(index+1)%data.plots.length; renderPlots(); }};
}}
document.querySelectorAll(".tab").forEach(function(button) {{ button.onclick=function() {{ document.querySelectorAll(".tab,.panel").forEach(function(item) {{ item.classList.remove("active"); }}); button.classList.add("active"); active=button.dataset.panel; document.getElementById(active).classList.add("active"); render(active); }}; }});
document.getElementById("filter").oninput=function() {{ render(active); }}; render(active);
</script></body></html>"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(page, encoding="utf-8")
