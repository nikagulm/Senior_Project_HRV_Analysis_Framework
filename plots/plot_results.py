import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons, Button, TextBox

RESULTS = "results"
TOLERANCE = 0.15  #s same as the evaluator
WINDOW = 10.0 #s of ECG shown at once

#records shown in the record explorer and saved as figures (any other record can still be typed in)
#100 easy, 108 many false beats for Pan-Tompkins, 119 Hilbert misses beats, 208 search-back matters, 118e06 noise
SHOWN = ["100", "108", "119", "208", "118e06"]
EXAMPLE = "208" #the record saved as the ECG figure with --save

COLOURS = {
    "Reference": "black",
    "Pan-Tompkins": "#1f4e9c", "PT fixed threshold": "#4f7fd0", "PT no search-back": "#6c95dc", "PT moving average": "#8fb0e8",
    "Hamilton": "#c0561b", "Ham fixed threshold": "#e07b39", "Ham no search-back": "#e8935a", "Ham moving average": "#efae80",
    "Hilbert": "#2e7d32",
}
METRICS = [("SDNN", "SDNN"), ("RMSSD", "RMSSD"), ("pNN50", "pNN50"), ("LF", "LF"), ("HF", "HF"), ("LF_HF", "LF/HF")]

results = pd.read_csv(f"{RESULTS}/results.csv", dtype={"record": str})
info = pd.read_csv(f"{RESULTS}/records.csv", dtype={"record": str}).set_index("record")
peaks = pd.read_csv(f"{RESULTS}/peaks.csv", dtype={"record": str})
nn = pd.read_csv(f"{RESULTS}/nn.csv", dtype={"record": str})

records = list(dict.fromkeys(results["record"]))
#noise records are named like 118e24, the clean MIT-BIH records are just numbers
noisy = [r for r in records if "e" in r]
clean = [r for r in records if "e" not in r]
detectors = [d for d in dict.fromkeys(results["detector"]) if d != "Reference"]
save = "--save" in sys.argv
if save:
    os.makedirs(f"{RESULTS}/plots", exist_ok=True)

#% difference of every HRV value from the reference of the same record
ref = results[results.detector == "Reference"].set_index("record")
det = results[results.detector != "Reference"].copy()
for m, _ in METRICS:
    det[f"d_{m}"] = (det[m] - det.record.map(ref[m])) / det.record.map(ref[m]) * 100


def match(found, reference):
    #returns (hits, false positives, missed reference beats), same rule as the evaluator
    used = np.zeros(len(reference), bool)
    hits, extra = [], []
    for t in found:
        i = np.searchsorted(reference, t)
        best = None
        for j in (i - 1, i):
            if 0 <= j < len(reference) and not used[j] and abs(reference[j] - t) <= TOLERANCE:
                if best is None or abs(reference[j] - t) < abs(reference[best] - t):
                    best = j
        if best is None:
            extra.append(t)
        else:
            used[best] = True
            hits.append(t)
    return np.array(hits), np.array(extra), reference[~used]


#window 1. record explorer
fig = plt.figure("Record explorer", figsize=(15, 9.5))
axECG = fig.add_axes([0.24, 0.52, 0.68, 0.42])
axHRV = fig.add_axes([0.24, 0.07, 0.68, 0.33])
axRecord = fig.add_axes([0.06, 0.86, 0.09, 0.05])
axRecPrev = fig.add_axes([0.02, 0.79, 0.08, 0.05])
axRecNext = fig.add_axes([0.11, 0.79, 0.08, 0.05])
axShow = fig.add_axes([0.02, 0.17, 0.17, 0.48])
axPrev = fig.add_axes([0.02, 0.07, 0.08, 0.05])
axNext = fig.add_axes([0.11, 0.07, 0.08, 0.05])
axShow.set_title("Detectors", fontsize=10)

state = {"record": SHOWN[0] if SHOWN[0] in records else records[0], "start": 0.0, "errors": np.array([])}
shown = {d: d in ("Pan-Tompkins", "Hamilton", "Hilbert") for d in detectors}
cache = {}


def load(record):
    if record not in cache:
        fs = int(info.loc[record, "fs"])
        signal = np.fromfile(f"{RESULTS}/ecg_{record}.bin", dtype=np.float32)
        refPeaks = np.sort(peaks[(peaks.record == record) & (peaks.detector == "Reference")].time.values)
        found = {}
        for d in detectors:
            p = np.sort(peaks[(peaks.record == record) & (peaks.detector == d)].time.values)
            found[d] = match(p, refPeaks)
        cache[record] = (fs, signal, refPeaks, found)
    return cache[record]


def error_times():
    _, _, _, found = load(state["record"])
    times = []
    for d in detectors:
        if shown[d]:
            _, extra, missed = found[d]
            times.extend(extra)
            times.extend(missed)
    return np.unique(np.round(np.array(times), 1))


def draw_ecg():
    fs, signal, refPeaks, found = load(state["record"])
    start = state["start"]
    a, b = int(start * fs), int((start + WINDOW) * fs)
    t = np.arange(a, min(b, len(signal))) / fs
    y = signal[a:b]

    axECG.clear()
    axECG.plot(t, y, color="#555555", linewidth=0.9)
    for r in refPeaks[(refPeaks >= start) & (refPeaks < start + WINDOW)]:
        axECG.axvline(r, color="#bbbbbb", linewidth=0.8, zorder=0)

    top = np.percentile(y, 99.5) if len(y) else 1
    span = max(top - np.percentile(y, 0.5), 0.5) if len(y) else 1
    rows = [d for d in detectors if shown[d]]
    for k, d in enumerate(rows):
        hits, extra, missed = found[d]
        level = top + span * (0.25 + 0.22 * k)
        inside = lambda v: v[(v >= start) & (v < start + WINDOW)]
        axECG.plot(inside(hits), [level] * len(inside(hits)), "v", color=COLOURS[d], markersize=7)
        axECG.plot(inside(extra), [level] * len(inside(extra)), "x", color="red", markersize=9, markeredgewidth=2.5)
        axECG.plot(inside(missed), [level] * len(inside(missed)), "o", markerfacecolor="none", markeredgecolor="red",
                   markersize=9, markeredgewidth=2)
        axECG.text(start + WINDOW * 1.003, level, d, color=COLOURS[d], va="center", fontsize=8, clip_on=False)
    axECG.set_xlim(start, start + WINDOW)
    axECG.set_ylim(np.percentile(y, 0.5) - 0.1 * span if len(y) else -1, top + span * (0.4 + 0.22 * len(rows)))
    axECG.set_xlabel("time (s)")
    axECG.set_ylabel("mV")
    errors = state["errors"]
    axECG.set_title(f"Record {state['record']}:   ▼ correct beat    ✕ false beat    ○ missed beat"
                    f"    (grey lines = annotated beats, {len(errors)} error locations)", fontsize=10)

    fig.canvas.draw_idle()


def draw_record():
    record = state["record"]
    state["errors"] = error_times()

    #HRV difference from the reference, the main result for this record
    axHRV.clear()
    rows = det[(det.record == record) & det.detector.isin([d for d in detectors if shown[d]])]
    n = max(len(rows), 1)
    width = 0.8 / n
    for k, (_, row) in enumerate(rows.iterrows()):
        values = [row[f"d_{m}"] for m, _ in METRICS]
        x = np.arange(len(METRICS)) + (k - (n - 1) / 2) * width
        bars = axHRV.bar(x, values, width, color=COLOURS[row.detector], label=row.detector)
        for bx, v in zip(x, values):
            axHRV.text(bx, v, f"{v:+.0f}%" if abs(v) >= 10 else f"{v:+.1f}%", ha="center", va="bottom" if v >= 0 else "top", fontsize=7)
    axHRV.axhline(0, color="black", linewidth=1)
    biggest = rows[[f"d_{m}" for m, _ in METRICS]].abs().values.max() if len(rows) else 0
    if biggest > 50:
        axHRV.set_yscale("symlog", linthresh=10)  #log scale only when errors are large
    else:
        axHRV.margins(y=0.25)
    axHRV.set_xticks(range(len(METRICS)))
    axHRV.set_xticklabels([label for _, label in METRICS])
    axHRV.set_ylabel("difference from reference (%)")
    ref_row = ref.loc[record]
    axHRV.set_title(f"HRV error per detector   (reference: SDNN {ref_row.SDNN:.1f} ms, RMSSD {ref_row.RMSSD:.1f} ms, "
                    f"LF {ref_row.LF:.0f} ms², HF {ref_row.HF:.0f} ms², LF/HF {ref_row.LF_HF:.2f})", fontsize=10)
    axHRV.legend(fontsize=8, ncol=3)
    axHRV.grid(True, axis="y", alpha=0.3)
    draw_ecg()


def on_record(label):
    label = label.strip()
    if label not in records:
        return
    state["record"] = label
    state["start"] = 300.0
    if recordBox.text != label:
        recordBox.set_val(label) #calls on_record again
        return
    draw_record()


def step_record(direction):
    current = state["record"]
    i = SHOWN.index(current) + direction if current in SHOWN else 0
    on_record(SHOWN[i % len(SHOWN)])


def on_show(label):
    shown[label] = not shown[label]
    draw_record()


def jump(direction):
    errors = state["errors"]
    centre = state["start"] + WINDOW / 2
    target = errors[errors > centre + 0.5] if direction > 0 else errors[errors < centre - 0.5][::-1]
    if len(target):
        state["start"] = max(target[0] - WINDOW / 2, 0)
        draw_ecg()


SHOWN = [r for r in SHOWN if r in records] or records[:1]
recordBox = TextBox(axRecord, "Record ", initial=SHOWN[0])
recordBox.on_submit(on_record)
recPrevButton = Button(axRecPrev, "◀ record")
recNextButton = Button(axRecNext, "record ▶")
recPrevButton.on_clicked(lambda e: step_record(-1))
recNextButton.on_clicked(lambda e: step_record(+1))
showButtons = CheckButtons(axShow, detectors, [shown[d] for d in detectors])
for label, d in zip(showButtons.labels, detectors):
    label.set_color(COLOURS[d])
showButtons.on_clicked(on_show)
prevButton = Button(axPrev, "◀ prev error")
nextButton = Button(axNext, "next error ▶")
prevButton.on_clicked(lambda e: jump(-1))
nextButton.on_clicked(lambda e: jump(+1))


def save_beat_figure(record, detectors_to_show, seconds=10.0):
    #the window starts after the first 5 minutes at the first detection error
    fs, signal, refPeaks, found = load(record)
    errors = []
    for d in detectors_to_show:
        _, extra, missed = found[d]
        errors.extend(extra)
        errors.extend(missed)
    errors = np.sort(np.array(errors))
    later = errors[errors >= 300.0]
    start = max(later[0] - 2.0, 300.0) if len(later) else 0.0
    a, b = int(start * fs), int((start + seconds) * fs)
    t = np.arange(a, min(b, len(signal))) / fs
    y = signal[a:b]
    f, ax = plt.subplots(figsize=(11, 4.2))
    ax.plot(t, y, color="#555555", linewidth=1.0)
    top = np.percentile(y, 99.5)
    for k, d in enumerate(detectors_to_show):
        hits, extra, missed = found[d]
        level = top + 0.45 + 0.4 * k
        inside = lambda v: v[(v >= start) & (v < start + seconds)]
        ax.plot(inside(hits), [level] * len(inside(hits)), "v", color=COLOURS[d], markersize=9)
        ax.plot(inside(extra), [level] * len(inside(extra)), "x", color="red", markersize=10, markeredgewidth=2.5)
        ax.plot(inside(missed), [level] * len(inside(missed)), "o", markerfacecolor="none", markeredgecolor="red",
                markersize=10, markeredgewidth=2)
        ax.text(start - 0.15, level, d, color=COLOURS[d], ha="right", va="center", fontsize=11, fontweight="bold")
    ax.plot([], [], "v", color="grey", markersize=9, label="detected beat")
    ax.plot([], [], "o", markerfacecolor="none", markeredgecolor="red", markersize=10, markeredgewidth=2, label="missed beat")
    ax.plot([], [], "x", color="red", markersize=10, markeredgewidth=2.5, label="false beat")
    ax.legend(loc="lower right", fontsize=10, ncol=3)
    ax.set_xlim(start, start + seconds)
    ax.set_ylim(np.percentile(y, 0.5) - 0.75, top + 0.45 + 0.4 * len(detectors_to_show))
    ax.set_yticks([])
    ax.set_xlabel("time (s)", fontsize=11)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    f.tight_layout()
    f.savefig(f"{RESULTS}/plots/beats_{record}.png", dpi=200)
    plt.close(f)


if save:
    #one ECG example 
    save_beat_figure(EXAMPLE if EXAMPLE in info.index else SHOWN[0], ["Pan-Tompkins", "Hamilton", "Hilbert"])
on_record(SHOWN[0])
draw_record()


#window 2. noise experiment
#records 118 and 119 with more and more noise. "clean" is the original record
LEVELS = [("clean", ""), ("24", "e24"), ("18", "e18"), ("12", "e12"), ("6", "e06"), ("0", "e00"), ("-6", "e_6")]
if noisy:
    fig3, axes = plt.subplots(3, 2, figsize=(14, 10), sharex=True, num="Noise experiment")
    for col, base in enumerate(["118", "119"]):
        names = [base + suffix for _, suffix in LEVELS]
        for row, (column, label) in enumerate([("Se", "Se (%)"), ("PPV", "PPV (%)"), ("d_RMSSD", "|ΔRMSSD| (%)")]):
            ax = axes[row, col]
            for d in detectors:
                values = []
                for name in names:
                    v = det[(det.record == name) & (det.detector == d)][column]
                    values.append(abs(v.iloc[0]) if len(v) else np.nan)
                main = d in ("Pan-Tompkins", "Hamilton", "Hilbert")
                ax.plot(range(len(LEVELS)), values, "o-" if main else "o--", color=COLOURS[d], linewidth=2.2 if main else 1.1, markersize=5 if main else 3, label=d)
            ax.set_ylabel(label)
            ax.grid(True, alpha=0.3)
            if column == "d_RMSSD":
                ax.set_yscale("symlog", linthresh=10)
            if row == 0:
                ax.set_title(f"Record {base} with electrode-motion noise", fontsize=11)
        axes[2, col].set_xticks(range(len(LEVELS)))
        axes[2, col].set_xticklabels([label for label, _ in LEVELS])
        axes[2, col].set_xlabel("signal-to-noise ratio (dB)")
    axes[0, 0].legend(fontsize=8, ncol=2)
    fig3.tight_layout()
    if save:
        fig3.savefig(f"{RESULTS}/plots/noise.png", dpi=150)


#window 3. FFT vs Goertzel
timing = f"{RESULTS}/fft_goertzel.csv"
if os.path.exists(timing):
    t = pd.read_csv(timing)
    fig4, ax = plt.subplots(figsize=(9, 5.5), num="FFT vs Goertzel")
    ax.plot(t.frequencies, t.fft_ms, "o-", color="#1f4e9c", linewidth=2, label="FFT (whole spectrum)")
    ax.plot(t.frequencies, t.goertzel_ms, "s-", color="#c0561b", linewidth=2, label="Goertzel (one recurrence per frequency)")
    #where the two lines cross, found between the last point where Goertzel is faster and the first where it is not
    faster = t.goertzel_ms < t.fft_ms
    for i in range(1, len(t)):
        if faster.iloc[i - 1] and not faster.iloc[i]:
            x0, x1 = np.log(t.frequencies.iloc[i - 1]), np.log(t.frequencies.iloc[i])
            d0 = np.log(t.goertzel_ms.iloc[i - 1] / t.fft_ms.iloc[i - 1])
            d1 = np.log(t.goertzel_ms.iloc[i] / t.fft_ms.iloc[i])
            k = np.exp(x0 + (x1 - x0) * d0 / (d0 - d1))
            ax.axvline(k, color="grey", linestyle=":", linewidth=1.5)
            ax.text(k * 1.08, t.goertzel_ms.min() * 1.5, f"lines cross at about {k:.0f} frequencies", color="grey", fontsize=9)
            break
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("number of frequencies needed")
    ax.set_ylabel("time per 5-minute segment (ms)")
    ax.set_title("FFT vs Goertzel on one 5-minute HRV segment", fontsize=11)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)
    fig4.tight_layout()
    if save:
        fig4.savefig(f"{RESULTS}/plots/fft_goertzel.png", dpi=150)

if save:
    print(f"Plots saved to {RESULTS}/plots/")

plt.show()