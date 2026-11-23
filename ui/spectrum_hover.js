const {species, traces, mixture, resolvingPower, intensityLabel, unit} = __HOVER_METADATA__;
const plot = document.getElementById("{plot_id}");
const minimumIsotopeShare = 0.001;
// Suppress signals below one millionth of the relevant plotted peak.
const minimumPeakFraction = 1e-6;

const tooltip = document.createElement("div");
tooltip.setAttribute("role", "tooltip");
tooltip.hidden = true;
Object.assign(tooltip.style, {
    position: "fixed",
    zIndex: 10,
    maxWidth: "min(440px, calc(100vw - 36px))",
    maxHeight: "calc(100vh - 36px)",
    overflow: "auto",
    overscrollBehavior: "contain",
    background: "#fff",
    color: "#242424",
    border: "1px solid #b8b8b8",
    borderRadius: "4px",
    padding: "8px 10px",
    font: "13px/1.45 system-ui, sans-serif",
});
document.body.append(tooltip);

function contributions(index, mass) {
    const entries = species[index].isotopologues.map(([center, probability, label]) => {
        const sigma = (center / resolvingPower) / 2.355;
        const offset = (mass - center) / sigma;
        const gaussian = (1 / (sigma * Math.sqrt(2 * Math.PI))) * Math.exp(-0.5 * offset ** 2);
        return {label, signal: probability * gaussian};
    });
    const signal = entries.reduce((sum, entry) => sum + entry.signal, 0);
    return {entries, signal};
}

function addLine(text, bold = false) {
    const line = document.createElement("div");
    if (bold) line.style.fontWeight = "bold";
    line.textContent = text;
    tooltip.append(line);
}

function showIsotopologues({entries, signal}) {
    addLine("Isotopologues (% of this species’ local signal)");
    const visible = entries.filter(entry => entry.signal / signal >= minimumIsotopeShare);
    visible.sort((a, b) => b.signal - a.signal);
    for (const entry of visible) {
        addLine(`${entry.label} · ${(100 * entry.signal / signal).toFixed(2)}%`);
    }
    const omitted = entries.length - visible.length;
    if (omitted) {
        addLine(`${omitted} below 0.1% omitted`);
    }
}

function describePoint(point) {
    const trace = traces[point.curveNumber];
    const cutoff = trace.peak * minimumPeakFraction;
    if (!(point.y > cutoff)) return false;

    const isTotal = trace.species === null;
    addLine(isTotal ? "Hypothetical total" : species[trace.species].formula, true);
    addLine(`m/z ${point.x.toFixed(6)} Th`);
    addLine(`${intensityLabel} ${point.y.toPrecision(5)}${unit}`);

    if (isTotal) {
        const components = [];
        for (const [index, weight] of mixture.weights.entries()) {
            if (weight === 0) continue;
            const local = contributions(index, point.x);
            const intensity = local.signal * weight * mixture.scale;
            if (intensity > cutoff) components.push({index, local, intensity});
        }
        components.sort((a, b) => b.intensity - a.intensity);
        for (const {index, local, intensity} of components) {
            const share = 100 * intensity / point.y;
            addLine(`${species[index].formula} · ${share.toPrecision(4)}% of local mixture signal`, true);
            showIsotopologues(local);
        }
    } else {
        showIsotopologues(contributions(trace.species, point.x));
    }
    return true;
}

let hideTimer;
function hideTooltip() {
    tooltip.hidden = true;
}

plot.on("plotly_hover", ({points, event}) => {
    clearTimeout(hideTimer);
    tooltip.replaceChildren();
    if (!describePoint(points[0])) {
        hideTooltip();
        return;
    }

    tooltip.hidden = false;
    tooltip.scrollTop = 0;
    const anchor = points[0].bbox;
    const x = event?.clientX ?? anchor.x1;
    const y = event?.clientY ?? anchor.y1;
    const left = x + 8 + tooltip.offsetWidth <= window.innerWidth
        ? x + 8 : x - tooltip.offsetWidth - 8;
    tooltip.style.left = `${Math.max(0, left)}px`;
    tooltip.style.top = `${Math.max(0, Math.min(y, window.innerHeight - tooltip.offsetHeight))}px`;
});

plot.on("plotly_unhover", () => {
    hideTimer = setTimeout(hideTooltip, 180);
});
tooltip.addEventListener("mouseenter", () => clearTimeout(hideTimer));
tooltip.addEventListener("mouseleave", hideTooltip);
plot.on("plotly_relayout", hideTooltip);
plot.on("plotly_restyle", hideTooltip);
