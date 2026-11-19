"""
Modal analysis of a cantilever beam: Gmsh + CalculiX + Python
==============================================================
A fully scripted 3D modal analysis:
  1. Hexahedral mesh with the Gmsh Python API
  2. Eigenfrequency analysis in CalculiX (run from the command line)
  3. Automatic mode identification (vertical bending, lateral bending, torsion, axial)
  4. Verification of every frequency against analytical theory
  5. Mode shape comparison with Euler-Bernoulli (MAC)
  6. Effective modal mass (participation)
  7. Mesh convergence of the frequencies
  8. Design study: beam length vs. resonance with a motor
  9. 3D mode shape images and an animated GIF

Units: mm, N, tonne, s  ->  stress in MPa, frequency in Hz
"""
import os
import subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib import animation, cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import gmsh

# ============================ SETTINGS ============================
CCX = r"C:\Users\enzoz\Downloads\PrePoMax\PrePoMax v3.0.0\Solver\ccx_dynamic.exe"

L, B, H = 200.0, 20.0, 10.0       # length (x), width (y), height (z) [mm]
E, NU = 210000.0, 0.3             # steel [MPa]
RHO = 7.85e-9                     # steel density [tonne/mm^3]
N_MODES = 10
MESH = (40, 4, 4)                 # elements in x, y, z
OUT = "results"
# ==================================================================

G = E / (2 * (1 + NU))
BETA_L = [1.87510, 4.69409, 7.85476, 10.99554, 14.13717]  # cantilever roots
os.makedirs(OUT, exist_ok=True)


# ------------------------------------------------------------------
# 1. MESH
# ------------------------------------------------------------------
def build_mesh(length, width, height, nx, ny, nz):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("beam")
    p = gmsh.model.geo.addPoint(0, 0, 0)
    line = gmsh.model.geo.extrude([(0, p)], length, 0, 0, [nx])
    surf = gmsh.model.geo.extrude([line[1]], 0, width, 0, [ny], recombine=True)
    gmsh.model.geo.extrude([surf[1]], 0, 0, height, [nz], recombine=True)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    nodes = {int(t): np.array(coords[3 * i:3 * i + 3]) for i, t in enumerate(tags)}
    _, etags, enodes = gmsh.model.mesh.getElements(3)
    elements = []
    for i, t in enumerate(etags[0]):
        conn = [int(n) for n in enodes[0][8 * i:8 * i + 8]]
        p0, p1, p3, p4 = (nodes[conn[j]] for j in (0, 1, 3, 4))
        if np.dot(np.cross(p1 - p0, p3 - p0), p4 - p0) < 0:   # positive Jacobian
            conn = [conn[0], conn[3], conn[2], conn[1], conn[4], conn[7], conn[6], conn[5]]
        elements.append((int(t), conn))
    gmsh.finalize()
    return nodes, elements


# ------------------------------------------------------------------
# 2. CALCULIX INPUT AND RUN
# ------------------------------------------------------------------
def write_inp(nodes, elements, path, n_modes):
    fix = [n for n, c in nodes.items() if abs(c[0]) < 1e-6]
    with open(path, "w") as f:
        f.write("*NODE, NSET=NALL\n")
        for n, c in nodes.items():
            f.write(f"{n}, {c[0]:.6f}, {c[1]:.6f}, {c[2]:.6f}\n")
        f.write("*ELEMENT, TYPE=C3D8I, ELSET=EALL\n")
        for e, conn in elements:
            f.write(f"{e}, " + ", ".join(map(str, conn)) + "\n")
        f.write("*NSET, NSET=FIX\n")
        for i in range(0, len(fix), 10):
            f.write(", ".join(map(str, fix[i:i + 10])) + "\n")
        f.write(f"*MATERIAL, NAME=STEEL\n*ELASTIC\n{E}, {NU}\n*DENSITY\n{RHO}\n")
        f.write("*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL\n")
        f.write("*BOUNDARY\nFIX, 1, 3, 0.0\n")
        f.write(f"*STEP\n*FREQUENCY\n{n_modes}\n*NODE FILE\nU\n*END STEP\n")


def run_ccx(job):
    res = subprocess.run([CCX, job], capture_output=True, text=True)
    if res.returncode != 0 or not os.path.exists(job + ".frd"):
        print(res.stdout[-2000:])
        raise RuntimeError("CalculiX failed")


def solve(length, width, height, mesh, job, n_modes=N_MODES):
    nodes, elements = build_mesh(length, width, height, *mesh)
    write_inp(nodes, elements, job + ".inp", n_modes)
    run_ccx(job)
    freqs, eff_mass, total_mass = read_dat(job + ".dat")
    shapes = read_frd(job + ".frd")
    return nodes, elements, freqs, eff_mass, total_mass, shapes


# ------------------------------------------------------------------
# 3. READ RESULTS
# ------------------------------------------------------------------
def read_dat(path):
    """Frequencies, effective modal mass (x, y, z) and total mass."""
    freqs, eff, total = [], [], None
    section = None
    with open(path) as f:
        for line in f:
            if "E I G E N V A L U E   O U T P U T" in line:
                section = "eig"
                continue
            if "P A R T I C I P A T I O N" in line:
                section = None
                continue
            if "E F F E C T I V E   M O D A L   M A S S" in line:
                section = "eff"
                continue
            if "T O T A L   E F F E C T I V E   M A S S" in line:
                section = "tot"
                continue
            parts = line.split()
            if section == "eig" and len(parts) == 5 and parts[0].isdigit():
                freqs.append(float(parts[3]))
            elif section == "eff" and len(parts) == 7 and parts[0].isdigit():
                eff.append([float(v) for v in parts[1:4]])
            elif section == "tot" and len(parts) == 6:
                total = [float(v) for v in parts[0:3]]
                section = None
    return np.array(freqs), np.array(eff), np.array(total)


def read_frd(path):
    """Mode shapes: list of {node: (ux, uy, uz)}, one per mode (fixed-width .frd format)."""
    shapes, current, reading = [], None, False
    with open(path) as f:
        for line in f:
            if line.startswith(" -4  DISP"):
                current, reading = {}, True
                continue
            if reading and line.startswith(" -1"):
                node = int(line[3:13])
                current[node] = np.array([float(line[13:25]), float(line[25:37]), float(line[37:49])])
            elif reading and line.startswith(" -3"):
                shapes.append(current)
                reading = False
    return shapes


# ------------------------------------------------------------------
# 4. MODE IDENTIFICATION
#    Each cross-section is split into rigid motions: mean ux, uy, uz and
#    twist about the centroid. The dominant one defines the mode type.
# ------------------------------------------------------------------
def section_motion(nodes, shape, width, height):
    xs = sorted({round(c[0], 6) for c in nodes.values()})
    yc, zc = width / 2, height / 2
    data = {x: [] for x in xs}
    for n, c in nodes.items():
        data[round(c[0], 6)].append((c, shape[n]))
    rows = []
    for x in xs:
        cs = np.array([d[0] for d in data[x]])
        us = np.array([d[1] for d in data[x]])
        dy, dz = cs[:, 1] - yc, cs[:, 2] - zc
        twist = np.sum(dy * us[:, 2] - dz * us[:, 1]) / np.sum(dy ** 2 + dz ** 2)
        rows.append([x, us[:, 0].mean(), us[:, 1].mean(), us[:, 2].mean(), twist])
    return np.array(rows)   # columns: x, ux, uy, uz, twist


def classify(nodes, shape, width, height):
    s = section_motion(nodes, shape, width, height)
    r_g = np.sqrt((width ** 2 + height ** 2) / 12)      # radius of gyration of the section
    energy = {
        "Axial": np.sum(s[:, 1] ** 2),
        "Lateral bending": np.sum(s[:, 2] ** 2),
        "Vertical bending": np.sum(s[:, 3] ** 2),
        "Torsion": np.sum((s[:, 4] * r_g) ** 2),
    }
    return max(energy, key=energy.get), s


# ------------------------------------------------------------------
# 5. ANALYTICAL THEORY
# ------------------------------------------------------------------
def theory(kind, n, length, width, height):
    A = width * height
    if kind == "Vertical bending":
        I = width * height ** 3 / 12
        return BETA_L[n - 1] ** 2 / (2 * np.pi * length ** 2) * np.sqrt(E * I / (RHO * A))
    if kind == "Lateral bending":
        I = height * width ** 3 / 12
        return BETA_L[n - 1] ** 2 / (2 * np.pi * length ** 2) * np.sqrt(E * I / (RHO * A))
    if kind == "Torsion":
        b, t = max(width, height), min(width, height)
        J = b * t ** 3 * (1 / 3 - 0.21 * t / b * (1 - t ** 4 / (12 * b ** 4)))   # Saint-Venant
        Ip = width * height * (width ** 2 + height ** 2) / 12
        return (2 * n - 1) / (4 * length) * np.sqrt(G * J / (RHO * Ip))
    if kind == "Axial":
        return (2 * n - 1) / (4 * length) * np.sqrt(E / RHO)


def eb_shape(n, x, length):
    """Euler-Bernoulli cantilever mode shape, normalized to 1 at the tip."""
    bl = BETA_L[n - 1]
    b = bl / length
    sigma = (np.cosh(bl) + np.cos(bl)) / (np.sinh(bl) + np.sin(bl))
    phi = np.cosh(b * x) - np.cos(b * x) - sigma * (np.sinh(b * x) - np.sin(b * x))
    return phi / phi[-1]


def mac(a, b):
    """Modal Assurance Criterion: 1 = identical shapes, 0 = unrelated."""
    return np.dot(a, b) ** 2 / (np.dot(a, a) * np.dot(b, b))


# ------------------------------------------------------------------
# 6. 3D PLOTS
# ------------------------------------------------------------------
FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


def boundary_faces(elements):
    count = {}
    for _, conn in elements:
        for f in FACES:
            quad = tuple(conn[i] for i in f)
            key = tuple(sorted(quad))
            count.setdefault(key, [0, quad])
            count[key][0] += 1
    return [v[1] for v in count.values() if v[0] == 1]


def draw_mode(ax, nodes, faces, shape, title, scale_frac=0.18, phase=1.0):
    umax = max(np.linalg.norm(u) for u in shape.values())
    scale = scale_frac * L / umax * phase
    mag = {n: np.linalg.norm(u) / umax for n, u in shape.items()}
    polys, colors = [], []
    for q in faces:
        polys.append([nodes[n] + scale * shape[n] for n in q])
        colors.append(cm.turbo(np.mean([mag[n] for n in q])))
    ax.add_collection3d(Poly3DCollection(polys, facecolors=colors, edgecolors=(0, 0, 0, 0.15), linewidths=0.2))
    # undeformed outline
    for y, z in [(0, 0), (B, 0), (0, H), (B, H)]:
        ax.plot([0, L], [y, y], [z, z], color="gray", lw=0.5, alpha=0.5)
    ax.set_xlim(0, L)
    ax.set_ylim(-0.25 * L + B / 2, 0.25 * L + B / 2)
    ax.set_zlim(-0.25 * L + H / 2, 0.25 * L + H / 2)
    ax.set_box_aspect((2, 1, 1))
    ax.view_init(elev=22, azim=-60)
    ax.set_axis_off()
    ax.set_title(title, fontsize=10)


# ==================================================================
# MAIN
# ==================================================================
def main():
    report = []

    def log(text=""):
        print(text)
        report.append(text)

    # ---------------- Baseline analysis ----------------
    nodes, elements, freqs, eff, total, shapes = solve(L, B, H, MESH, "modal_base")
    counters = {}
    modes = []
    for i, (f, shape) in enumerate(zip(freqs, shapes)):
        kind, sec = classify(nodes, shape, B, H)
        counters[kind] = counters.get(kind, 0) + 1
        n = counters[kind]
        f_th = theory(kind, n, L, B, H)
        modes.append(dict(mode=i + 1, f=f, kind=kind, n=n, f_th=f_th, sec=sec,
                          meff=eff[i] / total * 100))

    log(f"Steel cantilever beam {L:.0f} x {B:.0f} x {H:.0f} mm | mesh {MESH} | "
        f"{len(elements)} C3D8I elements | mass {RHO * L * B * H * 1e3:.4f} kg\n")
    log(f"{'Mode':>4} {'Type':<18} {'n':>2} {'FEA [Hz]':>10} {'Theory [Hz]':>12} {'Diff':>8} "
        f"{'Meff X':>7} {'Meff Y':>7} {'Meff Z':>7}")
    for m in modes:
        d = (m["f"] - m["f_th"]) / m["f_th"] * 100
        log(f"{m['mode']:>4} {m['kind']:<18} {m['n']:>2} {m['f']:>10.1f} {m['f_th']:>12.1f} {d:>+7.2f}% "
            f"{m['meff'][0]:>6.1f}% {m['meff'][1]:>6.1f}% {m['meff'][2]:>6.1f}%")
    cum = np.sum([m["meff"] for m in modes], axis=0)
    log(f"{'':>4} {'Cumulative effective mass':<45} {cum[0]:>6.1f}% {cum[1]:>6.1f}% {cum[2]:>6.1f}%")

    # ---------------- Mode shape vs Euler-Bernoulli ----------------
    log("\nMode shape check (vertical bending vs Euler-Bernoulli):")
    vert = [m for m in modes if m["kind"] == "Vertical bending"][:3]
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    colors = ["tab:blue", "tab:red", "tab:green"]
    for m, col in zip(vert, colors):
        x, uz = m["sec"][:, 0], m["sec"][:, 3]
        uz = uz / uz[-1]
        phi = eb_shape(m["n"], x, L)
        m["mac"] = mac(uz, phi)
        log(f"  Bending mode {m['n']}: MAC = {m['mac']:.5f}")
        xf = np.linspace(0, L, 300)
        ax.plot(xf, eb_shape(m["n"], xf, L), color=col, lw=1.5, label=f"Theory, mode {m['n']}")
        ax.plot(x, uz, "o", color=col, ms=4, label=f"FEA, mode {m['n']} (MAC = {m['mac']:.4f})")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("Position along the beam x [mm]")
    ax.set_ylabel("Normalized vertical displacement")
    ax.set_title("Vertical bending mode shapes: FEA vs Euler-Bernoulli")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(f"{OUT}/mode_shapes_vs_theory.png", dpi=150)
    plt.close(fig)

    # ---------------- 3D mode shapes ----------------
    faces = boundary_faces(elements)
    fig = plt.figure(figsize=(14, 8))
    for k in range(6):
        ax = fig.add_subplot(2, 3, k + 1, projection="3d")
        m = modes[k]
        draw_mode(ax, nodes, faces, shapes[k], f"Mode {m['mode']}: {m['kind']} {m['n']}\n{m['f']:.1f} Hz")
    fig.suptitle("First six mode shapes (color = normalized displacement)", fontsize=13)
    fig.tight_layout()
    fig.savefig(f"{OUT}/mode_shapes_3d.png", dpi=140)
    plt.close(fig)

    # ---------------- Animated GIF ----------------
    fig = plt.figure(figsize=(10, 6.5))
    axes = [fig.add_subplot(2, 2, k + 1, projection="3d") for k in range(4)]
    frames = 24

    def update(t):
        phase = np.sin(2 * np.pi * t / frames)
        for k, ax in enumerate(axes):
            ax.cla()
            m = modes[k]
            draw_mode(ax, nodes, faces, shapes[k], f"{m['kind']} {m['n']}: {m['f']:.1f} Hz", phase=phase)
        return []

    anim = animation.FuncAnimation(fig, update, frames=frames)
    anim.save(f"{OUT}/modes_animation.gif", writer=animation.PillowWriter(fps=12))
    plt.close(fig)

    # ---------------- Mesh convergence ----------------
    log("\nMesh convergence (first four modes):")
    meshes = [(10, 1, 1), (20, 2, 2), (40, 4, 4), (80, 8, 8)]
    conv = []
    for msh in meshes:
        nd, el, fr, *_ = solve(L, B, H, msh, f"modal_{msh[0]}x{msh[1]}x{msh[2]}", n_modes=4)
        conv.append((len(el), fr[:4]))
        log(f"  {len(el):>5} elements: " + "  ".join(f"{v:8.2f}" for v in fr[:4]) + "  Hz")
    f_fine, f_med = conv[-1][1], conv[-2][1]
    log("  Change medium -> fine: " + "  ".join(f"{(a - b) / b * 100:+.3f}%" for a, b in zip(f_fine, f_med)))

    fig, ax = plt.subplots(figsize=(7, 4.3))
    n_el = [c[0] for c in conv]
    for j in range(4):
        vals = [c[1][j] / conv[-1][1][j] for c in conv]
        ax.semilogx(n_el, vals, "o-", label=f"Mode {j + 1}")
    ax.axhline(1, color="gray", ls="--")
    ax.set_xlabel("Number of elements")
    ax.set_ylabel("Frequency / finest-mesh frequency")
    ax.set_title("Mesh convergence of the natural frequencies")
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{OUT}/mesh_convergence.png", dpi=150)
    plt.close(fig)

    # ---------------- Design study: resonance ----------------
    rpm = 1800
    f_exc = [rpm / 60, 2 * rpm / 60]          # 1x and 2x motor harmonics [Hz]
    band = 0.20                                # keep frequencies 20% away from excitation
    log(f"\nDesign study: motor at {rpm} rpm -> excitation at {f_exc[0]:.0f} Hz (1x) and {f_exc[1]:.0f} Hz (2x)")
    log(f"  Rule: the first vertical and lateral bending frequencies must stay at least {band * 100:.0f}% away from an excitation")
    lengths = [200, 250, 300, 350, 400, 450, 500, 550, 600]
    study = []
    for Ln in lengths:
        nx = int(round(40 * Ln / 200))
        _, _, fr, *_ = solve(Ln, B, H, (nx, 4, 4), f"modal_L{Ln}", n_modes=2)
        f1, f2 = fr[0], fr[1]          # 1st vertical and 1st lateral bending
        safe = all(abs(f - fe) / fe > band for f in (f1, f2) for fe in f_exc)
        study.append((Ln, f1, f2, safe))
        log(f"  L = {Ln:3d} mm: f1 = {f1:7.1f} Hz (theory {theory('Vertical bending', 1, Ln, B, H):7.1f}), "
            f"f2 = {f2:7.1f} Hz  ->  {'OK' if safe else 'RESONANCE RISK'}")

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    Ls = np.array([s[0] for s in study])
    Lf = np.linspace(180, 620, 300)
    ax.plot(Ls, [s[1] for s in study], "o", color="tab:blue", label="FEA, 1st vertical bending")
    ax.plot(Lf, [theory("Vertical bending", 1, v, B, H) for v in Lf], "-", color="tab:blue", alpha=0.5,
            label="Euler-Bernoulli, vertical")
    ax.plot(Ls, [s[2] for s in study], "s", color="tab:green", label="FEA, 1st lateral bending")
    ax.plot(Lf, [theory("Lateral bending", 1, v, B, H) for v in Lf], "-", color="tab:green", alpha=0.5,
            label="Euler-Bernoulli, lateral")
    for fe, lab in zip(f_exc, ["1x", "2x"]):
        ax.axhspan(fe * (1 - band), fe * (1 + band), color="tab:red", alpha=0.15)
        ax.axhline(fe, color="tab:red", lw=1)
        ax.text(190, fe * (1 + band) * 1.02, f"{lab} motor = {fe:.0f} Hz", color="tab:red", va="bottom", fontsize=8)
    for s in study:
        for f in (s[1], s[2]):
            if any(abs(f - fe) / fe <= band for fe in f_exc):
                ax.plot(s[0], f, "x", color="tab:red", ms=10, mew=2)
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_yticks([20, 30, 40, 60, 100, 150, 200, 300, 400])
    ax.set_xlabel("Beam length L [mm]")
    ax.set_ylabel("First natural frequency [Hz]")
    ax.set_title(f"Resonance check with a {rpm} rpm motor (red bands = ±{band * 100:.0f}%)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{OUT}/resonance_study.png", dpi=150)
    plt.close(fig)

    with open(f"{OUT}/report.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(report))
    print(f"\nFigures and report saved in the '{OUT}' folder.")


if __name__ == "__main__":
    main()
