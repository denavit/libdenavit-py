from libdenavit import opensees as ops
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import matplotlib.cm as color
from matplotlib.colors import Normalize
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter
import numpy as np
import os

def get_node_coords():
    node_coords = dict()
    node_tags = ops.getNodeTags()
    for i in node_tags:
        node_coords[i] = ops.nodeCoord(i)
    return node_coords

def get_node_coords_and_disp():
    node_coords = dict()
    node_disp = dict()
    node_tags = ops.getNodeTags()
    for i in node_tags:
        node_coords[i] = ops.nodeCoord(i)
        node_disp[i] = ops.nodeDisp(i)
    return (node_coords,node_disp)
    
def get_element_nodes():
    element_nodes = dict()
    ele_tags = ops.getEleTags()
    for i in ele_tags:
        element_nodes[i] = ops.eleNodes(i)
    return element_nodes

def plot_undeformed_2d(axis_equal=False):
    node_coords = get_node_coords()
    element_nodes = get_element_nodes()
    fig = plt.figure()
    for i in element_nodes:
        coordi = node_coords[element_nodes[i][0]]
        coordj = node_coords[element_nodes[i][1]]
        xplt = [coordi[0],coordj[0]]
        yplt = [coordi[1],coordj[1]]
        plt.plot(xplt,yplt,'ko-')
    if axis_equal:
        plt.gca().axis('equal')
    plt.show()
    return
    
def plot_deformed_2d(scale_factor=1.0,show_undeformed=True,axis_equal=False):
    (node_coords,node_disp) = get_node_coords_and_disp()
    element_nodes = get_element_nodes()
    fig = plt.figure()
    if show_undeformed:
        for i in element_nodes:
            coordi = node_coords[element_nodes[i][0]]
            coordj = node_coords[element_nodes[i][1]]
            xplt = [coordi[0],coordj[0]]
            yplt = [coordi[1],coordj[1]]
            plt.plot(xplt,yplt,'-',color='lightgrey')
    for i in element_nodes:
        coordi = node_coords[element_nodes[i][0]]
        coordj = node_coords[element_nodes[i][1]]
        dispi = node_disp[element_nodes[i][0]]
        dispj = node_disp[element_nodes[i][1]]
        xplt = [coordi[0]+scale_factor*dispi[0],coordj[0]+scale_factor*dispj[0]]
        yplt = [coordi[1]+scale_factor*dispi[1],coordj[1]+scale_factor*dispj[1]]
        plt.plot(xplt,yplt,'-',color='k')
    if axis_equal:
        plt.gca().axis('equal')
    plt.show()
    plt.closefig()
    return

def plot_PMM_Interaction_values(PMM_values, axis_equal=False, cmap_name='coolwarm', show=True):
    PMM_dict = dict(PMM_values)
    node_coords = get_node_coords()
    element_nodes = get_element_nodes()

    fig, ax = plt.subplots(figsize=(7, 6))

    PMM_all = list(PMM_dict.values())
    vmin, vmax = min(PMM_all), max(PMM_all)
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.get_cmap(cmap_name)

    for ele_tag, nodes in element_nodes.items():
        coord_i = node_coords[nodes[0]]
        coord_j = node_coords[nodes[1]]
        xplt = [coord_i[0], coord_j[0]]
        yplt = [coord_i[1], coord_j[1]]

        val = PMM_dict.get(ele_tag)
        color = cmap(norm(val)) if val is not None else 'gray'
        ax.plot(xplt, yplt, color=color, lw=3, marker='o', markersize=4)

        xm, ym = (coord_i[0] + coord_j[0]) / 2, (coord_i[1] + coord_j[1]) / 2
        if val is not None:
            ax.text(xm, ym + 0.02, f"{val:.2f}", color='black', fontsize=8,
                    ha='center', va='bottom')

    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax, shrink=0.8)
    cbar.set_label("P–M–M Interaction Ratio")

    if axis_equal:
        ax.axis('equal')

    ax.set_xlabel("X-coordinate")
    ax.set_ylabel("Y-coordinate")
    ax.set_title("Undeformed 2D Frame with P–M–M Interaction Colormap")
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax

def plot_sfd(scale=0.002, axis_equal=False, annotate=True, show=True):
    """
    Plot Shear Force Diagram (SFD) for the current OpenSees 2D frame model.
    Styled consistently with PMM interaction plots (uniform fonts and layout).

    Parameters
    ----------
    scale : float
        Scaling factor for shear diagram size relative to element length.
    axis_equal : bool
        If True, keeps equal aspect ratio for X and Y axes.
    annotate : bool
        If True, shows global max/min shear in the figure corner.
    show : bool
        If True, displays the plot interactively.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes
    """

    # --- Get model geometry ---
    node_coords = get_node_coords()
    element_nodes = get_element_nodes()

    # --- Initialize figure ---
    fig, ax = plt.subplots(figsize=(7, 6))

    all_V = []  # store all shear values

    # --- Loop through each element ---
    for ele_tag, nodes in element_nodes.items():
        ni, nj = nodes
        xi, yi = node_coords[ni]
        xj, yj = node_coords[nj]

        L = np.sqrt((xj - xi)**2 + (yj - yi)**2)
        dx, dy = (xj - xi) / L, (yj - yi) / L
        n_perp = np.array([-dy, dx])

        try:
            forces = ops.eleResponse(ele_tag, 'localForce')
        except:
            print(f"[warn] Skipping element {ele_tag} (no localForce data).")
            continue

        Vi, Vj = forces[1], -forces[4]
        all_V.extend([Vi, Vj])

        # Base element line
        ax.plot([xi, xj], [yi, yj], 'k-', lw=1.5, zorder=1)

        # Shear diagram lines
        shear_coords = np.array([
            [xi, yi] + scale * Vi * n_perp,
            [xj, yj] + scale * Vj * n_perp
        ])

        ax.plot([xi, shear_coords[0, 0]], [yi, shear_coords[0, 1]], 'b--', lw=1)
        ax.plot([xj, shear_coords[1, 0]], [yj, shear_coords[1, 1]], 'b--', lw=1)
        ax.plot(
            [shear_coords[0, 0], shear_coords[1, 0]],
            [shear_coords[0, 1], shear_coords[1, 1]],
            'b-', lw=2
        )

    # --- Annotate global extrema ---
    if annotate and all_V:
        Vmax, Vmin = max(all_V), min(all_V)
        ax.text(0.02, 0.97, f"Max V = {Vmax:.2f}", transform=ax.transAxes,
                fontsize=9, fontfamily='Times New Roman', color='blue', ha='left', va='top')
        ax.text(0.02, 0.93, f"Min V = {Vmin:.2f}", transform=ax.transAxes,
                fontsize=9, fontfamily='Times New Roman', color='blue', ha='left', va='top')

    # --- Axis formatting ---
    if axis_equal:
        ax.axis('equal')

    ax.set_xlabel("X-coordinate", fontsize=11, fontfamily='Times New Roman')
    ax.set_ylabel("Y-coordinate", fontsize=11, fontfamily='Times New Roman')
    ax.set_title("Undeformed 2D Frame – Shear Force Diagram",
                 fontsize=13, fontfamily='Times New Roman', fontweight='normal')

    # --- Tick parameters for consistency ---
    ax.tick_params(axis='both', labelsize=9)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontfamily('Times New Roman')

    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax

def plot_bmd(scale=0.001, axis_equal=True, annotate=True, show=True, save_path=None):
    """
    Plot Bending Moment Diagram (BMD) for the current OpenSees 2D frame model.

    Parameters
    ----------
    scale : float
        Scaling factor for moment diagram size relative to element length.
    axis_equal : bool
        If True, keeps equal aspect ratio for X and Y axes.
    annotate : bool
        If True, shows global max/min moment in the figure corner.
    show : bool
        If True, displays the plot interactively.
    save_path : str or None
        File path to save the figure (e.g., 'Frame_1/BMD.png'). If None, figure is not saved.
    """

    node_coords = get_node_coords()
    element_nodes = get_element_nodes()

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_title("Bending Moment Diagram (BMD)", fontfamily='Times New Roman', fontsize=14)
    ax.set_xlabel("X", fontfamily='Times New Roman', fontsize=12)
    ax.set_ylabel("Y", fontfamily='Times New Roman', fontsize=12)

    all_M = []  # store all moment values

    for ele_tag, nodes in element_nodes.items():
        ni, nj = nodes
        xi, yi = node_coords[ni]
        xj, yj = node_coords[nj]

        L = np.sqrt((xj - xi)**2 + (yj - yi)**2)
        dx, dy = (xj - xi) / L, (yj - yi) / L
        n_perp = np.array([-dy, dx])

        try:
            forces = ops.eleResponse(ele_tag, 'localForce')
        except:
            print(f"[warn] Skipping element {ele_tag} (no localForce data).")
            continue

        Mi, Mj = forces[2], -forces[5]
        all_M.extend([Mi, Mj])

        # Base element line
        ax.plot([xi, xj], [yi, yj], 'k-', lw=1.2)

        # Moment diagram
        moment_coords = np.array([
            [xi, yi] + scale * Mi * n_perp,
            [xj, yj] + scale * Mj * n_perp
        ])

        ax.plot([xi, moment_coords[0, 0]], [yi, moment_coords[0, 1]], 'r--', lw=1)
        ax.plot([xj, moment_coords[1, 0]], [yj, moment_coords[1, 1]], 'r--', lw=1)
        ax.plot(
            [moment_coords[0, 0], moment_coords[1, 0]],
            [moment_coords[0, 1], moment_coords[1, 1]],
            'r-', lw=2
        )

    # Global extrema annotation
    if annotate and all_M:
        Mmax, Mmin = max(all_M), min(all_M)
        ax.text(0.02, 0.97, f"Max M = {Mmax:.2f}", transform=ax.transAxes,
                fontsize=10, fontfamily='Times New Roman', color='red', ha='left', va='top')
        ax.text(0.02, 0.92, f"Min M = {Mmin:.2f}", transform=ax.transAxes,
                fontsize=10, fontfamily='Times New Roman', color='red', ha='left', va='top')

    if axis_equal:
        ax.axis('equal')

    ax.legend(handles=[plt.Line2D([0], [0], color='r', lw=2, label='Moment')],
              loc='best', fontsize=9, frameon=False)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"BMD saved to: {save_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)

def plot_afd(scale=0.05, axis_equal=True, annotate=True, show=True, save_path=None):
    """
    Plot Axial Force Diagram (AFD) for the current OpenSees 2D frame model.
    (Plotted perpendicular to the member, consistent with SFD and BMD.)

    Parameters
    ----------
    scale : float
        Scaling factor for axial force diagram relative to element length.
    axis_equal : bool
        If True, keeps equal aspect ratio for X and Y axes.
    annotate : bool
        If True, shows global max/min axial forces in the figure corner.
    show : bool
        If True, displays the plot interactively.
    save_path : str or None
        File path to save the figure (e.g., 'Frame_1/AFD.png').
    """

    # --- Get geometry ---
    node_coords = get_node_coords()
    element_nodes = get_element_nodes()

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_title("Axial Force Diagram (AFD)", fontfamily='Times New Roman', fontsize=14)
    ax.set_xlabel("X", fontfamily='Times New Roman', fontsize=12)
    ax.set_ylabel("Y", fontfamily='Times New Roman', fontsize=12)

    all_P = []  # store all axial forces

    # --- Loop over all elements ---
    for ele_tag, nodes in element_nodes.items():
        ni, nj = nodes
        xi, yi = node_coords[ni]
        xj, yj = node_coords[nj]

        # element vector and length
        L = np.sqrt((xj - xi)**2 + (yj - yi)**2)
        dx, dy = (xj - xi) / L, (yj - yi) / L

        # direction perpendicular to element (same as SFD/BMD)
        n_perp = np.array([-dy, dx])

        try:
            forces = ops.eleResponse(ele_tag, 'localForce')
        except:
            print(f"[warn] Skipping element {ele_tag} (no localForce data).")
            continue

        Pi, Pj = -forces[0], forces[3]  # axial forces (tension +ve)
        all_P.extend([Pi, Pj])

        # Base element line
        ax.plot([xi, xj], [yi, yj], 'k-', lw=1.2)

        # Axial force diagram (perpendicular projection)
        axial_coords = np.array([
            [xi, yi] + scale * Pi * n_perp,
            [xj, yj] + scale * Pj * n_perp
        ])

        # Vertical projection lines
        ax.plot([xi, axial_coords[0, 0]], [yi, axial_coords[0, 1]], 'g--', lw=1)
        ax.plot([xj, axial_coords[1, 0]], [yj, axial_coords[1, 1]], 'g--', lw=1)

        # Connect the tips (diagram body)
        ax.plot(
            [axial_coords[0, 0], axial_coords[1, 0]],
            [axial_coords[0, 1], axial_coords[1, 1]],
            'g-', lw=2
        )

    # --- Global extrema annotation ---
    if annotate and all_P:
        Pmax, Pmin = max(all_P), min(all_P)
        ax.text(0.02, 0.97, f"Max P = {Pmax:.2f}", transform=ax.transAxes,
                fontsize=10, fontfamily='Times New Roman', color='green', ha='left', va='top')
        ax.text(0.02, 0.92, f"Min P = {Pmin:.2f}", transform=ax.transAxes,
                fontsize=10, fontfamily='Times New Roman', color='green', ha='left', va='top')

    # --- Layout ---
    if axis_equal:
        ax.axis('equal')

    ax.legend(handles=[plt.Line2D([0], [0], color='g', lw=2, label='Axial')],
              loc='best', fontsize=9, frameon=False)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"AFD saved to: {save_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)


def animate_PMM_evolution(PMM_steps,title, cmap_name='coolwarm',
                          interval=10, save_path=None,
                          fps=30, dpi=150):
    """
    Animate the evolution of P–M–M interaction values for the current OpenSees model.
    If save_path is given, uses extension to choose writer:
        .mp4/.mov/.m4v -> FFMpegWriter
        .gif           -> PillowWriter
        otherwise      -> saves a single PNG (first frame)
    """
    # --- Get coordinates and connectivity from model ---
    node_coords = get_node_coords()
    element_nodes = get_element_nodes()

    # --- Normalize color scale across all steps ---
    all_vals = [v for step in PMM_steps for (_, v) in step]
    vmin, vmax = min(all_vals), max(all_vals)
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = color.get_cmap(cmap_name)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.axis('equal')
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    # --- Initialize plot lines and labels ---
    lines, labels = {}, {}
    initial_PMM = dict(PMM_steps[0])
    for ele_tag, nodes in element_nodes.items():
        xi, yi = node_coords[nodes[0]]
        xj, yj = node_coords[nodes[1]]
        (line,) = ax.plot(
            [xi, xj], [yi, yj],
            color=cmap(norm(initial_PMM.get(ele_tag, 0))),
            lw=3, marker='o', markersize=4,
        )
        lines[ele_tag] = line

        xm, ym = (xi + xj) / 2, (yi + yj) / 2
        val = initial_PMM.get(ele_tag, 0)
        labels[ele_tag] = ax.text(
            xm, ym + 0.02, f"{val:.2f}",
            fontsize=8, color='black',
            ha='center', va='bottom'
        )

    sm = color.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax, shrink=0.8)
    cbar.set_label("P–M Interaction Ratio")

    def update(frame_idx):
        PMM_dict = dict(PMM_steps[frame_idx])
        for ele_tag in lines:
            val = PMM_dict.get(ele_tag, 0)
            c = cmap(norm(val))
            lines[ele_tag].set_color(c)
            labels[ele_tag].set_text(f"{val:.2f}")
        ax.set_title(f"{title} (Step {frame_idx + 1})")
        # returning artists is helpful for blitting, but blit=False so this is optional
        return list(lines.values()) + list(labels.values())

    ani = FuncAnimation(
        fig, update,
        frames=len(PMM_steps),
        interval=interval,
        blit=False,
        repeat=True,
    )

    if save_path:
        ext = os.path.splitext(save_path)[1].lower()

        if ext in [".mp4", ".mov", ".m4v"]:
            # Requires ffmpeg installed and on PATH
            writer = FFMpegWriter(fps=fps, bitrate=1800)
            ani.save(save_path, writer=writer, dpi=dpi)
            print(f"MP4 animation saved to {save_path}")

        elif ext == ".gif":
            # Uses Pillow; should work out of the box if pillow is installed
            writer = PillowWriter(fps=fps)
            ani.save(save_path, writer=writer, dpi=dpi)
            print(f"GIF animation saved to {save_path}")

        else:
            # Fallback: just save first frame as PNG
            # (what you’re effectively seeing now)
            fig.savefig(save_path, dpi=dpi)
            print(f"Unknown extension '{ext}', saved first frame as image at {save_path}")
    else:
        plt.show()

    return ani
