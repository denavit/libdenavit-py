from libdenavit import opensees as ops
import matplotlib.pyplot as plt
from matplotlib import cm,colormaps
from matplotlib.colors import Normalize
import matplotlib.cm as color
from matplotlib.colors import Normalize
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter
import os
import numpy as np

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

def plot_element_values_2d(element_values, element_value_name='Element Value', axis_equal=False, show=True):

# Plots a scalar value associated with each element in a 2D OpenSees model.
# element_values must be a dictionary or an iterable of (element_tag, value) pairs.
# Each element tag must correspond to an element tag in the OpenSees model.
# Each element must have a single numerical value (e.g., P-M ratio, axial force, or utilization ratio).
# Example 1: element_values = {1: 0.25, 2: 0.50, 3: 0.85}
# Example 2: element_values = [(1, 0.25), (2, 0.50), (3, 0.85)]
# Elements without a supplied value are plotted in gray.
# element_value_name specifies the figure title and colorbar label.

    element_values = dict(element_values)
    node_coords = get_node_coords()
    element_nodes = get_element_nodes()
    fig, ax = plt.subplots()
    values = list(element_values.values())
    if not values:
        raise ValueError("No element values were provided.")
    vmin, vmax = min(values), max(values)
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = colormaps.get_cmap('coolwarm')
    for i in element_nodes:
        coordi = node_coords[element_nodes[i][0]]
        coordj = node_coords[element_nodes[i][1]]
        xplt = [coordi[0], coordj[0]]
        yplt = [coordi[1], coordj[1]]
        value = element_values.get(i)
        color = cmap(norm(value)) if value is not None else 'gray'
        ax.plot(xplt, yplt, 'o-', color=color)
        if value is not None:
            xm = (coordi[0] + coordj[0]) / 2
            ym = (coordi[1] + coordj[1]) / 2
            ax.text(xm, ym + 0.02, f"{value:.2f}")
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label(element_value_name)
    if axis_equal:
        ax.axis('equal')
    ax.set_xlabel("X-coordinate")
    ax.set_ylabel("Y-coordinate")
    ax.set_title(f"{element_value_name} Colormap")
    fig.tight_layout()
    if show:
        plt.show()
    return fig, ax

def plot_internal_force_diagram_2d(force_type='M', scale=None, axis_equal=False, annotate=True, show=True, save_path=None):
    force_options = {
        'P': {'indices': (0, 3), 'signs': (-1, 1), 'color': 'green', 'fill_color': 'lightgreen', 'scale': 0.05, 'name': 'Axial Force'},
        'V': {'indices': (1, 4), 'signs': (1, -1), 'color': 'blue', 'fill_color': 'lightblue', 'scale': 0.002, 'name': 'Shear Force'},
        'M': {'indices': (2, 5), 'signs': (1, -1), 'color': 'red', 'fill_color': 'lightcoral', 'scale': 0.001, 'name': 'Bending Moment'}
    }
    force_type = force_type.upper()
    if force_type not in force_options:
        raise ValueError("force_type must be 'P', 'V', or 'M'.")
    options = force_options[force_type]
    if scale is None:
        scale = options['scale']
    node_coords = get_node_coords()
    element_nodes = get_element_nodes()
    fig, ax = plt.subplots()
    all_values = []
    for i in element_nodes:
        coordi = np.array(node_coords[element_nodes[i][0]])
        coordj = np.array(node_coords[element_nodes[i][1]])
        L = np.linalg.norm(coordj - coordi)
        if L == 0:
            continue
        n_perp = np.array([-(coordj[1] - coordi[1]), coordj[0] - coordi[0]]) / L
        try:
            forces = ops.eleResponse(i, 'localForce')
            if forces is None or len(forces) < 6:
                print(f"Skipping element {i}: invalid localForce response.")
                continue
        except Exception as error:
            print(f"Skipping element {i}: {error}")
            continue
        index_i, index_j = options['indices']
        sign_i, sign_j = options['signs']
        value_i = sign_i * forces[index_i]
        value_j = sign_j * forces[index_j]
        all_values.extend([value_i, value_j])
        diagram_i = coordi + scale * value_i * n_perp
        diagram_j = coordj + scale * value_j * n_perp
        ax.fill([coordi[0], coordj[0], diagram_j[0], diagram_i[0]], [coordi[1], coordj[1], diagram_j[1], diagram_i[1]], color=options['fill_color'], alpha=0.35)
        ax.plot([coordi[0], coordj[0]], [coordi[1], coordj[1]], 'k-')
        ax.plot([coordi[0], diagram_i[0]], [coordi[1], diagram_i[1]], color=options['color'], linestyle='--')
        ax.plot([coordj[0], diagram_j[0]], [coordj[1], diagram_j[1]], color=options['color'], linestyle='--')
        ax.plot([diagram_i[0], diagram_j[0]], [diagram_i[1], diagram_j[1]], color=options['color'])
    if annotate and all_values:
        ax.text(0.02, 0.97, f"Max = {max(all_values):.2f}\nMin = {min(all_values):.2f}", transform=ax.transAxes, color=options['color'], va='top')
    if axis_equal:
        ax.axis('equal')
    ax.set_xlabel("X-coordinate")
    ax.set_ylabel("Y-coordinate")
    ax.set_title(f"Undeformed 2D Frame – {options['name']} Diagram")
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig, ax

def animate_element_values_2d(element_values_steps, element_value_name='Element Value', interval=10, axis_equal=False, show=True, save_path=None, fps=10):
# Animates the evolution of scalar values associated with elements in a 2D OpenSees model.
# element_values_steps must be a list containing element values at successive analysis steps. For element values look 'plot_element_values_2d'
# Each step can be a dictionary {element_tag: value} or an iterable of (element_tag, value) pairs.
# Each element tag must match an element tag in the OpenSees model.
# Each element must have a single numerical value per step (e.g., P-M ratio, axial force, or utilization ratio).
# Example 1: element_values_steps = [{1: 0.25, 2: 0.50}, {1: 0.40, 2: 0.65}, {1: 0.60, 2: 0.85}]
# Example 2: element_values_steps = [[(1, 0.25), (2, 0.50)], [(1, 0.40), (2, 0.65)], [(1, 0.60), (2, 0.85)]]
# Elements without a supplied value at a particular step are plotted in gray.
# element_value_name specifies the animation title and colorbar label.
# The color scale remains constant across all analysis steps for consistent comparison.
# interval specifies the delay between frames in milliseconds; fps specifies the saved video's frame rate.
# save_path supports MP4, MOV, M4V, and GIF formats.
    node_coords = get_node_coords()
    element_nodes = get_element_nodes()
    if not element_values_steps:
        raise ValueError("No element values were provided.")
    all_values = [value for step in element_values_steps for value in dict(step).values() if value is not None]
    if not all_values:
        raise ValueError("No element values were found.")
    norm = Normalize(vmin=min(all_values), vmax=max(all_values))
    cmap = colormaps.get_cmap('coolwarm')
    fig, ax = plt.subplots()
    lines = {}
    labels = {}
    initial_values = dict(element_values_steps[0])
    for i in element_nodes:
        coordi = node_coords[element_nodes[i][0]]
        coordj = node_coords[element_nodes[i][1]]
        xplt = [coordi[0], coordj[0]]
        yplt = [coordi[1], coordj[1]]
        value = initial_values.get(i)
        color = cmap(norm(value)) if value is not None else 'gray'
        line, = ax.plot(xplt, yplt, 'o-', color=color)
        lines[i] = line
        xm = (coordi[0] + coordj[0]) / 2
        ym = (coordi[1] + coordj[1]) / 2
        labels[i] = ax.text(xm, ym + 0.02, f"{value:.2f}" if value is not None else "")
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label(element_value_name)
    if axis_equal:
        ax.axis('equal')
    ax.set_xlabel("X-coordinate")
    ax.set_ylabel("Y-coordinate")
    def update(frame_index):
        values = dict(element_values_steps[frame_index])
        for i in lines:
            value = values.get(i)
            lines[i].set_color(cmap(norm(value)) if value is not None else 'gray')
            labels[i].set_text(f"{value:.2f}" if value is not None else "")
        ax.set_title(f"{element_value_name} (Step {frame_index + 1})")
        return list(lines.values()) + list(labels.values())
    update(0)
    fig.tight_layout()
    ani = FuncAnimation(fig, update, frames=len(element_values_steps), interval=interval, blit=False)
    if save_path:
        ext = os.path.splitext(save_path)[1].lower()
        if ext in ('.mp4', '.mov', '.m4v'):
            ani.save(save_path, writer=FFMpegWriter(fps=fps))
        elif ext == '.gif':
            ani.save(save_path, writer=PillowWriter(fps=fps))
        else:
            raise ValueError("save_path must end with .mp4, .mov, .m4v, or .gif.")
        print(f"Animation saved to: {save_path}")
    if show:
        plt.show()
    else:
        plt.close(fig)
    return ani

