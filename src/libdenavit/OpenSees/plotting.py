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
    return

def plot_element_values_2d(element_values, element_value_name='Element Value', axis_equal=False, show=True):
    # Plots a scalar value associated with each element in a 2D OpenSees model.
    # element_values must be a dictionary or an iterable of (element_tag, value) pairs.
    # Each element tag must correspond to an element tag in the OpenSees model.
    # Each element must have a single numerical value (e.g., P-M ratio, axial force, or utilization ratio).
    # Example 1: element_values = {1: 0.25, 2: 0.50, 3: 0.85}
    # Example 2: element_values = [(1, 0.25), (2, 0.50), (3, 0.85)]
    # Elements without a supplied value are plotted in gray.
    # element_value_name specifies the plot title and colorbar label.
    element_values = dict(element_values)
    node_coords = get_node_coords()
    element_nodes = get_element_nodes()
    values = [value for value in element_values.values() if value is not None]
    if not values:
        raise ValueError("No valid element values were provided.")
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
        plt.plot(xplt, yplt, 'o-', color=color)
        if value is not None:
            xm = (coordi[0] + coordj[0]) / 2
            ym = (coordi[1] + coordj[1]) / 2
            plt.text(xm, ym + 0.02, f"{value:.2f}")
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    plt.colorbar(sm, ax=plt.gca()).set_label(element_value_name)
    if axis_equal:
        plt.axis('equal')
    plt.xlabel("X-coordinate")
    plt.ylabel("Y-coordinate")
    plt.title(f"{element_value_name} Colormap")
    plt.tight_layout()
    if show:
        plt.show()
    return

def plot_internal_force_diagram_2d(force_type='M', scale=None, axis_equal=False, annotate=False, show=True, save_path=None, show_values=True):
    force_options = {
        'P': {'indices': (0, 3), 'signs': (-1, 1), 'color': 'green', 'fill_color': 'lightgreen', 'scale': 0.05, 'name': 'Axial Force'},
        'V': {'indices': (1, 4), 'signs': (1, -1), 'color': 'blue', 'fill_color': 'lightblue', 'scale': 0.002, 'name': 'Shear Force'},
        'M': {'indices': (2, 5), 'signs': (-1, 1), 'color': 'red', 'fill_color': 'lightcoral', 'scale': 0.001, 'name': 'Bending Moment'}
    }
    force_type = force_type.upper()
    if force_type not in force_options:
        raise ValueError("force_type must be 'P', 'V', or 'M'.")
    options = force_options[force_type]
    if scale is None:
        scale = options['scale']
    node_coords = get_node_coords()
    element_nodes = get_element_nodes()
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
        if force_type=='M':
            diagram_i = coordi - scale * value_i * n_perp
            diagram_j = coordj - scale * value_j * n_perp
        plt.fill([coordi[0], coordj[0], diagram_j[0], diagram_i[0]], [coordi[1], coordj[1], diagram_j[1], diagram_i[1]], color=options['fill_color'], alpha=0.15)
        plt.plot([coordi[0], coordj[0]], [coordi[1], coordj[1]], 'k-', linewidth=2.5)
        plt.plot([coordi[0], diagram_i[0]], [coordi[1], diagram_i[1]], color=options['color'], linestyle='--')
        plt.plot([coordj[0], diagram_j[0]], [coordj[1], diagram_j[1]], color=options['color'], linestyle='--')
        plt.plot([diagram_i[0], diagram_j[0]], [diagram_i[1], diagram_j[1]], color=options['color'])
        if show_values:
            x_offset = 0.5
            y_offset = 0.5
            plt.text(diagram_i[0] + x_offset, diagram_i[1] + y_offset, f"{value_i:.2f}", color=options['color'])
            plt.text(diagram_j[0] + x_offset, diagram_j[1] + y_offset, f"{value_j:.2f}", color=options['color'])
    if axis_equal:
        plt.axis('equal')
    plt.xlabel("X-coordinate")
    plt.ylabel("Y-coordinate")
    plt.title(f"Undeformed 2D Frame – {options['name']} Diagram")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    return

def animate_element_values_2d(element_values_steps, element_value_name='Element Value', interval=10, axis_equal=False, show=True, save_path=None, fps=10):
    # Animates scalar values associated with elements in a 2D OpenSees model.
    # Each step can be a dictionary {element_tag: value} or a list of (element_tag, value) pairs.
    # Example: element_values_steps = [{1: 0.25, 2: 0.50}, {1: 0.40, 2: 0.65}, {1: 0.60, 2: 0.85}]
    # Element tags must match OpenSees element tags, and each element has one numerical value per step.
    # Missing element values are displayed in gray.
    # The color scale remains constant across all analysis steps.
    # interval controls frame delay in milliseconds; fps controls saved video frame rate.
    # save_path supports MP4, MOV, M4V, and GIF.
    node_coords = get_node_coords()
    element_nodes = get_element_nodes()
    if not element_values_steps:
        raise ValueError("No element values were provided.")
    all_values = [value for step in element_values_steps for value in dict(step).values() if value is not None]
    if not all_values:
        raise ValueError("No element values were found.")
    norm = Normalize(vmin=min(all_values), vmax=max(all_values))
    cmap = colormaps.get_cmap('coolwarm')
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
        line, = plt.plot(xplt, yplt, 'o-', color=color)
        lines[i] = line
        xm = (coordi[0] + coordj[0]) / 2
        ym = (coordi[1] + coordj[1]) / 2
        labels[i] = plt.text(xm, ym + 0.02, f"{value:.2f}" if value is not None else "")
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    plt.colorbar(sm, ax=plt.gca()).set_label(element_value_name)
    if axis_equal:
        plt.axis('equal')
    plt.xlabel("X-coordinate")
    plt.ylabel("Y-coordinate")
    def update(frame_index):
        values = dict(element_values_steps[frame_index])
        for i in lines:
            value = values.get(i)
            lines[i].set_color(cmap(norm(value)) if value is not None else 'gray')
            labels[i].set_text(f"{value:.2f}" if value is not None else "")
        plt.title(f"{element_value_name} (Step {frame_index + 1})")
        return list(lines.values()) + list(labels.values())
    update(0)
    plt.tight_layout()
    ani = FuncAnimation(plt.gcf(), update, frames=len(element_values_steps), interval=interval, blit=False)
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
        plt.close(ani._fig)
    return ani

