import openseespy.opensees as ops
import matplotlib.pyplot as plt

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
    return