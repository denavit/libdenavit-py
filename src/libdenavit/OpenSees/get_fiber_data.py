import openseespy.opensees as ops
import json
import os


def get_fiber_data(section_tag: str, filename=None, keep_json=False, plot_fibers=False):
    
    if filename is None:
        from datetime import datetime
        now = datetime.now()
        filename = 'data_fiber_' + now.strftime('%Y%m%d_%H%M%S%f') + '.json'
    
    if os.path.exists(filename):
        ans = input('The file already exists. Enter "Y" to overwrite it or "N" to cancel: ')
        if ans == 'N':
            raise Exception('Cancelled by user to prevent overwriting the json file.')

    ops.printModel('-JSON', '-file', f'{filename}')
    output_data = json.load(open(str(filename)))
    sections_list = output_data['StructuralAnalysisModel']['properties']['sections']
    x = []
    y = []
    A = []
    m = []

    for ind, value in enumerate(sections_list):
        if value["name"] == str(section_tag):
            section = value

    for j in section['fibers']:
        y.append(j['coord'][0])
        x.append(j['coord'][1])
        A.append(j['area'])
        m.append(j['material'])

    if not keep_json:
        os.remove(f'{filename}')

    if plot_fibers:
        import matplotlib.pyplot as plt
        m_int = list(map(int,m))
        plt.scatter(x,y,A,m_int)
        plt.gca().axis('equal')
        plt.show()

    return x, y, A, m
