import openseespy.opensees as ops
from libdenavit import area_of_circular_segment, centroid_of_circular_segment


def circ_patch_2d(mat_tag, num_fibers_half_circle, D, Di=0, yc=0):
    # Height of each fiber (strip)
    dy = D / (2 * num_fibers_half_circle)

    # Loop over half circle starting from center
    for i in range(num_fibers_half_circle):
        y1 = D/2 * i / num_fibers_half_circle

        A1 = area_of_circular_segment(D / 2, y1)
        yc1 = centroid_of_circular_segment(D / 2, y1)
        A2 = area_of_circular_segment(D / 2, y1 + dy)
        yc2 = centroid_of_circular_segment(D / 2, y1 + dy)
        if y1 >= Di/2:
            iA = A1 - A2
            iy = (yc1 * A1 - yc2 * A2) / iA
        elif (y1 + dy) >= Di/2:
            A3 = area_of_circular_segment(Di / 2, y1)
            yc3 = centroid_of_circular_segment(Di / 2, y1)
            iA = A1 - A2 - A3
            iy = (yc1 * A1 - yc2 * A2 - yc3 * A3) / iA 
        else:
            A3 = area_of_circular_segment(Di / 2, y1)
            yc3 = centroid_of_circular_segment(Di / 2, y1)
            A4 = area_of_circular_segment(Di / 2, y1 + dy)
            yc4 = centroid_of_circular_segment(Di / 2, y1 + dy)
            iA = A1 - A2 - A3 + A4
            iy = (yc1 * A1 - yc2 * A2 - yc3 * A3 + yc4 * A4) / iA 

        ops.fiber(yc + iy, 0, iA, mat_tag)  # Fiber on top half of circle
        ops.fiber(yc - iy, 0, iA, mat_tag)  # Fiber on bottom half of circle

    return


def run_example():
    D = 10
    Di = 5
    
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)

    ops.node(0, 0, 0, 0)
    ops.node(1, 0, 0, 0)

    ops.fix(1, 1, 0, 1, 1, 0, 0)
    ops.uniaxialMaterial('Elastic', 1, 1000)
    ops.uniaxialMaterial('Elastic', 2, 2000)

    ops.section("Fiber", 1, '-GJ', 10e6)
    circ_patch_2d(1, 100, D, Di = Di)
    ops.element('zeroLengthSection', 1, 0, 1, 1)

    from libdenavit.OpenSees import get_fiber_data
    x, y, A, m = get_fiber_data("1")
    print(f'x: {x}\n', f'y: {y}\n', f'A: {A}\n', f'm: {m}\n')

    print(sum(A))
    import math
    print(math.pi/4*(D*D-Di*Di))

    import numpy as np
    y = np.array(y)
    A = np.array(A)
    I = np.sum(A*y*y)
    print(I)
    print(math.pi/64*(D**4-Di**4))    

if __name__ == "__main__":
    run_example()

