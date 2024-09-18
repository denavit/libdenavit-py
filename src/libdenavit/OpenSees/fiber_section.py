import openseespy.opensees as ops
from libdenavit import area_of_circular_segment, centroid_of_circular_segment
from math import ceil, tan, acos


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


def obround_patch_2d(mat_tag, num_fibers, D, a, axis, yc=0):
    if axis == 'y':
       
        # Number of fibers in each region
        nfa = ceil(a / (D + a) * num_fibers)
        nfd_half = ceil((D/2) / (D + a) * num_fibers)
        
        # Height of each fiber (strip) in each region
        dya  = a / nfa
        dyd = D / (2 * nfd_half)

        iAa = dya * D
        for i in range(nfa):
            iya = a/2 - (i*2 + 1)/2 * dya
            ops.fiber(yc+iya, 0, iAa, mat_tag)

        for i in range(nfd_half):
            y1 = D/2 * i/nfd_half
            A1 = area_of_circular_segment(D / 2, y1)
            yc1 = centroid_of_circular_segment(D / 2, y1)
            A2 = area_of_circular_segment(D / 2, y1 + dyd)
            yc2 = centroid_of_circular_segment(D / 2, y1 + dyd)
            iA = A1 - A2
            iy = (yc1 * A1 - yc2 * A2) / iA
            ops.fiber(yc+iy+a/2, 0, iA, mat_tag)
            ops.fiber(yc-iy-a/2, 0, iA, mat_tag)

    elif axis == 'x':
    
        # Number of fibers in each region
        num_fib_half = ceil(num_fibers/2)
        
        # Height of each fiber (strip) in each region
        dy = D / (2*num_fib_half)

        # Loop over half circle starting from center
        for i in range(num_fib_half):
            y1 = D/2 * i / num_fib_half

            A1 = area_of_circular_segment(D / 2, y1)
            yc1 = centroid_of_circular_segment(D / 2, y1)
            A2 = area_of_circular_segment(D / 2, y1 + dy)
            yc2 = centroid_of_circular_segment(D / 2, y1 + dy)
            iA = A1 - A2
            iy = (yc1 * A1 - yc2 * A2) / iA

            iA2 = a * dy
            ops.fiber(yc + iy, 0, iA+iA2, mat_tag)
            ops.fiber(yc - iy, 0, iA+iA2, mat_tag)
    
    else:
        raise ValueError(f'obround_patch_2d not implemented for {axis = })')


def obround_patch_2d_confined(mat_tag_cover, mat_tag_core, num_fibers, D, a, Dc, axis, yc=0):
    if axis == 'x':
        
        # Number of fibers in each region
        nf_mid  = ceil(Dc/2 * num_fibers / (D+a))
        nf_center  = ceil(a/2 * num_fibers / (D+a))
        nf_top = ceil((D - Dc)/2 * num_fibers / (D+a))
        
        # Height of each fiber (strip) in each region
        dy_mid  = (Dc/2) / nf_mid
        dy_center  = (a/2) / nf_center
        dy_top = ((D - Dc)/2) / nf_top

        for i in range(nf_center):
            y1 = a / 2 * i / nf_center
            A1 = area_of_circular_segment(Dc / 2, y1)
            yc1 = centroid_of_circular_segment(Dc / 2, y1)
            A2 = area_of_circular_segment(Dc / 2, y1 + dy_center)
            yc2 = centroid_of_circular_segment(Dc / 2, y1 + dy_center)
            iA_conf = A1 - A2
            iy_conf = -(yc1 * A1 - yc2 * A2) / iA_conf + a/2

            A3 = D * dy_center
            iA_unconf = A3 - iA_conf
            iy_unconf = -(((y1 + dy_center/2) - a/2) * A3 + iy_conf * iA_conf) / iA_unconf

            ops.fiber(yc + iy_conf, 0, iA_conf, mat_tag_core)
            ops.fiber(yc - iy_conf, 0, iA_conf, mat_tag_core)
            ops.fiber(yc + iy_unconf, 0, iA_unconf, mat_tag_core)
            ops.fiber(yc - iy_unconf, 0, iA_unconf, mat_tag_cover)

        for i in range(nf_mid):
            y1 = Dc / 2 * i / nf_mid
            A1 = area_of_circular_segment(Dc / 2, y1)
            yc1 = centroid_of_circular_segment(Dc / 2, y1)
            A2 = area_of_circular_segment(Dc / 2, y1 + dy_mid)
            yc2 = centroid_of_circular_segment(Dc / 2, y1 + dy_mid)
            iA_conf = A1 - A2
            iy_conf = (yc1 * A1 - yc2 * A2) / iA_conf

            A3 = area_of_circular_segment(D / 2, y1)
            yc3 = centroid_of_circular_segment(D / 2, y1)
            A4 = area_of_circular_segment(D / 2, y1 + dy_mid)
            yc4 = centroid_of_circular_segment(D / 2, y1 + dy_mid)
            iA_unconf = A3 - A4 - iA_conf
            iy_unconf = (yc3 * A3 - yc4 * A4 - iy_conf * iA_conf) / iA_unconf

            ops.fiber(yc + iy_conf + a/2, 0, iA_conf, mat_tag_core)
            ops.fiber(yc - iy_conf - a/2, 0, iA_conf, mat_tag_core)
            ops.fiber(yc + iy_unconf + a/2, 0, iA_unconf, mat_tag_core)
            ops.fiber(yc - iy_unconf - a/2, 0, iA_unconf, mat_tag_core)

        for i in range(nf_top):
            y1 = (D/2 - Dc/2) * i / nf_top
            A1 = area_of_circular_segment(D/2, y1 + Dc/2)
            yc1 = centroid_of_circular_segment(D/2, y1 + Dc/2)
            A2 = area_of_circular_segment(D/2, y1 + Dc/2 + dy_top)
            yc2 = centroid_of_circular_segment(D/2, y1 + Dc/2 + dy_top)
            iA = A1 - A2
            iy = (yc1 * A1 - yc2 * A2) / iA

            ops.fiber(yc + iy + a/2, 0, iA, mat_tag_cover)
            ops.fiber(yc - iy - a/2, 0, iA, mat_tag_cover)

    elif axis == 'y':
        
        # Find the height of the rectangle inside circle
        y_center = a / 2 * tan(acos((a / 2) / (Dc / 2)))
        nf = ceil(Dc / 2 * num_fibers / D)

        # Number of fibers in each region
        nf_center = ceil(y_center / (Dc / 2) * nf)
        nf_mid = ceil((Dc / 2 - y_center) / (Dc / 2) * nf)
        nf_top = ceil((D / 2 - Dc / 2) / (D / 2) * nf)

        # Height of each fiber (strip) in each region
        dy_center = y_center / nf_center
        dy_mid = (Dc / 2 - y_center) / nf_mid
        dy_top = (D / 2 - Dc / 2) / nf_top

        for i in range(nf_center):
            y1 = y_center * i / nf_center
            A1 = area_of_circular_segment(Dc / 2, y1)
            yc1 = centroid_of_circular_segment(Dc / 2, y1)
            A2 = area_of_circular_segment(Dc / 2, y1 + dy_center)
            yc2 = centroid_of_circular_segment(Dc / 2, y1 + dy_center)

            A3 = A1 - A2
            yc3 = (yc1 * A1 - yc2 * A2) / A3

            A4 = dy_center * a
            yc4 = y1 + dy_center / 2

            iA_1 = A3 + A4
            iy_1 = (yc3 * A3 + yc4 * A4) / iA_1

            A5 = area_of_circular_segment(D / 2, y1)
            yc5 = centroid_of_circular_segment(D / 2, y1)
            A6 = area_of_circular_segment(D / 2, y1 + dy_center)
            yc6 = centroid_of_circular_segment(D / 2, y1 + dy_center)

            iA_2 = A5 - A6 - A3
            iy_2 = (yc5 * A5 - yc6 * A6 - yc3 * A3) / iA_2

            ops.fiber(yc + iy_1, 0, iA_1, mat_tag_core)
            ops.fiber(yc - iy_1, 0, iA_1, mat_tag_core)
            ops.fiber(yc + iy_2, 0, iA_2, mat_tag_cover)
            ops.fiber(yc - iy_2, 0, iA_2, mat_tag_cover)

        for i in range(nf_mid):
            y1 = y_center + (Dc / 2 - y_center) * i / nf_mid
            A1 = area_of_circular_segment(Dc / 2, y1)
            yc1 = centroid_of_circular_segment(Dc / 2, y1)
            A2 = area_of_circular_segment(Dc / 2, y1 + dy_mid)
            yc2 = centroid_of_circular_segment(Dc / 2, y1 + dy_mid)
            iA_1 = (A1 - A2) * 2
            iy_1 = (yc1 * A1 * 2 - yc2 * A2 * 2) / iA_1

            A3 = area_of_circular_segment(D / 2, y1)
            yc3 = centroid_of_circular_segment(D / 2, y1)
            A4 = area_of_circular_segment(D / 2, y1 + dy_mid)
            yc4 = centroid_of_circular_segment(D / 2, y1 + dy_mid)

            A5 = a * dy_mid
            yc5 = y1 + dy_mid / 2

            iA_2 = A5 + A3 - A4 - iA_1
            iy_2 = (yc5 * A5 + yc3 * A3 - yc4 * A4 - iy_1 * iA_1) / iA_2

            ops.fiber(yc + iy_1, 0, iA_1, mat_tag_core)
            ops.fiber(yc - iy_1, 0, iA_1, mat_tag_core)
            ops.fiber(yc + iy_2, 0, iA_2, mat_tag_cover)
            ops.fiber(yc - iy_2, 0, iA_2, mat_tag_cover)

        for i in range(nf_top):
            y1 = Dc/2 + (D - Dc)/2 * i / nf_top
            A1 = area_of_circular_segment(D / 2, y1)
            yc1 = centroid_of_circular_segment(D / 2, y1)
            A2 = area_of_circular_segment(D / 2, y1 + dy_top)
            yc2 = centroid_of_circular_segment(D / 2, y1 + dy_top)

            A3 = a * dy_top
            yc3 = y1 + dy_top / 2

            iA = A3 + A1 - A2
            iy = (yc3 * A3 + yc1 * A1 - yc2 * A2) / iA

            ops.fiber(yc + iy, 0, iA, mat_tag_cover)
            ops.fiber(yc - iy, 0, iA, mat_tag_cover)

    else:
        raise ValueError(f'obround_patch_2d_confined not implemented for {axis = })')


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
    circ_patch_2d(1, 100, D, Di=Di)
    ops.element('zeroLengthSection', 1, 0, 1, 1)

    from libdenavit.OpenSees import get_fiber_data
    x, y, A, m = get_fiber_data("1")
    print(f'x: {x}\n', f'y: {y}\n', f'A: {A}\n', f'm: {m}\n')

    print(sum(A))
    import math
    print(math.pi / 4 * (D * D - Di * Di))

    import numpy as np
    y = np.array(y)
    A = np.array(A)
    I = np.sum(A * y * y)
    print(I)
    print(math.pi / 64 * (D ** 4 - Di ** 4))


def run_example_2(num_fibers):
    D = 4
    Dc = 3.5
    a = 3
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)

    ops.node(0, 0, 0, 0)
    ops.node(1, 0, 0, 0)

    ops.fix(1, 1, 0, 1, 1, 0, 0)
    ops.uniaxialMaterial('Elastic', 1, 1000)
    ops.uniaxialMaterial('Elastic', 2, 2000)

    ops.section("Fiber", 1, '-GJ', 10e6)
    obround_patch_2d_confined(1, 2, num_fibers, D, a, Dc, "x")
    ops.element('zeroLengthSection', 1, 0, 1, 1)

    from libdenavit.OpenSees import get_fiber_data
    x, y, A, m = get_fiber_data("1")
    #print(f'x: {x}\n', f'y: {y}\n', f'A: {A}\n', f'm: {m}\n')
    print(f"A  = {sum(A) :.6f} vs 24.566 from CAD")

    import numpy as np
    y = np.array(y)
    A = np.array(A)
    Ix = np.sum(A * y * y)
    print(f"{Ix = :.6f} vs 81.8407 from CAD")


def run_example_3(num_fibers):
    D = 4
    Dc = 3.5
    a = 3
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)

    ops.node(0, 0, 0, 0)
    ops.node(1, 0, 0, 0)

    ops.fix(1, 1, 0, 1, 1, 0, 0)
    ops.uniaxialMaterial('Elastic', 1, 1000)
    ops.uniaxialMaterial('Elastic', 2, 2000)

    ops.section("Fiber", 1, '-GJ', 10e6)
    obround_patch_2d_confined(1, 2, num_fibers, D, a, Dc, "y")
    ops.element('zeroLengthSection', 1, 0, 1, 1)

    from libdenavit.OpenSees import get_fiber_data
    x, y, A, m = get_fiber_data("1")
    #print(f'x: {x}\n', f'y: {y}\n', f'A: {A}\n', f'm: {m}\n')
    print(f"A  = {sum(A) :.6f} vs 24.566 from CAD")

    import numpy as np
    y = np.array(y)
    A = np.array(A)
    Iy = np.sum(A * y * y)
    print(f"{Iy = :.6f} vs 28.5700 from CAD")


if __name__ == "__main__":
    # print("\nExample 1: ")
    # run_example()
    print("\nExample 2: ")
    run_example_2(50)
    print("\nExample 3: ")
    run_example_3(50)
