import openseespy.opensees as ops
from libdenavit import area_of_circular_segment, centroid_of_circular_segment
from libdenavit.OpenSees import get_fiber_data


def circ_patch_2d(mat_tag, num_fibers_half_circle, D, yc=0):
    # Height of each fiber (strip)
    dy = D / (2 * num_fibers_half_circle)

    # Loop over half circle starting from center
    for i in range(num_fibers_half_circle):
        y1 = D/2 * i / num_fibers_half_circle

        A1 = area_of_circular_segment(D / 2, y1)
        yc1 = centroid_of_circular_segment(D / 2, y1)
        A2 = area_of_circular_segment(D / 2, y1 + dy)
        yc2 = centroid_of_circular_segment(D / 2, y1 + dy)
        iA = A2 - A1
        iy = (yc2 * A2 - yc1 * A1) / iA

        ops.fiber(yc + iy, 0, iA, mat_tag)  # Fiber on top half of circle
        ops.fiber(yc - iy, 0, iA, mat_tag)  # Fiber on bottom half of circle

    return


if __name__ == "__main__":
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)

    ops.node(0, 0, 0, 0)
    ops.node(1, 0, 0, 0)

    ops.fix(1, 1, 0, 1, 1, 0, 0)
    ops.uniaxialMaterial('Elastic', 1, 1000)
    ops.uniaxialMaterial('Elastic', 2, 2000)

    ops.section("Fiber", 1, '-GJ', 10e6)
    circ_patch_2d(1, 10, 10)
    ops.element('zeroLengthSection', 1, 0, 1, 1)

    x, y, A, m = get_fiber_data("1", keep_json=True)
    print(f'x: {x}\n', f'y: {y}\n', f'A: {A}\n', f'm: {m}\n')

    print(sum(A))

