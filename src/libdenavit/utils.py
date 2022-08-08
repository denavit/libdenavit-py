import math


def find_limit_point_in_list(the_list, val):

    # Quick return if the first entry of the list is the target value
    if the_list[0] == val:
        return 0, 0.0

    # Create a list with True where the corresponding value in the_list 
    # is on the other side of the target value and False otherwise 
    if the_list[0] > val:
        tf_list = [x<=val for x in the_list]
    else:
        tf_list = [x>=val for x in the_list]

    # Find the first instance of True
    try:
        ind = tf_list.index(True) - 1
    except ValueError:
        # If list does not pass target value, return None
        return None,None
    except Exception as e:
        raise e
    
    # Determine where target value is between list items
    if val == the_list[ind]:
        x = 0.0
    else:
        x = (val - the_list[ind]) / (the_list[ind+1] - the_list[ind])

    return ind, x


def interpolate_list(the_list, ind, x=0):

    if (ind is None) or (x is None):
        raise ValueError(f'One or both of ind and x are None ({ind = }) ({x = })')

    if x == 0:
        return the_list[ind]
        
    return the_list[ind] + (the_list[ind+1] - the_list[ind]) * x


def find_intersection_between_two_lines(Ax,Ay,Bx,By,Cx,Cy,Dx,Dy):
    # This function finds the intersection of lines AB and CD. The lines are defined by two points but the intesection
    # may occur and outside of the two points. If the lines are parallel, Ix and Iy are retured as empty list.

    if (Ax == Bx) and (Cx == Dx):
        Ix = []
        Iy = []

    elif (Ax == Bx):
        CDm = (Dy - Cy) / (Dx - Cx)
        Ix = Ax
        Iy = CDm * (Ix - Cx) + Cy
    elif (Cx == Dx):
        ABm = (By - Ay) / (Bx - Ax)
        Ix = Cx
        Iy = ABm * (Ix - Ax) + Ay
    else:
        ABm = (By - Ay) / (Bx - Ax)
        CDm = (Dy - Cy) / (Dx - Cx)
        if (ABm == CDm):
            Ix = []
            Iy = []
        else:
            Ix = (ABm * Ax - CDm * Cx - Ay + Cy) / (ABm - CDm)
            Iy = ABm * (Ix - Ax) + Ay

    return Ix, Iy


def area_of_circular_segment(radius, height):
    # radius is the radius of the circle
    # height is the distance from the center of the circle to the lower boundary of the segment
    # https://mathworld.wolfram.com/CircularSegment.html

    if radius <= 0:
        raise ValueError(f'radius must be greater than zero, {radius = }')
    if height < 0:
        raise ValueError(f'height must be greater than or equal to zero, {height = }')
    if height > radius:
        raise ValueError(f'height must be less than or equal to radius, {radius = }, {height = }')
    if height == radius:
        return 0.

    theta = 2 * math.acos(height / radius)
    area = 0.5 * radius ** 2 * (theta - math.sin(theta))
    return area


def centroid_of_circular_segment(radius, height):
    # radius is the radius of the circle
    # height is the distance from the center of the circle to the lower boundary of the segment
    # https://mathworld.wolfram.com/CircularSegment.html

    if radius <= 0:
        raise ValueError(f'radius must be greater than zero, {radius = }')
    if height < 0:
        raise ValueError(f'height must be greater than or equal to zero, {height = }')
    if height > radius:
        raise ValueError(f'height must be less than or equal to  radius, {radius = }, {height = }')
    if height == radius:
        return radius

    theta = 2 * math.acos(height / radius)
    centroid = (4 * radius * (math.sin(0.5 * theta)) ** 3) / (3 * (theta - math.sin(theta)))

    return centroid