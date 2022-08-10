import numpy as np
import math
from libdenavit import find_limit_point_in_list, interpolate_list, find_intersection_between_two_lines
import operator
from shapely.geometry import LineString

def cart2pol(x, y):
    assert type(x) == type(y), "x and y must be of the same type"
    if type(x) in [int, float]:
        x = [x]
        y = [y]

    assert len(x) == len(y), "x and y must be of the same length"

    rho = []
    phi = []
    for i, _ in enumerate(x):
        rho.append(np.sqrt(x[i]**2 + y[i]**2))
        phi.append(np.arctan2(y[i], x[i]))

    return rho, phi


class InteractionDiagram2d():

    def __init__(self, idx: list, idy: list):
        _, q = cart2pol(idx, idy)
        q = [i%(2*np.pi) for i in q]
        ind = np.argsort(q)
        q = np.sort(q)
        idx = [idx[i] for i in ind]
        idy = [idy[i] for i in ind]

        self.idx = idx
        self.idy = idy
        self.q = q


    def radial_distance(self, angles, degrees=False):
        if degrees:
            angles = np.deg2rad(angles)

        angles = [i%(2*np.pi) for i in angles]
        d = [None] * len(angles)

        for i in range(len(angles)):
            ind = None

            for j, k in enumerate(self.q):
                if k >= angles[i]:
                    ind = j
                    break

            if ind is None:
                d[i] = None

            elif ind == 0:
                if self.q[0] == angles[i]:
                    d[i] = math.hypot(self.idx[0], self.idy[0])
                else:
                    d[i] = 0
            else:
                Ax = self.idx[ind]
                Ay = self.idy[ind]
                Bx = self.idx[ind-1]
                By = self.idy[ind-1]
                Ix, Iy = find_intersection_between_two_lines(Ax, Ay, Bx, By, 0, 0, np.cos(angles[i]), np.sin(angles[i]))
                try:
                    if len(Ix) == 0:
                        raise ValueError('Bad intersection diagram')
                except:
                    d[i] = math.hypot(Ix, Iy)

        return d


    def compare_two(self, test_id, angles):
        d_base = self.radial_distance(angles)
        d_test = test_id.radial_distance(angles)

        errors = [(i-j) / i for i, j in zip(d_base, d_test)]
        return errors


    def check_points(self, pointsX, pointsY):
        d_pts, q_pts = cart2pol(pointsX, pointsY)
        d = self.radial_distance(q_pts)
        max_ratio = np.max(np.divide(d_pts, d))
        return max_ratio


    def find_intersection(self, pathX, pathY):

        line_1 = LineString(np.column_stack((self.idx, self.idy)))
        line_2 = LineString(np.column_stack((pathX, pathY)))
        intersection = line_1.intersection(line_2)

        # @todo get ind and x from as output
        return intersection.x, intersection.y,


    def find_x_given_y(self, Y, signX):
        if signX.lower() in ['+', 'positive', 'pos']:
            peakX = 1.1 * np.max(self.idx)
        elif signX.lower() in ['-', 'negative', 'neg']:
            peakX = 1.1 * np.min(self.idx)
        else:
            raise ValueError('signX must be positive or negative')

        npts = 10
        pathX = np.linspace(0, peakX, npts)
        pathY = Y * np.ones(npts)
        X, _ = self.find_intersection(pathX, pathY)
        return X


    def find_y_given_x(self, X, signY):
        if signY.lower() in ['+', 'positive', 'pos']:
            peakY = 1.1 * np.max(self.idy)
        elif signY.lower() in ['-', 'negative', 'neg']:
            peakY = 1.1 * np.min(self.idy)
        else:
            raise ValueError('signY must be positive or negative')

        npts = 1000
        pathX = X * np.ones(npts)
        pathY = np.linspace(0, peakY, npts)
        _, Y = self.find_intersection(pathX, pathY)
        return Y


    def plot(self, *args):
        plt.plot(self.idx, self.idy, *args)


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    a1 = [0.9, 1.0, 0.8, 0.0]
    b1 = [0.0, 0.4, 0.9, 1.0]
    c1 = InteractionDiagram2D(a1, b1)

    a2 = [0.8, 1.0, 0.8, 0.0]
    b2 = [0.0, 0.4, 1.0, 1.1]
    c2 = InteractionDiagram2D(a2, b2)

    angles = np.linspace(0, np.pi/2, 50)
    distance = c1.radial_distance(angles)

    print("Compare errors= \n",
          c1.compare_two(c2, angles))
    print("Check point (0.9, 0.0): \n",
          c1.check_points(0.9, 0.0))
    print("Find Intersection: \n",
          c1.find_intersection(a2, b2))
    print("Find X given Y: \n",
          c1.find_x_given_y(0.9, '+'))
    print("Find Y given X: \n",
          c1.find_y_given_x(0.9, '+'))

    plt.plot(distance*np.cos(angles), distance*np.sin(angles), 'bo', label='Angle Points')
    plt.plot(a1, b1, '-r', label='Interaction Diagram 1')
    plt.plot(a2, b2, '-g', label='Interaction Diagram 2')
    plt.plot(*c1.find_intersection(a2, b2), "ro")
    plt.legend()

    plt.show()

