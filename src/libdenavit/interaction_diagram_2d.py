import numpy as np
import math
from libdenavit import find_limit_point_in_list, interpolate_list, find_intersection_between_two_lines
import operator
from shapely.geometry import LineString, Polygon
import matplotlib.pyplot as plt

def cart2pol(x, y):
    x = np.asarray(x)
    y = np.asarray(y)

    rho = np.hypot(x, y)  # Computes sqrt(x**2 + y**2) efficiently
    phi = np.arctan2(y, x)  # Computes arctan2(y, x) element-wise

    return rho, phi


class InteractionDiagram2d():
    def __init__(self, idx: list, idy: list, is_closed=False):
        idx = idx.tolist() if isinstance(idx, np.ndarray) else idx
        idy = idy.tolist() if isinstance(idy, np.ndarray) else idy

        _, q = cart2pol(idx, idy)
        q = [i%(2*np.pi) for i in q]
        ind = np.argsort(q)
        q = np.sort(q)
        idx = [idx[i] for i in ind]
        idy = [idy[i] for i in ind]

        self.idx = idx
        self.idy = idy
        self.q = q
        self.is_closed = is_closed


    def radial_distance(self, angles, degrees=False):
        if degrees:
            angles = np.deg2rad(angles)

        if type(angles) in [float, int, np.float64]:
            angles = [angles]

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

        if len(d) == 1:
            return d[0]

        return d


    def compare_two(self, test_id, angles, degrees=False):
        d_base = self.radial_distance(angles, degrees)
        d_test = test_id.radial_distance(angles, degrees)
        d_base = np.array(d_base)
        d_test = np.array(d_test)
        errors = [(i-j) / i for i, j in zip(d_base, d_test)]
        return errors


    def check_points(self, pointsX, pointsY):
        d_pts, q_pts = cart2pol(pointsX, pointsY)
        d = self.radial_distance(q_pts)
        max_ratio = np.max(np.divide(d_pts, d))
        return max_ratio


    def find_intersection(self, pathX, pathY):
        if self.is_closed:
            polin = Polygon(np.column_stack((self.idx, self.idy)))
            line_1 = LineString(list(polin.exterior.coords))
        else:
            line_1 = LineString(np.column_stack((self.idx, self.idy)))
        
        line_2 = LineString(np.column_stack((pathX, pathY)))
        
        if len(pathX) == 2:  # If it's just a simple line segment
            # Calculate direction and extend the line
            dx = pathX[-1] - pathX[0]
            dy = pathY[-1] - pathY[0]
            length = np.sqrt(dx**2 + dy**2)
            
            if length > 0:
                # Normalize direction
                dx_norm = dx / length
                dy_norm = dy / length
                
                # Extend the line by a large factor
                max_boundary_extent = max(np.max(np.abs(self.idx)), np.max(np.abs(self.idy))) * 2
                
                extended_x = pathX[0] + dx_norm * max_boundary_extent
                extended_y = pathY[0] + dy_norm * max_boundary_extent
                
                line_2 = LineString([(pathX[0], pathY[0]), (extended_x, extended_y)])
        
        intersection = line_1.intersection(line_2)
        
        # Handle empty intersection
        if intersection.is_empty:
            raise Exception('No intersection found between the two paths')
        
        # Handle Point
        if intersection.geom_type == 'Point':
            return intersection.x, intersection.y
        
        # Handle MultiPoint
        elif intersection.geom_type == 'MultiPoint':
            x = [pt.x for pt in intersection.geoms]
            y = [pt.y for pt in intersection.geoms]
            return x, y
        
        # Handle LineString (overlapping segments)
        elif intersection.geom_type == 'LineString':
            coords = list(intersection.coords)
            x = [coord[0] for coord in coords]
            y = [coord[1] for coord in coords]
            return x, y
        
        # Handle MultiLineString or GeometryCollection
        elif intersection.geom_type in ['MultiLineString', 'GeometryCollection']:
            x = []
            y = []
            for geom in intersection.geoms:
                if geom.geom_type == 'Point':
                    x.append(geom.x)
                    y.append(geom.y)
                elif geom.geom_type == 'LineString':
                    coords = list(geom.coords)
                    x.extend([coord[0] for coord in coords])
                    y.extend([coord[1] for coord in coords])
            return x, y
        
        else:
            raise Exception(f'Unexpected intersection type: {intersection.geom_type}')


    def find_x_given_y(self, Y, signX):
        X_list = []
        if type(Y) in [int, float, np.float64]:
            Y_list = [Y]
        for y in Y_list:
            if signX.lower() in ['+', 'positive', 'pos']:
                peakX = 1.1 * np.max(self.idx)
            elif signX.lower() in ['-', 'negative', 'neg']:
                peakX = 1.1 * np.min(self.idx)
            else:
                raise ValueError('signX must be positive or negative')

            npts = 10
            pathX = np.linspace(0, peakX, npts)
            pathY = y * np.ones(npts)
            X_list.append(self.find_intersection(pathX, pathY)[0])

        if type(Y) in [int, float, np.float64]:
            return X_list[0]
        else:
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

    def is_outside(self, x, y):
        if x <= 0 or y <= 0:
            raise ValueError("Coordinates must be positive (x > 0, y > 0)")

        # compute polar coords of the point
        r_p, theta = cart2pol(x, y)

        # find boundary radius at angle theta
        r_b = self.radial_distance(theta)

        # 4) compare
        return r_p > r_b

    def plot(self, *args, **kwargs):
        if len(args) == 0 and len(kwargs) == 0:
            args = ['-o']

        if self.is_closed:
            plt.plot(self.idx+[self.idx[0]], self.idy+[self.idy[0]], *args, **kwargs)
        else:
            plt.plot(self.idx, self.idy, *args, **kwargs)


if __name__ == "__main__":

    # Interaction Diagram 1
    a1 = [0.9, 1.1, 0.8, 0.0]
    b1 = [0.0, 0.4, 0.9, 1.0]
    c1 = InteractionDiagram2d(a1, b1)

    # Interaction Diagram 2
    a2 = [0.8, 1.0, 0.8, 0.0]
    b2 = [0.0, 0.4, 1.0, 1.1]
    c2 = InteractionDiagram2d(a2, b2)

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

    c1.plot('-r', label='Interaction Diagram 1')
    plt.plot(distance*np.cos(angles), distance*np.sin(angles), 'bo', label='Angle Points')

    c2.plot('-g', label='Interaction Diagram 2')
    plt.plot(*c1.find_intersection(a2, b2), "k*", label='Intersection Point')

    plt.legend()

    plt.show()

