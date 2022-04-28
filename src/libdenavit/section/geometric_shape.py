from math import sqrt


class GeometricShape:

    def I(self, axis):
        if axis.lower() == 'x':
            return self.Ix
        elif axis.lower() == 'y':
            return self.Iy
        else:
            raise Exception(f'Unknown axis: {axis}')

    def S(self, axis):
        if axis.lower() == 'x':
            return self.Sx
        elif axis.lower() == 'y':
            return self.Sy
        else:
            raise Exception(f'Unknown axis: {axis}')

    def Z(self, axis):
        if axis.lower() == 'x':
            return self.Zx
        elif axis.lower() == 'y':
            return self.Zy
        else:
            raise Exception(f'Unknown axis: {axis}')

    def r(self, axis):
        if axis.lower() == 'x':
            return self.rx
        elif axis.lower() == 'y':
            return self.ry
        else:
            raise Exception(f'Unknown axis: {axis}')

    @property
    def rx(self):
        return sqrt(self.Ix / self.A)

    @property
    def ry(self):
        return sqrt(self.Iy / self.A)
