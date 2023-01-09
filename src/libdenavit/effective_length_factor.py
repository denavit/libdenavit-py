from math import inf, sqrt, tan, pi
from scipy.optimize import fsolve


# @todo - write sidesway_inhibited_effective_length_factor


def sidesway_uninhibited_effective_length_factor(GA, GB):
    
    # Check basic cases for a quick return
    if GA == inf and GB == inf:
        return inf
    elif (GA == inf and GB == 0) or (GA == 0 and GB == inf):
        return 2
    elif GA == 0 and GB == 0:
        return 1

    # Change infinity to large value to avoid trouble in fsolve
    if GA == inf:
        GA = 1e10
    if GB == inf:
        GB = 1e10

    # Make a good initial guess of the effective length factor
    # Equation 5.14
    # Unified Design of Steel Structures, 3rd Ed.
    # Geschwindner, Liu, and Carter
    Kguess = sqrt((1.6 * GA * GB + 4 * (GA + GB) + 7.5) / (GA + GB + 7.5))

    # Define functions to solve
    # Equation C-A-7-2 of the Commentary on the 2016 AISC Specification
    def fcn(K, GA, GB):
        pK = pi / K
        num = (GA * GB * pK ** 2 - 36) / (6 * (GA + GB)) - pK / tan(pK)
        return num

    def fcn_with_one_zero(K, G):
        pK = pi / K
        num = -6 / G - pK / tan(pK)
        return num

    # Use fsolve to determine the effective length factor
    if GA == 0:
        K = fsolve(fcn_with_one_zero, Kguess, args=GB)
    elif GB == 0:
        K = fsolve(fcn_with_one_zero, Kguess, args=GA)
    else:
        K = fsolve(fcn, Kguess, args=(GA, GB))

    # Basic check for there has been trouble
    if K[0] < 1:
        raise ValueError("K factor is less than 1. This is not valid.")

    return K[0]


def run_examples():
    target_K = 3.0
    GA = inf
    GB = 6/(pi/target_K)/tan(pi/target_K)
    K = sidesway_uninhibited_effective_length_factor(GA, GB)
    print(f'Target K = {target_K}')
    print(f'Calculated K = {K}')


if __name__ == "__main__":
    run_examples()
