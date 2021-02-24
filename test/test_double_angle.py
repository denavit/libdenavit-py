from libdenavit.section import DoubleAngle
from libdenavit.section.database import double_angle_database

properties = ['A', 'y_bar', 'yp', 'Ix', 'Zx', 'Sx', 'rx', 'Iy', 'Zy', 'Sy', 'ry']

for prop in properties:
    print('\n=== Checking %s ===' % prop)

    max_error_upper = 0.
    max_error_lower = 0.

    for key, iDoubleAngle in double_angle_database.items(): 
        
        # Get property from Python class
        s = DoubleAngle.from_name(key)
        X_calc = getattr(s, prop)

        # Get property from database
        if prop == 'y_bar':
            X_database = iDoubleAngle['y']
        else:
            X_database = iDoubleAngle[prop] 
        
        # Compare
        percent_diff = 100*(X_calc-X_database)/X_database
        
        if abs(percent_diff) > 4:
            print('%s --- %.4f / %.4f --- %.4f%%' % (key,X_calc,X_database,percent_diff))
        
        if percent_diff > max_error_upper:
            max_error_upper = percent_diff
        if percent_diff < max_error_lower:
            max_error_lower = percent_diff

    print('Error Summary:')
    print('Upper limit: %.4f%%' % max_error_upper)
    print('Lower limit: %.4f%%' % max_error_lower)
