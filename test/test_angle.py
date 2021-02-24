from libdenavit.section import Angle
from libdenavit.section.database import angle_database

properties = ['A', 'x_bar', 'y_bar', 'xp', 'yp', 'Ix', 'Zx', 'Sx', 'rx', 'Iy', 'Zy', 'Sy', 'ry', 'Iz', 'rz', 'J', 'Cw', 'tan_alpha']

for prop in properties:
    print('\n=== Checking %s ===' % prop)

    max_error_upper = 0.
    max_error_lower = 0.

    for key, iAngle in angle_database.items(): 
        
        # Get property from Python class
        s = Angle.from_name(key)
        X_calc = getattr(s, prop)
        
        # Get proproety from database
        if prop == 'x_bar':
            X_database = iAngle['x']
        elif prop == 'y_bar':
            X_database = iAngle['y']
        else:
            X_database = iAngle[prop]    
        
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

# Alternate way to do formatting (f-strings) -- requires Python 3.6+
#
# print(f'{key} --- {X_calc:.4f} / {X_database:.4f} --- {percent_diff:.4f}%')
#
# print(f'Upper limit: {max_error_upper:.4f}%')
# print(f'Lower limit: {max_error_lower:.4f}%')
