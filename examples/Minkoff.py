from libdenavit import Minkoff

# Set basic attributes
attrs = {
    'd': 28,
    'L': 50*12,
    'bt': 2.0,
    'tt': 3/16,
    'bb': 1.75,
    'tb': 5/32,
    'separation': 0.8125,
    'k': 1,
    'w': 0,
    'print_results': False,
}


print('Compare calculated values to tabulated values in Minkoff (1975)')

print('Load at Shear Center')
attrs['yp'] = 'Shear Center'
P = Minkoff(**attrs)
print(f' Calculated, P = {P:.3f} lbs')
print(f' Table 2,    P = 424 lbs')

print('Load at Top Chord')
attrs['yp'] = 'Top Chord Centroid'
P = Minkoff(**attrs)
print(f' Calculated, P = {P:.3f} lbs')
print(f' Table 2,    P = 301 lbs')

print('Load at Bottom Chord')
attrs['yp'] = 'Bottom Chord Centroid'
P = Minkoff(**attrs)
print(f' Calculated, P = {P:.3f} lbs')
print(f' Table 2,    P = 722 lbs')
 
    