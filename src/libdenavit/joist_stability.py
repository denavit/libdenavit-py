from libdenavit.section import Angle
from math import pi,sqrt

def Minkoff(**attrs):
    '''
    Returns the critical distributed load or point load of a joist according to the Minkoff equation.
    
    References:
        Minkoff, R. M. (1975). “Stability of Joists During Erection.” M.S. Thesis, Washington 
            University, St. Louis, Missouri.
        SJI. (2020). Standard Specifications for K-Series, LH-Series, and DLH-Series Open Web 
            Steel Joists, and for Joist Girders. ANSI/SJI 100-2020, Steel Joist Institute, 
            Florence, South Carolina.
    
    Parameters:
        d (float): joist depth
        L (float): unbraced length of joist
        tt (float): top chord angle thickness
        bt (float): top chord angle width (optional if yt, At, and Iyt are defined)
        tb (float): bottom chord angle thickness
        bb (float): bottom chord angle width (optional if yb, Ab, and Iyb are defined)        
        separation (float): distance between top and bottom chord angles 
            (optional if Iyt and Iyb are defined)
        yt (float): 
        yb (float):
        At (float): cross-sectional area of top chord (optional, default computed from tt and bt)
        Ab (float): cross-sectional area of bottom chord (optional, default computed from tb and bb)
        Iyt (float): moment of inerita of the top chord about the minor axis of the joist 
            (optional, default computed from tt, bt, and separation)
        Iyb (float): moment of inerita of the bottom chord about the minor axis of the joist 
            (optional, default computed from tb, bb, and separation)
        E (float): modulus of elasticity (optional, default 29,000,000 psi)
        poissons_ratio (float): Poisson's ratio (optional, default 0.3)
        G (float): shear modulus (optional, default computed from E and poissons_ratio)
        k (float): effective length factor (optional, default 0.85)
        yp (float or str): vertical position of load above top of top chord 
            (optional, default load applied at centroid)
            'Shear Center' = load applied at shear center
            'Top Chord Centroid' = load applied at top chord centroid
            'Bottom Chord Centroid' = load applied at bottom chord centroid
        P (float): point load (downward positive) (optional, default 300 lbs)
            if not None, then the function will return the critical distributed load
        w (float): distributed load (downward positive) (optional, default None)
            if not None, then the function will retrun the critical point load
        print_results (bool): prints additonal results
    '''
    # Joist Depth and Length
    d = attrs['d']
    L = attrs['L']

    # Get top chord properties
    tt = attrs['tt']
    bt = attrs.get('bt')
    separation = attrs.get('separation')
    
    if 'yt' in attrs:
        yt = attrs['yt']
    else:
        if bt is None:
            raise ValueError('yt or bt required as input')
        else:
            yt = Angle(bt,bt,tt).y_bar
            
    if 'At' in attrs:
        At = attrs['At']
    else:
        if bt is None:
            raise ValueError('At or bt required as input')
        else:
            At = 2*Angle(bt,bt,tt).A
            
    if 'Iyt' in attrs:
        Iyt = attrs['Iyt']
    else:
        if (bt is None) or (separation is None):
            raise ValueError('Iyt or bt and separation required as input')
        else:
            chord = Angle(bt,bt,tt)
            Iyt = 2*(chord.Iy + chord.A*(chord.y_bar+separation/2)**2)

    # Get bottom chord properties
    tb = attrs['tb']
    bb = attrs.get('bb')
    
    if 'yb' in attrs:
        yt = attrs['yb']
    else:
        if bb is None:
            raise ValueError('yb or bb required as input')
        else:
            yb = Angle(bb,bb,tb).y_bar
            
    if 'Ab' in attrs:
        Ab = attrs['Ab']
    else:
        if bb is None:
            raise ValueError('Ab or bb required as input')
        else:
            Ab = 2*Angle(bb,bb,tb).A
            
    if 'Iyb' in attrs:
        Iyb = attrs['Iyb']
    else:
        if (bb is None) or (separation is None):
            raise ValueError('Iyb or bb and separation required as input')
        else:
            chord = Angle(bb,bb,tb)
            Iyb = 2*(chord.Iy + chord.A*(chord.y_bar+separation/2)**2)

    # Material properties
    E = attrs.get('E',29000000)
    poissons_ratio = attrs.get('poissons_ratio',0.3)
    
    if 'G' in attrs:
        G = attrs['G']
    else:
        G = E/(2*(1+poissons_ratio))
       
    # Effective Length Factor
    k = attrs.get('k',0.85)
    
    # Position of load
    yp = attrs.get('yp')
    
    # Loading
    P = attrs.get('P')
    w = attrs.get('w')
    if (P is None) and (w is None):
        P = 300
        
    # Options
    print_results = attrs.get('print_results',False)
    
    
    # Calculations
    de = d - yt - yb 
    y = Ab*de/(At+Ab)

    Ix = At*y**2 + Ab*(de-y)**2
    Iy = Iyt + Iyb

    J = 1/3*(At*tt**2 + Ab*tb**2)
    Cw = de**2*Iyb*Iyt/Iy

    yo = -y + Iyb*de/Iy   # Location of shear center w.r.t. centroid (positive number means shear center is below the centroid)
    betax = (Ab*(de-y)**3 - At*y**3)/Ix - 2*yo

    # Vertical location of load, P, from shear center (positive when load is above the shear center)
    if yp is None:
        # Locate point load at joist center of gravity
        ae = yo
    elif yp == 'Shear Center':
        ae = 0
    elif yp == 'Top Chord Centroid':
        ae = yo + y
    elif yp == 'Bottom Chord Centroid':
        ae = yo + y - de
    else:
        # Locate point load at yp above the top chord (negative value is below the top chord)
        ae = yo + y + yt + yp

    if print_results:
        print(f'{tt = }')
        print(f'{yt = }')
        print(f'{At = }')
        print(f'{Iyt = }')

        print(f'{tb = }')
        print(f'{yb = }')
        print(f'{Ab = }')
        print(f'{Iyb = }')

        print(f'{d  = }')
        print(f'{L  = }')
        print(f'{de = }')
        print(f'{y  = }')
        print(f'{Ix = }')
        print(f'{Iy = }')    
        print(f'{J  = }')
        print(f'{Cw = }')
        print(f'{yo = }')
        print(f'{betax = }')
        print(f'{ae = }')


    if (P is not None) and (w is not None):
        raise ValueError('Both P and w are defined, not clear what to solve for')
    elif P is not None:
        a = ((pi**2 + 3)/24)**2
        b = P*((pi**2+3)/12)*((pi**2+4)/16) - (pi**4*E*Iy)/(2*(k*L)**3)*(betax*((pi**2-3)/24)-yo/2)
        c = P**2*((pi**2+4)/16)**2 - (pi**4*E*Iy)/(2*(k*L)**3)*(P*(betax*((pi**2-4)/16)-ae)+(pi**4*E*Cw)/(2*(k*L)**3)+(pi**2*G*J)/(2*k*L))
        W = (-b + sqrt(b**2-4*a*c))/(2*a)
        w = W/L
        
        if print_results:
            print(f'{a = }')
            print(f'{b = }')
            print(f'{c = }')
        
        return w
    elif w is not None:
        a = ((pi**2 + 4)/16)**2
        b = (w*L)*((pi**2+3)/12)*((pi**2+4)/16) - (pi**4*E*Iy)/(2*(k*L)**3)*(betax*((pi**2-4)/16)-ae)
        c = (w*L)**2*((pi**2 + 3)/24)**2 - (w*L)*(pi**4*E*Iy)/(2*(k*L)**3)*(betax*((pi**2-3)/24)-yo/2) - (pi**4*E*Iy)/(2*(k*L)**3)*((pi**4*E*Cw)/(2*(k*L)**3)+(pi**2*G*J)/(2*k*L))
        P = (-b + sqrt(b**2-4*a*c))/(2*a)
        
        if print_results:
            print(f'{a = }')
            print(f'{b = }')
            print(f'{c = }')
        
        return P
    else:
        raise ValueError('Should not have arrived here')
    