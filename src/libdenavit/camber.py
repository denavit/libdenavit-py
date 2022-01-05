def camber(xi,L,c,camber_type='CircularArc'):
    '''
    Returns the camber at a position along the length of a beam. 
    
    x = Position along the length of the beam where the camber is to be computed
    L = Length of the beam
    c = Maximum camber (occurs at mid-span, x = L/2) 
    '''
    
    # Check for quick return
    # also avoids divide by zero errors later
    if c == 0:
        return 0.
    
    # Compute camber    
    if camber_type.lower() == 'circulararc':
        r = (c**2 + (L/2)**2)/(2*c)
        z = sqrt(r**2 - (x - L/2)**2) + (c-r)
    else:
        return Exception(f'Unknown camber type {camber_type}')

    return z