import numpy as np

def camber(x,L,c,camber_type='CircularArc'):
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
        sign_c = np.sign(c)
        c = abs(c)
        r = (c**2 + (L/2)**2)/(2*c)
        z = sign_c * (np.sqrt(r**2 - np.power(x - L/2,2)) + (c-r))
    else:
        return Exception(f'Unknown camber type {camber_type}')

    return z
    
def run_example():
    L = 100
    c = 1
    
    x = np.linspace(0,L,100)
    y = camber(x,L,c)
    
    x1 = 0.5*L
    y1 = camber(x1,L,c) 
    
    x2 = 0.10*L
    y2 = camber(x2,L,c)
    
    import matplotlib.pyplot as plt
    plt.plot(x,y,'k-')
    plt.plot(x1,y1,'ro')
    plt.plot(x2,y2,'bo')
    plt.show()
    
if __name__ == "__main__":
    run_example()