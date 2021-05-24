
def available_strength(Rn,strength_type,phi,Omega):
    if strength_type.lower() == 'nominal':
        return Rn
    elif strength_type.lower() == 'design' or strength_type.lower() == 'lrfd':
        return phi*Rn
    elif strength_type.lower() == 'allowable' or strength_type.lower() == 'asd':
        return Rn/Omega
    else:
        raise Exception('Unknown strength_type: %s' % strength_type)  
