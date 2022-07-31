
def find_limit_point_in_list(the_list, val):

    # Quick return if the first entry of the list is the target value
    if the_list[0] == val:
        return 0, 0.0

    # Create a list with True where the corresponding value in the_list 
    # is on the other side of the target value and False otherwise 
    if the_list[0] > val:
        tf_list = [x<=val for x in the_list]
    else:
        tf_list = [x>=val for x in the_list]

    # Find the first instance of True
    try:
        ind = tf_list.index(True) - 1
    except ValueError:
        # If list does not pass target value, return None
        return None,None
    except Exception as e:
        raise e
    
    # Determine where target value is between list items
    if val == the_list[ind]:
        x = 0.0
    else:
        x = (val - the_list[ind]) / (the_list[ind+1] - the_list[ind])

    return ind, x


def interpolate_list(the_list, ind, x=0):

    if (ind is None) or (x is None):
        raise ValueError(f'One or both of ind and x are None ({ind = }) ({x = })')

    if x == 0:
        return the_list[ind]
        
    return the_list[ind] + (the_list[ind+1] - the_list[ind]) * x
