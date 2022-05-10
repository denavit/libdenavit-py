import bisect

class GoalSeekMonotonic:
    '''
    A class for perforing a goal seek type analysis. This class only works if 
    the relationship between the input and output is monotonically increasing 
    or monotonically increasing.
    '''
    
    def __init__(self,target_output,tolerance,tolerance_mode='BothSides',starting_input=0,basic_input_increment=1):
        self.target_output          = target_output
        self.tolerance              = tolerance
        self.tolerance_mode         = tolerance_mode
        self.starting_input         = starting_input
        self.basic_input_increment  = basic_input_increment
        
        self._data_input  = []
        self._data_output = []
        self._is_increasing = None
        
    def add_data(self,new_input,new_output):
        ind = bisect.bisect_right(self._data_input,new_input)
        self._data_input.insert(ind,new_input)
        self._data_output.insert(ind,new_output)
        
        # Check if output is monotonic and set self._is_increasing
        is_increasing = self._is_increasing
        for i in range(1,len(self._data_output)):
            if is_increasing is None:
                if self._data_output[i-1] != self._data_output[i]:
                    is_increasing = (self._data_output[i-1] < self._data_output[i])
            else:
                if is_increasing != (self._data_output[i-1] < self._data_output[i]):
                    raise Exception('Output is not monotonic')
        self._is_increasing = is_increasing
        
        #print(f'   Input  = {self._data_input}')
        #print(f'   Output = {self._data_output}')
        
    def check_tolerance(self,output):
        if self.tolerance_mode == 'BothSides':
            pass_tf = abs(output-self.target_output) <= self.tolerance
        elif self.tolerance_mode == 'UnderOnly':
            pass_tf = -self.tolerance <= output-self.target_output <= 0
        elif self.tolerance_mode == 'OverOnly':
            pass_tf = 0 <= output-self.target_output <= self.tolerance
        else:
            raise Exception(f'Unknown tolerance mode: {self.tolerance_mode}')
        return pass_tf
        
    def add_and_check(self,new_input,new_output):
        self.add_data(new_input,new_output)
        return self.check_tolerance(new_output)
        
    def suggest_input(self):       
        l = len(self._data_input)
        
        # No output yet
        if l == 0:
            suggested_input = self.starting_input
            return suggested_input
              
        # Only one output so far, use basic increment
        if l == 1:
            suggested_input = self.starting_input + self.basic_input_increment
            return suggested_input
              
        # All output too high so far
        if self.target_output < self._data_output[0]:
            
            if self._data_output[0] == self._data_output[1]:
                # Output not changing, cannot use linear extrapolation, use basic increment
                if self._is_increasing or (self._is_increasing is None):
                    suggested_input = self._data_input[0]-(self.basic_input_increment*l)
                    return suggested_input
                else:
                    suggested_input = self._data_input[0]+(self.basic_input_increment*l)  
                    return suggested_input                  
            
            else:
                # Data shows a change, use linear extrapolation
                in1  = self._data_input[0]
                in2 = self._data_input[1]
                out1  = self._data_output[0]
                out2 = self._data_output[1]
                suggested_input =  in1 + (in2-in1)*(self.target_output-out1)/(out2-out1)            
                return suggested_input
               
               
        # All output too high so far        
        if self.target_output > self._data_output[-1]:
        
            if self._data_output[-1] == self._data_output[-2]:
                # Output not changing, cannot use linear extrapolation, use basic increment
                if self._is_increasing or (self._is_increasing is None):
                    suggested_input = self._data_input[-1]+(self.basic_input_increment*l)
                    return suggested_input
                else:
                    suggested_input = self._data_input[-1]-(self.basic_input_increment*l)
                    return suggested_input
                
            else:
                # Data shows a change, use linear extrapolation
                in1  = self._data_input[-1]
                in2 = self._data_input[-2]
                out1  = self._data_output[-1]
                out2 = self._data_output[-2]
                suggested_input = in1 + (in2-in1)*(self.target_output-out1)/(out2-out1)
                return suggested_input


        # Determine suggested input using Newton's method using secant slope 
        # between two nearest outputs.  
        ind = bisect.bisect_right(self._data_output,self.target_output)      
        dist_hi = abs(self._data_output[ind] - self.target_output)
        dist_lo = abs(self._data_output[ind-1] - self.target_output)
        
        if dist_hi < dist_lo:
            ind1 = ind
            if ind == (l-1):
                # hi is the highest in the output, second closest will be lo
                ind2 = ind-1
            else:
                dist_hi2 = abs(self._data_output[ind+1] - self.target_output)
                if dist_hi2 < dist_lo:
                    ind2 = ind+1
                else:
                    ind2 = ind-1
        else:
            ind1 = ind-1
            if ind == 1:
                # lo is the lowest in the output, second closest will be hi
                ind2 = ind
            else:
                dist_lo2 = abs(self._data_output[ind-2] - self.target_output)
                if dist_lo2 < dist_hi:
                    ind2 = ind-2
                else:
                    ind2 = ind
                    
        in1  = self._data_input[ind1]
        in2  = self._data_input[ind2]
        out1 = self._data_output[ind1]
        out2 = self._data_output[ind2]
        suggested_input = in1 + (in2-in1)*(self.target_output-out1)/(out2-out1)
        
        if suggested_input >= self._data_input[ind]:
            suggested_input = 0.5*(self._data_input[ind-1]+self._data_input[ind])
        if suggested_input <= self._data_input[ind-1]:
            suggested_input = 0.5*(self._data_input[ind-1]+self._data_input[ind])
        
        return suggested_input
    
def run_example():
    A = GoalSeekMonotonic(1.0,0.001,tolerance_mode='UnderOnly')
    
    print(f'Suggest {A.suggest_input()}')
    print(A.add_and_check(0.0,0.5))
    print(f'Suggest {A.suggest_input()}')
    print(A.add_and_check(1.0,0.7))
    print(f'Suggest {A.suggest_input()}')
    print(A.add_and_check(2.0,0.9))
    print(f'Suggest {A.suggest_input()}')   
    print(A.add_and_check(3.0,1.7))
    print(f'Suggest {A.suggest_input()}') 
    print(A.add_and_check(2.125,1.05))
    print(f'Suggest {A.suggest_input()}') 
    print(A.add_and_check(2.083,0.999))
    print(f'Suggest {A.suggest_input()}') 
    print(A.add_and_check(2.084,0.9995))
    print(f'Suggest {A.suggest_input()}') 
    
if __name__ == "__main__":
    run_example()