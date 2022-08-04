from math import pi, sin
from libdenavit import find_limit_point_in_list, interpolate_list
from libdenavit.OpenSees import AnalysisResults
import openseespy.opensees as ops
import numpy as np


class NonSwayColumn2d:
    def __init__(self, section, length, et, eb, dxo=0.0, n_elem=8):
        # Physical parameters
        self.section = section
        self.length = length
        self.et = et
        self.eb = eb
        self.dxo = dxo
        
        # General options
        self.include_initial_geometric_imperfections = True
        
        # OpenSees analysis options
        self.ops_n_elem = n_elem
        self.ops_element_type = "mixedBeamColumn"
        self.ops_geom_transf_type = "Corotational"
    
    @property
    def ops_mid_node(self):
        if self.ops_n_elem % 2 == 0:
            return self.ops_n_elem / 2
        else:
            raise ValueError(f'Number of elements should be even {self.ops_n_elem = }')
    
    def build_ops_model(self, section_id, section_args, section_kwargs):
        ops.wipe()
        ops.model('basic', '-ndm', 2, '-ndf', 3)
        
        for index in range(self.ops_n_elem + 1):
            if self.include_initial_geometric_imperfections:
                x = sin(index / self.ops_n_elem * pi) * self.dxo
            else:
                x = 0.
            y = index / self.ops_n_elem * self.length
            ops.node(index, x, y)
            ops.mass(index, 1, 1, 1)
        
        ops.fix(0, 1, 1, 0)
        ops.fix(self.ops_n_elem, 1, 0, 0)
        
        ops.geomTransf(self.ops_geom_transf_type, 100)
        
        if type(self.section).__name__ == "RC":
            self.section.build_ops_fiber_section(section_id, *section_args, **section_kwargs)
        
        ops.beamIntegration("Lobatto", 1, 1, 3)
        
        for index in range(self.ops_n_elem):
            ops.element(self.ops_element_type, index, index, index + 1, 100, 1)
    
    def run_ops_analysis(self, analysis_type, section_args, section_kwargs, e=1.0, P=0,
                         perc_drop=0.05, maximum_abs_disp_limit_ratio=0.1, num_steps_vertical=10, disp_incr_factor=0.00005):
        """ Run an OpenSees analysis of the column
        
        Parameters
        ----------
        analysis_type : str
            The type of analysis to run, options are
                - 'proportional_limit_point'
                - 'nonproportional_limit_point'
                - 'proportional_target_force' (not yet implemented)
                - 'nonproportional_target_force' (not yet implemented)
                - 'proportional_target_disp' (not yet implemented)
                - 'nonproportional_target_disp' (not yet implemented)
        section_args : list
            Non-keyworded arguments for the section's build_ops_fiber_section 
        section_kwargs : dict
            Keyworded arguments for the section's build_ops_fiber_section
        
        Loading Notes
        -------------
        - The vertical load applied to column is P = LFV
        - The moment applied to bottom of column is M = LFH*eb
        - The moment applied to top of column is M = -LFH*et
        - For proportional analyses, LFV and LFH are increased simultaneously 
          with a ratio of LFH/LFV = e (P is ignored)
        - For non-proportional analyses, LFV is increased to P first then held 
          constant, then LFH is increased (e is ignored)
          
        """
        
        self.build_ops_model(1, section_args, section_kwargs)
        
        # Initilize analysis results
        results = AnalysisResults()
        results.applied_axial_load = []
        results.applied_moment_top = []
        results.applied_moment_bot = []
        results.maximum_abs_moment = []
        results.maximum_abs_disp = []
        results.lowest_eigenvalue = []

        # Define function to find limit point
        def find_limit_point():
            ind,x = find_limit_point_in_list(results.lowest_eigenvalue,0)
            if ind is None:
                # @todo - if eigenvalue does not reach zero, we can define the limit point based on deformation or strain
                results.applied_axial_load_at_limit_point = None
                results.applied_moment_top_at_limit_point = None
                results.applied_moment_bot_at_limit_point = None
                results.maximum_abs_moment_at_limit_point = None
                results.maximum_abs_disp_at_limit_point   = None
            else:
                results.applied_axial_load_at_limit_point = interpolate_list(results.applied_axial_load,ind,x)
                results.applied_moment_top_at_limit_point = interpolate_list(results.applied_moment_top,ind,x)
                results.applied_moment_bot_at_limit_point = interpolate_list(results.applied_moment_bot,ind,x)
                results.maximum_abs_moment_at_limit_point = interpolate_list(results.maximum_abs_moment,ind,x)
                results.maximum_abs_disp_at_limit_point   = interpolate_list(results.maximum_abs_disp,ind,x)
        
        # Run analysis
        if analysis_type.lower() == 'proportional_limit_point':
            
            # time = LFV
            ops.timeSeries('Linear', 100)
            ops.pattern('Plain', 200, 100)
            ops.load(self.ops_n_elem, 0, -1, self.et * e)
            ops.load(0, 0, 0, -self.eb * e)
            ops.constraints('Plain')
            ops.numberer('RCM')
            ops.system('UmfPack')
            ops.test('NormUnbalance', 1e-2, 10)
            ops.algorithm('Newton')
            
            # @todo - we may eventually need more sophisticated selection of dof to control
            if self.et * e == 0. and self.eb * e == 0.:
                # Axial only analysis
                dU = -self.length * disp_incr_factor
                ops.integrator('DisplacementControl', self.ops_n_elem, 2, dU)
            else:
                dU = self.length * disp_incr_factor
                ops.integrator('DisplacementControl', self.ops_mid_node, 1, dU)
            
            ops.analysis('Static')
            
            # Define recorder
            def record():
                time = ops.getTime()
                results.applied_axial_load.append(time)
                results.applied_moment_top.append(self.et * e * time)
                results.applied_moment_bot.append(-self.eb * e * time)
                results.maximum_abs_moment.append(self.ops_get_maximum_abs_moment())
                results.maximum_abs_disp.append(self.ops_get_maximum_abs_disp())
                results.lowest_eigenvalue.append(ops.eigen(1)[0])

            record()
            
            maximum_applied_axial_load = 0.
            while True:
                ok = ops.analyze(1)
                
                if ok != 0:
                    results.exit_message = 'Analysis Failed'
                    break
                
                record()

                # Check for drop in applied load
                current_applied_axial_load = results.applied_axial_load[-1]
                maximum_applied_axial_load = max(maximum_applied_axial_load, current_applied_axial_load)
                if current_applied_axial_load < (1 - perc_drop) * maximum_applied_axial_load:
                    results.exit_message = 'Limit Point Reached'
                    break

                # Check for lowest eigenvalue less than zero
                if results.lowest_eigenvalue[-1] < 0:
                    results.exit_message = 'Limit Point Reached'
                    break
                
                # Check for maximum displacement
                if results.maximum_abs_disp[-1] > maximum_abs_disp_limit_ratio * self.length:
                    results.exit_message = 'Deformation Limit Reached'
                    break

            find_limit_point()
            return results
        
        elif analysis_type.lower() == 'nonproportional_limit_point':
            # region Run vertical load (time = LFV)
            ops.timeSeries('Linear', 100)
            ops.pattern('Plain', 200, 100)
            ops.load(self.ops_n_elem, 0, -1, 0)
            ops.constraints('Plain')
            ops.numberer('RCM')
            ops.system('UmfPack')
            ops.test('NormUnbalance', 1e-2, 10)
            ops.algorithm('Newton')
            ops.integrator('LoadControl', P / num_steps_vertical)
            ops.analysis('Static')
            
            # Define recorder
            def record():
                time = ops.getTime()
                results.applied_axial_load.append(time)
                results.applied_moment_top.append(0)
                results.applied_moment_bot.append(0)
                results.maximum_abs_moment.append(self.ops_get_maximum_abs_moment())
                results.maximum_abs_disp.append(self.ops_get_maximum_abs_disp())
                results.lowest_eigenvalue.append(ops.eigen(1)[0])
            
            record()
            
            for i in range(num_steps_vertical):
                ok = ops.analyze(1)
                
                if ok != 0:
                    results.exit_message = 'Analysis Failed In Vertical Loading'
                    return results
                
                record()
                
                if results.maximum_abs_disp[-1] > maximum_abs_disp_limit_ratio * self.length:
                    results.exit_message = 'Deformation Limit Reached In Vertical Loading'
                    return results
            
            # endregion
            
            # region Run lateral load (time = LFH)
            ops.loadConst('-time', 0.0)
            
            ops.timeSeries('Linear', 101)
            ops.pattern('Plain', 201, 101)
            ops.load(self.ops_n_elem, 0, 0, self.et)
            ops.load(0, 0, 0, -self.eb)
            
            # @todo - we may eventually need more sophisticated selection of dof to control
            dU = self.length * disp_incr_factor
            ops.integrator('DisplacementControl', self.ops_mid_node, 1, dU)
            
            ops.analysis('Static')
            
            # Define recorder
            def record():
                time = ops.getTime()
                results.applied_axial_load.append(P)
                results.applied_moment_top.append(self.et * time)
                results.applied_moment_bot.append(-self.eb * time)
                results.maximum_abs_moment.append(self.ops_get_maximum_abs_moment())
                results.maximum_abs_disp.append(self.ops_get_maximum_abs_disp())
                results.lowest_eigenvalue.append(ops.eigen(1)[0])

            record()
            
            maximum_time = 0
            while True:
                ok = ops.analyze(1)
                
                if ok != 0:
                    results.exit_message = 'Analysis Failed'
                    break
                
                record()

                # Check for drop in applied load (time = the horzontal load factor)
                current_time = ops.getTime()
                maximum_time = max(maximum_time, current_time)
                if current_time < (1 - perc_drop) * maximum_time:
                    results.exit_message = 'Limit Point Reached'
                    break
                    
                # Check for lowest eigenvalue less than zero
                if results.lowest_eigenvalue[-1] < 0:
                    results.exit_message = 'Limit Point Reached'
                    break
                
                # Check for maximum displacement
                if results.maximum_abs_disp[-1] > maximum_abs_disp_limit_ratio * self.length:
                    results.exit_message = 'Deformation Limit Reached'
                    break

            find_limit_point()
            return results
        
        else:
            raise ValueError(f'Analysis type {analysis_type} not implemented')

    def run_ops_interaction(self, section_args, section_kwargs, num_points=10):
        
        # Run one axial load only analyis to determine maximum axial strength
        results = self.run_ops_analysis('proportional_limit_point', section_args, section_kwargs, e=0)
        P = [results.applied_axial_load_at_limit_point]
        M1 = [0]
        M2 = [results.maximum_abs_moment_at_limit_point]

        # Loop axial linearly spaced axial loads witn non-proportional analyses
        for i in range(1,num_points):
            iP = P[0] * (num_points-1-i) / (num_points-1)
            if iP == 0:
                iP = 0.001*P[0] # @todo - change this to run a cross sectional analysis
            results = self.run_ops_analysis('nonproportional_limit_point', section_args, section_kwargs, P=iP)
            P.append(iP)
            M1.append(max(results.applied_moment_top))
            M2.append(max(results.maximum_abs_moment))

        return P, M1, M2

    def run_ops_interaction_proportional(self, section_args, section_kwargs, e_list, **kwargs):
        P  = []
        M1 = []
        M2 = []
        
        for e in e_list:
            results = self.run_ops_analysis('proportional_limit_point', section_args, section_kwargs, e=e, **kwargs)
            P.append(results.applied_axial_load_at_limit_point)
            M1.append(results.applied_moment_top_at_limit_point)
            M2.append(results.maximum_abs_moment_at_limit_point)

        return P, M1, M2

    def ops_get_maximum_abs_moment(self):
        # This code assumed (but does not check) that moment at j-end of 
        # one element equals the moment at the i-end of the next element.
        moment = [abs(ops.eleForce(0, 3))]
        for i in range(self.ops_n_elem):
            moment.append(abs(ops.eleForce(i, 6)))
        
        return max(moment)
    
    def ops_get_maximum_abs_disp(self):
        disp = []
        for i in range(self.ops_n_elem + 1):
            disp.append(abs(ops.nodeDisp(i, 1)))
        
        return max(disp)
