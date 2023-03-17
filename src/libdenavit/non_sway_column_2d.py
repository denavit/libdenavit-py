from math import pi, sin
from libdenavit import find_limit_point_in_list, interpolate_list, InteractionDiagram2d, CrossSection2d
from libdenavit.OpenSees import AnalysisResults
import openseespy.opensees as ops
import numpy as np


class NonSwayColumn2d:
    def __init__(self, section, length, et, eb, dxo=0.0, n_elem=6, axis=None):
        # Physical parameters
        self.section = section
        self.length = length
        self.et = et
        self.eb = eb
        self.dxo = dxo
        self.axis = axis
        
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
            self.section.build_ops_fiber_section(section_id, *section_args, **section_kwargs, axis=self.axis)
        else:
            raise ValueError(f'Unknown cross section type {type(self.section).__name__}')

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
        - Compression is positive
        - The vertical load applied to column is P = LFV
        - The moment applied to bottom of column is M = LFH*eb
        - The moment applied to top of column is M = -LFH*et
        - For proportional analyses, LFV and LFH are increased simultaneously
          with a ratio of LFH/LFV = e (P is ignored)
        - For non-proportional analyses, LFV is increased to P first then held
          constant, then LFH is increased (e is ignored)
          
        """
        eigenvalue_as_limit = True
        load_drop_as_limit = True
        deformation_as_limit = True
        strain_as_limit = True

        self.build_ops_model(1, section_args, section_kwargs)
        
        # Initilize analysis results
        results = AnalysisResults()
        results.applied_axial_load = []
        results.applied_moment_top = []
        results.applied_moment_bot = []
        results.maximum_abs_moment = []
        results.maximum_abs_disp = []
        results.lowest_eigenvalue = []
        results.maximum_concrete_compression_strain = []
        results.maximum_steel_strain = []
        # Define function to find limit point
        def find_limit_point():
            ind,x = find_limit_point_in_list(results.lowest_eigenvalue,0)
            if ind is None:
                ind, x = find_limit_point_in_list(results.maximum_concrete_compression_strain, -0.005)
            if ind is None:
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
                dU = self.length * disp_incr_factor
                ops.integrator('DisplacementControl', self.ops_mid_node, 1, dU)
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
                results.maximum_concrete_compression_strain.append(self.ops_get_concrete_compression_strain())
                results.maximum_steel_strain.append(self.get_maximum_steel_strain())

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

                if load_drop_as_limit:
                    if current_applied_axial_load < (1 - perc_drop) * maximum_applied_axial_load:
                        results.exit_message = 'Load Drop Limit Reached'
                        break

                # Check for lowest eigenvalue less than zero
                if eigenvalue_as_limit:
                    if results.lowest_eigenvalue[-1] < 0:
                        results.exit_message = 'Eigenvalue Limit Reached'
                        break

                # Check for maximum displacement
                if deformation_as_limit:
                    if results.maximum_abs_disp[-1] > maximum_abs_disp_limit_ratio * self.length:
                        results.exit_message = 'Deformation Limit Reached'
                        break

                # Check for strain in extreme compressive fiber
                if strain_as_limit:
                    if results.maximum_concrete_compression_strain[-1] < -0.005:
                        results.exit_message = 'Extreme Compressive Fiber Strain Limit Reached'
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
                results.maximum_concrete_compression_strain.append(self.ops_get_concrete_compression_strain())
                results.maximum_steel_strain.append(self.get_maximum_steel_strain())
            
            record()
            
            for i in range(num_steps_vertical):
                ok = ops.analyze(1)
                
                if ok != 0:
                    results.exit_message = 'Analysis Failed In Vertical Loading'
                    return results
                
                record()
                if deformation_as_limit:
                    if results.maximum_abs_disp[-1] > maximum_abs_disp_limit_ratio * self.length:
                        results.exit_message = 'Deformation Limit Reached In Vertical Loading'
                        return results

                # Check for lowest eigenvalue less than zero
                if eigenvalue_as_limit:
                    if results.lowest_eigenvalue[-1] < 0:
                        results.exit_message = 'Eigenvalue Limit Reached In Vertical Loading'
                        return results

                # Check for strain in extreme compressive fiber
                if strain_as_limit:
                    if results.maximum_concrete_compression_strain[-1] < -0.005:
                        results.exit_message = 'Extreme Compressive Fiber Strain Limit Reached In Vertical Loading'
                        return results
            
            # endregion
            
            # region Run lateral load (time = LFH)
            ops.loadConst('-time', 0.0)
            
            ops.timeSeries('Linear', 101)
            ops.pattern('Plain', 201, 101)
            ops.load(self.ops_n_elem, 0, 0, self.et)
            ops.load(0, 0, 0, -self.eb)

            # @todo - we may eventually need more sophisticated selection of dof to control
            #dU = self.length * disp_incr_factor
            #ops.integrator('DisplacementControl', self.ops_mid_node, 1, dU)
            
            if self.et == 0:
                dU = np.sign(-self.eb)*disp_incr_factor/10
                ops.integrator('DisplacementControl', 0, 3, dU)
            else:
                dU = np.sign(self.et)*disp_incr_factor/10
                ops.integrator('DisplacementControl', self.ops_n_elem, 3, dU)
            
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
                results.maximum_concrete_compression_strain.append(self.ops_get_concrete_compression_strain())
                results.maximum_steel_strain.append(self.get_maximum_steel_strain())

            record()
            
            maximum_time = 0
            while True:
                ok = ops.analyze(1)

                if ok != 0:
                    print('Trying ModifiedNewton')
                    ops.algorithm('ModifiedNewton')
                    ok = ops.analyze(1)

                if ok != 0:
                    print('Trying KrylovNewton')
                    ops.algorithm('KrylovNewton')
                    ok = ops.analyze(1)

                if ok != 0:
                    print('Trying KrylovNewton and Greater Tolerance')
                    ops.algorithm('KrylovNewton')
                    ops.test('NormUnbalance', 1e-1, 10)
                    ok = ops.analyze(1)
                
                if ok == 0:
                    ops.algorithm('Newton')
                    ops.test('NormUnbalance', 1e-2, 10)
                else:
                    results.exit_message = 'Analysis Failed'
                    break
                
                record()

                # Check for drop in applied load (time = the horzontal load factor)
                current_time = ops.getTime()
                maximum_time = max(maximum_time, current_time)
                if load_drop_as_limit:
                    if current_time < (1 - perc_drop) * maximum_time:
                        results.exit_message = 'Load Drop Limit Reached'
                        break
                    
                # Check for lowest eigenvalue less than zero
                if eigenvalue_as_limit:
                    if results.lowest_eigenvalue[-1] < 0:
                        results.exit_message = 'Eigenvalue Limit Reached'
                        break
                
                # Check for maximum displacement
                if deformation_as_limit:
                    if results.maximum_abs_disp[-1] > maximum_abs_disp_limit_ratio * self.length:
                        results.exit_message = 'Deformation Limit Reached'
                        break

                # Check for strain in extreme compressive fiber
                if strain_as_limit:
                    if results.maximum_concrete_compression_strain[-1] < -0.005:
                        results.exit_message = 'Extreme Compressive Fiber Strain Limit Reached'
                        break

            find_limit_point()
            return results
        
        else:
            raise ValueError(f'Analysis type {analysis_type} not implemented')

    def run_ops_interaction(self, section_args, section_kwargs, num_points=10, prop_disp_incr_factor=1e-7,
                            nonprop_disp_incr_factor=1e-4, section_load_factor=1e-1):

        # Run one axial load only analyis to determine maximum axial strength
        results = self.run_ops_analysis('proportional_limit_point', section_args, section_kwargs, e=0,
                                        disp_incr_factor=prop_disp_incr_factor)
        P = [results.applied_axial_load_at_limit_point]
        M1 = [0]
        M2 = [results.maximum_abs_moment_at_limit_point]
        exit_message = [results.exit_message]
        if P is None or P == [None]:
            raise ValueError('Analysis failed at axial only loading')

        # Loop axial linearly spaced axial loads witn non-proportional analyses
        for i in range(1,num_points):
            iP = P[0] * (num_points-1-i) / (num_points-1)
            if iP == 0:
                cross_section = CrossSection2d(self.section, self.axis)
                results = cross_section.run_ops_analysis('nonproportional_limit_point', section_args, section_kwargs, P=0, load_incr_factor=section_load_factor)
                P.append(iP)
                M1.append(results.maximum_abs_moment_at_limit_point)
                M2.append(results.maximum_abs_moment_at_limit_point)

            else:
                results = self.run_ops_analysis('nonproportional_limit_point', section_args, section_kwargs, P=iP, disp_incr_factor=nonprop_disp_incr_factor)
                P.append(iP)
                M1.append(max(results.applied_moment_top))
                M2.append(max(results.maximum_abs_moment))
                exit_message.append(results.exit_message)
            plot_at_step=False
            if plot_at_step:
                import matplotlib.pyplot as plt
                plt.figure()
                plt.title(f'Axial Load = {iP}, column: D= {self.section.depth("x")}, L={self.length}')
                plt.plot(results.mid_node_disp, results.applied_moment_top, '-ro')
                plt.plot(results.mid_node_disp, results.maximum_abs_moment, '-bo')
                plt.legend(['M1', 'M2'])

                plt.figure()
                plt.plot(results.mid_node_disp, results.lowest_eigenvalue)
                plt.show()

        return {'P': list(np.array(P)), 'M1': M1, 'M2': M2, 'exit_message': exit_message} # @todo what is going on with P?

    def run_ops_interaction_proportional(self, section_args, section_kwargs, e_list, **kwargs):
        P  = []
        M1 = []
        M2 = []
        
        for e in e_list:
            results = self.run_ops_analysis('proportional_limit_point', section_args, section_kwargs, e=e, **kwargs)
            P.append(results.applied_axial_load_at_limit_point)
            M1.append(results.applied_moment_top_at_limit_point)
            M2.append(results.maximum_abs_moment_at_limit_point)

        return {'P': list(np.array(P)), 'M1': M1, 'M2': M2}

    def run_AASHTO_interaction(self, EI_type, num_points=10, section_factored=True, Pc_factor=0.75, beta_dns=0, minimum_eccentricity=False):
    
        # beta_dns is the ratio of the maximum factored sustained axial load divided by
        # the total factored axial load associated with the same load combination
        # default is zero (i.e., short term loading)
        
        # Note that this function uses 
        #   M1 to mean applied first-order moment 
        #   M2 to mean internal second-order moment
        # this notation is differnt than what is used in AASHTO.
        
        # Parameters
        k = 1  # Effective length factor (always one for this non-sway column)
        EIeff = self.section.EIeff(self.axis, EI_type, beta_dns)
        Pc = pi**2 * EIeff / (k * self.length)**2
        h = self.section.depth(self.axis)
        
        # Get cross-sectional interaction diagram
        P_id, M_id, _ = self.section.section_interaction_2d(self.axis, 100, factored=section_factored)
        id2d = InteractionDiagram2d(M_id, P_id, is_closed=True)

        # Run one axial load only analysis to determine maximum axial strength
        if minimum_eccentricity:
            P_path  = np.linspace(0, max(1.001*min(P_id),-0.999*Pc_factor*Pc), 1000)
            M2_path = np.zeros_like(P_path)
            for i,P in enumerate(P_path):
                delta = max(self.Cm/(1 - (-P)/(Pc_factor*Pc)), 1.0)
                if self.section.units.lower() == "us":
                    M1_min = -P*(0.6 + 0.03 * h)  # ACI 6.6.4.5.4
                elif self.section.units.lower() == "si":
                    M1_min = -P * (0.015 + 0.03 * h)
                else:
                    raise ValueError("The unit system defined in the section is not supported")
                M2_path[i] = delta*M1_min

            iM2, iP = id2d.find_intersection(M2_path, P_path)

            P_list  = [iP]
            M1_list = [0]
            M2_list = [iM2]

        else:
            M1_min = 0
            buckling_load = -Pc_factor*Pc
            if buckling_load < min(P_id):
                P_list = [min(P_id)]
                M1_list = [0]
                M2_list = [0]
            else:
                P_list  = [buckling_load]
                M1_list = [0]
                M2_list = [id2d.find_x_given_y(buckling_load, 'pos')]            

        # Loop axial linearly spaced axial loads witn non-proportional analyses
        for i in range(1,num_points):
            iP = 0.999*P_list[0] * (num_points-i-1) / (num_points-1)
            iM2 = id2d.find_x_given_y(iP, 'pos')
            delta = max(self.Cm / (1 - (-iP) / (Pc_factor * Pc)), 1.0)
            iM1 = iM2/delta
            P_list.append(iP)
            M1_list.append(iM1)
            M2_list.append(iM2)

        results = {'P':list(-1*np.array(P_list)),'M1':M1_list,'M2':M2_list}
        return results

    def ops_get_concrete_compression_strain(self):
        strain = []
        for i in range(self.ops_n_elem):
            axial_strain = ops.nodeDisp(i, 2)
            curvature = ops.nodeDisp(self.ops_mid_node, 3)
            strain.append(self.section.maximum_concrete_compression_strain(axial_strain, curvature, self.axis))
        return min(strain)

    def get_maximum_steel_strain(self):
        strain = []
        for i in range(self.ops_n_elem):
            axial_strain = ops.nodeDisp(i, 2)
            curvature = ops.nodeDisp(i, 3)
            strain.append(self.section.maximum_steel_strain(axial_strain, curvature, self.axis))
        return max(strain)

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

    @property
    def Cm(self):
        Cm = 0.6 + 0.4 * min([self.et, self.eb], key=abs) / max([self.et, self.eb], key=abs)
        return Cm

    def calculated_EI(self, P_list, M1_list, P_CS = None, M_CS = None, section_factored=True, Pc_factor=0.75):
        P_list = np.array(P_list)
        EIgross = self.section.EIgross(self.axis)
    
        if (M_CS is None) or (P_CS is None):
            P_CS, M_CS, _ = self.section.section_interaction_2d(self.axis, 100, factored=section_factored, only_compressive=True)
        id2d = InteractionDiagram2d(M_CS, P_CS, is_closed=False)
        
        EI_list = []

        for P, M1 in zip(P_list, M1_list):
            if P < min(P_CS) or P == max(P_CS):
                EI_list.append(float("nan"))
                continue
                       
            M2 = id2d.find_x_given_y(P, 'pos')
            
            if M1 >= M2:
                EI_list.append(float("nan"))
                continue
            
            delta = M2/M1
            Pc = delta * P / (Pc_factor * (delta - self.Cm))
            k = 1  # Effective length factor (always one for this non-sway column)
            EI = Pc * (k*self.length/pi)**2
            EI_list.append(EI)

        results = {'P':np.array(P_list),'EI':np.array(EI_list),'EIgross':EIgross}
        return results