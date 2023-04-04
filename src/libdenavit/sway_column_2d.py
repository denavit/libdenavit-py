from math import inf, pi, sin
import matplotlib.pyplot as plt
from libdenavit import find_limit_point_in_list, interpolate_list, InteractionDiagram2d
from libdenavit import sidesway_uninhibited_effective_length_factor, CrossSection2d
from libdenavit.OpenSees import AnalysisResults
import openseespy.opensees as ops
import numpy as np


class SwayColumn2d:
    _K = None
    def __init__(self, section, length, k_bot, k_top, gamma, dxo=0.0, Dxo=0.0, n_elem=6, axis=None):
        # Physical parameters
        # Note that the rotational spring stiffnesses (k_top and b_bot) 
        # can be defined from G using k = (6*EI_col)/(G*L)
        self.section = section
        self.length = length
        self.dxo = dxo
        self.Dxo = Dxo
        self.k_top = k_top
        self.k_bot = k_bot
        self.gamma = gamma
        self.axis = axis

        # General options
        self.include_initial_geometric_imperfections = True

        # OpenSees analysis options
        self.ops_n_elem = n_elem
        self.ops_element_type = "mixedBeamColumn"
        self.ops_geom_transf_type = "Corotational"
        self.ops_integration_points = 3

    @property
    def lever_arm(self):
        if self.k_bot == 0 and self.k_top > 0:
            return self.length
        elif self.k_bot > 0 and self.k_top == 0:  
            return self.length
        elif self.k_bot == self.k_top:
            return self.length/2
        else:
            raise ValueError(f'lever_arm not implemented for k_bot = {self.k_bot} and k_top = {self.k_top}')

    def build_ops_model(self, section_id, section_args, section_kwargs):
        ops.wipe()
        ops.model('basic', '-ndm', 2, '-ndf', 3)
        
        # Check if node numbering does not interfere with end stiffnesses
        if self.ops_n_elem >= 1000:
            raise ValueError(f'To have more than 1000 elements, the node numbering scheme needs to be changed')

        # Defining nodes
        for index in range(self.ops_n_elem + 1):
            if self.include_initial_geometric_imperfections:
                x = sin(index / self.ops_n_elem * pi) * self.dxo + index / self.ops_n_elem * self.Dxo
            else:
                x = 0.
            y = index / self.ops_n_elem * self.length
            ops.node(index, x, y)
            ops.mass(index, 1, 1, 1)

        # Define end fixities
        if self.k_bot == inf:
            ops.fix(0, 1, 1, 1)
        elif self.k_bot == 0:
            ops.fix(0, 1, 1, 0)
        elif type(self.k_bot) in [int, float]:
            ops.fix(0, 1, 1, 0)
            ops.node(1000, 0, 0)
            ops.fix(1000, 1, 1, 1)
            ops.uniaxialMaterial('Elastic', 1000, self.k_bot)
            ops.element("zeroLength", 1000, 1000, 0, '-mat', 1000, '-dir', 6)
        else:
            raise ValueError(f'k values are not supported')

        if self.k_top == inf:
            ops.fix(self.ops_n_elem, 0, 0, 1)
        elif type(self.k_bot) in [int, float]:
            if self.include_initial_geometric_imperfections:
                ops.node(1001, self.Dxo, self.length)
            else:
                ops.node(1001, 0, self.length)
            ops.fix(1001, 1, 1, 1)
            ops.uniaxialMaterial('Elastic', 1001, self.k_top)
            ops.element("zeroLength", 1001, 1001, self.ops_n_elem, '-mat', 1001, '-dir', 6)
        else:
            raise ValueError(f'k values are not supported')

        # Define leaning column
        if self.gamma != 0:
            ops.node(1002, 0, 0)
            if self.include_initial_geometric_imperfections:
                ops.node(1003, self.Dxo, self.length)
            else:
                ops.node(1003, 0, self.length)
            ops.fix(1002, 1, 1, 1)
            ops.fix(1003, 0, 0, 1)
            ops.equalDOF(self.ops_n_elem, 1003, 1)
            ops.uniaxialMaterial('Elastic', 1002, 1)
            ops.element('corotTruss', 1002, 1002, 1003, 10e8, 1002)  # @todo - not sure about this part

        ops.geomTransf(self.ops_geom_transf_type, 100)

        if type(self.section).__name__ == "RC":
            self.section.build_ops_fiber_section(section_id, *section_args, **section_kwargs)
        else:
            raise ValueError(f'Unknown cross section type {type(self.section).__name__}')

        ops.beamIntegration("Lobatto", 1, 1, self.ops_integration_points)

        for index in range(self.ops_n_elem):
            ops.element(self.ops_element_type, index, index, index + 1, 100, 1)

    def run_ops_analysis(self, analysis_type, section_args, section_kwargs, e=1.0, P=0, num_steps_vertical=10,
                         disp_incr_factor=0.00005):
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
        section_args : list or tuple
            Non-keyworded arguments for the section's build_ops_fiber_section
        section_kwargs : dict
            Keyworded arguments for the section's build_ops_fiber_section
        
        Loading Notes
        -------------
        - The vertical load applied to column is P = LFV
        - The Horizontal load applied at the top of the column is H that is dependent on the boundary conditions
        - For proportional analyses, LFV and LFH are increased simultaneously
          with a ratio (P is ignored)
        - For non-proportional analyses, LFV is increased to P first then held
          constant, then LFH is increased (e is ignored)
        """

        eigenvalue_limit = 0
        percent_load_drop_limit = 0.05
        deformation_limit = 0.1 * self.length
        concrete_strain_limit = -0.01
        steel_strain_limit = 0.05
        print_limit_point = True

        self.build_ops_model(1, section_args, section_kwargs)

        # region Initialize analysis results
        results = AnalysisResults()
        results.applied_axial_load = []
        results.applied_horizonal_load = []
        results.maximum_abs_moment = []
        results.maximum_abs_disp = []
        results.lowest_eigenvalue = []
        results.moment_at_top = []
        results.moment_at_bottom = []
        results.maximum_concrete_compression_strain = []
        results.maximum_steel_strain = []
        # endregion

        # Define a function to find limit point
        def find_limit_point():
            if print_limit_point:
                print(results.exit_message)
            if 'Analysis Failed' in results.exit_message:
                results.applied_axial_load_at_limit_point = np.nan
                results.applied_moment_top_at_limit_point = np.nan
                results.applied_moment_bot_at_limit_point = np.nan
                results.maximum_abs_moment_at_limit_point = np.nan
                results.maximum_abs_disp_at_limit_point   = np.nan
                return
            elif 'Eigenvalue Limit' in results.exit_message:
                ind,x = find_limit_point_in_list(results.lowest_eigenvalue, eigenvalue_limit)
            elif 'Extreme Compressive Fiber Strain Limit Reached' in results.exit_message:
                ind, x = find_limit_point_in_list(results.maximum_concrete_compression_strain, concrete_strain_limit)
            elif 'Extreme Steel Fiber Strain Limit Reached' in results.exit_message:
                ind, x = find_limit_point_in_list(results.maximum_steel_strain, steel_strain_limit)
            elif 'Deformation Limit Reached' in results.exit_message:
                ind, x = find_limit_point_in_list(results.maximum_abs_disp, deformation_limit)
            elif 'Load Drop Limit Reached' in results.exit_message:
                ind, x = find_limit_point_in_list(results.applied_axial_load, load_drop_limit)

            results.applied_axial_load_at_limit_point = interpolate_list(results.applied_axial_load, ind, x)
            results.applied_horizontal_load_at_limit_point = interpolate_list(results.applied_horizonal_load, ind, x)
            results.maximum_abs_moment_at_limit_point = interpolate_list(results.maximum_abs_moment, ind, x)
            results.maximum_abs_disp_at_limit_point   = interpolate_list(results.maximum_abs_disp, ind, x)

        # Run analysis
        if analysis_type.lower() == 'proportional_limit_point':

            # time = LFV
            ops.timeSeries('Linear', 100)
            ops.pattern('Plain', 200, 100)

            ops.load(self.ops_n_elem, e / self.lever_arm, -1, 0)
            if self.gamma != 0:
                ops.load(1003, 0, -self.gamma, 0)

            ops.constraints('Plain')
            ops.numberer('RCM')
            ops.system('UmfPack')
            ops.test('NormUnbalance', 1e-2, 10)
            ops.algorithm('Newton')

            # Axial only analysis
            dU = self.length * disp_incr_factor
            ops.integrator('DisplacementControl', self.ops_n_elem, 1, dU)

            ops.analysis('Static')

            # Define recorder
            def record():
                time = ops.getTime()
                results.applied_axial_load.append(time)
                results.applied_horizonal_load.append(time * e / self.lever_arm)
                results.maximum_abs_moment.append(self.ops_get_maximum_abs_moment())
                results.maximum_abs_disp.append(self.ops_get_maximum_abs_disp()[0])
                results.lowest_eigenvalue.append(ops.eigen(1)[0])
                results.moment_at_top.append(ops.eleForce(self.ops_n_elem - 1, 6))
                results.moment_at_bottom.append(ops.eleForce(0, 3))
                results.maximum_concrete_compression_strain.append(self.ops_get_section_strains()[0])
                results.maximum_steel_strain.append(self.ops_get_section_strains()[1])

            record()

            maximum_applied_axial_load = 0.
            while True:
                ok = ops.analyze(1)

                if ok != 0:
                    dU = self.length * disp_incr_factor/100
                    ops.integrator('DisplacementControl', self.ops_n_elem, 1, dU)
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
                    dU = self.length * disp_incr_factor
                    ops.integrator('DisplacementControl', self.ops_n_elem, 1, dU)
                    ops.algorithm('Newton')
                    ops.test('NormUnbalance', 1e-2, 10)

                else:
                    results.exit_message = 'Analysis Failed'
                    print('Analysis Failed')
                    break

                record()

                # Check for drop in applied load
                if percent_load_drop_limit is not None:
                    current_applied_axial_load = results.applied_axial_load[-1]
                    maximum_applied_axial_load = max(maximum_applied_axial_load, current_applied_axial_load)
                    load_drop_limit = (1 - percent_load_drop_limit) * maximum_applied_axial_load

                    if abs(current_applied_axial_load) < abs(load_drop_limit):
                        results.exit_message = 'Load Drop Limit Reached'
                        break

                # Check for lowest eigenvalue less than zero
                if eigenvalue_limit is not None:
                    if results.lowest_eigenvalue[-1] < eigenvalue_limit:
                        results.exit_message = 'Eigenvalue Limit Reached'
                        break

                # Check for maximum displacement
                if deformation_limit is not None:
                    if results.maximum_abs_disp[-1] > deformation_limit:
                        results.exit_message = 'Deformation Limit Reached'
                        break

                # Check for strain in extreme compressive concrete fiber
                if concrete_strain_limit is not None:
                    if results.maximum_concrete_compression_strain[-1] < concrete_strain_limit:
                        results.exit_message = 'Extreme Compressive Fiber Strain Limit Reached'
                        break

                # Check for strain in extreme steel fiber
                if steel_strain_limit is not None:
                    if results.maximum_steel_strain[-1] > steel_strain_limit:
                        results.exit_message = 'Extreme Steel Fiber Strain Limit Reached'
                        break

            find_limit_point()
            return results

        elif analysis_type.lower() == 'nonproportional_limit_point':
            # region Run vertical load (time = LFV)
            ops.timeSeries('Linear', 100)
            ops.pattern('Plain', 200, 100)
            ops.load(self.ops_n_elem, 0, -1, 0)
            if self.gamma != 0:
                ops.load(1003, 0, -self.gamma, 0)
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
                results.applied_horizonal_load.append(0.)
                results.maximum_abs_moment.append(self.ops_get_maximum_abs_moment())
                results.maximum_abs_disp.append(self.ops_get_maximum_abs_disp()[0])
                results.lowest_eigenvalue.append(ops.eigen(1)[0])
                results.moment_at_top.append(ops.eleForce(self.ops_n_elem - 1, 6))
                results.moment_at_bottom.append(ops.eleForce(0, 3))
                results.maximum_concrete_compression_strain.append(self.ops_get_section_strains()[0])
                results.maximum_steel_strain.append(self.ops_get_section_strains()[1])

            record()

            for i in range(num_steps_vertical):
                ok = ops.analyze(1)

                if ok != 0:
                    results.exit_message = 'Analysis Failed In Vertical Loading'
                    print('Analysis Failed In Vertical Loading')
                    return results

                record()

                if deformation_limit is not None:
                    if results.maximum_abs_disp[-1] > deformation_limit:
                        results.exit_message = 'Deformation Limit Reached In Vertical Loading'
                        return results

                # Check for lowest eigenvalue less than zero
                if eigenvalue_limit is not None:
                    if results.lowest_eigenvalue[-1] < eigenvalue_limit:
                        results.exit_message = 'Eigenvalue Limit Reached In Vertical Loading'
                        return results

                # Check for strain in extreme compressive concrete fiber
                if concrete_strain_limit is not None:
                    if results.maximum_concrete_compression_strain[-1] < concrete_strain_limit:
                        results.exit_message = 'Extreme Compressive Fiber Strain Limit Reached In Vertical Loading'
                        return results

                # Check for strain in extreme steel fiber
                if steel_strain_limit is not None:
                    if results.maximum_steel_strain[-1] > steel_strain_limit:
                        results.exit_message = 'Extreme Steel Fiber Strain Limit Reached In Vertical Loading'
                        return results
            # endregion

            # region Run lateral load (time = LFH)
            ops.loadConst('-time', 0.0)
            
            ops.timeSeries('Linear', 101)
            ops.pattern('Plain', 201, 101)
            ops.load(self.ops_n_elem, 1, 0, 0)

            dU = self.length * disp_incr_factor
            ops.integrator('DisplacementControl', self.ops_n_elem, 1, dU)

            ops.analysis('Static')

            # Define recorder
            def record():
                results.applied_axial_load.append(P)
                results.applied_horizonal_load.append(ops.getTime())
                results.maximum_abs_moment.append(self.ops_get_maximum_abs_moment())
                results.maximum_abs_disp.append(self.ops_get_maximum_abs_disp()[0])
                results.lowest_eigenvalue.append(ops.eigen(1)[0])
                results.moment_at_top.append(ops.eleForce(self.ops_n_elem - 1, 6))
                results.moment_at_bottom.append(ops.eleForce(0, 3))
                results.maximum_concrete_compression_strain.append(self.ops_get_section_strains()[0])
                results.maximum_steel_strain.append(self.ops_get_section_strains()[1])

            record()

            maximum_time = 0
            while True:
                ok = ops.analyze(1)

                if ok != 0:
                    dU = disp_incr_factor
                    ops.integrator('DisplacementControl', self.ops_n_elem, 3, dU)
                    ok = ops.analyze(1)

                if ok != 0:
                    dU = self.length * disp_incr_factor/1000
                    ops.integrator('DisplacementControl', self.ops_n_elem, 1, dU)
                    ok = ops.analyze(1)

                if ok != 0:
                    dU = disp_incr_factor/1000
                    ops.integrator('DisplacementControl', self.ops_n_elem, 3, dU)
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
                    dU = self.length * disp_incr_factor
                    ops.integrator('DisplacementControl', self.ops_n_elem, 1, dU)
                    ops.algorithm('Newton')
                    ops.test('NormUnbalance', 1e-2, 10)
                else:
                    results.exit_message = 'Analysis Failed'
                    print('Analysis Failed')
                    break

                record()

                # Check for drop in applied load (time = the horizontal load factor)
                if percent_load_drop_limit is not None:
                    current_time = ops.getTime()
                    maximum_time = max(maximum_time, current_time)
                    load_drop_limit = (1 - percent_load_drop_limit) * maximum_time
                    if current_time < load_drop_limit:
                        results.exit_message = 'Load Drop Limit Reached'
                        break

                # Check for lowest eigenvalue less than zero
                if eigenvalue_limit is not None:
                    if results.lowest_eigenvalue[-1] < eigenvalue_limit:
                        results.exit_message = 'Eigenvalue Limit Reached'
                        break

                # Check for maximum displacement
                if deformation_limit is not None:
                    if results.maximum_abs_disp[-1] > deformation_limit:
                        results.exit_message = 'Deformation Limit Reached'
                        break

                # Check for strain in extreme compressive fiber
                if concrete_strain_limit is not None:
                    if results.maximum_concrete_compression_strain[-1] < concrete_strain_limit:
                        results.exit_message = 'Extreme Compressive Fiber Strain Limit Reached'
                        break

                # Check for strain in extreme steel fiber
                if steel_strain_limit is not None:
                    if results.maximum_steel_strain[-1] > steel_strain_limit:
                        results.exit_message = 'Extreme Steel Fiber Strain Limit Reached'
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
        P_applied = results.applied_axial_load
        print(f'{P= }')
        M1 = [0]
        M2 = [results.maximum_abs_moment_at_limit_point]
        exit_message = [results.exit_message]
        if P is np.nan or P == [np.nan]:
            raise ValueError('Analysis failed at axial only loading')

        # Loop axial linearly spaced axial loads with non-proportional analyses
        for i in range(1, num_points):
            iP = P[0] * (num_points - 1 - i) / (num_points - 1)
            if iP == 0:
                cross_section = CrossSection2d(self.section, self.axis)
                results = cross_section.run_ops_analysis('nonproportional_limit_point', section_args, section_kwargs,
                                                         P=0, load_incr_factor=section_load_factor)
                P.append(iP)
                M1.append(results.maximum_abs_moment_at_limit_point)
                M2.append(results.maximum_abs_moment_at_limit_point)
                exit_message.append(results.exit_message)
            else:
                results = self.run_ops_analysis('nonproportional_limit_point', section_args, section_kwargs, P=abs(iP),
                                                disp_incr_factor=nonprop_disp_incr_factor)
                P.append(iP)
                M1.append(abs(results.applied_horizontal_load_at_limit_point * self.lever_arm))
                M2.append(results.maximum_abs_moment_at_limit_point)
                exit_message.append(results.exit_message)

            plot_at_step = False
            if plot_at_step:
                import matplotlib.pyplot as plt
                plt.figure()
                plt.title(f'Axial Load = {iP}, column: D= {self.section.depth("x")}, L={self.length}')
                plt.plot(results.maximum_abs_disp, results.applied_moment_top, '-ro')
                plt.plot(results.maximum_abs_disp, results.maximum_abs_moment, '-bo')
                plt.legend(['M1', 'M2'])

                plt.figure()
                plt.plot(results.maximum_abs_disp, results.lowest_eigenvalue)
                plt.show()

        results = {'P': P, 'M1': M1, 'M2': M2, 'exit_message': exit_message}
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
            max_disp = max(disp)
        return max_disp, disp.index(max_disp)

    @property
    def Cm(self):
        if self.k_top == 0 or self.k_bot == 0:
            Cm = 0.6
        elif self.k_top == self.k_bot:
            Cm = 0.2
        else:
            raise ValueError('Cm not implemented for unequal top and bottom stiffness')
        return Cm

    def run_AASHTO_interaction(self, EI_type, num_points=10, section_factored=True, Pc_factor=0.75, beta_dns=0,
                               minimum_eccentricity=False):

        # beta_dns is the ratio of the maximum factored sustained axial load divided by
        # the total factored axial load associated with the same load combination
        # default is zero (i.e., short term loading)

        # Note that this function uses
        #   M1 to mean applied first-order moment
        #   M2 to mean internal second-order moment
        # this notation is different from what is used in AASHTO.

        # Parameters
        EIeff = self.section.EIeff(self.axis, EI_type, beta_dns)
        if self._K is None:
            k_s = self.effective_length_factor(EIeff)
        else:
            k_s = self._K
        Pc_s = pi ** 2 * EIeff / (k_s * self.length) ** 2
        h = self.section.depth(self.axis)

        # Get cross-sectional interaction diagram
        P_id, M_id, _ = self.section.section_interaction_2d(self.axis, 100, factored=section_factored, only_compressive=True)
        id2d = InteractionDiagram2d(M_id, P_id, is_closed=False)

        # Run one axial load only analysis to determine maximum axial strength
        if minimum_eccentricity:
            P_path = np.linspace(0, max(1.001 * min(P_id), -0.999 * Pc_factor * Pc), 1000)
            M2_path = np.zeros_like(P_path)
            for i, P in enumerate(P_path):
                delta = max(self.Cm / (1 - (-P) / (Pc_factor * Pc)), 1.0)
                if self.section.units.lower() == "us":
                    M1_min = -P * (0.6 + 0.03 * h)  # ACI 6.6.4.5.4
                elif self.section.units.lower() == "si":
                    M1_min = -P * (0.015 + 0.03 * h)
                else:
                    raise ValueError("The unit system defined in the section is not supported")
                M2_path[i] = delta * M1_min

            iM2, iP = id2d.find_intersection(M2_path, P_path)

            P_list = [iP]
            M1_list = [0]
            M2_list = [iM2]

        else:
            buckling_load = Pc_factor * Pc_s
            if buckling_load > max(P_id):
                P_list = [min(P_id)]
                M1_list = [0]
                M2_list = [0]
            else:
                P_list = [buckling_load]
                M1_list = [0]
                M2_list = [id2d.find_x_given_y(buckling_load, 'pos')]

        # Loop axial linearly spaced axial loads with non-proportional analyses
        for i in range(1, num_points):
            iP = P_list[0] * (num_points - 1 - i) / (num_points - 1)
            iM2 = id2d.find_x_given_y(iP, 'pos')
            # delta_ns = max(self.Cm / (1 - (iP) / (Pc_factor * Pc)), 1.0)
            # iM1_ns = iM2 / delta_ns

            delta_s = 1 / (1 - iP / (Pc_factor * Pc_s))
            # if delta_s < 1:
            #     delta_s = 1
            iM1_s = iM2 / delta_s
            P_list.append(iP)
            M1_list.append(iM1_s)
            M2_list.append(iM2)
        results = {'P': P_list, 'M1': M1_list, 'M2': M2_list}
        return results


    def effective_length_factor(self, EI):
        if self.k_bot == 0:
            G_bot = inf
        elif self.k_bot == inf:
            G_bot = 0
        else:
            G_bot = 6 * EI / (self.k_bot * self.length)

        if self.k_top == 0:
            G_top = inf
        elif self.k_top == inf:
            G_top = 0
        else:
            G_top = 6 * EI / (self.k_top * self.length)

        K = sidesway_uninhibited_effective_length_factor(G_bot, G_top)

        return K

    def calculated_EI(self, P_list, M1_list, P_CS=None, M_CS=None, section_factored=True, Pc_factor=0.75):
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

            if M1 == 0:
                EI_list.append(float("nan"))
                continue
            else:
                delta_s = M2 / M1
            Pc = (P) / (Pc_factor * (1-1/delta_s))
            print(f"###P = {P}###")
            print(f'{delta_s= }')
            print(f'{Pc= }')

            if self._K is None:
                k = self.effective_length_factor(EIeff)
            else:
                k = self._K
            EI = Pc * (k * self.length / pi) ** 2
            EI_list.append(EI)

        results = {'P': np.array(P_list), 'EI': np.array(EI_list), 'EIgross': EIgross}
        return results
