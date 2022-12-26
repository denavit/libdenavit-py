from math import inf, pi, sin
from libdenavit import find_limit_point_in_list, interpolate_list, InteractionDiagram2d
from libdenavit.OpenSees import AnalysisResults
import openseespy.opensees as ops
import numpy as np


class SwayColumn2d:
    def __init__(self, section, length, k_bot, k_top, gamma, dxo=0.0, Dxo=0.0, n_elem=6):
        # Physical parameters
        self.section = section
        self.length = length
        self.dxo = dxo
        self.Dxo = Dxo
        self.k_top = k_top
        self.k_bot = k_bot
        self.gamma = gamma

        # General options
        self.include_initial_geometric_imperfections = True

        # OpenSees analysis options
        self.ops_n_elem = n_elem
        self.ops_element_type = "mixedBeamColumn"
        self.ops_geom_transf_type = "Corotational"

    @property
    def lever_arm(self):
        if self.k_bot == 0 and self.k_top > 0:
            return self.length
        elif self.k_bot > 0 and self.k_top == 0:  
            return self.length
        elif self.k_bot == self.k_top:
            return self.length/2
        else:
            raise ValueError(f'lever_arm not implemented for k_bot = {k_bot = } and k_top = {self.k_top}')

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
        # endregion

        def find_limit_point():
            # Define a function to find limit point
            ind, x = find_limit_point_in_list(results.lowest_eigenvalue, 0)
            if ind is None:
                # @todo - if eigenvalue does not reach zero, we can define the limit point based on deformation or strain
                results.applied_axial_load_at_limit_point = None
                results.applied_horizontal_load_at_limit_point = None
                results.maximum_abs_moment_at_limit_point = None
                results.maximum_abs_disp_at_limit_point = None
            else:
                results.applied_axial_load_at_limit_point = interpolate_list(results.applied_axial_load, ind, x)
                results.applied_horizontal_load_at_limit_point = interpolate_list(results.applied_horizonal_load, ind, x)
                results.maximum_abs_moment_at_limit_point = interpolate_list(results.maximum_abs_moment, ind, x)
                results.maximum_abs_disp_at_limit_point = interpolate_list(results.maximum_abs_disp, ind, x)

        # Run analysis
        if analysis_type.lower() == 'proportional_limit_point':

            # time = LFV
            ops.timeSeries('Linear', 100)
            ops.pattern('Plain', 200, 100)
               
            ops.load(self.ops_n_elem, e/self.lever_arm, -1, 0)               
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
                results.applied_horizonal_load.append(time*e/self.lever_arm)
                results.maximum_abs_moment.append(self.ops_get_maximum_abs_moment())
                results.maximum_abs_disp.append(self.ops_get_maximum_abs_disp()[0])
                results.lowest_eigenvalue.append(ops.eigen(1)[0])
                results.moment_at_top.append(ops.eleForce(self.ops_n_elem-1, 6))
                results.moment_at_bottom.append(ops.eleForce(0, 3))

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
            if self.gamma != 0:
                ops.load(1003, 0, -self.gamma, 0)
            ops.constraints('Plain')
            ops.numberer('RCM')
            ops.system('UmfPack')
            ops.test('NormUnbalance', 1e-2, 10, 1)
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
                results.moment_at_top.append(ops.eleForce(self.ops_n_elem-1, 6))
                results.moment_at_bottom.append(ops.eleForce(0, 3))

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
                results.moment_at_top.append(ops.eleForce(self.ops_n_elem-1, 6))
                results.moment_at_bottom.append(ops.eleForce(0, 3))

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
                    ops.test('NormUnbalance', 1e-1, 10, 1)
                    ok = ops.analyze(1)

                if ok == 0:
                    ops.algorithm('Newton')
                    ops.test('NormUnbalance', 1e-2, 10, 1)
                else:
                    results.exit_message = 'Analysis Failed'
                    break

                record()

                # Check for drop in applied load (time = the horizontal load factor)
                current_time = ops.getTime()
                maximum_time = max(maximum_time, current_time)
                if current_time < (1 - perc_drop) * maximum_time:
                    results.exit_message = 'Load Drop Limit Point Reached'
                    break

                # Check for lowest eigenvalue less than zero
                if results.lowest_eigenvalue[-1] < 0:
                    results.exit_message = 'Eigenvalue Limit Point Reached'
                    break

                # Check for maximum displacement
                if results.maximum_abs_disp[-1] > maximum_abs_disp_limit_ratio * self.length:
                    results.exit_message = 'Deformation Limit Reached'
                    break

            find_limit_point()
            return results

        else:
            raise ValueError(f'Analysis type {analysis_type} not implemented')

    def run_ops_interaction(self, section_args, section_kwargs, num_points=10, prop_disp_incr_factor=1e-7,nonprop_disp_incr_factor=1e-4):

        # Run one axial load only analyis to determine maximum axial strength
        results = self.run_ops_analysis('proportional_limit_point', section_args, section_kwargs, e=0,
                                        disp_incr_factor=prop_disp_incr_factor)
        P = [results.applied_axial_load_at_limit_point]
        M1 = [0]
        M2 = [results.maximum_abs_moment_at_limit_point]
        exit_message = [results.exit_message]
        if P is None:
            raise ValueError('Analysis failed at axial only loading')

        # Loop axial linearly spaced axial loads with non-proportional analyses
        for i in range(1, num_points):
            iP = P[0] * (num_points - 1 - i) / (num_points - 1)
            if iP == 0:
                iP = 0.001 * P[0]  # @todo - change this to run a cross sectional analysis
            results = self.run_ops_analysis('nonproportional_limit_point', section_args, section_kwargs, P=iP, disp_incr_factor=nonprop_disp_incr_factor)
            P.append(results.applied_axial_load_at_limit_point)
            M1.append(results.applied_horizontal_load_at_limit_point)
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

        # Store results in AnalysisResults object
        results = AnalysisResults()
        results.P = np.array(P)
        results.M1 = np.array(M1)
        results.M2 = np.array(M2)
        results.exit_message = exit_message

        return results

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


if __name__ == "__main__":
    import math
    from libdenavit.section import RC, Circle, ReinfCirc

    # Input Properties
    D = 48
    num_bars = 24
    fc = 4
    fy = 60
    rhosr = 0.02
    Ab = rhosr * math.pi / 4 * D ** 2 / num_bars
    length = 15 * D
    axis = 'x'
    k_top = 2000000
    k_bot = 2000000

    # Define RC section object
    conc_cross_section = Circle(D)
    reinforcement = ReinfCirc(0.5 * D - 3, num_bars, Ab)
    section = RC(conc_cross_section, reinforcement, fc, fy, 'US')

    col = SwayColumn2d(section, length, k_bot, k_top, 0)
    section_args = (1, "ElasticPP", "Concrete04_no_confinement", 20, 20)
    section_kwargs = {'axis': axis}

    # result = col.run_ops_analysis('proportional_limit_point', section_args, section_kwargs)
    result = col.run_ops_analysis('nonproportional_limit_point', section_args, section_kwargs, P=40e2)
    print(f"P: {result.applied_axial_load}")
    print(f"M_top: {result.moment_at_top}")
    print(f"M_bottom: {result.moment_at_bottom}")
    print(f"Max abs disp: {result.maximum_abs_disp}")
    print(result.exit_message)
    # endregion

