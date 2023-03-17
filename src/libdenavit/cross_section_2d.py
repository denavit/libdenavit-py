from math import pi, sin
from libdenavit import find_limit_point_in_list, interpolate_list
from libdenavit.OpenSees import AnalysisResults
import openseespy.opensees as ops


class CrossSection2d:
    def __init__(self, section, axis=None):
        self.section = section
        self.ops_element_type = "zeroLengthSection"
        self.axis = axis

    def build_ops_model(self, section_id, section_args, section_kwargs):
        ops.wipe()
        ops.model('basic', '-ndm', 2, '-ndf', 3)

        if type(self.section).__name__ == "RC":
            self.section.build_ops_fiber_section(section_id, *section_args, **section_kwargs)

        ops.node(1, 0, 0)
        ops.node(2, 0, 0)

        ops.fix(1, 1, 1, 1)
        ops.fix(2, 0, 1, 0)

        ops.mass(2, 1, 1, 1)

        if self.ops_element_type == "zeroLengthSection":
            # tag ndI ndJ  secTag
            ops.element(self.ops_element_type, 1, 1, 2, section_id)
        else:
            raise ValueError(f"ops_element_type {self.ops_element_type} not recognized")

    def run_ops_analysis(self, analysis_type, section_args, section_kwargs, e=0, P=0,
                         perc_drop=0.05, num_steps_vertical=20, load_incr_factor=1e-5):

        self.build_ops_model(1, section_args, section_kwargs)

        # Initialize analysis results
        results = AnalysisResults()
        results.applied_axial_load = []
        results.maximum_abs_moment = []
        results.lowest_eigenvalue = []
        results.extreme_comp_strain = []
        results.maximum_concrete_compression_strain = []
        results.maximum_steel_strain = []

        # Define function to find limit point
        def find_limit_point():
            ind, x = find_limit_point_in_list(results.lowest_eigenvalue, 0)
            if ind is None:
                ind, x = find_limit_point_in_list(results.maximum_concrete_compression_strain, -0.005)
            if ind is None:
                results.applied_axial_load_at_limit_point = None
                results.maximum_abs_moment_at_limit_point = None
            else:
                results.applied_axial_load_at_limit_point = interpolate_list(results.applied_axial_load, ind, x)
                results.maximum_abs_moment_at_limit_point = interpolate_list(results.maximum_abs_moment, ind, x)

        # Run analysis
        if analysis_type.lower() == 'proportional_limit_point':
            # time = LFV
            ops.timeSeries('Linear', 1)
            ops.pattern('Plain', 1, 1)
            ops.load(2, -1, 0, e)
            ops.integrator('LoadControl', load_incr_factor)
            ops.system('SparseGeneral', '-piv')
            ops.test('NormUnbalance', 1e-2, 10)
            ops.numberer('Plain')
            ops.constraints('Plain')
            ops.algorithm('Newton')
            ops.analysis('Static')
            ops.analyze(1)

            # Define recorder
            def record():
                time = ops.getTime()
                results.applied_axial_load.append(time)
                results.maximum_abs_moment.append(0)
                results.lowest_eigenvalue.append(ops.eigen("-fullGenLapack", 1)[0])
                axial_strain = ops.nodeDisp(2, 1)
                curvatureX = ops.nodeDisp(2, 3)
                results.maximum_concrete_compression_strain.append(
                    self.section.maximum_concrete_compression_strain(axial_strain, curvatureX=curvatureX))
                results.maximum_steel_strain.append(
                    self.section.maximum_tensile_steel_strain(axial_strain, curvatureX=curvatureX))

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
                    results.exit_message = 'Load Drop Limit Point is Reached'
                    break

                # Check for lowest eigenvalue less than zero
                if results.lowest_eigenvalue[-1] < 0:
                    results.exit_message = 'Eigen Value Limit Point is Reached'
                    break

                # Check for strain in extreme compressive fiber
                if results.maximum_concrete_compression_strain[-1] < -0.005:
                    results.exit_message = 'Extreme Compressive Fiber Limit Point is Reached'
                    break

            find_limit_point()
            return results

        elif analysis_type.lower() == 'nonproportional_limit_point':
            # region Run vertical load (time = LFV)
            ops.timeSeries('Linear', 100)
            ops.pattern('Plain', 200, 100)
            ops.load(2, -1, 0, 0)
            ops.constraints('Plain')
            ops.numberer('RCM')
            ops.system('UmfPack')
            ops.test('NormUnbalance', 1, 10)
            ops.algorithm('Newton')
            ops.integrator('LoadControl', P / num_steps_vertical)
            ops.analysis('Static')

            # Define recorder
            def record():
                time = ops.getTime()
                results.applied_axial_load.append(time)
                results.maximum_abs_moment.append(0)
                results.lowest_eigenvalue.append(ops.eigen("-fullGenLapack", 1)[0])
                axial_strain = ops.nodeDisp(2, 1)
                curvatureX = ops.nodeDisp(2, 3)
                results.maximum_concrete_compression_strain.append(
                    self.section.maximum_concrete_compression_strain(axial_strain, curvatureX=curvatureX))
                results.maximum_steel_strain.append(
                    self.section.maximum_tensile_steel_strain(axial_strain, curvatureX=curvatureX))

            record()
            for i in range(num_steps_vertical):
                ok = ops.analyze(1)

                if ok != 0:
                    results.exit_message = 'Analysis Failed In Vertical Loading'
                    return results

                record()

            # endregion

            # Run lateral load (time = LFH)
            ops.loadConst('-time', 0.0)

            ops.timeSeries('Linear', 101)
            ops.pattern('Plain', 201, 101)
            ops.load(2, 0, 0, 1)
            ops.integrator('LoadControl', load_incr_factor)
            ops.analysis('Static')

            # Define recorder
            def record():
                time = ops.getTime()
                results.applied_axial_load.append(P)
                results.maximum_abs_moment.append(abs(ops.eleForce(1, 3)))
                results.lowest_eigenvalue.append(ops.eigen("-fullGenLapack", 1)[0])
                axial_strain = ops.nodeDisp(2, 1)
                curvatureX = ops.nodeDisp(2, 3)
                results.maximum_concrete_compression_strain.append(
                    self.section.maximum_concrete_compression_strain(axial_strain, curvatureX))
                results.maximum_steel_strain.append(
                    self.section.maximum_tensile_steel_strain(axial_strain, curvatureX))

            record()

            maximum_time = 0
            while True:
                ok = ops.analyze(1)

                if ok != 0:
                    results.exit_message = 'Analysis Failed'
                    break

                record()

                # Check for drop in applied load (time = the horizontal load factor)
                current_time = ops.getTime()
                maximum_time = max(maximum_time, current_time)
                if current_time < (1 - perc_drop) * maximum_time:
                    results.exit_message = 'Load Drop Limit Point is Reached'
                    break

                # Check for lowest eigenvalue less than zero
                if results.lowest_eigenvalue[-1] < 0:
                    results.exit_message = 'Eigen Value Limit Point is Reached'
                    break
                if results.maximum_concrete_compression_strain[-1] < -0.005:
                    results.exit_message = 'Extreme Compressive Fiber Limit Point is Reached'
                    break

            find_limit_point()
            return results

        else:
            raise ValueError(f'Analysis type {analysis_type} not implemented')

    def run_ops_interaction(self, section_args, section_kwargs, num_points=10, prop_load_incr_factor=1e-2,
                            nonprop_load_incr_factor=1e-2):

        # Run one axial load only analysis to determine maximum axial strength
        print("Running cross-section axial only analysis...")
        results = self.run_ops_analysis('proportional_limit_point', section_args, section_kwargs,
                                        load_incr_factor=prop_load_incr_factor)
        print("Axial only analysis is completed.")
        P = [max(results.applied_axial_load)]
        M = [0]
        if P in [None, [None]]:
            raise ValueError('Analysis failed at axial only loading')

        # Loop axial linearly spaced axial loads with non-proportional analyses
        print("Running cross-section non-proportional analysis...")
        for i in range(1, num_points):
            iP = P[0] * (num_points - 1 - i) / (num_points - 1)
            results = self.run_ops_analysis('nonproportional_limit_point', section_args, section_kwargs, P=iP,
                                            load_incr_factor=nonprop_load_incr_factor)
            P.append(iP)
            M.append(max(results.maximum_abs_moment))
        print("Non-proportional analysis is completed.")
        return {'P': P, "M": M}

    def run_AASHTO_interaction(self, axis, section_factored=True):
        P_id, M_id, _ = self.section.section_interaction_2d(axis, 100, factored=section_factored)
        return {'P': P_id,'M': M_id}
