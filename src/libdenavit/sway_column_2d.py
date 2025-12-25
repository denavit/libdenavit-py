import warnings
from math import inf, pi, sin
import matplotlib.pyplot as plt
from libdenavit import opensees as ops
from libdenavit import find_limit_point_in_list, interpolate_list, InteractionDiagram2d
from libdenavit import sidesway_uninhibited_effective_length_factor, CrossSection2d
from libdenavit.OpenSees import AnalysisResults
import numpy as np
from scipy.optimize import fsolve
from libdenavit.analysis_helpers import try_analysis_options, ops_get_section_strains, ops_get_maximum_abs_moment, ops_get_maximum_abs_disp, check_analysis_limits
from libdenavit.column_2d import Column2d


class SwayColumn2d(Column2d):   
    def __init__(self, section, length, k_bot, k_top, gamma, **kwargs):
        """
            Represents a sway 2D column

            This class defines a sway column element with physical parameters such as section properties,
            length, and boundary conditions. It also allows customization of analysis options.

            Parameters:
                section: The section object representing the cross-sectional properties.
                length: The length of the entire column.
                k_bot: The rotational spring stiffness at the bottom of the column.
                k_top: The rotational spring stiffness at the top of the column.
                gamma: The gravity load ratio for leaning column effect.
                kwargs: Additional keyword arguments for customization.
                          Dxo (float, optional): Lateral displacement at the top. Default is 0.0.
                          effective_length_factor_override (float or None, optional): Override for effective length factor. Default is None.
                          axis (str, optional): Axis. Default is None.
                          n_elem (int, optional): Number of elements for OpenSees analysis. Default is 8.
                          element_type (str, optional): Type of OpenSees element. Default is 'mixedBeamColumn'.
                          ops_geom_transf_type (str, optional): OpenSees geometric transformation type. Default is 'Corotational'.
                          ops_integration_points (int, optional): Number of integration points for OpenSees analysis. Default is 3.
        """

        super().__init__(section, length, **kwargs)
        
        # Physical parameters
        # Note that the rotational spring stiffnesses (k_top and k_bot) 
        # can be defined from G using k = (6*EI_col)/(G*L)
        self.k_top = k_top
        self.k_bot = k_bot
        self.gamma = gamma
        
        # Specific defaults for sway analysis
        self.Dxo = kwargs.get('Dxo', 0.0)
        self.effective_length_factor_override = kwargs.get('effective_length_factor_override', None)


    def _initialize_results(self):
        """Initialize analysis results object with required attributes."""
        # Get base attributes and add Sway-specific ones
        results = super()._initialize_results()
        # Add Sway-specific attributes
        for attr in ['applied_horizontal_load', 'moment_at_top', 'moment_at_bottom']:
            setattr(results, attr, [])
        return results


    def _set_limit_point_values(self, results, ind, x):
        """Override to include horizontal load values."""
        super()._set_limit_point_values(results, ind, x)
        results.applied_horizontal_load_at_limit_point = interpolate_list(results.applied_horizontal_load, ind, x)


    def _extract_analysis_config(self, **kwargs):
        """Override to adjust default values for sway analysis."""
        config = super()._extract_analysis_config(**kwargs)
        # Adjust defaults for sway analysis if not explicitly provided
        if 'num_steps_vertical' not in kwargs:
            config['num_steps_vertical'] = 10
        if 'disp_incr_factor' not in kwargs:
            config['disp_incr_factor'] = 5e-05
        return config


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


    def _initialize_results(self):
        """Initialize analysis results object with required attributes."""
        # Get base attributes and add Sway-specific ones
        results = super()._initialize_results()
        # Add Sway-specific attributes
        for attr in ['applied_horizontal_load', 'moment_at_top', 'moment_at_bottom']:
            setattr(results, attr, [])
        return results


    def build_ops_model(self, section_id, section_args, section_kwargs, **kwargs):
        ops.wipe()
        ops.model('basic', '-ndm', 2, '-ndf', 3)
        
        # Check if node numbering does not interfere with end stiffnesses
        if self.ops_n_elem >= 1000:
            raise ValueError(f'To have more than 1000 elements, the node numbering scheme needs to be changed')

        # Defining nodes
        for index in range(self.ops_n_elem + 1):
            if isinstance(self.dxo, (int, float)) and isinstance(self.Dxo, (int, float)):
                x = sin(index / self.ops_n_elem * pi) * self.dxo + index / self.ops_n_elem * self.Dxo
            else:
                raise ValueError(f'Unknown value of dxo ({self.dxo}) or Dxo ({self.Dxo})')

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
            ops.node(1001, self.Dxo, self.length)
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
            self.section.build_ops_fiber_section(section_id, *section_args, **section_kwargs, axis=self.axis)
        else:
            raise ValueError(f'Unknown cross section type {type(self.section).__name__}')

        ops.beamIntegration("Lobatto", 1, 1, self.ops_integration_points)

        for index in range(self.ops_n_elem):
            ops.element(self.ops_element_type, index, index, index + 1, 100, 1)


    def _set_limit_point_values(self, results, ind, x):
        """Override to include horizontal load values."""
        super()._set_limit_point_values(results, ind, x)
        results.applied_horizontal_load_at_limit_point = interpolate_list(results.applied_horizontal_load, ind, x)


    def run_ops_interaction(self, **kwargs):
        # Parse keyword arguments
        section_id = kwargs.get('section_id', 1)
        section_args = kwargs.get('section_args', ())
        section_kwargs = kwargs.get('section_kwargs', {})
        num_points = kwargs.get('num_points', 10)
        prop_disp_incr_factor = kwargs.get('prop_disp_incr_factor', 1e-7)
        nonprop_disp_incr_factor = kwargs.get('nonprop_disp_incr_factor', 1e-4)
        section_load_factor = kwargs.get('section_load_factor', 1e-1)
        plot_load_deformation = kwargs.get('plot_load_deformation', False)
        full_results = kwargs.get('full_results', False)

        if plot_load_deformation:
            fig_at_step, ax_at_step = plt.subplots(2, 1, figsize=(10, 6), gridspec_kw={'height_ratios': [3, 1]})

        # Run one axial load only analyis to determine maximum axial strength
        results = self.run_ops_analysis('proportional_limit_point', e=0, section_id=section_id,
                                        section_args=section_args, section_kwargs=section_kwargs,
                                        disp_incr_factor=prop_disp_incr_factor, deformation_limit=None)
        P = [results.applied_axial_load_at_limit_point]
        M1 = [0]
        M2 = [results.maximum_abs_moment_at_limit_point]
        if full_results:
            M1t_path = [results.moment_at_top]
            M1b_path = [results.moment_at_bottom]
            M2_path = [results.maximum_abs_moment]
            applied_h_path = [results.applied_horizontal_load]
            abs_disp_path = [results.maximum_abs_disp]


        exit_message = [results.exit_message]
        if P is np.nan or P == [np.nan]:
            raise ValueError('Analysis failed at axial only loading')

        # Loop axial linearly spaced axial loads with non-proportional analyses
        for i in range(1, num_points):
            iP = P[0] * (num_points - 1 - i) / (num_points - 1)
            if iP == 0:
                cross_section = CrossSection2d(self.section, self.axis)
                results = cross_section.run_ops_analysis('nonproportional_limit_point',
                                                         section_id=section_id, section_args=section_args,
                                                         section_kwargs=section_kwargs,P=0,
                                                         load_incr_factor=section_load_factor)
                P.append(iP)
                M1.append(results.maximum_abs_moment_at_limit_point)
                M2.append(results.maximum_abs_moment_at_limit_point)

                if full_results:
                    M1t_path.append(results.maximum_abs_moment)
                    M1b_path.append(results.maximum_abs_moment)
                    M2_path.append(results.maximum_abs_moment)
                    applied_h_path.append([0]*len(results.maximum_abs_moment))
                    abs_disp_path.append([0]*len(results.maximum_abs_moment))

                exit_message.append(results.exit_message)
            else:
                results = self.run_ops_analysis('nonproportional_limit_point', section_id=section_id,
                                                section_args=section_args, section_kwargs=section_kwargs, P=abs(iP),
                                                disp_incr_factor=nonprop_disp_incr_factor, deformation_limit=None)
                P.append(iP)
                M1.append(abs(results.applied_horizontal_load_at_limit_point * self.lever_arm))
                M2.append(results.maximum_abs_moment_at_limit_point)

                if full_results:
                    M1t_path.append(results.moment_at_top)
                    M1b_path.append(results.moment_at_bottom)
                    M2_path.append(results.maximum_abs_moment)
                    applied_h_path.append(results.applied_horizontal_load)
                    abs_disp_path.append(results.maximum_abs_disp)

                exit_message.append(results.exit_message)

            if plot_load_deformation:
                if iP==0:
                    continue
                ax_at_step[0].plot(results.maximum_abs_disp, np.array(results.applied_horizontal_load) * self.lever_arm, '-o', label=f'{iP:,.0f}', markersize=5)
                ax_at_step[0].legend()

                ax_at_step[1].plot(results.maximum_abs_disp, results.lowest_eigenvalue, label=f'{iP:,.0f}', markersize=5)
        if plot_load_deformation:
            ax_at_step[0].set_xlabel('Displacement (in)')
            ax_at_step[0].set_ylabel('Applied Moment (kips)')

            ax_at_step[1].set_xlabel('Displacement (in)')
            ax_at_step[1].set_ylabel('Eigenvalue')
            ax_at_step[1].set_ylim(-100,)

            fig_at_step.tight_layout()
            plt.show()

        if full_results:
            return {'P': np.array(P), 'M1': np.array(M1), 'M2': np.array(M2), 'exit_message': exit_message,
                    'M1t_path': np.array(M1t_path), 'M1b_path': np.array(M1b_path),'M2_path': np.array(M2_path),
                    'applied_h_path': np.array(applied_h_path), 'abs_disp_path': np.array(abs_disp_path)}
        else:
            return {'P': np.array(P), 'M1': np.array(M1), 'M2': np.array(M2), 'exit_message': exit_message}

    
    def _run_sway_proportional_limit_point(self, config, results):
        """Run proportional limit point analysis for sway column."""
        
        def update_dU(disp_incr_factor, div_factor=1, analysis_type='proportional_limit_point'):
            if analysis_type == "proportional_limit_point":
                dU = self.length * disp_incr_factor / div_factor
                ops.integrator('DisplacementControl', self.ops_n_elem, 1, dU)
            elif analysis_type == "nonproportional_limit_point":
                dU = self.length * disp_incr_factor / div_factor
                ops.integrator('DisplacementControl', self.ops_n_elem, 1, dU)

        def reset_analysis_options(disp_incr_factor):
            update_dU(disp_incr_factor)
            ops.algorithm('Newton')
            ops.test('NormUnbalance', 1e-3, 10)

        # time = LFV
        ops.timeSeries('Linear', 100)
        ops.pattern('Plain', 200, 100)

        ops.load(self.ops_n_elem, config['e'] / self.lever_arm, -1, 0)
        if self.gamma != 0:
            ops.load(1003, 0, -self.gamma, 0)

        ops.constraints('Plain')
        ops.numberer('RCM')
        ops.system('UmfPack')
        ops.test('NormUnbalance', 1e-3, 10)
        ops.algorithm('Newton')

        # Axial only analysis
        dU = self.length * config['disp_incr_factor']
        ops.integrator('DisplacementControl', self.ops_n_elem, 1, dU)

        ops.analysis('Static')

        # Define recorder
        def record():
            time = ops.getTime()
            section_strains = ops_get_section_strains(self)

            results.applied_axial_load.append(time)
            results.applied_horizontal_load.append(time * config['e'] / self.lever_arm)
            results.maximum_abs_moment.append(ops_get_maximum_abs_moment(self))
            results.maximum_abs_disp.append(ops_get_maximum_abs_disp(self))
            results.lowest_eigenvalue.append(ops.eigen('-fullGenLapack', 1)[0])
            results.moment_at_top.append(ops.eleForce(self.ops_n_elem - 1, 6))
            results.moment_at_bottom.append(ops.eleForce(0, 3))
            results.maximum_concrete_compression_strain.append(section_strains[0])
            results.maximum_steel_strain.append(section_strains[1])
            if self.axis == 'x':
                results.curvature.append(section_strains[2])
            elif self.axis == 'y':
                results.curvature.append(section_strains[3])

        record()

        maximum_applied_axial_load = 0.
        disp_incr_factor = config['disp_incr_factor']
        
        while True:
            ok = ops.analyze(1)

            if ok != 0 and config['try_smaller_steps']:
                for div_factor in [1e1, 1e2, 1e3, 1e4, 1e5, 1e6]:
                    update_dU(disp_incr_factor, div_factor, analysis_type='proportional_limit_point')
                    ok = ops.analyze(1)
                    if ok == 0 and div_factor in [1e3, 1e4, 1e5, 1e6]:
                        disp_incr_factor /= 10
                        break
                    elif ok == 0:
                        break
                    else:
                        ok = try_analysis_options()
                        if ok == 0 and div_factor in [1e3, 1e4, 1e5, 1e6]:
                            disp_incr_factor /= 10
                            break
                        elif ok == 0:
                            break

            if ok != 0 and not config['try_smaller_steps']:
                ok = try_analysis_options()

            if ok == 0:
                reset_analysis_options(disp_incr_factor)

            elif ok != 0:
                results.exit_message = 'Analysis Failed'
                warnings.warn('Analysis Failed')
                break

            record()

            # Check for drop in applied load
            if config['percent_load_drop_limit'] is not None:
                current_applied_axial_load = results.applied_axial_load[-1]
                maximum_applied_axial_load = max(maximum_applied_axial_load, current_applied_axial_load)
                load_drop_limit = (1 - config['percent_load_drop_limit']) * maximum_applied_axial_load

                if abs(current_applied_axial_load) < abs(load_drop_limit):
                    results.exit_message = 'Load Drop Limit Reached'
                    break

            # Check for limits
            exit_message = check_analysis_limits(results, eigenvalue_limit=config['eigenvalue_limit'], 
                                   deformation_limit=config['deformation_limit'],
                                   concrete_strain_limit=config['concrete_strain_limit'], 
                                   steel_strain_limit=config['steel_strain_limit'],
                                   section_type=type(self.section).__name__)
            if exit_message:
                results.exit_message = exit_message
                break
        
        return results


    def _run_proportional_analysis(self, config, results):
        """Run proportional analysis for sway column."""
        return self._run_sway_proportional_limit_point(config, results)


    def _run_sway_nonproportional_limit_point(self, config, results):
        """Run nonproportional limit point analysis for sway column."""
        
        def update_dU(disp_incr_factor, div_factor=1, analysis_type='nonproportional_limit_point'):
            if analysis_type == "proportional_limit_point":
                dU = self.length * disp_incr_factor / div_factor
                ops.integrator('DisplacementControl', self.ops_n_elem, 1, dU)
            elif analysis_type == "nonproportional_limit_point":
                dU = self.length * disp_incr_factor / div_factor
                ops.integrator('DisplacementControl', self.ops_n_elem, 1, dU)

        def reset_analysis_options(disp_incr_factor):
            update_dU(disp_incr_factor)
            ops.algorithm('Newton')
            ops.test('NormUnbalance', 1e-3, 10)

        # region Run vertical load (time = LFV)
        ops.timeSeries('Linear', 100)
        ops.pattern('Plain', 200, 100)
        ops.load(self.ops_n_elem, 0, -1, 0)
        if self.gamma != 0:
            ops.load(1003, 0, -self.gamma, 0)
        ops.constraints('Plain')
        ops.numberer('RCM')
        ops.system('UmfPack')
        ops.test('NormUnbalance', 1e-3, 10)
        ops.algorithm('Newton')
        ops.integrator('LoadControl', config['P'] / config['num_steps_vertical'])
        ops.analysis('Static')

        # Define recorder
        def record():
            time = ops.getTime()
            section_strains = ops_get_section_strains(self)

            results.applied_axial_load.append(time)
            results.applied_horizontal_load.append(0.)
            results.maximum_abs_moment.append(ops_get_maximum_abs_moment(self))
            results.maximum_abs_disp.append(ops_get_maximum_abs_disp(self))
            results.lowest_eigenvalue.append(ops.eigen('-fullGenLapack', 1)[0])
            results.moment_at_top.append(ops.eleForce(self.ops_n_elem - 1, 6))
            results.moment_at_bottom.append(ops.eleForce(0, 3))
            results.maximum_concrete_compression_strain.append(section_strains[0])
            results.maximum_steel_strain.append(section_strains[1])
            if self.axis == 'x':
                results.curvature.append(section_strains[2])
            elif self.axis == 'y':
                results.curvature.append(section_strains[3])

        record()

        for i in range(config['num_steps_vertical']):
            ok = ops.analyze(1)

            if ok != 0:
                results.exit_message = 'Analysis Failed In Vertical Loading'
                warnings.warn('Analysis Failed In Vertical Loading')
                return results

            record()

            # Check for limits
            exit_message = check_analysis_limits(results, eigenvalue_limit=config['eigenvalue_limit'], 
                                            deformation_limit=config['deformation_limit'],
                                            concrete_strain_limit=config['concrete_strain_limit'], 
                                            steel_strain_limit=config['steel_strain_limit'],
                                            section_type=type(self.section).__name__)
            if exit_message:
                results.exit_message = exit_message
                break
            
        # endregion

        # region Run lateral load (time = LFH)
        ops.loadConst('-time', 0.0)
        ops.timeSeries('Linear', 101)
        ops.pattern('Plain', 201, 101)
        ops.load(self.ops_n_elem, 1, 0, 0)
        dU = self.length * config['disp_incr_factor']
        ops.integrator('DisplacementControl', self.ops_n_elem, 1, dU)
        ops.analysis('Static')

        # Define recorder
        def record():
            results.applied_axial_load.append(config['P'])
            section_strains = ops_get_section_strains(self)

            results.applied_horizontal_load.append(ops.getTime())
            results.maximum_abs_moment.append(ops_get_maximum_abs_moment(self))
            results.maximum_abs_disp.append(ops_get_maximum_abs_disp(self))
            results.lowest_eigenvalue.append(ops.eigen('-fullGenLapack', 1)[0])
            results.moment_at_top.append(ops.eleForce(self.ops_n_elem - 1, 6))
            results.moment_at_bottom.append(ops.eleForce(0, 3))
            results.maximum_concrete_compression_strain.append(section_strains[0])
            results.maximum_steel_strain.append(section_strains[1])
            if self.axis == 'x':
                results.curvature.append(section_strains[2])
            elif self.axis == 'y':
                results.curvature.append(section_strains[3])

        record()

        maximum_time = 0
        disp_incr_factor = config['disp_incr_factor']
        
        while True:
            ok = ops.analyze(1)

            if ok != 0 and config['try_smaller_steps']:
                if ok != 0 and config['try_smaller_steps']:
                    for div_factor in [1e1, 1e2, 1e3, 1e4, 1e5, 1e6]:
                        update_dU(disp_incr_factor, div_factor, analysis_type='nonproportional_limit_point')
                        ok = ops.analyze(1)
                        if ok == 0 and div_factor in [1e3, 1e4, 1e5, 1e6]:
                            disp_incr_factor /= 10
                            break
                        elif ok == 0:
                            break
                        else:
                            ok = try_analysis_options()
                            if ok == 0 and div_factor in [1e3, 1e4, 1e5, 1e6]:
                                disp_incr_factor /= 10
                                break
                            elif ok == 0:
                                break

            if ok != 0 and not config['try_smaller_steps']:
                ok = try_analysis_options()

            if ok == 0:
                reset_analysis_options(disp_incr_factor)

            elif ok != 0:
                results.exit_message = 'Analysis Failed'
                warnings.warn('Analysis Failed')
                break

            record()

            # Check for drop in applied load (time = the horizontal load factor)
            if config['percent_load_drop_limit'] is not None:
                current_time = ops.getTime()
                maximum_time = max(maximum_time, current_time)
                load_drop_limit = (1 - config['percent_load_drop_limit']) * maximum_time
                if current_time < load_drop_limit:
                    results.exit_message = 'Load Drop Limit Reached'
                    break

            # Check for limits
            exit_message = check_analysis_limits(results, eigenvalue_limit=config['eigenvalue_limit'], 
                                            deformation_limit=config['deformation_limit'],
                                            concrete_strain_limit=config['concrete_strain_limit'], 
                                            steel_strain_limit=config['steel_strain_limit'],
                                            section_type=type(self.section).__name__)
            if exit_message:
                results.exit_message = exit_message
                break
            
        return results


    def _run_nonproportional_analysis(self, config, results):
        """Run nonproportional analysis for sway column."""
        return self._run_sway_nonproportional_limit_point(config, results)


    @property
    def Cm(self):
        if self.k_top == 0 or self.k_bot == 0:
            Cm = 0.6
        elif self.k_top == self.k_bot:
            Cm = 0.2
        else:
            raise ValueError('Cm not implemented for unequal top and bottom stiffness')
        return Cm


    def run_AASHTO_interaction(self, EI_type, **kwargs):
        # beta_dns is the ratio of the maximum factored sustained axial load divided by
        # the total factored axial load associated with the same load combination
        # default is zero (i.e., short term loading)

        # Note that this function uses
        #   M1 to mean applied first-order moment
        #   M2 to mean internal second-order moment
        # this notation is different from what is used in AASHTO.
        num_points = kwargs.get('num_points', 10)
        section_factored = kwargs.get('section_factored', True)
        Pc_factor = kwargs.get('Pc_factor', 0.75)
        beta_dns = kwargs.get('beta_dns', 0)
        minimum_eccentricity = kwargs.get('minimum_eccentricity', False)

        # Get cross-sectional interaction diagram
        P_id, M_id, _ = self.section.section_interaction_2d(self.axis, 100, factored=section_factored,
                                                            only_compressive=True)
        id2d = InteractionDiagram2d(M_id, P_id, is_closed=False)

        if minimum_eccentricity:
            raise NotImplementedError('Minimum eccentricity not implemented')

        else:
            EIeff = self.section.EIeff(self.axis, EI_type, beta_dns, P=max(P_id), M=0, col=self)
            k = self.effective_length_factor(EIeff)
            Pc = pi ** 2 * EIeff / (k * self.length) ** 2
            buckling_load = Pc_factor * Pc

        P_list, M1_list, M2_list = [], [], []

        if buckling_load > max(P_id):
            # Buckling does not happend since the maximum axial strength is less than the lower bound buckling load
            P_list.append(max(P_id))
            M1_list.append(0)
            M2_list.append(0)

        else:
            # Buckling happens
            if EI_type.lower() in ['aci-a', 'aci-b']:
                P_list.append(buckling_load)
                M1_list.append(0)
                M2_list.append(id2d.find_x_given_y(buckling_load, 'pos'))

            else:
                def error(P):
                    P = P[0]
                    EIeff = self.section.EIeff(self.axis, EI_type, beta_dns, P=P, M=0, col=self)
                    k = self.effective_length_factor(EIeff)
                    Pc = pi ** 2 * EIeff / (k * self.length) ** 2
                    return P - Pc_factor * Pc

                # Find P such that error = zero
                Pguess = 0.9 * max(P_id)
                solution, _, ier, _ = fsolve(error, Pguess, full_output=True)
                if ier != 1:
                    Pguess = 0.1 * max(P_id)
                    solution, _, ier, _ = fsolve(error, Pguess, full_output=True)
                    if ier != 1:
                        raise Exception("Buckling load calculation did not converge.")
                buckling_load = solution[0]

                # Buckling
                error = []
                max_M_section = max(M_id)
                M2_trials = np.arange(0, max_M_section, max_M_section/1000)
                for M2 in M2_trials:
                    EIeff = self.section.EIeff(self.axis, EI_type, beta_dns, P=buckling_load, M=M2, col=self)
                    k = self.effective_length_factor(EIeff)
                    Pc = pi ** 2 * EIeff / (k * self.length) ** 2
                    error.append(buckling_load - Pc_factor * Pc)

                M2 = M2_trials[error.index(min(error, key=abs))]

                P_list.append(buckling_load)
                M1_list.append(0)
                M2_list.append(M2)

        # Loop axial linearly spaced axial loads witn non-proportional analyses
        for i in range(1, num_points):
            iP = 0.999 * P_list[0] * (num_points - i - 1) / (num_points - 1)
            iM2_section = id2d.find_x_given_y(iP, 'pos')

            EIeff = self.section.EIeff(self.axis, EI_type, beta_dns, P=iP, M=iM2_section, col=self)
            k = self.effective_length_factor(EIeff)
            Pc = pi ** 2 * EIeff / (k * self.length) ** 2

            iM1_list = [0]
            iM2_list = np.arange(0, iM2_section, iM2_section / 1000)

            for iM2 in iM2_list:
                EIeff = self.section.EIeff(self.axis, EI_type, beta_dns, P=iP, M=iM2, col=self)
                k = self.effective_length_factor(EIeff)
                Pc = pi ** 2 * EIeff / (k * self.length) ** 2
                if Pc_factor * Pc < iP:
                    break
                delta_s = max(1 / (1 - (iP) / (Pc_factor * Pc)), 1.0)
                iM1_list.append(iM2 / delta_s)

            iM1 = max(iM1_list)
            iM2 = iM2_list[iM1_list.index(iM1) -1]
            P_list.append(iP)
            M1_list.append(iM1)
            M2_list.append(iM2)
        results = {'P': np.array(P_list), 'M1': np.array(M1_list), 'M2': np.array(M2_list)}
        return results


    def effective_length_factor(self, EI):
        if self.effective_length_factor_override is not None:
            return self.effective_length_factor_override
    
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


    def calculated_EI_ops(self, P_list, M1_list, M2_ops, Pc_factor=1, G_bot=None, G_top=None) -> dict:
        """
            Back-calculate the effective flexural stiffness (EI) based on OpenSees results.

            Parameters:
                P_list (array-like): Array of axial loads.
                M1_list (array-like): Array of applied first-order moments.
                M2_list (array-like): Array of internal second-order moments.
                Pc_factor (float, optional): The factor to use in calculating the critical buckling load.
                                            Default is 1.

            Returns:
            dict: A dictionary containing back-calculated EI values for operational load conditions:
                - 'P': Array of axial loads
                - 'M1': Array of applied first-order moments
                - 'EI_ops': Array of back-calculated effective flexural stiffness values
                - 'EIgross': Gross flexural stiffness of the section
        """

        P_list = np.array(P_list)
        M1_list = np.array(M1_list)
        M2_ops = np.array(M2_ops)
        EIgross = self.section.EIgross(self.axis)

        M2_list = []
        EI_list_ops = []
        for P, M1 in zip(P_list, M1_list):
            M2 = np.interp(P, np.flip(P_list), np.flip(M2_ops))
            M2_list.append(M2)

            if M1 > M2:
                EI_list_ops.append(float("nan"))
                continue

            delta = M2 / M1
            Pc = P / (Pc_factor * (1 - 1 / delta))
            k = sidesway_uninhibited_effective_length_factor(G_bot, G_top)
            EI = Pc * (k * self.length / pi) ** 2
            if EI>EIgross:
                EI = EIgross
            EI_list_ops.append(EI)

        return {"P": np.array(P_list), "M1": np.array(M1_list), "M2":np.array(M2_list), "EI_ops": np.array(EI_list_ops), "EIgross": EIgross}


    def calculated_EI_design(self, P_list, M1_list, P_design, M2_design, section_factored=False, Pc_factor=1, G_bot=None, G_top=None) -> dict:
        """
            Back-calculate the effective flexural stiffness (EI) based on OpenSees and AASHTO values.

            Parameters:
                P_list (array-like): Array of axial loads.
                M1_list (array-like): Array of applied first-order moments.
                Pc_factor (float, optional): The factor to use in calculating the critical buckling load.
                                            Default is 1.

            Returns:
            dict: A dictionary containing back-calculated EI values for operational load conditions:
                - 'P': Array of axial loads
                - 'M1': Array of applied first-order moments
                - 'EI_AASHTO': Array of back-calculated effective flexural stiffness values
                - 'EIgross': Gross flexural stiffness of the section
        """
        P_list = np.array(P_list)
        M1_list = np.array(M1_list)
        P_design = np.array(P_design)
        M2_design = np.array(M2_design)

        EIgross = self.section.EIgross(self.axis)
        M2_list = []
        EI_list_AASHTO = []

        for P, M1 in zip(P_list, M1_list):
            M2 = np.interp(P, np.flip(P_design), np.flip(M2_design))
            M2_list.append(M2)

            if P < min(P_design) or P > max(P_design) or M1 > M2:
                EI_list_AASHTO.append(float("nan"))
                continue

            delta = M2 / M1
            Pc = P / (Pc_factor * (1-1/delta))
            k = sidesway_uninhibited_effective_length_factor(G_bot, G_top)
            EI = Pc * (k * self.length / pi) ** 2
            if EI>EIgross:
                EI = EIgross
            EI_list_AASHTO.append(EI)
        return {"P": np.array(P_list), "M1": np.array(M1_list), "M2":np.array(M2_list), "EI_AASHTO": np.array(EI_list_AASHTO),
                "EIgross": EIgross}


    def _extract_analysis_config(self, **kwargs):
        """Override to adjust default values for sway analysis."""
        config = super()._extract_analysis_config(**kwargs)
        # Adjust defaults for sway analysis if not explicitly provided
        if 'num_steps_vertical' not in kwargs:
            config['num_steps_vertical'] = 10
        if 'disp_incr_factor' not in kwargs:
            config['disp_incr_factor'] = 5e-05
        return config
    
    
    def _find_limit_point(self, results, config, analysis_type):
        """Find and set limit point values with sway-specific logic."""
        if 'Analysis Failed' in results.exit_message:
            if analysis_type.lower() == 'proportional_limit_point':
                ind, x = find_limit_point_in_list(results.applied_axial_load, max(results.applied_axial_load))
            elif analysis_type.lower() == 'nonproportional_limit_point':
                ind, x = find_limit_point_in_list(results.applied_horizontal_load, max(results.applied_horizontal_load))
            warnings.warn('Analysis failed')
        elif 'Eigenvalue Limit' in results.exit_message:
            ind, x = find_limit_point_in_list(results.lowest_eigenvalue, config['eigenvalue_limit'])
        elif 'Extreme Compressive Concrete Fiber Strain Limit Reached' in results.exit_message:
            ind, x = find_limit_point_in_list(results.maximum_concrete_compression_strain, config['concrete_strain_limit'])
        elif 'Extreme Steel Fiber Strain Limit Reached' in results.exit_message:
            ind, x = find_limit_point_in_list(results.maximum_steel_strain, config['steel_strain_limit'])
        elif 'Deformation Limit Reached' in results.exit_message:
            ind, x = find_limit_point_in_list(results.maximum_abs_disp, config['deformation_limit'])
        elif 'Load Drop Limit Reached' in results.exit_message:
            if analysis_type.lower() == 'proportional_limit_point':
                ind, x = find_limit_point_in_list(results.applied_axial_load, max(results.applied_axial_load))
            elif analysis_type.lower() == 'nonproportional_limit_point':
                ind, x = find_limit_point_in_list(results.applied_horizontal_load, max(results.applied_horizontal_load))
        else:
            raise Exception('Unknown limit point')

        self._set_limit_point_values(results, ind, x)
        
    
