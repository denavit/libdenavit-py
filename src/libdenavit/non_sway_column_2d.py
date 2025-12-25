from math import pi, sin, log10
from libdenavit import opensees as ops
from libdenavit import find_limit_point_in_list, interpolate_list, InteractionDiagram2d, CrossSection2d
from libdenavit.OpenSees import AnalysisResults
import numpy as np
import matplotlib.pyplot as plt
import warnings
from scipy.optimize import fsolve
import io
import sys
from libdenavit.analysis_helpers import try_analysis_options, ops_get_section_strains, ops_get_maximum_abs_moment, ops_get_maximum_abs_disp, check_analysis_limits
from libdenavit.column_2d import Column2d


class NonSwayColumn2d(Column2d):
    def __init__(self, section, length, et, eb, **kwargs):
        """
            Represents a non-sway 2D column

            This class defines a non-sway column element with physical parameters such as section properties,
            length, and material properties. It also allows customization of analysis options.

            Parameters:
                section: The section object representing the cross-sectional properties.
                length: The length of the entire column.
                et: The eccentricity at the top of the column.
                eb: The eccentricity at the bottom of the column.
                kwargs: Additional keyword arguments for customization.
                          dxo (float or None, optional): Initial geometric imperfection.
                                                         Default is 0.0. If None, then no imperfection is included.
                          axis (str, optional): Axis. Default is None.
                          n_elem (int, optional): Number of elements for OpenSees analysis. Default is 6.
                          element_type (str, optional): Type of OpenSees element. Default is 'mixedBeamColumn'.
                          ops_geom_transf_type (str, optional): OpenSees geometric transformation type. Default is 'Corotational'.
                          ops_integration_points (int, optional): Number of integration points for OpenSees analysis. Default is 3.
        """

        super().__init__(section, length, **kwargs)
        
        # Physical parameters
        self.et = et
        self.eb = eb
        
        # Specific defaults for non-sway analysis
        self.creep = kwargs.get('creep', False)
        self.P_sus = kwargs.get('P_sus', 0.0)
        self.t_sus = kwargs.get('t_sus', 10000)


    def _get_default_deformation_limit(self):
        """Override default deformation limit for non-sway columns."""
        return 0.1 * self.length


    def _initialize_results(self):
        """Initialize analysis results object with required attributes."""
        # Get base attributes and add NonSway-specific ones
        results = super()._initialize_results()
        # Add NonSway-specific attributes
        for attr in ['applied_moment_top', 'applied_moment_bot']:
            setattr(results, attr, [])
        return results


    def _set_limit_point_values(self, results, ind, x):
        """Override to include moment values."""
        super()._set_limit_point_values(results, ind, x)
        results.applied_moment_top_at_limit_point = interpolate_list(results.applied_moment_top, ind, x)
        results.applied_moment_bot_at_limit_point = interpolate_list(results.applied_moment_bot, ind, x)


    def _run_proportional_analysis(self, config, results):
        """Run proportional analysis for non-sway column."""
        if self.creep:
            return self._run_ops_proportional_with_creep(config, results)
        else:
            return self._run_ops_proportional_no_creep(config, results)


    def _run_nonproportional_analysis(self, config, results):
        """Run nonproportional analysis for non-sway column."""
        return self._run_ops_nonproportional_limit_point(config, results)


    def run_ops_interaction(self, **kwargs):
    
        # Parse keyword arguments
        section_id = kwargs.get('section_id', 1)
        section_args = kwargs.get('section_args', [])
        section_kwargs = kwargs.get('section_kwargs', {})
        num_points = kwargs.get('num_points', 10)
        prop_disp_incr_factor = kwargs.get('prop_disp_incr_factor', 1e-6)
        nonprop_disp_incr_factor = kwargs.get('nonprop_disp_incr_factor', 1e-5)
        section_load_factor = kwargs.get('section_load_factor', 1e-1)
        plot_load_deformation = kwargs.get('plot_load_deformation', False)
        full_results = kwargs.get('full_results', False)

        if plot_load_deformation:
            fig_at_step, ax_at_step = plt.subplots(2, 1, figsize=(10, 6), gridspec_kw={'height_ratios': [3, 1]})

        # Run one axial load only analyis to determine maximum axial strength
        results = self.run_ops_analysis('proportional_limit_point', e=0, section_id=section_id,
                                        section_args=section_args, section_kwargs=section_kwargs,
                                        disp_incr_factor=prop_disp_incr_factor)
        P = [results.applied_axial_load_at_limit_point]
        M1 = [0]
        M2 = [results.maximum_abs_moment_at_limit_point]
        if full_results:
            M1t_path = [results.applied_moment_top]
            M1b_path = [results.applied_moment_bot]
            M2_path = [results.maximum_abs_moment]
            disp_path = [results.maximum_abs_disp]


        exit_message = [results.exit_message]
        if P is np.nan or P == [np.nan]:
            raise ValueError('Analysis failed at axial only loading')

        # Loop axial linearly spaced axial loads witn non-proportional analyses
        for i in range(1,num_points):
            iP = P[0] * (num_points-1-i) / (num_points-1)
            if iP == 0:
                cross_section = CrossSection2d(self.section, self.axis)
                results = cross_section.run_ops_analysis('nonproportional_limit_point', P=0,
                                                         section_id=section_id, section_args=section_args,
                                                         load_incr_factor=section_load_factor)
                P.append(iP)
                M1.append(results.maximum_abs_moment_at_limit_point)
                M2.append(results.maximum_abs_moment_at_limit_point)
                if full_results:
                    M1t_path.append(results.maximum_abs_moment)
                    M1b_path.append(results.maximum_abs_moment)
                    M2_path.append(results.maximum_abs_moment)
                    disp_path.append([0]*len(results.maximum_abs_moment))

                exit_message.append(results.exit_message)
            else:
                results = self.run_ops_analysis('nonproportional_limit_point', P=iP,
                                                section_id=section_id, section_args=section_args,
                                                disp_incr_factor=nonprop_disp_incr_factor)
                P.append(iP)
                M1.append(results.applied_moment_top_at_limit_point)
                M2.append(results.maximum_abs_moment_at_limit_point)
                if full_results:
                    M1t_path.append(results.applied_moment_top)
                    M1b_path.append(results.applied_moment_bot)
                    M2_path.append(results.maximum_abs_moment)
                    disp_path.append(results.maximum_abs_disp)

                exit_message.append(results.exit_message)

            if plot_load_deformation:
                if iP==0:
                    print(f'{results.maximum_abs_moment=}')
                else:
                    ax_at_step[0].plot(results.maximum_abs_disp, results.applied_moment_top, '-o', label=f'{iP:,.0f}', markersize=5)
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
                    'M1t_path': M1t_path, 'M1b_path': M1b_path, 'M2_path': M2_path, 'disp_path': disp_path}
        else:
            return {'P': np.array(P), 'M1': np.array(M1), 'M2': np.array(M2), 'exit_message': exit_message}


    def run_ops_interaction_proportional(self, e_list, **kwargs):
        results = [self.run_ops_analysis('proportional_limit_point', e=e, **kwargs) for e
                   in e_list]
        P = [result.applied_axial_load_at_limit_point for result in results]
        M1 = [result.applied_moment_top_at_limit_point for result in results]
        M2 = [result.maximum_abs_moment_at_limit_point for result in results]

        return {'P': np.array(P), 'M1': np.array(M1), 'M2': np.array(M2)}


    def run_AASHTO_interaction(self, EI_type, **kwargs):
        """
        Perform AASHTO LRFD-based interaction analysis for the column.

            Parameters:
                EI_type (str): The type of effective flexural stiffness of member to use in the analysis.
                num_points (int, optional): The number of points to use in the interaction diagram. Default is 10.
                section_factored (bool, optional): Whether to use factored section properties. Default is True.
                Pc_factor (float, optional): The factor to use in calculating the buckling load. Default is 0.75.
                betadns (float, optional): The ratio of the maximum factored sustained axial load to the total factored axial load
                                        for the same load combination. Default is 0 (short-term loading).
                minimum_eccentricity (bool, optional): Whether to consider minimum eccentricity in the analysis. Default is False.

            Note:
              This function uses the notation:
                - M1 to represent the applied first-order moment.
                - M2 to represent the internal second-order moment.
                - This notation differs from the notation used in AASHTO.

            Returns:
            dict: A dictionary containing interaction diagram data:
                - 'P': Array of axial loads
                - 'M1': Array of applied first-order moments
                - 'M2': Array of internal second-order moments
        """
        num_points = kwargs.get('num_points', 10)
        section_factored = kwargs.get('section_factored', True)
        Pc_factor = kwargs.get('Pc_factor', 0.75)
        beta_dns = kwargs.get('beta_dns', 0)
        minimum_eccentricity = kwargs.get('minimum_eccentricity', False)

        # Get cross-sectional interaction diagram
        P_id, M_id, _ = self.section.section_interaction_2d(self.axis, 100, factored=section_factored,
                                                            only_compressive=True)
        id2d = InteractionDiagram2d(M_id, P_id, is_closed=False)

        k = 1  # Effective length factor (always one for this non-sway column)

        # Run one axial load only analysis to determine maximum axial strength
        if minimum_eccentricity:
            raise NotImplementedError('Minimum eccentricity not implemented')

        # Compute buckling load based on maximum axial strength (this should be a lower bound)
        else:
            EIeff = self.section.EIeff(self.axis, EI_type, beta_dns, P=max(P_id), M=0, col=self)
            Pc = pi ** 2 * EIeff / (k * self.length) ** 2
            buckling_load = Pc_factor * Pc

        P_list, M1_list, M2_list, EIeff_list = [], [], [], []

        if buckling_load > max(P_id):
            # Buckling does not happen since the maximum axial strength is less than the lower bound buckling load
            P_list.append(max(P_id))
            M1_list.append(0)
            M2_list.append(0)
            EIeff_list.append(EIeff)

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
                M2_trials = np.arange(0, max_M_section, max_M_section / 1000)
                for M2 in M2_trials:
                    EIeff = self.section.EIeff(self.axis, EI_type, beta_dns, P=buckling_load, M=M2, col=self)
                    Pc = pi ** 2 * EIeff / (k * self.length) ** 2
                    error.append(buckling_load - Pc_factor * Pc)

                M2 = M2_trials[error.index(min(error))]

                P_list.append(buckling_load)
                M1_list.append(0)
                M2_list.append(M2)
                EIeff_list.append(self.section.EIeff(self.axis, EI_type, beta_dns, P=buckling_load, M=M2, col=self))

        # Loop axial linearly spaced axial loads with non-proportional analyses
        for i in range(1, num_points):
            iP = 0.999 * P_list[0] * (num_points - i - 1) / (num_points - 1)
            if EI_type.lower() in ['aci-a', 'aci-b']:
                iM2 = id2d.find_x_given_y(iP, 'pos')
                EIeff = self.section.EIeff(self.axis, EI_type, beta_dns)
                k = 1
                Pc = pi ** 2 * EIeff / (k * self.length) ** 2
                delta = max(self.Cm / (1 - (iP) / (Pc_factor * Pc)), 1.0)
                iM1 = iM2 / delta
            else:
                iM2_section = id2d.find_x_given_y(iP, 'pos')
                k = 1  # Effective length factor (always one for this non-sway column)

                iM1_list = [0]
                iM2_list = np.arange(0, iM2_section, iM2_section/1000)
                for iM2 in iM2_list:
                    EIeff = self.section.EIeff(self.axis, EI_type, beta_dns, P=iP, M=iM2, col=self)
                    Pc = pi ** 2 * EIeff / (k * self.length) ** 2
                    if Pc_factor * Pc < iP:
                        break

                    delta = max(self.Cm / (1 - (iP) / (Pc_factor * Pc)), 1.0)

                    iM1_list.append(iM2 / delta)

                iM1 = max(iM1_list)
                iM2 = iM2_list[iM1_list.index(iM1)-1]
            P_list.append(iP)
            M1_list.append(iM1)
            M2_list.append(iM2)
            EIeff_list.append(self.section.EIeff(self.axis, EI_type, beta_dns, P=iP, M=iM2, col=self))
        results = {'P':np.array(P_list),'M1':np.array(M1_list),'M2':np.array(M2_list), 'EIeff':np.array(EIeff_list)}
        return results

    @property
    def Cm(self) -> float:
        return 0.6 + 0.4 * min(self.et, self.eb, key=abs) / max(self.et, self.eb, key=abs)


    def calculated_EI_ops(self, P_list, M1_list, M2_ops, Pc_factor=1) -> dict:
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
            Pc = P / (1-self.Cm/delta) / Pc_factor
            k = 1  # Effective length factor (always one for this non-sway column)
            EI = Pc * (k * self.length / pi) ** 2
            if EI>EIgross:
                EI = EIgross
            EI_list_ops.append(EI)

        return {"P":np.array(P_list), "M1":np.array(M1_list), "M2":np.array(M2_list), "Calculated EI":np.array(EI_list_ops),
                "EIgross":EIgross}


    def calculated_EI_design(self, P_list, M1_list, P_design, M2_design, section_factored=False, Pc_factor=1) -> dict:
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

            if M1>M2:
                EI_list_AASHTO.append(EIgross)
                continue

            delta = M2 / M1
            Pc = P / (1 - self.Cm / delta) / Pc_factor
            k = 1  # Effective length factor (always one for this non-sway column)
            EI = Pc * (k * self.length / pi) ** 2
            if EI>EIgross:
                EI = EIgross
            EI_list_AASHTO.append(EI)
            
        return {"P": np.array(P_list), "M1": np.array(M1_list), "M2":np.array(M2_list), "Calculated EI": np.array(EI_list_AASHTO),
                "EIgross": EIgross}
        
    
    @property
    def ops_mid_node(self):
        if self.ops_n_elem % 2 == 0:
            return self.ops_n_elem // 2
        raise ValueError(f'Number of elements should be even {self.ops_n_elem = }')


    def build_ops_model(self, section_id, section_args, section_kwargs, **kwargs):
        """
           Build the OpenSees finite element model for the non-sway 2D column.

           This method constructs the finite element model in OpenSees for the non-sway 2D column element.

           Parameters:
               section_id: An integer id for the section
               section_args: Positional arguments for building the section using OpenSees via section.build_ops_fiber_section().
                             (For RC sections the args are: section_id, start_material_id, steel_mat_type, conc_mat_type, nfy, nfx)
               section_kwargs: Keword arguments for building the section using OpenSees via section.build_ops_fiber_section().
                               (For RC sections, no kwargs are necessary).
               kwargs: Additional keyword arguments.
                         start_node_fixity (tuple, optional): Fixity conditions at the start node. Default is (1, 1, 0).
                         end_node_fixity (tuple, optional): Fixity conditions at the end node. Default is (1, 0, 0).

           Returns:
               None
       """

        # region Extract kwargs
        creep_props_dict = kwargs.get('creep_props_dict', dict())
        shrinkage_props_dict = kwargs.get('shrinkage_props_dict', dict())
        start_node_fixity = kwargs.get('start_node_fixity', (1, 1, 0))
        end_node_fixity = kwargs.get('end_node_fixity', (1, 0, 0))
        # endregion

        # region Build OpenSees model
        ops.wipe()
        ops.model('basic', '-ndm', 2, '-ndf', 3)
        # endregion

        # region Define Nodes and Fixities and Geometric Transformation
        for index in range(self.ops_n_elem + 1):
            if isinstance(self.dxo, (int, float)):
                x = sin(index / self.ops_n_elem * pi) * self.dxo
            elif self.dxo == None:
                x = 0.
            else:
                raise ValueError(f'Unknown value of dxo ({self.dxo})')
            y = index / self.ops_n_elem * self.length
            ops.node(index, x, y)
            ops.mass(index, 1, 1, 1)

        ops.fix(0, *start_node_fixity)
        ops.fix(self.ops_n_elem, *end_node_fixity)

        ops.geomTransf(self.ops_geom_transf_type, 100)
        # endregion and

        # region Define Fiber Section
        if type(self.section).__name__ == "RC":
            self.section.build_ops_fiber_section(section_id, *section_args, **section_kwargs, axis=self.axis,
                                                 creep=self.creep, creep_props_dict=creep_props_dict,
                                                 shrinkage_props_dict=shrinkage_props_dict)
        elif type(self.section).__name__ == "CCFT":
            self.section.build_ops_fiber_section(section_id, *section_args, **section_kwargs, creep=self.creep, axis=self.axis)
        elif type(self.section).__name__ == "I_shape":
            self.section.build_ops_fiber_section(section_id, *section_args, **section_kwargs, axis=self.axis)
        else:
            raise ValueError(f'Unknown cross section type {type(self.section).__name__}')
        # endregion

        ops.beamIntegration("Lobatto", 1, 1, self.ops_integration_points)

        for index in range(self.ops_n_elem):
            if self.ops_element_type == 'elasticBeamColumn':
                raise NotImplementedError('Elastic beam column element is not implemented for this class.')
            else:
                ops.element(self.ops_element_type, index, index, index + 1, 100, 1)


    def _set_limit_point_values(self, results, ind, x):
        """Override to include moment values."""
        super()._set_limit_point_values(results, ind, x)
        results.applied_moment_top_at_limit_point = interpolate_list(results.applied_moment_top, ind, x)
        results.applied_moment_bot_at_limit_point = interpolate_list(results.applied_moment_bot, ind, x)

    
    def _initialize_results(self):
        """Initialize analysis results object with required attributes."""
        # Get base attributes and add NonSway-specific ones
        results = super()._initialize_results()
        # Add NonSway-specific attributes
        for attr in ['applied_moment_top', 'applied_moment_bot']:
            setattr(results, attr, [])
        return results


    def _run_ops_proportional_no_creep(self, config, results):
        """Run proportional limit point analysis without creep."""
        # time = LFV
        ops.timeSeries('Linear', 100)
        ops.pattern('Plain', 200, 100)

        sgn_et = int(np.sign(self.et))
        sgn_eb = int(np.sign(self.eb))
        
        if sgn_et != sgn_eb and (sgn_eb != 0 and sgn_et != 0):
            ecc_sign = sgn_et if abs(self.et) >= abs(self.eb) else sgn_eb
        else:
            ecc_sign = sgn_et
        
        ops.constraints('Plain')
        ops.numberer('RCM')
        ops.system('UmfPack')
        ops.test('NormUnbalance', 1e-3, 10)
        ops.algorithm('Newton')
        

        if config['e'] == 0.0:
            # axial-only: vertical load only, no moments
            ops.load(self.ops_n_elem, 0, -1, 0)  # reference vertical load
            ops.load(0, 0, 0, 0.0)                 # no moment
            dF = 1.0 / max(1, config['num_steps_vertical'])
            ops.integrator('LoadControl', dF)
        else:
            ops.load(self.ops_n_elem, 0, -1, self.et * config['e'] * ecc_sign)
            ops.load(0, 0, 0, -self.eb * config['e'] * ecc_sign)
            
            if sgn_et != sgn_eb and (sgn_eb != 0 and sgn_et != 0):
                if max(self.et, self.eb, key=abs) == self.et:
                    dof = 3 * self.ops_n_elem // 4
                else:
                    dof = 1 * self.ops_n_elem // 4
                dU = self.length * config['disp_incr_factor'] / 2
                ops.integrator('DisplacementControl', dof, 1, dU)
            else:
                dU = self.length * config['disp_incr_factor']
                ops.integrator('DisplacementControl', self.ops_mid_node, 1, dU)

        ops.analysis('Static', '-noWarnings')

        # Define recorder
        def record():
            time = ops.getTime()
            section_strains = ops_get_section_strains(self)

            results.applied_axial_load.append(time)
            results.applied_moment_top.append(self.et * config['e'] * time * ecc_sign)
            results.applied_moment_bot.append(-self.eb * config['e'] * time * ecc_sign)
            results.maximum_abs_moment.append(ops_get_maximum_abs_moment(self))
            results.maximum_abs_disp.append(ops_get_maximum_abs_disp(self))
            results.lowest_eigenvalue.append(ops.eigen('-fullGenLapack', 1)[0])
            if type(self.section).__name__ == "RC":
                results.maximum_concrete_compression_strain.append(section_strains[0])
                results.maximum_steel_strain.append(section_strains[1])
            elif type(self.section).__name__ == "I_shape":
                results.maximum_compression_strain.append(section_strains[0])
                results.maximum_tensile_strain.append(section_strains[1])

            if self.axis == 'x':
                results.curvature.append(section_strains[2])
            elif self.axis == 'y':
                results.curvature.append(section_strains[3])
            else:
                raise ValueError(f'The value of axis ({self.axis}) is not supported.')

        def update_dU(disp_incr_factor, div_factor=1):
            sgn_et = int(np.sign(self.et))
            sgn_eb = int(np.sign(self.eb))
            if sgn_et != sgn_eb and (sgn_eb != 0 and sgn_et != 0):
                if max(self.et, self.eb, key=abs) == self.et:
                    dof = 3 * self.ops_n_elem // 4
                else:
                    dof = 1 * self.ops_n_elem // 4
                dU = self.length * disp_incr_factor / 2 / div_factor
                ops.integrator('DisplacementControl', dof, 1, dU)
            else:
                dU = self.length * disp_incr_factor / div_factor
                ops.integrator('DisplacementControl', self.ops_mid_node, 1, dU)

        def reset_analysis_options(disp_incr_factor):
            update_dU(disp_incr_factor)
            ops.algorithm('Newton')
            ops.test('NormUnbalance', 1e-3, 10)

        record()
        

        maximum_applied_axial_load = 0.
        disp_incr_factor = config['disp_incr_factor']
        
        
        while True:
            ok = ops.analyze(1)

            if ok != 0 and config['try_smaller_steps']:
                for div_factor in [1e1, 1e2, 1e3, 1e4, 1e5, 1e6]:
                    update_dU(disp_incr_factor, div_factor)
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

            if ok == 0 and config['try_smaller_steps']:
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
                if current_applied_axial_load < (1 - config['percent_load_drop_limit']) * maximum_applied_axial_load:
                    results.exit_message = 'Load Drop Limit Reached'
                    break

            # Check for limits
            exit_message = check_analysis_limits(
                    results,
                    eigenvalue_limit=config['eigenvalue_limit'],
                    deformation_limit=config['deformation_limit'],
                    concrete_strain_limit=config['concrete_strain_limit'],
                    steel_strain_limit=config['steel_strain_limit'],
                    section_type=type(self.section).__name__,
                )
            if exit_message:
                results.exit_message = exit_message
                break

        return results


    def _run_ops_proportional_with_creep(self, config, results):
        """Run proportional limit point analysis with creep."""
        # region Determine the sign of the eccentricity
        sgn_et = int(np.sign(self.et))
        sgn_eb = int(np.sign(self.eb))
        if sgn_et != sgn_eb:
            if max(self.et, self.eb, key=abs) == self.et:
                ecc_sign = sgn_et
            else:
                ecc_sign = sgn_eb
        else:
            ecc_sign = sgn_et
        # endregion

        # region Define recorder
        def record(lam=0):
            section_strains = ops_get_section_strains(self)

            # Backup the original stderr
            original_stderr = sys.stderr
            try:
                # Redirect stderr to nowhere
                sys.stderr = io.StringIO()
                time = ops.getLoadFactor(200) + ops.getLoadFactor(2000)
            except:
                try:
                    time = ops.getLoadFactor(200)
                except:
                    time = 0
            finally:
                # Restore stderr
                sys.stderr = original_stderr

            results.applied_axial_load.append(time + lam)
            results.applied_moment_top.append(self.et * config['e'] * (time + lam) * ecc_sign)
            results.applied_moment_bot.append(-self.eb * config['e'] * (time + lam) * ecc_sign)
            results.maximum_abs_moment.append(ops_get_maximum_abs_moment(self))
            results.maximum_abs_disp.append(ops_get_maximum_abs_disp(self))
            results.lowest_eigenvalue.append(ops.eigen('-fullGenLapack', 1)[0])
            if type(self.section).__name__ == "RC":
                results.maximum_concrete_compression_strain.append(section_strains[0])
                results.maximum_steel_strain.append(section_strains[1])
            elif type(self.section).__name__ == "I_shape":
                results.maximum_compression_strain.append(section_strains[0])
                results.maximum_tensile_strain.append(section_strains[1])

            if self.axis == 'x':
                results.curvature.append(section_strains[2])
            elif self.axis == 'y':
                results.curvature.append(section_strains[3])
            else:
                raise ValueError(f'The value of axis ({self.axis}) is not supported.')
        # endregion

        # region Do one analysis with no load
        ops.setCreep(1)
        ops.setTime(self.section.Tcr)
        
        ops.timeSeries('Constant', 1)
        ops.pattern('Plain', 1, 1)

        ops.integrator('LoadControl', 0)
        ops.system('Umfpack')
        ops.test('NormUnbalance', 1e-3, 20)
        ops.algorithm('ModifiedNewton')
        ops.constraints('Plain')
        ops.numberer('RCM')
        ops.analysis('Static')
        
        ok = ops.analyze(1)
        record()
        
        # endregion

        # region Run the sustained load phase
        load_step_for_sustained = 1000
        
        ops.setCreep(0)
        ops.setTime(0)
        ops.wipeAnalysis()
        ops.system('Umfpack')
        ops.test('NormUnbalance', 1e-3, 20)
        ops.algorithm('ModifiedNewton')
        ops.constraints('Plain')
        ops.numberer('RCM')
        ops.integrator('LoadControl', self.P_sus/load_step_for_sustained)
        ops.analysis('Static','-noWarnings')
        
        t = self.section.Tcr
        tfinish = self.t_sus + self.section.Tcr

        ops.timeSeries('Linear', 100)
        ops.pattern('Plain', 200, 100, '-factor', 1)
        ops.load(self.ops_n_elem, 0, -1, self.et * config['e'] * ecc_sign)
        ops.load(0, 0, 0, -self.eb * config['e'] * ecc_sign)
        
        while t < tfinish:
            # Apply sustained load at Tcr
            if t == self.section.Tcr:
                for i in np.arange(1, load_step_for_sustained+1, 1):
                                        
                    ok = ops.analyze(1)
                    record()
                    
                    if ok < 0:
                        print(f'Analysis failed before full sustained load is reached, load: {i/load_step_for_sustained*self.P_sus}')
                        results.exit_message = 'Analysis failed before full sustained load is reached'
                        return results

                    if config['deformation_limit'] is not None:
                        if results.maximum_abs_disp[-1] > config['deformation_limit']:
                            print(f'Analysis failed before full sustained load is reached, load: {i/load_step_for_sustained*self.P_sus}')
                            results.exit_message = 'Analysis failed before full sustained load is reached'
                            return results

                ops.loadConst('-time', self.P_sus)
                ops.integrator('LoadControl', 0)
                                
                ops.setCreep(1)
                
            ops.setTime(t)
            ok = ops.analyze(1)
            # if t == self.section.Tcr and ok < 0:
                # ok = ops.analyze(1)
            
            record()
            
            if ok < 0:
                if t == self.section.Tcr:
                    print(f'Analysis failed on the first step of maintaining sustained load')
                    results.exit_message = 'Analysis failed on the first step of maintaining sustained load'
                    return results
                else:
                    print(f'Analysis failed while maintaining sustained load, time: {t}')
                    results.exit_message = f'Analysis failed while maintaining sustained load: {t}'
                return results

            if config['deformation_limit'] is not None:
                if results.maximum_abs_disp[-1] > config['deformation_limit']:
                    print(f'Analysis failed while maintaining sustained load: {t}')
                    results.exit_message = 'Deformation Limit Reached'
                    return results

            logt0 = log10(t)
            logt1 = logt0 + 0.01
            t1 = 10 ** logt1
            t = t1
            
        # endregion


        def update_dU(disp_incr_factor, div_factor=1):
            if config['e'] == 0.0:
                # axial-only: shrink load step, stay in LoadControl
                dF = (1.0 / max(1, config['num_steps_vertical'])) / div_factor
                ops.integrator('LoadControl', dF)
            else:
                # bending present: DisplacementControl as before
                sgn_et = int(np.sign(self.et))
                sgn_eb = int(np.sign(self.eb))
                if sgn_et != sgn_eb and (sgn_eb != 0 and sgn_et != 0):
                    dof = 3 * self.ops_n_elem // 4 if abs(self.et) >= abs(self.eb) else 1 * self.ops_n_elem // 4
                    dU = self.length * disp_incr_factor / 2 / div_factor
                    ops.integrator('DisplacementControl', dof, 1, dU)
                else:
                    dU = self.length * disp_incr_factor / div_factor
                    ops.integrator('DisplacementControl', self.ops_mid_node, 1, dU)

        def reset_analysis_options(disp_incr_factor):
            update_dU(disp_incr_factor)
            ops.algorithm('Newton')
            ops.test('NormUnbalance', 1e-3, 10)

        # region run final loading phase
        ops.setCreep(0)
        ops.setTime(self.section.Tcr+self.t_sus+1)
        ops.integrator('LoadControl', 0)
        ops.analyze(1)

        record()

        ops.timeSeries('Ramp', 1000, self.section.Tcr+self.t_sus+1, self.section.Tcr+self.t_sus+1, '-smooth', 0, '-factor', 100)
        ops.pattern('Plain', 2000, 1000)

        ops.wipeAnalysis()
        if config['e'] == 0:
            dF = self.P_sus / config['num_steps_vertical']
            ops.load(self.ops_n_elem, 0, -1, 0)
            ops.load(0, 0, 0, 0)
            ops.integrator('LoadControl', dF)
        elif sgn_et != sgn_eb:
            if max(self.et, self.eb, key=abs) == self.et:
                dof = 3 * self.ops_n_elem // 4
                ecc_sign = sgn_et
            else:
                dof = 1 * self.ops_n_elem // 4
                ecc_sign = sgn_eb
            dU = self.length * config['disp_incr_factor'] / 20
            ops.load(self.ops_n_elem, 0, -1, self.et * config['e'] * ecc_sign)
            ops.load(0, 0, 0, -self.eb * config['e'] * ecc_sign)
            ops.integrator('DisplacementControl', dof, 1, dU)
        else:
            ecc_sign = sgn_et
            dU = self.length * config['disp_incr_factor'] / 10
            ops.load(self.ops_n_elem, 0, -1, self.et * config['e'] * ecc_sign)
            ops.load(0, 0, 0, -self.eb * config['e'] * ecc_sign)
            ops.integrator('DisplacementControl', self.ops_mid_node, 1, dU)
            
            # dF = self.P_sus / num_steps_vertical
            # ops.load(self.ops_n_elem, 0, -1, 0)
            # ops.load(0, 0, 0, 0)
            # ops.integrator('LoadControl', dF)
        
        
        ops.constraints('Plain')
        ops.numberer('RCM')
        ops.system('UmfPack')
        ops.test('NormUnbalance', 1e-3, 20)
        ops.algorithm('ModifiedNewton')
        ops.analysis('Static', '-noWarnings')
        ops.analyze(1)

        maximum_applied_axial_load = 0.
        disp_incr_factor = config['disp_incr_factor']
        
        while True:
            ok = ops.analyze(1)
            
            # print(ops.getTime())
            # print(ops.getLoadFactor(200))
            # print(ops.getLoadFactor(2000))
            if ok != 0 and config['try_smaller_steps']:
                for div_factor in [1e1, 1e2, 1e3, 1e4, 1e5, 1e6]:
                    update_dU(disp_incr_factor, div_factor)
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

            if ok == 0 and config['try_smaller_steps']:
                reset_analysis_options(disp_incr_factor)

            elif ok != 0:
                results.exit_message = 'Analysis Failed'
                warnings.warn('Analysis Failed')
                break

            record()

            # region Check the exit conditions
            # Check for drop in applied load
            if config['percent_load_drop_limit'] is not None:
                current_applied_axial_load = results.applied_axial_load[-1]
                maximum_applied_axial_load = max(maximum_applied_axial_load, current_applied_axial_load)
                if current_applied_axial_load < (1 - config['percent_load_drop_limit']) * maximum_applied_axial_load:
                    results.exit_message = 'Load Drop Limit Reached'
                    break
            
            # Check for limits
            exit_message = check_analysis_limits(
                    results,
                    eigenvalue_limit=config['eigenvalue_limit'],
                    deformation_limit=config['deformation_limit'],
                    concrete_strain_limit=config['concrete_strain_limit'],
                    steel_strain_limit=config['steel_strain_limit'],
                    section_type=type(self.section).__name__,
                )
            if exit_message:
                results.exit_message = exit_message
                break

        # endregion
        return results

    
    def _run_ops_nonproportional_limit_point(self, config, results):
        """Run nonproportional limit point analysis."""
        def update_dU(disp_incr_factor, div_factor=1):
            sgn_et = int(np.sign(self.et))
            sgn_eb = int(np.sign(self.eb))
            if sgn_et != sgn_eb and (sgn_eb != 0 and sgn_et != 0):
                if max(self.et, self.eb, key=abs) == self.et:
                    dof = 3 * self.ops_n_elem // 4
                else:
                    dof = 1 * self.ops_n_elem // 4
                dU = self.length * disp_incr_factor / 2 / div_factor
                ops.integrator('DisplacementControl', dof, 1, dU)
            else:
                dU = self.length * disp_incr_factor / div_factor
                ops.integrator('DisplacementControl', self.ops_mid_node, 1, dU)

        def reset_analysis_options(disp_incr_factor):
            update_dU(disp_incr_factor)
            ops.algorithm('Newton')
            ops.test('NormUnbalance', 1e-3, 10)

        # region Run vertical load (time = LFV)
        ops.timeSeries('Linear', 100)
        ops.pattern('Plain', 200, 100)
        ops.load(self.ops_n_elem, 0, -1, 0)
        ops.constraints('Plain')
        ops.numberer('RCM')
        ops.system('UmfPack')
        ops.test('NormUnbalance', 1e-3, 10)
        ops.algorithm('Newton')
        ops.integrator('LoadControl', config['P'] / config['num_steps_vertical'])
        ops.analysis('Static', '-noWarnings')
        
        # Define recorder
        def record():
            time = ops.getTime()
            section_strains = ops_get_section_strains(self)

            results.applied_axial_load.append(time)
            results.applied_moment_top.append(0)
            results.applied_moment_bot.append(0)
            results.maximum_abs_moment.append(ops_get_maximum_abs_moment(self))
            results.maximum_abs_disp.append(ops_get_maximum_abs_disp(self))
            results.lowest_eigenvalue.append(ops.eigen('-fullGenLapack', 1)[0])
            if type(self.section).__name__ == "RC":
                results.maximum_concrete_compression_strain.append(section_strains[0])
                results.maximum_steel_strain.append(section_strains[1])
            elif type(self.section).__name__ == "I_shape":
                results.maximum_compression_strain.append(section_strains[0])
                results.maximum_tensile_strain.append(section_strains[1])

            if self.axis == 'x':
                results.curvature.append(section_strains[2])
            elif self.axis == 'y':
                results.curvature.append(section_strains[3])
            else:
                raise ValueError(f'The value of axis ({self.axis}) is not supported.')

        
        record()
        
        for i in range(config['num_steps_vertical']):
            ok = ops.analyze(1)
            
            if ok != 0:
                results.exit_message = 'Analysis Failed In Vertical Loading'
                warnings.warn('Analysis Failed In Vertical Loading')
                return results
            
            record()
            
            # Check for  limits
            exit_message = check_analysis_limits(
                    results,
                    eigenvalue_limit=config['eigenvalue_limit'],
                    deformation_limit=config['deformation_limit'],
                    concrete_strain_limit=config['concrete_strain_limit'],
                    steel_strain_limit=config['steel_strain_limit'],
                    section_type=type(self.section).__name__,
                )
            if exit_message:
                results.exit_message = exit_message
                break
            
        # endregion
        
        # region Run lateral load (time = LFH)
        ops.loadConst('-time', 0.0)
        ops.timeSeries('Linear', 101)
        ops.pattern('Plain', 201, 101)
        sgn_et = int(np.sign(self.et))
        sgn_eb = int(np.sign(self.eb))
        if sgn_et != sgn_eb and (sgn_eb != 0 and sgn_et != 0):
            if max(self.et, self.eb, key=abs) == self.et:
                dof = 3 * self.ops_n_elem // 4
                ecc_sign = sgn_et
            else:
                dof = 1 * self.ops_n_elem // 4
                ecc_sign = sgn_eb
            dU = self.length * config['disp_incr_factor'] / 2
            ops.load(self.ops_n_elem, 0, 0, self.et * config['e'] * ecc_sign)
            ops.load(0, 0, 0, -self.eb * config['e'] * ecc_sign)
            ops.integrator('DisplacementControl', dof, 1, dU)
        else:
            ecc_sign = sgn_et
            dU = self.length * config['disp_incr_factor']
            ops.load(self.ops_n_elem, 0, 0, self.et * config['e'] * ecc_sign)
            ops.load(0, 0, 0, -self.eb * config['e'] * ecc_sign)
            ops.integrator('DisplacementControl', self.ops_mid_node, 1, dU)
        
        ops.analysis('Static', '-noWarnings')
        
        # Define recorder
        def record():
            time = ops.getTime()
            section_strains = ops_get_section_strains(self)

            results.applied_axial_load.append(config['P'])
            results.applied_moment_top.append(self.et * time * ecc_sign)
            results.applied_moment_bot.append(-self.eb * time * ecc_sign)
            results.maximum_abs_moment.append(ops_get_maximum_abs_moment(self))
            results.maximum_abs_disp.append(ops_get_maximum_abs_disp(self))
            results.lowest_eigenvalue.append(ops.eigen('-fullGenLapack', 1)[0])
            if type(self.section).__name__ == "RC":
                results.maximum_concrete_compression_strain.append(section_strains[0])
                results.maximum_steel_strain.append(section_strains[1])
            elif type(self.section).__name__ == "I_shape":
                results.maximum_compression_strain.append(section_strains[0])
                results.maximum_tensile_strain.append(section_strains[1])

            if self.axis == 'x':
                results.curvature.append(section_strains[2])
            elif self.axis == 'y':
                results.curvature.append(section_strains[3])
            else:
                raise ValueError(f'The value of axis ({self.axis}) is not supported.')

        record()
        
        maximum_moment = 0
        disp_incr_factor = config['disp_incr_factor']

        while True:
            ok = ops.analyze(1)

            if ok != 0 and config['try_smaller_steps']:
                for div_factor in [1e1, 1e2, 1e3, 1e4, 1e5, 1e6]:
                    update_dU(disp_incr_factor, div_factor)
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

            # Check for drop in applied load (time = the horzontal load factor)
            if config['percent_load_drop_limit'] is not None:
                current_moment = results.maximum_abs_moment[-1]
                maximum_moment = max(current_moment, maximum_moment)
                if current_moment < (1 - config['percent_load_drop_limit']) * maximum_moment:
                    results.exit_message = 'Load Drop Limit Reached'
                    break
                
            # Check for  limits
            exit_message = check_analysis_limits(
                    results,
                    eigenvalue_limit=config['eigenvalue_limit'],
                    deformation_limit=config['deformation_limit'],
                    concrete_strain_limit=config['concrete_strain_limit'],
                    steel_strain_limit=config['steel_strain_limit'],
                    section_type=type(self.section).__name__,
                )
            if exit_message:
                results.exit_message = exit_message
                break
        
        return results
    


    def ops_get_section_strains(self):
        maximum_concrete_compression_strain = []
        maximum_tensile_steel_strain = []
        for i in range(self.ops_n_elem):
            for j in range(self.ops_integration_points):
                axial_strain, curvatureX, curvatureY = 0, 0, 0
                if self.axis == 'x':
                    axial_strain, curvatureX = ops.eleResponse(i,  # element tag
                                                   'section', j+1, # select integration point
                                                   'deformation')  # response type
                elif self.axis == 'y':
                    axial_strain, curvatureY = ops.eleResponse(i,  # element tag
                                                   'section', j+1, # select integration point
                                                   'deformation')  # response type
                else:
                    raise ValueError("The axis is not supported.")

                maximum_concrete_compression_strain.append(self.section.maximum_concrete_compression_strain(
                                                           axial_strain, curvatureX=curvatureX, curvatureY=curvatureY))
                maximum_tensile_steel_strain.append(self.section.maximum_tensile_steel_strain(
                                                           axial_strain, curvatureX=curvatureX, curvatureY=curvatureY))
        return min(maximum_concrete_compression_strain), max(maximum_tensile_steel_strain), curvatureX, curvatureY


    def ops_get_maximum_abs_moment(self) -> float:
        # This code assumed (but does not check) that moment at j-end of 
        # one element equals the moment at the i-end of the next element.
        return max(abs(ops.eleForce(i, 6)) for i in range(self.ops_n_elem))


    def ops_get_maximum_abs_disp(self) -> float:
        return max(abs(ops.nodeDisp(i, 1)) for i in range(self.ops_n_elem))


    @property
    def Cm(self) -> float:
        return 0.6 + 0.4 * min(self.et, self.eb, key=abs) / max(self.et, self.eb, key=abs)


    def calculated_EI_ops(self, P_list, M1_list, M2_ops, Pc_factor=1) -> dict:
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
            Pc = P / (1-self.Cm/delta) / Pc_factor
            k = 1  # Effective length factor (always one for this non-sway column)
            EI = Pc * (k * self.length / pi) ** 2
            if EI>EIgross:
                EI = EIgross
            EI_list_ops.append(EI)

        return {"P":np.array(P_list), "M1":np.array(M1_list), "M2":np.array(M2_list), "Calculated EI":np.array(EI_list_ops),
                "EIgross":EIgross}


    def calculated_EI_design(self, P_list, M1_list, P_design, M2_design, section_factored=False, Pc_factor=1) -> dict:
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

            if M1>M2:
                EI_list_AASHTO.append(EIgross)
                continue

            delta = M2 / M1
            Pc = P / (1 - self.Cm / delta) / Pc_factor
            k = 1  # Effective length factor (always one for this non-sway column)
            EI = Pc * (k * self.length / pi) ** 2
            if EI>EIgross:
                EI = EIgross
            EI_list_AASHTO.append(EI)
            
        return {"P": np.array(P_list), "M1": np.array(M1_list), "M2":np.array(M2_list), "Calculated EI": np.array(EI_list_AASHTO),
                "EIgross": EIgross}