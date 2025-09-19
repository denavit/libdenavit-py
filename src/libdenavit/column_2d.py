from abc import abstractmethod
from libdenavit.OpenSees import AnalysisResults
from libdenavit import interpolate_list, find_limit_point_in_list
import warnings

class Column2d:
    """
    A base class for 2D column models.
    """
    def __init__(self, section, length, **kwargs):
        """
        Initializes common properties for a column.

        Args:
            section: The cross-section object.
            length (float): The length of the column.
            **kwargs: Additional keyword arguments for customization.
        """
        self.section = section
        self.length = length

        # Common analysis and geometry parameters with default values
        self.axis = kwargs.get('axis', None)
        self.ops_n_elem = kwargs.get('ops_n_elem', 8)
        self.ops_element_type = kwargs.get('ops_element_type', 'mixedBeamColumn')
        self.ops_geom_transf_type = kwargs.get('ops_geom_transf_type', 'Corotational')
        self.ops_integration_points = kwargs.get('ops_integration_points', 3)
        self.dxo = kwargs.get('dxo', 0.0)
        self.include_initial_geometric_imperfections = kwargs.get('include_initial_geometric_imperfections', True)
        
    def _init_results(self, attrs: list) -> AnalysisResults:
        """Create an AnalysisResults with the given attribute names as empty lists."""
        results = AnalysisResults()
        for name in attrs:
            setattr(results, name, [])
        return results
    
    def _get_default_deformation_limit(self):
        """Get default deformation limit. Override in subclasses."""
        return 0.1 * self.length
    
    def _extract_analysis_config(self, **kwargs):
        """Extract and validate analysis configuration parameters."""
        config = {
            'section_id': kwargs.get('section_id', 1),
            'section_args': kwargs.get('section_args', []),
            'section_kwargs': kwargs.get('section_kwargs', {}),
            'e': kwargs.get('e', 1.0),
            'P': kwargs.get('P', 0),
            'num_steps_vertical': kwargs.get('num_steps_vertical', 1000),
            'disp_incr_factor': kwargs.get('disp_incr_factor', 1e-5),
            'eigenvalue_limit': kwargs.get('eigenvalue_limit', 0),
            'deformation_limit': kwargs.get('deformation_limit', 'default'),
            'concrete_strain_limit': kwargs.get('concrete_strain_limit', -0.01),
            'steel_strain_limit': kwargs.get('steel_strain_limit', 0.05),
            'percent_load_drop_limit': kwargs.get('percent_load_drop_limit', 0.05),
            'try_smaller_steps': kwargs.get('try_smaller_steps', True),
            'creep_props_dict': kwargs.get('creep_props_dict', dict()),
            'shrinkage_props_dict': kwargs.get('shrinkage_props_dict', dict()),
        }
        
        # Set default deformation limit using subclass override
        if config['deformation_limit'] == 'default':
            config['deformation_limit'] = self._get_default_deformation_limit()
            
        return config
    
    def _initialize_results(self):
        """Initialize analysis results object with required attributes. Override in subclasses."""
        # Base attributes that most column types will need
        attrs = [
            'applied_axial_load', 'maximum_abs_moment', 'maximum_abs_disp', 'lowest_eigenvalue',
            'maximum_concrete_compression_strain', 'maximum_steel_strain', 
            'maximum_compression_strain', 'maximum_tensile_strain', 'curvature'
        ]
        return self._init_results(attrs)
    
    def _set_limit_point_values(self, results, ind, x):
        """Set limit point values. Override in subclasses for additional values."""
        results.applied_axial_load_at_limit_point = interpolate_list(results.applied_axial_load, ind, x)
        results.maximum_abs_moment_at_limit_point = interpolate_list(results.maximum_abs_moment, ind, x)
        results.maximum_abs_disp_at_limit_point = interpolate_list(results.maximum_abs_disp, ind, x)
    
    def _find_limit_point(self, results, config, analysis_type=None):
        """Find and set limit point values."""
        if 'Analysis Failed' in results.exit_message:
            ind, x = find_limit_point_in_list(results.applied_axial_load, max(results.applied_axial_load))
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
            ind, x = find_limit_point_in_list(results.applied_axial_load, max(results.applied_axial_load))
        else:
            # No message at all → set one & fallback
            results.exit_message = 'Analysis Ended Without Explicit Limit'
            ind, x = find_limit_point_in_list(results.applied_axial_load, max(results.applied_axial_load))

        self._set_limit_point_values(results, ind, x)
    
    def run_ops_analysis(self, analysis_type, **kwargs):
        """Template method for running OpenSees analysis."""
        config = self._extract_analysis_config(**kwargs)
        
        self.build_ops_model(config['section_id'], config['section_args'], config['section_kwargs'], 
                            creep_props_dict=config['creep_props_dict'],
                            shrinkage_props_dict=config['shrinkage_props_dict'])
        
        results = self._initialize_results()
        
        # Delegate to specific analysis implementation
        if analysis_type.lower() == 'proportional_limit_point':
            results = self._run_proportional_analysis(config, results)
        elif analysis_type.lower() == 'nonproportional_limit_point':
            results = self._run_nonproportional_analysis(config, results)
        else:
            raise ValueError(f'Analysis type {analysis_type} not implemented')
        
        self._find_limit_point(results, config, analysis_type)
        return results

    @abstractmethod
    def _run_proportional_analysis(self, config, results):
        """Run proportional analysis. Must be implemented by subclasses."""
        pass

    @abstractmethod  
    def _run_nonproportional_analysis(self, config, results):
        """Run nonproportional analysis. Must be implemented by subclasses."""
        pass
    