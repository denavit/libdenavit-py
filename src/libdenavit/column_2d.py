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
    