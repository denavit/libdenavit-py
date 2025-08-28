from libdenavit.OpenSees import AnalysisResults



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
    
    def _init_results(self, attrs: list) -> AnalysisResults:
        """Create an AnalysisResults with the given attribute names as empty lists."""
        results = AnalysisResults()
        for name in attrs:
            setattr(results, name, [])
        return results    
    