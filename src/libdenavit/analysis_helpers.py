from libdenavit import opensees as ops


def try_analysis_options():
    """
    Tries different analysis algorithms and tolerances in OpenSees.

    This function attempts to converge the analysis using a sequence of
    different algorithms ('ModifiedNewton', 'KrylovNewton') and tolerances.

    Returns:
        int: The status of the analysis (0 if successful, -1 if failed).
    """
    options = [('ModifiedNewton', 1e-3),
               ('KrylovNewton', 1e-3),
               ('KrylovNewton', 1e-2)]

    for algorithm, tolerance in options:
        ops.algorithm(algorithm)
        ops.test('NormUnbalance', tolerance, 10)
        ok = ops.analyze(1)
        if ok == 0:
            break
    return ok


def ops_get_section_strains(column):
    """
    Get section strains for both RC and I_shape sections
    Returns: [maximum_compression_strain, maximum_tensile_strain, curvatureX, curvatureY]
    """
    if type(column.section).__name__ == "RC":
        # Original RC implementation
        maximum_concrete_compression_strain = []
        maximum_tensile_steel_strain = []
        max_kx = 0.0
        max_ky = 0.0
        for i in range(column.ops_n_elem):
            for j in range(column.ops_integration_points):
                axial_strain, curvatureX, curvatureY = 0, 0, 0
                if column.axis == 'x':
                    axial_strain, curvatureX = ops.eleResponse(i,  # element tag
                                                              'section', j+1,  # select integration point
                                                              'deformation')  # response type
                elif column.axis == 'y':
                    axial_strain, curvatureY = ops.eleResponse(i,  # element tag
                                                              'section', j+1,  # select integration point
                                                              'deformation')  # response type
                else:
                    raise ValueError("The axis is not supported.")
                
                maximum_concrete_compression_strain.append(column.section.maximum_concrete_compression_strain(
                                                          axial_strain, curvatureX=curvatureX, curvatureY=curvatureY))
                maximum_tensile_steel_strain.append(column.section.maximum_tensile_steel_strain(
                                                   axial_strain, curvatureX=curvatureX, curvatureY=curvatureY))
                max_kx = max(max_kx, abs(curvatureX))
                max_ky = max(max_ky, abs(curvatureY))

        return min(maximum_concrete_compression_strain), max(maximum_tensile_steel_strain), max_kx, max_ky

    elif type(column.section).__name__ == "I_shape":
        # I_shape implementation using section methods
        maximum_compression_strains = []
        maximum_tensile_strains = []
        max_kx = 0.0
        max_ky = 0.0

        for i in range(column.ops_n_elem):
            for j in range(column.ops_integration_points):
                axial_strain, curvatureX, curvatureY = 0, 0, 0
                if column.axis == 'x':
                    axial_strain, curvatureX = ops.eleResponse(i,  # element tag
                                                              'section', j+1,  # select integration point
                                                              'deformation')  # response type
                elif column.axis == 'y':
                    axial_strain, curvatureY = ops.eleResponse(i,  # element tag
                                                              'section', j+1,  # select integration point
                                                              'deformation')  # response type
                else:
                    raise ValueError("The axis is not supported.")
                
                # Use I_shape methods
                maximum_compression_strains.append(column.section.maximum_compression_strain(
                                                  axial_strain, curvatureX=curvatureX, curvatureY=curvatureY))
                maximum_tensile_strains.append(column.section.maximum_tensile_strain(
                                              axial_strain, curvatureX=curvatureX, curvatureY=curvatureY))
                max_kx = max(max_kx, abs(curvatureX))
                max_ky = max(max_ky, abs(curvatureY))

        return min(maximum_compression_strains), max(maximum_tensile_strains), max_kx, max_ky


def ops_get_maximum_abs_moment(column) -> float:
    """
    Gets the maximum absolute moment from the column elements in OpenSees.

    Args:
        column: The column object with analysis properties.

    Returns:
        float: The maximum absolute moment.
    """
    moment = [abs(ops.eleForce(0, 3))]  # Moment at the start of the first element
    for i in range(column.ops_n_elem):
        moment.append(abs(ops.eleForce(i, 6))) # Moment at the end of each element
    return max(moment)


def ops_get_maximum_abs_disp(column) -> float:
    """
    Gets the maximum absolute lateral displacement from the nodes in OpenSees.

    Args:
        column: The column object with analysis properties.

    Returns:
        float: The maximum absolute lateral displacement.
    """
    return max(abs(ops.nodeDisp(i, 1)) for i in range(column.ops_n_elem + 1))


def check_analysis_limits(results, **limits):
    """
    Checks if the analysis has reached a predefined limit.

    Args:
        results: The AnalysisResults object containing the current state.
        limits (dict): A dictionary of limit values.

    Returns:
        str or None: An exit message string if a limit is reached, otherwise None.
    """
    # Unpack limits from the dictionary
    eigenvalue_limit = limits.get('eigenvalue_limit')
    deformation_limit = limits.get('deformation_limit')
    concrete_strain_limit = limits.get('concrete_strain_limit')
    steel_strain_limit = limits.get('steel_strain_limit')
    
    # Optional; defaults to RC so existing callers keep working.
    section_type = limits.get('section_type', 'RC')

    if eigenvalue_limit is not None and results.lowest_eigenvalue[-1] < eigenvalue_limit:
        return 'Eigenvalue Limit Reached'

    if deformation_limit is not None and results.maximum_abs_disp[-1] > deformation_limit:
            return 'Deformation Limit Reached'
    
    if section_type == "RC":
        # For RC: Check concrete compression and steel tension
        if concrete_strain_limit is not None and results.maximum_concrete_compression_strain[-1] < concrete_strain_limit:
            return 'Concrete Compression Strain Limit Reached'
        
        if steel_strain_limit is not None and results.maximum_steel_strain[-1] > steel_strain_limit:
            return 'Steel Tensile Strain Limit Reached'
        
    elif section_type == "I_shape":
        # compression strains are negative; exceed limit if magnitude > steel_strain_limit
        if steel_strain_limit is not None and results.maximum_compression_strain[-1] < -steel_strain_limit:
            return 'Steel Compression Strain Limit Reached'
        if steel_strain_limit is not None and results.maximum_tensile_strain[-1] > steel_strain_limit:
            return 'Steel Tensile Strain Limit Reached'

    return None