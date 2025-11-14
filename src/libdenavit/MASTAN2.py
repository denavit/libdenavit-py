import numpy as np
from numpy.linalg import norm
from math import sin, cos
from scipy.io import savemat
from pathlib import Path

def el_webdir(nele, ends, beta_ang, coord, coord_end, defl_total):
	y=np.zeros(3)
	WEBDIR = np.zeros([nele, 3]) 
	for ele in range(nele):
		delt = np.subtract(coord[int(ends[ele][1])-1][:], coord[int(ends[ele][0])-1][:])                    
		x = delt/norm(delt,2)
		if abs(x[1]) < 0.9999:
			y[0] = 0.0
			y[1] = 1.0
			y[2] = 0.0
		else:
			y[0] = -1.0
			y[1] = 0.0
			y[2] = 0.0
		delt = np.subtract(coord_end[int(ends[ele][1])-1][:], coord_end[int(ends[ele][0])-1][:])
		x = delt/norm(delt)
		z = np.cross(x, y)
		z = z/norm(z)
		y = np.cross(z, x)
		y = y/norm(y)
		beta_ang_i = x[0]*defl_total[int(ends[ele][0])-1][3] + x[1]*defl_total[int(ends[ele][0])-1][4] + x[2]*defl_total[int(ends[ele][0])-1][5]
		beta_ang_j = x[0]*defl_total[int(ends[ele][1])-1][3] + x[1]*defl_total[int(ends[ele][1])-1][4] + x[2]*defl_total[int(ends[ele][1])-1][5]
		beta = beta_ang[ele] + (beta_ang_i + beta_ang_j)/2
		y[0] = y[0]*cos(beta) + z[0]*sin(beta)
		y[1] = y[1]*cos(beta) + z[1]*sin(beta)
		y[2] = y[2]*cos(beta) + z[2]*sin(beta)
		y = y/norm(y)
		WEBDIR[ele][:] = y
	return WEBDIR


def save_MASTAN2(**attrs):    

    ### Unpack input parameters and set defaults
    
    # Required 
    node_info       = np.array(attrs['node_info'])
    elem_info       = np.array(attrs['elem_info'])
    sect_info       = np.array(attrs['sect_info'])
    mat_info        = np.array(attrs['mat_info'])
    support_info    = np.array(attrs['support_info'])

    # Optional input (basic)
    model_title     = attrs.get('model_title', 'Model')
    sect_name       = attrs.get('sect_name', None)
    mat_name        = attrs.get('mat_name', None)
    nload_info      = attrs.get('nload_info', None)
    uniload_info    = attrs.get('uniload_info', None)    
    thermal_info    = attrs.get('thermal_info', None)
    ground_motion_data = attrs.get('ground_motion_data', None)
   
    # Optional input (advanced)
    analysis_info = attrs.get('analysis_info', np.array([]))
    periods_info = attrs.get('periods_info', np.array([]))   
    first_el_settings = attrs.get('first_el_settings', np.array([]))
    sec_el_settings = attrs.get('sec_el_settings', np.array([]))
    first_inel_settings = attrs.get('first_inel_settings', np.array([]))
    sec_inel_settings = attrs.get('sec_inel_settings', np.array([]))
    ebuck_settings = attrs.get('ebuck_settings', np.array([]))
    ibuck_settings = attrs.get('ibuck_settings', np.array([]))
    period_settings = attrs.get('period_settings', np.array([]))
    deflect_settings = attrs.get('deflect_settings', np.array([]))
    axial_settings = attrs.get('axial_settings', np.array([]))
    sheary_settings = attrs.get('sheary_settings', np.array([]))
    shearz_settings = attrs.get('shearz_settings', np.array([]))
    momx_settings = attrs.get('momx_settings', np.array([]))
    momy_settings = attrs.get('momy_settings', np.array([]))
    momz_settings = attrs.get('momz_settings', np.array([]))
    bimom_settings = attrs.get('bimom_settings', np.array([]))
    first_elastic_dynamic_settings = attrs.get('first_elastic_dynamic_settings', np.array([]))
    sec_elastic_dynamic_settings = attrs.get('sec_elastic_dynamic_settings', np.array([]))
    yldsurfval_settings = attrs.get('yldsurfval_settings', np.array([]))
    normincr_settings = attrs.get('normincr_settings', np.array([]))
    
    # Save Name
    save_name = attrs.get('save_name', 'Model')
    save_path = attrs.get('save_path', Path('.'))
    
    ### Set more advanced defaults
    if sect_name is None:
        sect_name = []
        for i in range(len(sect_info)):
            sect_name.append([f'Section {i}'])
    
    if mat_name is None:
        mat_name = []
        for i in range(len(mat_info)):
            mat_name.append([f'Material {i}'])
            
    if nload_info is None:
        nload_info = np.zeros([len(node_info),6])
        
    if uniload_info is None:
        uniload_info = np.zeros([len(elem_info),3])

    if thermal_info is None:
        thermal_info = np.zeros([len(elem_info),4])

    if ground_motion_data is None:
        # Set default ground motion data
        gm_record_val = [[0.0, 0.0],
                         [0.1, 0.4],
                         [0.2, -0.2],
                         [0.3, 0.0]]

        ground_motion_data = {
            'gm_dir': 1,
            'gm_accel_g': 386.08858,
            'gm_record_val_1': gm_record_val[:][0],
            'gm_record_val_2': gm_record_val[:][1]
        }
        
        ground_motion_data = [1, 386.08858, gm_record_val[:][0], gm_record_val[:][1]]
      
    ### Adjust input parameters types
    if type(sect_name) is str:
        sect_name = [sect_name]
       
    if type(mat_name) is str:
        mat_name = [mat_name]
               
    ### Expand elem_info
    webdir = el_webdir(np.size(elem_info, 0), elem_info[:, 0:2], elem_info[:, 5], node_info, node_info, np.zeros((np.size(node_info, 0), 6)))
    elem_info = np.concatenate((elem_info[:, 0:5], np.ones((np.size(elem_info, 0), 2)), 
						webdir, np.ones((np.size(elem_info, 0), 1)), 
						elem_info[:, 5:7], np.ones((np.size(elem_info, 0), 2)), 
						elem_info[:, 7:9], np.ones((np.size(elem_info, 0), 2)), 
						elem_info[:, 9:17]), 1)
	
    ### Expand nload_info and uniload_info and define fixity_info and settle_info 
    a = np.empty((np.size(nload_info, 0), 6))
    a[:] = np.nan
    nload_info = np.concatenate((nload_info, a), 1)
    fixity_info = np.array(support_info)
    fixity_info[~np.isnan(support_info)] = 0
    fixity_info = np.concatenate((fixity_info, a), 1)
    settle_info = np.array(support_info)
    settle_info[~np.isnan(settle_info) & (settle_info==0)] = np.nan
    settle_info = np.concatenate((settle_info, a), 1)
    b = np.empty((np.size(uniload_info, 0), 6))
    b[:] = np.nan
    uniload_info = np.concatenate((uniload_info, np.zeros((np.size(uniload_info, 0), 3)), b, thermal_info), 1)

    ### Save model data to .mat file
    model_data = dict(model_title = model_title,
                      node_info = node_info.astype('double'),
                      elem_info = elem_info,
                      sect_info = sect_info,
                      sect_name = np.array(sect_name, dtype=object),
                      mat_info = mat_info,
                      mat_name = np.array(mat_name, dtype=object),
                      fixity_info = fixity_info,
                      nload_info = nload_info,
                      uniload_info = uniload_info,
                      settle_info = settle_info,
                      analysis_info = analysis_info,
                      periods_info = periods_info,
                      first_el_settings = first_el_settings,
                      sec_el_settings = sec_el_settings,
                      first_inel_settings = first_inel_settings,
                      sec_inel_settings = sec_inel_settings,
                      ebuck_settings = ebuck_settings,
                      ibuck_settings = ibuck_settings,
                      period_settings = period_settings,
                      ground_motion_data = ground_motion_data,
                      first_elastic_dynamic_settings=first_elastic_dynamic_settings,
                      sec_elastic_dynamic_settings=sec_elastic_dynamic_settings,
                      deflect_settings = deflect_settings,
                      axial_settings = axial_settings,
                      sheary_settings = sheary_settings,
                      shearz_settings = shearz_settings,
                      momx_settings = momx_settings,
                      momy_settings = momy_settings,
                      momz_settings = momz_settings,
                      bimom_settings = bimom_settings,
                      yldsurfval_settings = yldsurfval_settings,
                      normincr_settings = normincr_settings,
                     )    

    savemat(Path(save_path,save_name + '.mat'), model_data)
                                 

if __name__ == "__main__":
    from math import nan, inf

    ##########################
    ### Example Structure  ###
    ##########################
    
    # Notes:
    # - Input arrays can be defined as numpy arrays or lists of lists. 
    
    
    # model_title = Brief text string describing section
    model_title = 'Example_Truss'


    ### Node information ###
    #  Requires three arrays to be constructed: 
    #  node_info, support_info, nload_info
    
    # node_info = [x coordinate, y coordinate, z coordinate]
    node_info = [[  0,   0, 0],     # Node 1's XYZ Coordinates
                 [180, 180, 0],     # Node 2's XYZ Coordinates
                 [360,   0, 0],     # Node 3's XYZ Coordinates
                 [540, 180, 0],     # Node 4's XYZ Coordinates
                 [720,   0, 0],     # Node 5's XYZ Coordinates
                ]       

    # support_info = [support condition of nodes' six d.o.f. (free=NaN, fixed=0, val for prescribed displacement)]
    support_info = [[  0,   0,   0, nan, nan, nan],
                    [nan, nan, nan, nan, nan, nan],
                    [nan, nan, nan, nan, nan, nan],
                    [nan, nan, nan, nan, nan, nan],
                    [  0,   0,   0, nan, nan, nan],
                   ]
    
    # nload_info = [concentrated force or moment on nodes' six d.o.f.]
    nload_info = [[0,    0, 0, 0, 0, 0],
                  [0,    0, 0, 0, 0, 0],
                  [0, -100, 0, 0, 0, 0],
                  [0,    0, 0, 0, 0, 0],
                  [0,    0, 0, 0, 0, 0],
                 ]


    ### Element information ###
    #  Requires three arrays to be constructed:
    #  elem_info, uniload_info, thermal_info

    # elem_info = [Node i, Node j, Section #, Material #, Beta Angle (rads), ...
    #				Flexure condition Node i (rigid=0,pin=1,spring=2), Flexure condition Node j (rigid=0,pin=1,spring=2), ...
    #				Warping condition Node i (fixed=0,free=1,cont=2), Warping condition Node j (fixed=0,free=1,cont=2), ...
    #				Major-axis spring stiffness node i (val = 0 (pin) to inf (rigid)),...
    #				Minor-axis spring stiffness node i (val = 0 (pin) to inf (rigid)),...
    #				Major-axis spring stiffness node j (val = 0 (pin) to inf (rigid)),...
    #				Minor-axis spring stiffness node j (val = 0 (pin) to inf (rigid)),...
    #				Major-axis spring moment capacity Mpz node i (val = value to inf (unlimited)),...
    #				Minor-axis spring moment capacity Mpz node i (val = value to inf (unlimited)),...
    #				Major-axis spring moment capacity Mpz node j (val = value to inf (unlimited)),...
    #				Minor-axis spring moment capacity Mpz node j (val = value to inf (unlimited))]
    elem_info = [[1, 2, 1, 1, 0, 0, 0, 0, 0, inf, inf, inf, inf, inf, inf, inf, inf],
                 [1, 3, 1, 1, 0, 0, 0, 0, 0, inf, inf, inf, inf, inf, inf, inf, inf],
                 [2, 3, 1, 1, 0, 0, 0, 0, 0, inf, inf, inf, inf, inf, inf, inf, inf],
                 [2, 4, 1, 1, 0, 0, 0, 0, 0, inf, inf, inf, inf, inf, inf, inf, inf],
                 [3, 4, 1, 1, 0, 0, 0, 0, 0, inf, inf, inf, inf, inf, inf, inf, inf],
                 [3, 5, 1, 1, 0, 0, 0, 0, 0, inf, inf, inf, inf, inf, inf, inf, inf],
                 [4, 5, 1, 1, 0, 0, 0, 0, 0, inf, inf, inf, inf, inf, inf, inf, inf],
                ]

    # uniload_info = [uniformly distributed loads along elements' local axes (x,y,z)]
    uniload_info = [[0,  0, 0],
                    [0, -1, 0],
                    [0,  0, 0],
                    [0,  0, 0],
                    [0,  0, 0],
                    [0, -1, 0],
                    [0,  0, 0]
                   ]

    # thermal_info = [thermal effects along elements' local axes, including
    #					   thermal coef, dT(centroid), Tgradient(y'), Tgradient(z')]
    thermal_info = [[0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                   ]


    ### Section information ###
    #  Requires two arrays to be constructed:
    #  sect_info and sect_name
    
    # sect_info = [Area (1), Moment of inertia Izz (2), Moment of inertia Iyy (3), Torsion constant J (4),
    #				Warping coefficient Cw (5), Plastic section modulus Zzz (6), Plastic section modulus Zyy (7), ...
    #				Shear area Ayy (8), Shear area Azz (9), Yield surface factor for P (10), Yield surface factor for Mz (11),
    #               Yield surface factor for My (12), Is symmetric? 1=yes/0=no (13), Shear center Ysc (14), Shear center Zsc (15),
    #               Wagner beta_y (16), Wagner beta_z (17), Wagner beta_w (18),Phi in radians (19), Product of inertia Iyz (20)]
    sect_info = [4.44, 48, 3.41, 0.137, 51.8, 13.6, 2.67, 1.98695, 2.1105, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]

    
    # section_name = Brief text string describing section
    section_name = 'W8X15'

    
    ### Material information ###
    #  Requires two arrays to be constructed:
    #  mat_info and mat_name

    # mat_info = [Modulus of elasticity E, Poisson Ratio v, Yield strength Fy, Weight density Wt. Dens.]
    mat_info = [29000, 0.3, 50, 0.0002835]

    # mat_name = Brief text string describing material
    material_name = 'Steel'
                     

    save_MASTAN2(
        model_title = model_title,
        node_info = node_info,
        support_info = support_info,
        nload_info = nload_info,
        elem_info = elem_info,
        uniload_info = uniload_info,
        thermal_info = thermal_info,
        sect_info = sect_info,
        sect_name = section_name,
        mat_info = mat_info,
        mat_name = material_name,
       )