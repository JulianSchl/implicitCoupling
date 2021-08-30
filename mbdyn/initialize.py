import logging
import pathlib
import numpy as np
import decimal
import pathlib

def __init__(self, case_name='membrane.msh', config_file_name='../precice-config.xml', 
			 participant_name=None, mesh_name=None, inter_mesh=False,
			 init_data=False, connectivity=False):

	
	### initial values for preCICE coupling
	
	self.displacement_relative_delta = []
	self.displacement_absolute_delta = []
	self.displacement_absolute_delta_str = []
	self.rotation_absolute_delta = []
	self.rotation_absolute_delta_str = []
	self.mbdyn_mesh = []
	self.mbdyn_mesh_str = []
	self.nodes_absolute = []
	self.nodes_absolute_str = []
	self.tvc_socket = []
	self.differ = []
	self.fluid_folder_name = 'openfoam'
	self.current_path = str(pathlib.Path(__file__).parent.absolute())
	self.configuration_file_name = self.current_path + "/../precice-config.xml"
	self.participant_name = "Solid"
	self.mesh_name = ['Solid-Mesh-Right','Solid-Mesh-Left']
	self.mesh_id = []
	self.patch = []
	self.patch2 = []
	self.patch_L = []
	self.transform = []
	self.transform_L = []
	self.vertices = []
	self.vertex_ids = []
	self.read_data = []
	self.force_tensor = []
	self.read_data_id = []
	self.write_data_id = []
	self.write_data_name = ['Displacement_Data_Right','Displacement_Data_Left']
	self.read_data_name = ['Force_Data_Right','Force_Data_Left']
	self.write_data = []
	self.write_data2 = []
	self.write_data_tmp = []
	self.centre_of_gravity = []
	self.exchange_socket = []
	self.displacement = []
	self.displacement_absolute = []
	self.rotation = []
	self.message = []

	if not (len(self.write_data_name) == len(self.read_data_name)):
		exit("You need to define the number of patches")
	self.patches = len(self.write_data_name)


	self.solver_process_index = 0
	self.solver_process_size = 1

	self.direction = True

	self.zeit = 10**-1
	self.omega = 2*np.pi/self.zeit
	self.last_rot_angle = 0
	self.iteration = 0

