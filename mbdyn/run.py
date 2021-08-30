import logging
#from matplotlib import pyplot as plt 
import precice
import time
import subprocess
from .csvreader import *	#csvreader.py
import pathlib
import numpy as np
import decimal
import socket
import errno
from scipy import interpolate

class MbdynAdapter:

	#initialize values
	
	from .initialize import __init__
	
	#run preCICE coupling
	def run(self):
		interface = precice.Interface(self.participant_name, self.configuration_file_name,
				                      self.solver_process_index, self.solver_process_size)

		for i in self.mesh_name:
			self.mesh_id.append(interface.get_mesh_id(str(i)))
		dimensions = interface.get_dimensions()
		
		for i in range(self.patches):
			XA, YA, ZA = np.hsplit(csvImport(self.current_path + '/../' + self.fluid_folder_name + '/patch' + str(i+1) + '.csv'),3)	#if mesh is in file
			self.vertices.append(np.array([XA.flatten(),YA.flatten(),ZA.flatten()]).T)

		
		#interpolate from beam to airfoil
		for i in range(self.patches):
			x = self.vertices[i][:,0]
			x = np.roll(x, -1)
			x = np.split(x,2)
			y = self.vertices[i][:,1]
			y = np.roll(y, -1)
			y = np.split(y,2)
			y_min = np.array((np.array(y[0]).min(axis=0),np.array(y[1]).min(axis=0))).max(axis=0)
			y_max = np.array((np.array(y[0]).max(axis=0),np.array(y[1]).max(axis=0))).min(axis=0)
			y_center = y_min + (y_max-y_min)/2
			num_of_points = 11
			
			#x = f(y)
			f = [interpolate.interp1d(y[0],x[0]),interpolate.interp1d(y[1],x[1])]
			ynew = [np.linspace(y_min,y_max,11),np.linspace(y_min,y_max,11)]
			xnew = [f[0](ynew[0]),f[1](ynew[1])]
			
			#x[0] rechts
			#x[1] links
			#von links nach rechts
			airfoil_mesh = np.concatenate(np.array((np.array((xnew[1],ynew[1],np.zeros((xnew[1].size)))) ,np.array((xnew[0],ynew[0],np.zeros((xnew[0].size)))))),axis=1)
			self.vertices[i] = airfoil_mesh.T

			self.vertex_ids.append(interface.set_mesh_vertices(self.mesh_id[i],self.vertices[i]))
			self.read_data_id.append(interface.get_data_id(self.read_data_name[i], self.mesh_id[i]))
			self.write_data_id.append(interface.get_data_id(self.write_data_name[i], self.mesh_id[i]))
		
		dt = interface.initialize()
		
	

		while interface.is_coupling_ongoing():
			
			if interface.is_action_required(precice.action_write_iteration_checkpoint()):
				print("DUMMY: Writing iteration checkpoint")
				interface.mark_action_fulfilled(precice.action_write_iteration_checkpoint())
			
			self.write_data.clear()
			for i in range(self.patches):
				if self.iteration % 2 == 0:
					iteration = self.iteration/2
				self.write_data.append(np.array((np.ones((22))*0+1e-6*iteration,np.ones((22))*0,np.ones((22))*0)).T)
				
			if interface.is_write_data_required(dt):
				for i in range(self.patches):
					interface.write_block_vector_data(self.write_data_id[i], self.vertex_ids[i], self.write_data[i])
			
			print("DUMMY: Advancing in time")
			dt = interface.advance(dt)
			if interface.is_action_required(precice.action_read_iteration_checkpoint()):
				print("DUMMY: Reading iteration checkpoint")
				interface.mark_action_fulfilled(precice.action_read_iteration_checkpoint())
				print("repeating timestep")
			else:
				print("precice is advancing in time")
				#MBDyn advance

			self.iteration += 1
		print("preCICE finalizing")
		interface.finalize()
		print("DUMMY: Closing python solver dummy...")
