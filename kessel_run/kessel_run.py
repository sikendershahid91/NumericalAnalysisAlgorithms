#!/usr/bin/python3

import scipy.io as sio 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import random
import sys
from timeit import default_timer as timer

class Shuttle:
	def __init__(self):
		self.mass = 1
		self.x_coordinate = random.uniform(-5, 5)
		self.y_coordinate = -10
		self.log_position = np.array([self.x_coordinate,self.y_coordinate]).reshape([1,2])
		self.trajectory = random.gauss(np.pi/2, np.pi/4)
		self.velocity = random.uniform(2,5)
		self.velocity_x = self.velocity * np.cos(self.trajectory)
		self.velocity_y = self.velocity * np.sin(self.trajectory)
		self.success = False
		self.time = 0.00

	def update_position(self, dt):
		# Update positions using Euler's finite difference method
		self.x_coordinate = self.x_coordinate + dt*self.velocity_x
		self.y_coordinate = self.y_coordinate + dt*self.velocity_y
		self.log_position = np.append(self.log_position,[[self.x_coordinate, self.y_coordinate]], axis = 0) 

	def update_velocity(self, dt, acceleration_x, acceleration_y):
		# Update Velocity using Euler's finite difference method
		self.velocity_x = self.velocity_x + dt * acceleration_x
		self.velocity_y = self.velocity_y + dt * acceleration_y


class Clusters:
	def __init__(self, x_np_array, y_np_array, m_np_array):
		self.x_coordinates = x_np_array
		self.y_coordinates = y_np_array
		self.masses = m_np_array
		self.force = 0

	def force_on_target(self, target):
		# Calculating the force on the target shuttle 
		distance_x = target.x_coordinate - self.x_coordinates
		distance_y = target.y_coordinate - self.y_coordinates
		distance_norm = np.sqrt((np.power(distance_x,2) + np.power(distance_y,2)))
		force_x = ((self.masses * distance_x)/np.power(abs(distance_norm),3)).sum()
		force_y = ((self.masses * distance_y)/np.power(abs(distance_norm),3)).sum()
		self.force = -(np.sqrt(np.power(force_x,2)+np.power(force_y,2)))
		return -force_x, -force_y


def main():
	# You may change the number of MC samples from the commandline 
	# $ ./kessel_run.py [T] [dt] [shuttle_num]
	if(len(sys.argv) == 2):
		T = sys.argv[1]
	T = sys.argv[1] if len(sys.argv) >= 2 else 1000
	dt = sys.argv[2] if len(sys.argv) == 3 else 0.01
	shuttle_num = sys.argv[3] if len(sys.argv) == 4 else 5

	# loading matlab_data and initilizing clusters and 2 shuttles
	matlab_data = sio.loadmat('cluster1.mat')
	clusters = Clusters(matlab_data['hX'], matlab_data['hY'], matlab_data['hM'])
	
	# Following list is a log for succesful runs
	success_shuttle = [] 
	program = True
	while program:
		shuttles = [Shuttle() for _ in range(int(shuttle_num))]
		# MC simulation
		for shuttle in shuttles:
			shuttle.update_position(dt)
			# starting a timer and will store it for comparison on graph
			start = timer()
			for t in range(1,int(T)-1):
				force_x, force_y = clusters.force_on_target(shuttle); 
				
				# Error check 1
				if(abs(clusters.force) >= 4):
					shuttle.success = False
					break
				# Euler methods 
				shuttle.update_velocity(dt,force_x,force_y) 
				shuttle.update_position(dt)

				# Error check 2
				if(abs(shuttle.x_coordinate) > 10):
					shuttle.success = False
					break

				# Success check
				if(shuttle.y_coordinate > 10):
					shuttle.success = True
					break

			end = timer()
			shuttle.time = end - start

		# save successful runs 
		for shuttle in shuttles:
			if shuttle.success:
				success_shuttle.append([shuttle.time,shuttle])
		success_shuttle.sort()

		# If need more runs to compare 
		if len(success_shuttle)>1:
			print("Error: only one successful shuttle, running another "+str(5)+" simulation")
			# Rerun if no two shuttle survived. 
			program = False

	# Plotting Shortest and Longest runs while displaying clusters and time taken to run 
	n_short = np.arange(len(success_shuttle[0][1].log_position[:,0]))
	n_long = np.arange(len(success_shuttle[-1][1].log_position[:,0]))
	plt.scatter(clusters.x_coordinates, clusters.y_coordinates)
	plt.scatter(success_shuttle[0][1].log_position[:,0], success_shuttle[0][1].log_position[:,1],marker='s', c=n_short, cmap=plt.cm.coolwarm)
	plt.scatter(success_shuttle[-1][1].log_position[:,0], success_shuttle[-1][1].log_position[:,1],marker='d',c=n_long, cmap=plt.cm.coolwarm)
	plt.title("Longest and Shortest path")
	plt.xlim([-10,10])
	plt.ylim([-10,10])
	plt.annotate("Time: {0:.4f} seconds".format(success_shuttle[0][1].time), xy=(success_shuttle[0][1].x_coordinate,success_shuttle[0][1].y_coordinate),
		xytext=(5,5),arrowprops=dict(arrowstyle="->"))
	plt.annotate("Time: {0:.4f} seconds".format(success_shuttle[-1][1].time), xy=(success_shuttle[-1][1].x_coordinate,success_shuttle[-1][1].y_coordinate),
		xytext=(5,3),arrowprops=dict(arrowstyle="->"))
	plt.show()
		
		
if __name__ == '__main__':
		main()	