import numpy as np 
import matplotlib.pyplot as plt 

"""
TODO:
~ Figure out element by element differences in 2D arrays for separation
~ Density function
~ Acceleration function
~ Updating every time step
~ Plotting 
"""

# Setting the random number generator seed 
np.random.seed(42)


def distance(pos):
	"""
	Caculate the distance of a point of the origin
	"""

	x = pos[0]
	y = pos[1]
	z = pos[2]

	return np.sqrt(x**2 + y**2 + z**2) 


def separation(ri, rj_list):
	"""
	Calculates the separation between a particle with all the other particles

	Parameters:
	ri 			: the distance from origin of the particle whose density is being calculated (float)
	rj_list 	: a 1D list of the distances from the origin for all the N particles (yes, this includes ri as well, but we remove the 0 value in the function) 
	"""

	# Calculating the separation of particle i with all the particles j from 1 to N (including the separation with itself which is zero)
	rsep = np.abs(rj_list - ri)

	# Removing the separation of particle i with itself and converting rsep to an array
	rsep = np.array([i for i in rsep if i != 0])

	return rsep


def WGauss(r, h):
	"""
	Returns the Gaussian kernel function
	"""

	W = 1/(h**3 * np.pi**(3/2)) * np.exp(-r**2/ h**2)

	return W


def GradWGauss(ri, rj, h):
	"""
	Returns the gradient of the kernel function
	"""

	gradw = -2/(h**5 * np.pi**(3/2)) * np.exp(-np.abs(ri - rj)**2/ h**2)
	dx = np.abs(ri[0] - rj[0])
	dy = np.abs(ri[1] - rj[1])
	dz = np.abs(ri[2] - rj[2])
	gradw_x = gradw * x
	gradw_y = gradw * y 
	gradw_z = gradw * z 

	return gradw_x, gradw_y, gradw_z


def density(ri, rj_list, m, h):
	"""
	Returns the density of the kernel for a particle

	Parameters:
	ri 			: the distance from origin of the particle whose density is being calculated (float)
	rj_list 	: a 1D list of the distances from the origin for all the N particles 
	m 			: mass of each particle (float)
	h 			: smoothing length (float)
	"""

	# Calculating the separations between the particles
	rsep = separation(ri, rj_list)

	# Calculating the density for the particle i
	rho_i = np.sum(m*WGauss(rsep, h))

	return rho_i


def pressure(rho, k, G):
	"""
	Calculates the pressure for a particle

	Parameters:
	rho 		: a 1D array of the densities of the N particles
	k 			: the normalization constant for the polytropic equation of state (float)
	G 			: the polytropic index (float)
	"""

	return k * rho**(1 + 1/G)


def acceleration(P, rho, m, h):
	"""
	Gets acceleration at each timestep 
	"""

	





def main():

	N = 3			# Total number of particles forming the star
	M = 2 			# Star mass
	R = 1			# Radius of the star
	m = M/N 		# Mass of a single particle

	t = 10			# Total simulation time
	dt = 0.1 		# Timestep value
	Nt = int(t/dt)	# Total number of timesteps

	h = 0.1 		# Smoothing length
	k = 0.1 		# Normalization constant for calculating the density
	G = 1.4 		# Polytropic constant
	nu = 0.5 		# Damping constant 

	# Randomly generating positions as (x, y, z) coordinates and initializing zero velocities (vx, vy, vz) for every particle
	positions = np.random.randn(N, 3)
	velocities = np.zeros((N, 3))
	print(positions, velocities)

	# Calculating the distances from the origin for each particle
	dists = []
	for coords in positions:
		dists.append(distance(coords))

	# Calculating the density for each particle 
	rhos = []
	for i in range(len(dists)):
		rhos.append(density(dists[i], dists, m, h))
	
	# Calculating the pressure for each particle
	Ps = pressure(np.array(rhos), k, G)

	# Calculating the acceleration for each particle


if __name__ == "__main__":

	main()