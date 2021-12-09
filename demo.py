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


def distance(x, y, z):
	"""
	Caculate the distance of a point of the origin
	"""

	return np.sqrt(x**2 + y**2 + z**2) 


def separation(pos):
	"""
	Calculate the separation between two points
	"""

	# a, b = np.triu_indices(n=4, m=3)

	# return np.abs(pos[a] - pos[b])


def WGauss(r, h):
	"""
	Returns the Gaussian kernel function
	"""

	W = 1/(h**3 * np/pi**(3/2)) * np.exp(-r**2/ h**2)

	return W


def GradWGauss(r, h, x, y, z):
	"""
	Returns the gradient of the kernel function
	"""

	gradw = -2/(h**5 * np.pi**(3/2)) * np.exp(-r**2/ h**2)
	gradw_x = gradw * x
	gradw_y = gradw * y 
	gradw_z = gradw * z 

	return gradw_x, gradw_y, gradw_z


def density(pos, m, h):
	"""
	Returns the density of the kernel (or something like that, it's late)
	"""

	r = separation(pos, pos)

	return r


def pressure(rho, k, G):
	"""
	Calculates the pressure
	"""

	return k * rho**(1 + 1/G)


def acceleration(P, ):
	"""
	Gets acceleration at each timestep 
	"""

	





def main():

	N = 5			# Total number of particles forming the star
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


	positions = np.random.randn(N, 3)
	velocities = np.zeros((N, 3))

	# print(positions, velocities)

	demo1 = np.array([[1, 2, 3], [1, 4, 5], [6, 8, 7], [1, 1, 1]])
	diffs = separation(demo1)
	print(diffs)
	# demo2 = np.array([[4, 5, 6], [2, 3, 4]])


if __name__ == "__main__":

	main()