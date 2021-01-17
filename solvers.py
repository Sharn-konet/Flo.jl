from numba.core.types.containers import Tuple
from numba.core.types.functions import Function
import numpy as np
import json

from numba import jitclass, vectorize, float64

from example_funcs import StrangeAttractors

# Plot in 3D
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

import plotly.express as px
import plotly.graph_objects as go

def read_solver(solver):
	with open("solvers.json", "r") as file:
		solvers = json.load(file)

	BT = ButcherTableau(**{att: np.array(eval(expression)) for att, expression in solvers['explicit_methods'][solver].items()}, name = solver)
	return BT

class ButcherTableau(object):
	def __init__(self, alpha: np.array, beta: np.array, gamma: np.array, name = "Custom Solver"):
		self.alpha = alpha
		self.beta = beta
		self.gamma = gamma
		self.embedded = alpha.ndim == 2

		assert((len(alpha) == len(beta)) or (alpha.shape[1] == beta.shape[0]))
		assert(gamma.shape[0] == gamma.shape[1])

	def __str__(self):
		rep = ""
		if self.title is not None: rep += "Model: " + self.title + "\n\n"

		# rep += "Butcher Tableau:\n===============\n\n"
		# for b in self.beta:
		# 	rep += str(b) + '|' + [str()]

	def __repr__(self):
		return self.title

# @jitclass([
# 	('alpha', float64[:]),
# 	('beta', float64[:]),
# 	('gamma', float64[:,:]),
# 	('positions', float64[:,:]),
# 	# ('behaviour', Function),
# 	# ('extra_args', Tuple),
# 	('clock', float),
# 	('ttl', float64[:])
# ])
class Swarm(object):
	def __init__(self, tableau: ButcherTableau, function: Function, *args, size: int = 1000, dimensions = 2, initial_t = 0):
		self.alpha, self.beta, self.gamma = tableau.alpha, tableau.beta, tableau.gamma
		self.tableau = tableau

		np.random.seed(1)

		self.positions = (np.random.rand(dimensions, size))*5

		# Currently particles dont reference the view properly
		self.particles = [Particle(self.positions[:,i], ttl = 100) for i in range(size)]

		self.behaviour = function
		self.extra_args = args

		self.error_history = np.empty(shape = (dimensions, size, 0))

		self.clock = initial_t

		self.ttl = np.empty(size)

	def swarm_step(self, step: float = 0.001) -> float:

		# Set initial stepsize
		h = step

		# Preallocate array
		derivative_evaluations = np.zeros([self.positions.shape[0], self.positions.shape[1], len(self.beta)])

		# Calculate derivative 
		for j in range(len(self.beta)): # j is index of particular derivative
			derivative_evaluations[:,:,j] = self.behaviour(self.clock + self.beta[j]*h, self.positions + h*np.matmul(derivative_evaluations, self.gamma[j,:]), *self.extra_args)

		# Calculate gradient at current tn for each embedded method
		f = np.matmul(derivative_evaluations, self.alpha.T)

		# Calculate estimates for y for each method in embedded method
		if self.tableau.embedded:
			y_estimates =  h*f + self.positions[:,:,np.newaxis] 
			# Get estimate of total error
			error = abs(np.diff(y_estimates, axis = 2))
			self.error_history = np.append(self.error_history, error, axis = 2)
			y_estimates = y_estimates[:,:,0]
		else:
			y_estimates = h*f + self.positions
				
		self.clock += h
		self.positions = y_estimates

class Particle(Swarm):
	def __init__(self, position: np.ndarray, ttl: float) -> None:
		self.position = position
		self.ttl = ttl

	def __repr__(self):
		return str(list(self.position))
		
	def calculate_step():
		pass

if __name__ == "__main__":
	BT = read_solver("Dormand-Prince")
	particles = Swarm(BT, StrangeAttractors.lorenz, size = 200, dimensions = 3)

	frames = []

	initial_scatter = go.Scatter3d(
		x=particles.positions[0,:],
		y=particles.positions[1,:],
		z=particles.positions[2,:],
		mode = 'markers'
	)

	frames.append(go.Frame(data = initial_scatter))

	for _ in range(1000):
		particles.swarm_step(0.001)
		frames.append(go.Frame(
			data = [
				go.Scatter3d(
					x=particles.positions[0,:],
					y=particles.positions[1,:],
					z=particles.positions[2,:],
					mode = 'markers',
					showlegend=False)
				]
		))


	iterations = np.array([*range(particles.error_history.shape[2])])
	swarm_size = np.array([*range(particles.error_history.shape[1])])
	iterations, swarm_size = np.meshgrid(iterations, swarm_size)

	fig = go.Figure(data=[go.Surface(z = particles.error_history[0,:,:])])

	fig.update_layout(title='Error Surface', autosize=True)

	fig.show()

	fig = go.Figure(
		data = initial_scatter,
		layout = go.Layout(
			title = go.layout.Title(text = "Lorenz Attractor"),
			scene = dict(
				xaxis=dict(range=[min([np.min(frame['data'][0]['x']) for frame in frames]), 
								  max([np.max(frame['data'][0]['x']) for frame in frames])], 
					autorange=False),
				yaxis=dict(range=[min([np.min(frame['data'][0]['y']) for frame in frames]), 
								  max([np.max(frame['data'][0]['y']) for frame in frames])], 
					autorange=False),
				zaxis=dict(range=[min([np.min(frame['data'][0]['z']) for frame in frames]), 
								  max([np.max(frame['data'][0]['z']) for frame in frames])], 
					autorange = False)
				),
			updatemenus=[dict(
				type="buttons",
				buttons=[dict(label="Play",
							method="animate",
							args=[None, {"frame": {"duration": 1/60, 
                                                    "redraw": True},
                                                    "fromcurrent": True, 
                                                    "transition": {"duration": 0.00}}])])]
		),
		frames = frames
	)

	fig.update_layout(scene_aspectmode='cube')

	fig.show()