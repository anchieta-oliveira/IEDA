#!/usr/bin/env python3

# Description
###############################################################################

""" Module to store and manipulate data from coordinates.

This module contains the Coordinates class, which is used to store and manipulate
data from coordinates. The class has methods to set and get data from coordinates.

Example
-------
	>>> coord = Coordinates()
	>>> coord.set_coordinates(1.0, 2.0, 3.0)
	>>> print(coord.coordinates)
	(1.0, 2.0, 3.0)

"""

# Imports
###############################################################################
import logging
import math

import numpy as np

from scipy.optimize import minimize

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# License
###############################################################################
'''
IEDA
Authors: José de Anchieta de Oliveira Filho
[The Federal University of Rio de Janeiro]
Carlos Chagas Filho Institute of Biophysics
Laboratory for Molecular Modeling and Dynamics
Av. Carlos Chagas Filho 373 - CCS - bloco G1-19,
Cidade Universitária - Rio de Janeiro, RJ, CEP: 21941-902
E-mail address: anchieta.oliveira@biof.ufrj.br
This project is licensed under Creative Commons license (CC-BY-4.0) (Ver qual)
'''

# Classes
###############################################################################
class Coordinates:
	""" Class to store and manipulate data from coordinates.
	
	Attributes
	----------
	x : float
		Coordinate x.
	y : float
		Coordinate y.
	z : float
		Coordinate z.
	"""

	def __init__(self, x: float = 0.0, y: float = 0.0, z: float = 0.0) -> None:
		''' Constructor for Coordinates class.

		Parameters
		----------
		x : float
			Coordinate x.
		y : float
			Coordinate y.
		z : float
			Coordinate z.
		'''

		self.x: float = x
		self.y: float = y
		self.z: float = z

	def get_x(self) -> float:
		''' Return the x coordinate.
		
		Returns
		-------
		float
			Coordinate x.
		'''

		return self.x

	def get_y(self) -> float:
		''' Return the y coordinate.

		Returns
		-------
		float
			Coordinate y.
		'''

		return self.y

	def get_z(self) -> float:
		''' Return the z coordinate.
		
		Returns
		-------
		float
			Coordinate z.
		'''

		return self.z
	
	def get_tuple(self):
		''' Return the coordinates as a tuple.

		Returns
		-------
		tuple[float, float, float]
			Coordinates as a tuple (x, y, z).
		'''

		return (self.x, self.y, self.z)
	
	def get_array(self) -> np.array:
		''' Return the coordinates as a numpy array.

		Returns
		-------
		np.array
			Coordinates as a numpy
		'''

		return np.array(self.get_tuple(), dtype = np.float64)

	def __change_x(self,x: float) -> None:
		''' Change the x coordinate.
		
		Parameters
		----------
		x : float
			New x coordinate.
		'''

		self.x = x
	
	def __change_y(self,y: float) -> None:
		''' Change the y coordinate.

		Parameters
		----------
		y : float
			New y coordinate.
		'''

		self.y = y 

	def __change_z(self, z: float) -> None:
		''' Change the z coordinate.

		Parameters
		----------
		z : float
			New z coordinate.
		'''

		self.z = z

	def convert_to(self, unit: str) -> None:
		''' Convert the coordinates to a different unit.

		Parameters
		----------
		unit : str
			Unit to convert the coordinates to. Can be 'angstrom' or 'bohr'.
		
		Raises
		------
		ValueError
			If the unit is not 'angstrom' or 'bohr'.
		'''

		# Check if the unit is valid
		if unit not in ['angstrom', 'bohr']:
			raise ValueError("Invalid unit. Choose 'angstrom' or 'bohr'.")

		if unit == "bohr":
			a_to_bohr = 1.8897259886
			self.x *= a_to_bohr
			self.y *= a_to_bohr
			self.z *= a_to_bohr

	def set_coordinates(self, x: float, y: float, z: float) -> None:
		''' Set the coordinates. '''

		self.__change_x(x)
		self.__change_y(y)
		self.__change_z(z)

	def measure_distance(self, x: float, y: float, z: float) -> float:
		''' Measure the distance between two points.

		Parameters
		----------
		x : float
			X coordinate of the second point.
		y : float
			Y coordinate of the second point.
		z : float
			Z coordinate of the second point.

		Returns
		-------
		float
			Distance between the two points.
		'''

		return math.sqrt((self.x - x) ** 2 + (self.y - y) ** 2 + (self.z - z) ** 2)

	def move_to(self, x: float, y: float, z: float) -> None:
		''' Move the coordinates to a new position.

		Parameters
		----------
		x : float
			New x coordinate.
		y : float
			New y coordinate.
		z : float
			New z coordinate.
		'''

		self.__change_x(x)
		self.__change_y(y)
		self.__change_z(z)

	def move_to_distance_from(self, distance: float, x: float, y:float, z:float) -> None:
		''' Move the coordinates to a new position at a certain distance from a point.

		Parameters
		----------
		distance : float
			Distance from the point.
		x : float
			X coordinate of the point.
		y : float
			Y coordinate of the point.
		z : float
			Z coordinate of the point.
		'''

		vector = (self.x - x, self.y - y, self.z- z)
		distance_atual = self.measure_distance(x, y, z)
		normalized_vector = (vector[0]/distance_atual, vector[1]/distance_atual, vector[2]/distance_atual)
		new_corrd = (x+ normalized_vector[0]*distance, y+ normalized_vector[1]*distance, z+ normalized_vector[2]*distance) 
		self.move_to(new_corrd[0], new_corrd[1], new_corrd[2])

	def move_atom_to_new_ref_atoms(self, ref_atoms: list, new_atoms: list) -> None:
		''' Move the atom to a new position based on reference atoms.

		Parameters
		----------
		ref_atoms : list
			List of reference atoms.
		new_atoms : list
			List of new atoms.
		'''

		coord = self.get_tuple()
		ref_atoms = np.array(ref_atoms)
		new_atoms = np.array(new_atoms)
		coord = np.array(coord)

		def rotation_matrix(angles: list) -> np.array:
			''' Calculate the rotation matrix from the angles.
			
			Parameters
			----------
			angles : list
				List of angles (alpha, beta, gamma).

			Returns
			-------
			np.array
				Rotation matrix.
			'''

			alpha, beta, gamma = angles
			Rx = np.array([[1, 0, 0], [0, np.cos(alpha), -np.sin(alpha)], [0, np.sin(alpha), np.cos(alpha)]], dtype=np.float64)
			Ry = np.array([[np.cos(beta), 0, np.sin(beta)], [0, 1, 0], [-np.sin(beta), 0, np.cos(beta)]], dtype=np.float64)
			Rz = np.array([[np.cos(gamma), -np.sin(gamma), 0], [np.sin(gamma), np.cos(gamma), 0], [0, 0, 1]], dtype=np.float64)
			
			return np.dot(Rz, np.dot(Ry, Rx))

        # Função para calcular o erro de ajuste entre os pontos de referência e os novos pontos
		def error_func(params: list) -> float:
			''' Calculate the error function for the optimization.
			
			Parameters
			----------
			params : list
				List of parameters (alpha, beta, gamma, tx, ty, tz).

			Returns
			-------
			float
				Error function value.
			'''

			R = rotation_matrix(params[:3])
			t = params[3:]
			new_atoms_transformed = np.dot(ref_atoms - coord, R.T) + coord + t
			return np.sum(np.square(new_atoms_transformed - new_atoms), dtype=np.float64)
		
        # Minimização do erro de ajuste para encontrar a melhor rotação e translação
		initial_params = np.zeros(6)
		result = minimize(error_func, initial_params, method='BFGS')

        # Obter os parâmetros encontrados
		R = rotation_matrix(result.x[:3])
		t = result.x[3:]

		new_coord = np.dot(coord - coord, R.T) + coord + t
		
		self.set_coordinates(x=new_coord[0], y=new_coord[1], z=new_coord[2])

	def rotate(self, axis: str = 'z', angle: float = 90.0) -> np.array:
		''' Rotate the coordinates around an axis.

		Parameters
		----------
		axis : str
			Axis to rotate the coordinates around. Can be 'x', 'y', or 'z'.
		angle : float
			Angle of rotation in degrees.

		Returns
		-------
		np.array
			Rotated coordinates.
		'''

		# Lowercase the axis
		axis = axis.lower()

		# Convert the angle to radians
		angle_radians = np.radians(angle)

		# Define the rotation matrix components
		if axis == 'x':
			rotation_matrix = np.array([[1, 0, 0],
										[0, np.cos(angle_radians), -np.sin(angle_radians)],
										[0, np.sin(angle_radians), np.cos(angle_radians)]])
		elif axis == 'y':
			rotation_matrix = np.array([[np.cos(angle_radians), 0, np.sin(angle_radians)],
										[0, 1, 0],
										[-np.sin(angle_radians), 0, np.cos(angle_radians)]])
		elif axis == 'z':
			rotation_matrix = np.array([[np.cos(angle_radians), -np.sin(angle_radians), 0],
										[np.sin(angle_radians), np.cos(angle_radians), 0],
										[0, 0, 1]])
		else:
			raise ValueError(f"Axis must be 'x', 'y', or 'z' but got '{axis}'.")

		# Apply the rotation matrix
		points = np.array([self.x, self.y, self.z], dtype=np.float64)
		rotated_points = np.dot(points, rotation_matrix.T)
		self.x, self.y, self.z = rotated_points

		return rotated_points

# Functions
###############################################################################
## Private ##

## Public ##

# Aliases
###############################################################################
