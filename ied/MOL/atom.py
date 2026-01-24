#!/usr/bin/env python3

# Description
###############################################################################

""" Module to store and manipulate data from atoms.

This module contains the Atom class, which is used to store and manipulate
data from atoms. The class has methods to set and get data from atoms.

Example
-------
	>>> atom = Atom()
	>>> atom.set_coordinates(1.0, 2.0, 3.0)
	>>> print(atom.coordinates)
	(1.0, 2.0, 3.0)

"""

# Imports
###############################################################################
import copy
import logging

from ied.MOL.element_data import data
from ied.MOL.coordinates import Coordinates

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
class Atom:
	""" Class to store and manipulate data from atoms.

	Attributes
	----------
	id : int
		Atom ID.
	name : str
		Atom name.
	element : str
		Atom element.
	index : int
		Atom index.
	coordinates : tuple
		Atom coordinates.
	bfactor : float
		Atom B-factor.
	occupancy : float
		Atom occupancy.
	chain : str
		Atom chain.
	altloc : str
		Atom alternate location.
	resname : str
		Atom residue name.
	resid : int
		Atom residue ID.
	charge : float
		Atom charge.
	segment : str
		Atom segment.
	mass : float
		Atom mass.
	atomic_number : int
		Atom atomic number.
	atomic_radii : float
		Atom atomic radii.
	symbol : str
		Atom symbol.

	"""
	
	def __init__(self, 	
				id: int = 0,
				name: str = "",
				element: str = "",
				atomic_number: int = 0,
				index: int = 0,
				coordinates: tuple = (.0, .0, .0),
				chain: str = "X",
				resid: int = 0,
				resname: str = "",
				altloc: str = "",
				charge: float = .0,
				bfactor: float = .0,
				occupancy: float = .0,
				resinid: float = -1,
				segment: str = "HKN", 
				mass: float = -1,
				atomic_radii: float = -1,
				symbol: str = "",
			):
		''' Constructor for Atom class.

		Parameters
		----------
		id : int, optional
			Atom ID. Default is 0.
		name : str, optional
			Atom name. Default is "".
		element : str, optional
			Atom element. Default is "".
		index : int, optional
			Atom index. Default is 0.
		coordinates : tuple, optional
			Atom coordinates. Default is (.0, .0, .0).
		bfactor : float, optional
			Atom B-factor. Default is .0.
		occupancy : float, optional
			Atom occupancy. Default is .0.
		chain : str, optional
			Atom chain. Default is "X".
		altloc : str, optional
			Atom alternate location. Default is "".
		resname : str, optional
			Atom residue name. Default is "".
		resid : int, optional
			Atom residue ID. Default is 0.
		charge : float, optional
			Atom charge. Default is .0.
		segment : str, optional
			Atom segment. Default is "HKN".
		mass : float, optional
			Atom mass. Default is -1.
		atomic_number : int, optional
			Atom atomic number. Default is 0.
		atomic_radii : float, optional
			Atom atomic radii. Default is -1.
		symbol : str, optional
			Atom symbol. Default is "".
		'''
		
		self.id: int = id
		self.name: str = name
		self.element: str = element
		self.index: int = index
		self.coordinates = Coordinates(coordinates[0], coordinates[1], coordinates[2])
		self.bfactor: float = bfactor
		self.occupancy: float = occupancy
		self.chain: str = chain
		self.altloc: str = altloc
		self.resname: str = resname
		self.resinid: float = resinid
		self.charge: float = charge
		self.segment: str = segment
		self.resid = resid
		self.atomic_number: int = atomic_number
		self.mass: float = mass
		self.electronic_configuration: str
		self.atomic_radii: float  = atomic_radii
		self.symbol: str = symbol
		self.guess_atomic_number()
		self.guess_atomic_radii()

	def set_index(self, index: int) -> None:
		''' Set the atom index.
		
		Parameters
		----------
		index : int
			Atom index.
		'''

		self.index = index

	def set_id(self, id: int) -> None:
		''' Set the atom ID.

		Parameters
		----------
		id : int
			Atom ID
		'''

		self.id = id

	def set_name(self, name: str) -> None:
		''' Set the atom name.
		
		Parameters
		----------
		name : str
			Atom name.
		'''

		self.name = name

	def set_coordinates(self, x: float = 0, y: float = 0, z: float = 0) -> None:
		''' Set the atom coordinates.
		
		Parameters
		----------
		x : float, optional
			X coordinate. Default is 0.
		y : float, optional
			Y coordinate. Default is 0.
		z : float, optional
			Z coordinate. Default is 0.
		'''

		self.coordinates.move_to(x=x, y=y, z=z)

	def set_bfactor(self, value: float) -> None:
		''' Set the atom B-factor.
		
		Parameters
		----------
		value : float
			Atom B-factor.
		'''

		self.bfactor = value
	
	def set_occupancy(self, value: float) -> None:
		''' Set the atom occupancy.
		
		Parameters
		----------
		value : float
			Atom occupancy.
		'''

		self.occupancy = value
	
	def set_chain(self, chain: str) -> None:
		''' Set the atom chain.
		
		Parameters
		----------
		chain : str
			Atom chain.
		'''

		self.chain = chain
	
	def set_element(self, element: str) -> None:
		''' Set the atom element.
		
		Parameters
		----------
		element : str
			Atom element.
		'''

		self.element = element
		
	def set_resid(self, resid: int) -> None:
		''' Set the atom residue ID.
		
		Parameters
		----------
		resid : int
			Atom residue ID.
		'''

		self.resid = resid

	def set_segment(self, segment: str) -> None:
		''' Set the atom segment.
		
		Parameters
		----------
		segment : str
			Atom segment.
		'''

		self.segment = segment

	def set_mass(self, mass: float):
		''' Set the atom mass.
		
		Parameters
		----------
		mass : float
			Atom mass.
		'''

		self.mass = mass
	
	def set_atomic_number(self, atomic_number: int):
		''' Set the atom atomic number.
		
		Parameters
		----------
		atomic_number : int
			Atom atomic number.
		'''

		self.atomic_number = atomic_number

	def guess_atomic_number(self) -> int:
		''' Guess the atom atomic number.

		Returns
		-------
		atomic_number : int
			Atom atomic number.
		'''

		try:
			self.atomic_number = data[self.symbol]['atomic_number']
		except:
			try:
				self.atomic_number = data[self.name[0]]['atomic_number']
			except:
				logging.error(f"Could not guess atomic number for atom '{self.symbol}'.")
				self.atomic_number = 0
		
		return self.atomic_number

	def guess_atomic_radii(self) -> float:
		''' Guess the atom atomic radii.
		
		Returns
		-------
		atomic_radii : float
			Atom atomic radii.
		'''

		try:
			self.atomic_radii = data[self.symbol]['atomic_radii']
		except:
			try:
				self.atomic_radii = data[self.name[0]]['atomic_radii']
			except:
				logging.error(f"Could not guess atomic radii for atom '{self.symbol}'. Setting to 1.")
				self.atomic_radii = 1

		return self.atomic_radii

	def __copy__(self):
		''' Copy the Atom object.	

		Returns
		-------
		Atom
			Copy of the Atom object.
		'''

		cls = self.__class__
		result = cls.__new__(cls)
		result.__dict__.update(self.__dict__)
		return result

	def __deepcopy__(self, memo: dict):
		''' Deep copy the Atom object.

		Parameters
		----------
		memo : dict
			Dictionary to store the objects.

		Returns
		-------
		Atom
			Deep copy of the Atom object.
		'''

		cls = self.__class__
		result = cls.__new__(cls)
		memo[id(self)] = result
		for k, v in self.__dict__.items():
			setattr(result, k, copy.deepcopy(v, memo))
		return result

# Functions
###############################################################################
## Private ##

## Public ##

# Aliases
###############################################################################
