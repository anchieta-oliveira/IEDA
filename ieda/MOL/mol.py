
# Imports
###############################################################################
import logging
import numpy as np
from scipy.spatial import distance
from abc import ABC, abstractmethod

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


class Mol(ABC):
	""" Class to store and manipulate Molecule.

	Attributes
	----------
	text : str
		PDB file content.
	description : str
		PDB file description.
	natoms : int
		Number of atoms.
	data : np.recarray
		Structured array with the data.
	atoms : list
		List of atoms.
	"""
	def __new__(cls, path: str = "", *args, **kwargs):
		if cls is Mol and path != "":
			ext = path.split('.')[-1].lower()
			if ext == "pdb":
				from ieda.MOL.PDB import PDB
				return super().__new__(PDB)
			elif ext == "xyz":
				from ieda.MOL.xyz import XYZ
				return super().__new__(XYZ)
			else:
				logging.error(f"Unsupported file extension: {ext}")
				raise ValueError(f"Unsupported file extension: {ext}")
			
		return super().__new__(cls)

	def __init__(self, path: str = "") -> None:
		''' Constructor for Mol class.
		
		Parameters
		----------
		path : str, optional
			Path to the Mol file. Default is an empty string..
		'''

		self.text:str = ""
		self.description:str = ""
		self.natoms:int = 0
		self.data:np.array = np.array([])
		self.atoms:list = []
		self.resids:list = []
		self.name = "None"

		if path != "":
			self.read(path)
			
	
	@abstractmethod
	def read(self, path: str) -> None:
		''' Abstract method to read a Mol file.
		
		Parameters
		----------
		path : str
			Path to the Mol file.
		'''
		pass
	
	@abstractmethod
	def write(self, path: str) -> None:
		''' Write the Mol file content to a file.

		Parameters
		----------
		path : str
			Path to save the PDB file.
		
		Raises
		------
		ValueError
			If the path is invalid or the file already exists.
		'''
		pass

	def make(self, data: np.recarray) -> None:
		''' Make the Mol object from a structured array.

		Parameters
		----------
		data : np.recarray
			Structured array containing the Mol data.
		'''
		self.natoms = self.data.inid.size 
		self.data = data

	def set_name(self, name: str) -> None:
		''' Set the PDB name.

		Parameters
		----------
		name : str
			PDB name.
		'''

		self.name = name


	def get_name(self) -> str:
		''' Get the PDB name.

		Returns
		-------
		str
			PDB name.
		'''
		logging.debug(f"Getting molecule name: {self.name}")

		return self.name
	

	def set_bfactor(self, value: float, ids: list) -> None:
		''' Set the B-factor for a list of atoms.

		Parameters
		----------
		value : float
			B-factor value.
		ids : list
			List of atom IDs.
		'''
		self.data.bfactors[ids] = value
		for i in ids:
			self.atoms[i].bfactor = value

		logging.debug(f"Set B-factor to {value} for atoms with IDs: {ids}")


	def set_occupancy(self, value: float, ids: list) -> None:
		''' Set the occupancy for a list of atoms.

		Parameters
		----------
		value : float
			Occupancy value.
		ids : list
			List of atom IDs.
		'''
		self.data.occupancys[ids] = value
		for i in ids:
			self.atoms[i].occupancy = value
		
		logging.debug(f"Set occupancy to {value} for atoms with IDs: {ids}")


	def get_distance_matrix(self) -> np.array:
		''' Get the distance matrix.

		Returns
		-------
		np.array
			Distance matrix.
		'''
		logging.debug("Calculating distance matrix.")
		return distance.cdist(self.data.coordinates, self.data.coordinates, 'euclidean')
		


	def get_center(self) -> np.array:
		''' Get the center of the Mol object.

		Returns
		-------
		np.array
			Center of the PDB object.
		'''

		logging.debug("Calculating center of the molecule.")
		center = np.mean(self.data.coordinates, axis = 0)
		
		logging.debug(f"Center coordinates: {center}")

		return center