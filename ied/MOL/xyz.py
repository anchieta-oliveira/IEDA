#!/usr/bin/env python3

# Description
###############################################################################

""" Module to store and manipulate data from XYZ files.

This module contains the XYZ class, which is used to store and manipulate
data from XYZ files. The class has methods to read, write, and show data from
XYZ files. The data is stored in a structured array, which can be accessed
using the get_text method.

Example
-------
	>>> xyz = XYZ()
	>>> xyz.read("file.xyz")
	>>> print(xyz.get_text())
	>>> xyz.show("vmd")
	>>> xyz.write("file.xyz")

"""

# Imports
###############################################################################
import os
import shutil
import logging
import tempfile
import subprocess
import numpy as np
from ied.MOL import Mol
from ied.MOL.atom import Atom

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
class XYZ(Mol):
	""" Class to store and manipulate data from XYZ files.
	
	Attributes
	----------
	name : str
		XYZ file name.
	description : str
		XYZ file description.
	atoms : list
		List of Atom objects.
	natoms : int
		Number of atoms.
	inid : np.array
		Array of atom IDs.
	elements : np.array
		Array of atom elements.
	xs : np.array
		Array of x coordinates.
	ys : np.array
		Array of y coordinates.
	zs : np.array
		Array of z coordinates.
	coordinates : np.array
		Array of coordinates.
	data : np.array
		Structured array to store the data.
	"""

	def __init__(self, path:str="") -> None:
		''' Constructor for XYZ class. '''

		self.description:str = ""
		self.atoms = []
		self.natoms:int = 0
		self.data:np.array = np.array([])

		if path != "":
			self.read(path)


	def read(self, path: str) -> None:
		''' Read data from an XYZ file.

		Parameters
		----------
		path : str
			Path to the XYZ file.
		
		Raises
		------
		FileNotFoundError
			If the file is not found.
		'''
		if os.path.isfile(path) == False:
			logging.error(f"File '{path}' not found.")
			raise FileNotFoundError(f"File '{path}' not found.")

		logging.debug(f"Reading XYZ file from path: {path}")

		self.name = path.split("/")[-1].split(".")[-2]
		inid = []
		xs = []
		ys = []
		zs = []
		coordinates = []
		elements = []
		
		with open(path, "r") as file:
			xyz_lines = file.readlines()
			logging.debug(f"Total lines read from file: {len(xyz_lines)}")
		
		self.natoms = int(xyz_lines[0])
		self.description = xyz_lines[1].strip()

		for i, line in enumerate(xyz_lines[2:]):
			inid.append(i)
			ls = line.split()
			elements.append(ls[0].strip())
			xs.append(float(ls[1]))
			ys.append(float(ls[2]))
			zs.append(float(ls[3]))
			coordinates.append(np.array([float(ls[1]), float(ls[2]), float(ls[3])]))

			self.atoms.append(Atom(index=i, name=elements[i], altloc="", resname="UNK", 
						   				resid=0, chain="Z", coordinates=(xs[i], ys[i], zs[i]), 
										occupancy=0.00, bfactor=0.00, segment="UNKK", charge=0.00, 
										element=elements[i], id=i, resinid=0
										))

		logging.debug("Updating structured data array.")

		dtype = [
			('inid', 'i4'),
			('xs', 'f8'),
			('ys', 'f8'),
			('zs', 'f8'),
			('elements', 'O'), 
		]
		
		struct = np.zeros(len(inid), dtype=dtype)

		struct['inid'] = np.array(inid, dtype=np.int32)
		struct['xs'] = np.array(xs, dtype=np.float32)
		struct['ys'] = np.array(ys, dtype=np.float32)
		struct['zs'] = np.array(zs, dtype=np.float32)
		struct['elements'] = np.array(elements, dtype=object)

		self.data = struct.view(np.recarray)


	def get_text(self) -> str:
		''' Return the text representation of the data.

		Returns
		-------
		str
			Text representation of the data.
		'''
		logging.debug("Generating text representation of the XYZ data.")

		txt = f"{self.natoms}\n{self.description}\n"

		for i in range(self.data.shape[0]):
			txt += f"{self.elements[i]:<2}\t{self.xs[i]:8.3f}{self.ys[i]:8.3f}{self.zs[i]:8.3f}\n"

		return txt
	
	def show(self, software: str = "vmd") -> None:
		''' Show the data using a visualization software.

		Parameters
		----------
		software : str, optional
			Visualization software. Default is "vmd".
		
		Raises
		------
		ValueError
			If the software is not available.
		'''
		logging.debug(f"Showing XYZ data using {software}.")

		# Check if the software is available
		if software not in ["vmd", "pymol"]:
			raise ValueError(f"Software '{software}' not available. Use 'vmd' or 'pymol'.")
		
		if shutil.which(software) is None:
			logging.error(f"{software} is available.")
			raise ValueError(f"{software} is not available. Please install it or check your PATH.")
		
		text_pdb = self.get_text()

		with tempfile.NamedTemporaryFile(mode='w+', delete=True, suffix=".pdb") as temp_file:
			# Write the text to a temporary file
			temp_file.write(text_pdb)
			temp_file.flush()

			# Get the name of the temporary file
			temp_file_name = temp_file.name 

			# Print the name of the temporary file (optional, only for visualization purposes)
			if software == "vmd":
				cmd = f"vmd -xyz {temp_file_name}"
			elif software == "pymol":
				cmd = f"pymol {temp_file_name}"
			else:
				cmd = ""
			
			if cmd:
				logging.debug(f"Running command: {cmd}")
				subprocess.run(cmd, shell = True)

	def write(self, path: str):
		''' Write data to an XYZ file.

		Parameters
		----------
		path : str
			Path to the XYZ file.
		
		Returns
		-------
		bool
			True if the data was written successfully, False otherwise.
		'''
		if os.path.isdir(path):
			raise ValueError(f"Invalid path: '{path}' is a directory.")
		else:
			logging.debug(f"Path '{path}' is valid.")
		
		if os.path.exists(path):
			os.rename(path, path + ".bak")
			logging.warning(f"File '{path}' already exists. Renamed to '{path}.bak'.")
		
		with open(path, "w") as file_pdb:
			file_pdb.write(self.get_text())

		logging.debug(f"XYZ data written to file: {path}")
		

	def move_center_to(self, center: tuple = (0, 0, 0)) -> None:
		''' Move the center of the molecule to a new position.

		Parameters
		----------
		center : tuple, optional
			New center position. Default is (0, 0, 0).
		'''
		logging.debug(f"Moving molecule center to: {center}")

		current_center = self.get_center()
		translation_vector = np.array(center) - np.array(current_center)
		self.coordinates += translation_vector
		self.xs, self.ys, self.zs = self.coordinates[:, 0], self.coordinates[:, 1], self.coordinates[:, 2]

		for i, at in enumerate(self.atoms):
			at.coordinates.x, at.coordinates.y, at.coordinates.z = self.coordinates[i] 
