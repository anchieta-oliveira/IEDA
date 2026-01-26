#!/usr/bin/env python3

# Description
###############################################################################

""" 
Module to store and manipulate data from PDB files.

This module contains the PDB class, which is used to store and manipulate
data from PDB files. The class has methods to read, write, and show data 
from PDB files.

Example:
--------

    >>> pdb = PDB()
    >>> pdb.read('example.pdb')
    >>> pdb.show()

"""

'''
PDB format example (Protein Data Bank):
-------------------------------------------------------
COLUMNS        DATA TYPE       FIELD         DEFINITION
-------------------------------------------------------

 1 -  6        Record name     "ATOM  "      Atom record.
 7 - 11        Integer         serial        Atom serial number.
13 - 16        Atom            name          Atom name.
17             Character       altLoc        Alternate location indicator.
18 - 20        Residue name    resName       Residue name.
22             Character       chainID       Chain identifier.
23 - 26        Integer         resSeq        Residue sequence number.
27             AChar           iCode         Code for insertion of residues.
31 - 38        Real(8.3)       x             X coordinate in Angstroms.
39 - 46        Real(8.3)       y             Y coordinate in Angstroms.
47 - 54        Real(8.3)       z             Z coordinate in Angstroms.
55 - 60        Real(6.2)       occupancy     Occupancy.
61 - 66        Real(6.2)       tempFactor    Temperature factor.
73 - 76        LString(4)      segID         Segment identifier (left-justified).
77 - 78        LString(2)      element       Element symbol (right-justified).
79 - 80        LString(2)      charge        Charge on the atom.
'''

# Imports
###############################################################################
import os
import copy
import logging
import tempfile
import subprocess
import numpy as np
from ieda.MOL.mol import Mol
from ieda.MOL.atom import Atom
from scipy.spatial.transform import Rotation as R

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
class PDB(Mol):
	""" Class to store and manipulate data from PDB files.

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

	def __init__(self, path: str = "") -> None:
		''' Constructor for PDB class.
		
		Parameters
		----------
		path : str, optional
			Path to the PDB file. Default is an empty string. If path is not provided the PDB object is created empty.
		'''
		super().__init__()
		self.text:str = ""
		self.description:str = ""
		self.natoms:int = 0
		self.data:np.array = np.array([])
		self.atoms:list = []
		self.resids:list = []
		self.name = "None"

		if path != "":
			self.read(path)
			

	def __to_float(self, s: str) -> float:
		''' Convert a string to a float.

		Parameters
		----------
		s : str
			String to be converted.

		Returns
		-------
		float
			Converted float.
		'''

		if any(not c.isdigit() for c in s) or '' == s:
			return .0
		else:
			return float(s)

	def remark_index(self) -> None:
		''' Remark the atom indexes. '''

		self.data.serials = np.array([i for i in range(self.natoms)])
		for i, at in enumerate(self.atoms):
			at.index = i+1

	def remark_inid(self) -> None:
		''' Remark the atom IDs. '''

		self.data.inid = np.array([i for i in range(self.natoms)])
		for i, at in enumerate(self.atoms):
			at.id = i+1

	def make(self, atoms: list=[], data: np.recarray=np.recarray) -> None:
		''' Create a PDB object.
		
		Parameters
		----------
		atoms : list, optional
			List of atoms. Default is an empty list.
		data : np.recarray, optional
			Structured array with the data. Default is an empty array.
		'''

		logging.debug("Creating PDB object from atoms or data.")

		if len(atoms) == 0 and data.size == 0:
			logging.warning("No atoms or data provided to create PDB object.")
			return
		
		if len(atoms) > 0:
			inid = []
			resinid = []
			serials = []
			names = []
			altlocs = []
			resnames = []
			chainids = []
			resseq = []
			xs = []
			ys = []
			zs = []
			occupancys = []
			bfactors = []
			segids = []
			elements = []
			charges = []
			coordinates = []
			icodes = []
			self.atoms = copy.deepcopy(atoms)

			for at in atoms:
				inid.append(at.id)
				resinid.append(at.resinid)
				serials.append(at.index)
				names.append(at.name)
				altlocs.append(at.altloc)
				resnames.append(at.resname)
				chainids.append(at.chain)
				resseq.append(at.resid)
				xs.append(at.coordinates.x)
				ys.append(at.coordinates.y)
				zs.append(at.coordinates.z)
				coordinates.append(np.array([at.coordinates.x, at.coordinates.y, at.coordinates.z]))
				occupancys.append(at.occupancy)
				bfactors.append(at.bfactor)
				segids.append(at.segment)
				elements.append(at.element)
				charges.append(at.charge)
				icodes.append('')
			
			inid = np.array(inid)
			resinid = np.array(resinid)
			serials = np.array(serials)
			names = np.array(names)
			altlocs = np.array(altlocs)
			resnames = np.array(resnames)
			chainids = np.array(chainids)
			resseq = np.array(resseq)
			xs = np.array(xs)
			ys = np.array(ys)
			zs = np.array(zs)
			occupancys = np.array(occupancys)
			bfactors = np.array(bfactors)
			segids = np.array(segids)
			elements = np.array(elements)
			charges = np.array(charges,  dtype=np.float32)
			coordinates = np.array(coordinates,  dtype=np.float32)
			icodes = np.array(icodes)
		
		
		elif data.size > 0:
			inid = data.inid
			resinid = data.resinid
			serials = data.serials
			names = data.names 
			altlocs = data.altlocs
			resnames = data.resnames
			chainids = data.chainids
			resseq = data.resseq
			icodes = data.icodes
			xs = data.xs
			ys = data.ys
			zs = data.zs
			occupancys = data.occupancys
			bfactors = data.bfactors
			segids = data.segids
			elements = data.elements 
			charges = data.charges
			coordinates = np.vstack((data['xs'],data['ys'],data['zs'])).T		
		
		dtype = [
			('inid', inid.dtype),
			('serials', serials.dtype),
			('names', names.dtype),
			('altlocs', altlocs.dtype),
			('resnames', resnames.dtype),
			('chainids', chainids.dtype),
			('resinid', resinid.dtype),
			('resseq', resseq.dtype),
			('icodes', icodes.dtype),
			('xs', xs.dtype),
			('ys', ys.dtype),
			('zs', zs.dtype),
			('coordinates', np.float32, 3),
			('occupancys', occupancys.dtype),
			('bfactors', bfactors.dtype),
			('segids', segids.dtype),
			('elements', elements.dtype),
			('charges', charges.dtype)
		]

		struct = np.zeros(inid.size, dtype=dtype)
		struct['inid'] = inid
		struct['resinid'] = resinid
		struct['serials'] = serials
		struct['names'] = names
		struct['altlocs'] = altlocs
		struct['resnames'] = resnames
		struct['chainids'] = chainids
		struct['resseq'] = resseq
		struct['icodes'] = icodes
		struct['xs'] = xs
		struct['ys'] = ys
		struct['zs'] = zs
		struct['coordinates'] = coordinates
		struct['occupancys'] = occupancys
		struct['bfactors'] = bfactors
		struct['segids'] = segids
		struct['elements'] = elements
		struct['charges'] = charges
		
		self.data = struct.view(np.recarray)
					

	def read(self, path: str) -> np.recarray:
		''' Read a PDB file.

		Parameters
		----------
		path : str
			Path to the PDB file.

		Returns
		-------
		np.recarray
			Structured array with the data.
		'''
		if os.path.isfile(path) == False:
			logging.error(f"File '{path}' not found.")
			raise FileNotFoundError(f"File '{path}' not found.")
		
		self.name = path.split("/")[-1].split(".")[-2]
		logging.debug(f"Reading PDB file from path: {path}")
			
		inid = []
		serials = []
		names = []
		altlocs = []
		resnames = []
		chainids = []
		resseq = []
		resinid = []
		xs = []
		ys = []
		zs = []
		occupancys = []
		bfactors = []
		segids = []
		elements = []
		charges = []
		coordinates = []
		icodes = []

		with open(path, "r") as file:
			pdb_lines = file.readlines()
			logging.debug(f"PDB file '{path}' successfully read. Total lines: {len(pdb_lines)}")

		try:
			i = 0
			rinid = 0

			for line in pdb_lines:
				if line.startswith('ATOM') or line.startswith('HETATM'):
					record = line[0:6].strip()
					inid.append(i)
					index = line[6:11].strip()
					if any(not c.isdigit() for c in index) or "*****" in index:
						index = 00000
					index = int(index)
					serials.append(index)
			
					atom = line[12:16].strip()
					names.append(atom)
					
					resname = line[17:21].strip()
					resnames.append(line[17:21].strip())
					
					altloc = line[16].strip()
					altlocs.append(altloc)
					
					chain = line[21].strip()
					chainids.append(chain)
					
					residue_index = int(line[22:26].strip())
					resseq.append(residue_index)
					if residue_index != resseq[i-1]:  
						rinid += 1

					resinid.append(rinid)

					icodes.append(line[26].strip())

					x = float(line[30:38].strip())
					xs.append(x)
					
					y = float(line[38:46].strip())
					ys.append(y)
					
					z = float(line[46:54].strip())
					zs.append(z)
					
					coordinates.append(np.array([x, y, z], dtype=np.float32))

					occupancy = float(line[54:60].strip())
					occupancys.append(occupancy)
					
					beta = float(line[60:66].strip())
					bfactors.append(beta)
					
					segID = line[72:76].strip()
					segids.append(segID)
					
					element = line[76:78].strip()
					elements.append(element)
					
					charge = self.__to_float(line[78:80].strip())
					charges.append(charge)
					
					self.atoms.append(Atom(index=index, name=atom, altloc=altloc, resname=resname, 
											resid=residue_index, chain=chain, coordinates=(x, y, z), 
											occupancy=occupancy, bfactor=beta, segment=segID, charge=charge, 
											element=element, id=i, resinid=rinid
											))
					i += 1
			
			inid = np.array(inid)
			self.natoms = inid.size
			resinid = np.array(resinid)
			serials = np.array(serials)
			names = np.array(names)
			altlocs = np.array(altlocs)
			resnames = np.array(resnames)
			chainids = np.array(chainids)
			resseq = np.array(resseq)
			xs = np.array(xs)
			ys = np.array(ys)
			zs = np.array(zs)
			occupancys = np.array(occupancys)
			bfactors = np.array(bfactors)
			segids = np.array(segids)
			elements = np.array(elements)
			charges = np.array(charges)
			coordinates = np.array(coordinates,  dtype=np.float32)
			icodes = np.array(icodes)

			dtype = [
				('inid', inid.dtype),
				('serials', serials.dtype),
				('names', names.dtype),
				('altlocs', altlocs.dtype),
				('resnames', resnames.dtype),
				('chainids', chainids.dtype),
				('resinid', resinid.dtype),
				('resseq', resseq.dtype),
				('icodes', icodes.dtype),
				('xs', xs.dtype),
				('ys', ys.dtype),
				('zs', zs.dtype),
				('coordinates', np.float32, 3),
				('occupancys', occupancys.dtype),
				('bfactors', bfactors.dtype),
				('segids', segids.dtype),
				('elements', elements.dtype),
				('charges', charges.dtype)
			]

			struct = np.zeros(inid.size, dtype=dtype)

			struct['inid'] = inid
			struct['resinid'] = resinid
			struct['serials'] = serials
			struct['names'] = names
			struct['altlocs'] = altlocs
			struct['resnames'] = resnames
			struct['chainids'] = chainids
			struct['resseq'] = resseq
			struct['icodes'] = icodes
			struct['xs'] = xs
			struct['ys'] = ys
			struct['zs'] = zs
			struct['coordinates'] = coordinates
			struct['occupancys'] = occupancys
			struct['bfactors'] = bfactors
			struct['segids'] = segids
			struct['elements'] = elements
			struct['charges'] = charges
			
			self.data = struct.view(np.recarray)
		
		except Exception as e:
			logging.error(f"Error reading PDB file '{path}': {e}")
			raise e

		return self.data

	def add_atoms(self, atoms: list) -> None:
		''' Add atoms to the PDB object.

		Parameters
		----------
		atoms : list
			List of atoms.
		'''
		logging.debug(f"Adding {len(atoms)} atoms to the PDB object.")

		new_rows = len(atoms)
		dtype = self.data.dtype 
		new_data = np.zeros(new_rows, dtype=dtype)

		for i, at in enumerate(atoms):
			new_data[i]['inid']      = at.id
			new_data[i]['resinid']   = at.resinid
			new_data[i]['serials']   = at.index
			new_data[i]['names']     = at.name
			new_data[i]['altlocs']   = at.altloc
			new_data[i]['resnames']  = at.resname
			new_data[i]['chainids']  = at.chain
			new_data[i]['resseq']    = at.resid
			new_data[i]['icodes']    = ''
			new_data[i]['xs']        = at.coordinates.x
			new_data[i]['ys']        = at.coordinates.y
			new_data[i]['zs']        = at.coordinates.z
			new_data[i]['occupancys'] = at.occupancy
			new_data[i]['bfactors']   = at.bfactor
			new_data[i]['segids']     = at.segment
			new_data[i]['elements']   = at.element
			new_data[i]['charges']    = at.charge

		combined = np.concatenate([self.data, new_data]).view(np.recarray)

		self.data = combined
		self.atoms += copy.deepcopy(atoms)
		

	def get_text(self) -> str:
		''' Get the PDB file content.
		
		Returns
		-------
		str
			PDB file content.
		'''
		logging.debug("Generating PDB file content as text.")
		txt = ""
		for i in range(self.data.shape[0]):
			txt += f"ATOM  {self.data.serials[i]:>5} {self.data.names[i]:<4}{self.data.altlocs[i]:<1}{self.data.resnames[i]:>3} {self.data.chainids[i]:>1}{self.data.resseq[i]:>4}{self.data.icodes[i]:>1}   {self.data.xs[i]:8.3f}{self.data.ys[i]:8.3f}{self.data.zs[i]:8.3f}{self.data.occupancys[i]:6.2f}{self.data.bfactors[i]:6.2f}{self.data.segids[i]:<4}{self.data.elements[i]:>2}{self.data.charges[i]:2}\n"
		
		return txt
	
	def show(self, software: str = "vmd") -> None:
		''' Show the PDB file content.

		Parameters
		----------
		software : str, optional
			Software to show the PDB file. Default is "vmd". Options are "vmd" and "pymol".

		Raises
		------
		ValueError
			If the software is invalid.
		'''

		logging.debug(f"Showing PDB file using {software} software.")

		# Check if the software is valid
		if software not in ["vmd", "pymol"]:
			raise ValueError(f"Invalid software. Choose 'vmd' or 'pymol'. Got: '{software}'")

		text_pdb = self.get_text()

		with tempfile.NamedTemporaryFile(mode='w+', delete=True, suffix=".pdb") as temp_file:
			# Escrever conteúdo no arquivo temporário
			temp_file.write(text_pdb)
			temp_file.flush()

			# Obter o nome do arquivo temporário
			temp_file_name = temp_file.name 

			# Imprimir o nome do arquivo temporário (opcional, apenas para fins de visualização)
			if software == "vmd":
				subprocess.run(f"vmd -pdb {temp_file_name}", shell=True)
			elif software == "pymol":
				subprocess.run(f"pymol {temp_file_name}", shell=True)


	def write(self, path: str) -> None:
		''' Write the PDB file content to a file.

		Parameters
		----------
		path : str
			Path to save the PDB file.
		
		Raises
		------
		ValueError
			If the path is invalid or the file already exists.
		'''

		logging.debug(f"Writing PDB file to path: {path}")

		# Check if the path is valid and if the file already exists
		if os.path.isdir(path):
			raise ValueError(f"Invalid path: '{path}' is a directory.")
		else:
			logging.debug(f"Path '{path}' is valid.")
		
		if os.path.exists(path):
			os.rename(path, path + ".bak")
			logging.warning(f"File '{path}' already exists. Renamed to '{path}.bak'.")
			
		with open(path, "w") as file_pdb:
			file_pdb.write(self.get_text())

		logging.debug(f"PDB file successfully written to '{path}'.")	

	def updata_atoms(self) -> None:
		''' Update the atoms list. '''
		logging.debug("Updating atoms list from structured data array.")
		logging.debug(f"It metod is deprecated and will be removed in future versions.")

		self.atoms = [Atom(
							id=e.inid, index=e.serials, name=e.names, altloc=e.altlocs, 
					  		resname=e.resnames, chain=e.chainids, resid=e.resseq, coordinates=(e.xs, e.ys, e.zs),
							occupancy=e.occupancys, bfactor=e.bfactors, segment=e.segids, element=e.elements, charge=e.charges
							) 
							for e in self.data]	
		
	def rotate(self, angle: float, axis: str) -> None:
		''' Rotate the PDB object.

		Parameters
		----------
		angle : float
			Angle to rotate the PDB object.
		axis : str
			Axis to rotate the PDB object. Options are 'x', 'y', and 'z'.
		'''
		logging.debug(f"Rotating PDB object by {angle} degrees around {axis}-axis.")

		# Criar a rotação usando scipy
		rotation = R.from_euler(axis, angle, degrees=True)  # 'degrees=True' se o ângulo estiver em graus

		# Aplicar a rotação a todas as coordenadas
		self.data.coordinates = rotation.apply(self.data.coordinates)
		self.data.xs, self.data.ys, self.data.zs = self.data.coordinates[:, 0], self.data.coordinates[:, 1], self.data.coordinates[:, 2]

		for i, at in enumerate(self.atoms):
			at.coordinates.x, at.coordinates.y, at.coordinates.z = self.coordinates[i] 


	def move_center_to(self, center: tuple = (0, 0, 0)) -> None:
		''' Move the center of the PDB object to a specific point.

		Parameters
		----------
		center : tuple, optional
			Center of the PDB object. Default is (0, 0, 0).
		'''
		logging.debug(f"Moving PDB object center to {center}.")

		current_center = self.get_center()
		translation_vector = np.array(center) - np.array(current_center)
		self.coordinates += translation_vector
		self.data.xs, self.data.ys, self.data.zs = self.data.coordinates[:, 0], self.data.coordinates[:, 1], self.data.coordinates[:, 2]

		for i, at in enumerate(self.atoms):
			at.coordinates.x, at.coordinates.y, at.coordinates.z = self.coordinates[i] 

	def __copy__(self):
		''' Copy the PDB object.

		Returns
		-------
		PDB
			Copied PDB object.
		'''
		logging.debug("Creating a shallow copy of the PDB object.")

		cls = self.__class__
		result = cls.__new__(cls)
		result.__dict__.update(self.__dict__)

		return result
		
	def __deepcopy__(self, memo: dict):
		''' Deep copy the PDB object.

		Paremeters
		----------
		memo : dict
			Memory dictionary.

		Returns
		-------
		PDB
			Deep copied PDB object.
		'''
		logging.debug("Creating a deep copy of the PDB object.")

		cls = self.__class__
		result = cls.__new__(cls)
		memo[id(self)] = result

		for k, v in self.__dict__.items():
			setattr(result, k, copy.deepcopy(v, memo))

		return result
