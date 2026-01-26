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
PSF Files:
---------------------------------------
https://www.ks.uiuc.edu/Training/TutorialsOverview/namd/namd-tutorial-unix-html/node23.html
'''

# Imports
###############################################################################
import logging
import numpy as np



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
class PSF:
	""" Class to store and manipulate data from PSF files.

	Attributes
	----------
	text : str
		PSF file content.

	"""

	def __init__(self, path: str = "") -> None:
		''' Constructor for PSF class.
		
		Parameters
		----------
		path : str, optional
			Path to the PSF file. Default is an empty string. If path is not provided the PSF object is created empty.
		'''

		self.text:str = ""
		self.data:np.array = np.array([])
		self.atoms:list = []
		self.resids:list = []
		self.name = "None"
		self.inid = np.array([], dtype=np.int32)
		self.resinid = np.array([], dtype=np.int32)
		self.serials = np.array([], dtype=np.int32)
		self.names = np.array([], dtype=np.str_)
		self.resnames = np.array([], dtype=np.str_)
		self.resseq = np.array([], dtype=np.int32)
		self.charges = np.array([], dtype=np.float32)
		if path != "":
			self.read(path)

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

		if path:
			self.name = path.split("/")[-1].split(".")[-2]
		else:
			self.name = ""
		inid = []
		serials = []
		names = []
		resnames = []
		resseq = []
		resinid = []
		segids = []
		charges = []

		try:
			with open(path, "r") as file:
				pdb_lines = file.readlines()

			i = 0
			rinid = 0
			ntitle = False; natom = False
			for line in pdb_lines:
				if '!NTITLE' in line: ntitle = True
				if '!NATOM' in line:  natom = True; ntitle = False; continue
				if '!NBOND:' in line: nbond = True; natom = False; continue
				
				if ntitle: pass
				if natom and line != "\n":
		
					lsp = line.split()
					inid.append(i)
					index = lsp[0].strip()
					if any(not c.isdigit() for c in index) or "*****" in index:
						index = 00000
					index = int(index)
					serials.append(index)

					segids.append(lsp[1])

					atom = lsp[4]
					names.append(atom)
					
					resname = lsp[3]
					resnames.append(resname)
										
					residue_index = int(lsp[2])
					resseq.append(residue_index)
					if residue_index != resseq[i-1]:  
						rinid += 1
					
					resinid.append(rinid)
					
					charge = float(lsp[6])
					charges.append(charge)
					i += 1
			
			self.inid = np.array(inid)
			self.resinid = np.array(resinid)
			self.serials = np.array(serials)
			self.names = np.array(names)
			self.resnames = np.array(resnames)
			self.resseq = np.array(resseq)
			self.segids = np.array(segids)
			self.charges = np.array(charges)
			self.update_data()	
		except FileNotFoundError:
			print(f"File '{path}' not found.")
			return None

		return self.data


	def update_data(self) -> None:
		''' Update the structured array with the data. '''

		dtype = [
			('inid', self.inid.dtype),
			('serials', self.serials.dtype),
			('names', self.names.dtype),
			('resnames', self.resnames.dtype),
			('resinid', self.resinid.dtype),
			('resseq', self.resseq.dtype),
			('segids', self.segids.dtype),
			('charges', self.charges.dtype)
		]

		self.data = np.recarray(self.inid.size, dtype=dtype)

		self.data['inid'] = self.inid
		self.data['resinid'] = self.resinid
		self.data['serials'] = self.serials
		self.data['names'] = self.names
		self.data['resnames'] = self.resnames
		self.data['resseq'] = self.resseq
		self.data['segids'] = self.segids
		self.data['charges'] = self.charges

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