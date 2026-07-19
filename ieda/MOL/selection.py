#!/usr/bin/env python3

# Description
###############################################################################

""" Module to store and manipulate data from selections.

This module contains the Selection class, which is used to store and manipulate
data from selections. The class has methods to select data from a PDB file.

Example
-------
	>>> selection = Selection('resid 1 to 10')
	>>> selection.result.save('selection.pdb')

"""

# Imports
###############################################################################
import inspect
import logging

import numpy as np
from ieda.MOL.mol import Mol
from ieda.MOL.PDB import PDB
from ieda.MOL.resid_data import *

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
This project is licensed under the MIT License.
'''

# Classes
###############################################################################

class Selection:
	""" Class to store and manipulate data from selections.

	Attributes
	----------
	selection : str
		Selection string.
	mol : None
		Molecule object.
	name : str
		Selection name.
	result : Mol
		Selected data.
	functions : dict
		Dictionary of selection functions.
	sel_data : np.recarray
		Selected data in structured array.
	"""

	def __init__(self, 
				selection: str,
				mol: Mol, 
				name: str = ""
		        ):
		''' Constructor for Selection class.
		
		Parameters
		----------
		selection : str
			Selection string.
		mol : None
			Molecule object.
		name : str
			Selection
		'''
		
		self.selection:str = selection
		self.mol = mol
		self.result:Mol
		self.name = name
		self.functions:dict
		self.sel_data:np.recarray
		self.__functions()
		self.__call_function()

	def __list_to_float(self, l: list) -> list:
		''' Convert a list of strings to a list of floats.
		
		Parameters
		----------
		l : list
			List of strings.

		Returns
		-------
		list
			List of floats.
		'''

		return list(map(lambda x: float(x), l))
	
	def __num_to(self, l: list) -> list:
		''' Convert a list of strings to a list of floats.

		Parameters
		----------
		l : list
			List of strings.

		Returns
		-------
		list
			List of floats.
		'''

		if "to" in list(l):
			return list(range(int(l[0]), int(l[-1]) + 1))
		return l

	def all(self, list: list, data: np.recarray) -> np.recarray:
		''' Select all data.
		
		Parameters
		----------
		list : list
			List of data.
		data : np.recarray
			Data to be selected.

		Returns
		-------
		np.recarray
			Selected data.
		'''

		return  data

	def select_resids(self, resseq: list, data: np.recarray) -> np.recarray:
		''' Select residues by residue ID.
		
		Parameters
		----------
		resseq : list
			List of residue IDs.
		data : np.recarray
			Data to be selected.
		
		Returns
		-------
		np.recarray
			Selected data.
		'''
		logging.debug(f"Selecting residues with resseq: {resseq}")
		return data[np.isin(data.resseq, self.__list_to_float(resseq))]

	def select_resnames(self, resnames: list, data: np.recarray) -> np.recarray:
		''' Select residues by residue name.

		Parameters
		----------
		resnames : list
			List of residue names.
		data : np.recarray
			Data to be selected.

		Returns
		-------
		np.recarray
			Selected data.
		'''
		logging.debug(f"Selecting residues with resnames: {resnames}")
		return data[np.isin(data.resnames, resnames)]

	def select_resinid(self, resinid: list, data: np.recarray) -> np.recarray:
		''' Select residues by residue ID.

		Parameters
		----------
		resinid : list
			List of residue IDs.
		data : np.recarray
			Data to be selected.

		Returns
		-------
		np.recarray
			Selected data.
		'''
		logging.debug(f"Selecting residues with resinid: {resinid}")
		resinid = self.__num_to(resinid)
		return data[np.isin(data.resinid, self.__list_to_float(resinid))]

	def select_atoms(self, serials: list, data: np.recarray) -> np.recarray:
		''' Select atoms by serial number.

		Parameters
		----------
		serials : list
			List of serial numbers.
		data : np.recarray
			Data to be selected.

		Returns
		-------
		np.recarray
			Selected data.
		'''
		logging.debug(f"Selecting atoms with serials: {serials}")
		serials = self.__num_to(serials)
		return data[np.isin(data.serials, self.__list_to_float(serials))]
	
	def select_atoms_ids(self, inids: list, data: np.recarray) -> np.recarray:
		''' Select atoms by atom ID.

		Parameters
		----------
		inids : list
			List of atom IDs.
		data : np.recarray
			Data to be selected.

		Returns
		-------
		np.recarray
			Selected data.
		'''
		logging.debug(f"Selecting atoms with inids: {inids}")
		inids = self.__num_to(inids)
		return data[np.isin(data.inid, self.__list_to_float(inids))]

	def select_name(self, names: list, data: np.recarray) -> np.recarray:
		''' Select atoms by atom name.

		Parameters
		----------
		names : list
			List of atom names.

		Returns
		-------
		np.recarray
			Selected data.
		'''
		logging.debug(f"Selecting atoms with names: {names}")
		try:
			dt = data[np.isin(data.names, names)]
		except Exception as e:
			logging.error(f"Error selecting atoms by name: {e}")
			dt = []	
		return dt

	def select_chain(self, chainids: list, data: np.recarray) -> np.recarray:
		''' Select chains by chain ID.
		
		Parameters
		----------
		chainids : list
			List of chain IDs.
		data : np.recarray
			Data to be selected.

		Returns
		-------
		np.recarray
			Selected data.
		'''
		logging.debug(f"Selecting chains with chainids: {chainids}")
		return data[np.isin(data.chainids, chainids)]

	def select_segname(self, segids: list, data: np.recarray) -> np.recarray:
		''' Select segments by segment ID.

		Parameters
		----------
		segids : list
			List of segment IDs.
		data : np.recarray
			Data to be selected.

		Returns
		-------
		np.recarray
			Selected data.
		'''
		logging.debug(f"Selecting segments with segids: {segids}")
		return data[np.isin(data.segids, segids)]

	def select_beta(self, beta: list, data: np.recarray) -> np.recarray:
		''' Select atoms by beta factor.

		Parameters
		----------
		beta : list
			List of beta factor values.
		data : np.recarray
			Data to be selected.

		Returns
		-------
		np.recarray
			Selected data.

		Raises
		------
		ValueError
			If the operator is invalid.
		'''

		if beta[0] == ">":
			return data[data.bfactors > float(beta[1])]
		elif beta[0] == "<":
			return data[data.bfactors < float(beta[1])]
		elif beta[0] == ">=":
			return data[data.bfactors >= float(beta[1])]
		elif beta[0] == "<=":
			return data[data.bfactors <= float(beta[1])]
		elif beta[0] == "=" or  beta[0] == "==":
			return data[data.bfactors == float(beta[1])]
		else:
			raise ValueError("Invalid operator for beta factor selection.")	
	
	def select_occupancy(self, occupancys: list, data: np.recarray) -> np.recarray:
		''' Select atoms by occupancy.

		Parameters
		----------
		occupancys : list
			List of occupancy values.
		data : np.recarray
			Data to be selected.

		Returns
		-------
		np.recarray
			Selected data.

		Raises
		------
		ValueError
			If the operator is invalid.
		'''
		logging.debug(f"Selecting atoms with occupancy condition: {occupancys}")
		if occupancys[0] == ">":
			return data[data.occupancys > float(occupancys[1])]
		elif occupancys[0] == "<":
			return data[data.occupancys < float(occupancys[1])]
		elif occupancys[0] == ">=":
			return data[data.occupancys >= float(occupancys[1])]
		elif occupancys[0] == "<=":
			return data[data.occupancys <= float(occupancys[1])]
		elif occupancys[0] == "=" or  occupancys[0] == "==":
			return data[data.occupancys == float(occupancys[1])]
		else:
			raise ValueError("Invalid operator for occupancy selection.")

	def select_within(self, args: list, data: np.recarray) -> np.recarray:
		''' Select atoms within a certain distance from a point.	

		Parameters
		----------
		args : list
			List of arguments.
		data : np.recarray

		Returns
		-------
		np.recarray
			Selected data.
		'''
		logging.debug(f"Selecting atoms within distance: {args}")
		dist = float(args[0])
		sels = self.__get_sels(selection=' '.join(args[2:]))
		d_tmp = self.__call_fun_in_sel(sels=sels, data=data)
		ids_sel = np.where(np.isin(data['inid'], d_tmp['inid']))[0]
		dmx = self.mol.get_distance_matrix()

		return data[np.any(dmx[ids_sel] <= dist, axis=0)]

	def select_not(self, args: list, data: np.recarray) -> np.recarray:
		''' Select atoms that are not in the selection.
		
		Parameters
		----------
		args : list
			List of arguments.
		data : np.recarray
			Data to be selected.

		Returns
		-------
		np.recarray
			Selected data.
		'''
		logging.debug(f"Selecting atoms not in selection: {args}")
		sels = self.__get_sels(selection=' '.join(args))
		d_tmp = self.__call_fun_in_sel(sels=sels, data=data)

		return data[np.isin(data.inid, d_tmp.inid, invert=True)]

	def select_same(self, args: list, data: np.recarray) -> np.recarray:
		''' Select atoms that are the same in the selection.
		
		Parameters
		----------
		args : list
			List of arguments.
		data : np.recarray
			Data to be selected.

		Returns
		-------
		np.recarray
			Selected data.
		'''
		logging.debug(f"Selecting atoms same in selection: {args}")
		sels = self.__get_sels(selection=' '.join(args[2:]))
		d_tmp = self.__call_fun_in_sel(sels=sels, data=data)

		return self.functions[args[0]](d_tmp[list(inspect.signature(self.functions[args[0]]).parameters.values())[0].name], data) 

	def select_protein(self, args: list, data: np.recarray) -> np.recarray:
		''' Select protein atoms.

		Parameters
		----------
		args : list
			List of arguments.
		data : np.recarray
			Data to be selected.

		Returns
		-------
		np.recarray
			Selected data.
		'''
		logging.debug(f"Selecting protein atoms. Names: {[*resids_data]}")
		return self.select_resnames(resnames=[*resids_data], data=data)

	def select_backbone(self, args: list, data: np.recarray) -> np.recarray:
		''' Select backbone atoms.

		Parameters
		----------
		args : list
			List of arguments.
		data : np.recarray
			Data to be selected.

		Returns
		-------
		np.recarray
			Selected data.
		'''
		logging.debug(f"Selecting backbone atoms. Names: ['N', 'HN', 'HA', 'CA', 'C', 'O']")
		result = self.select_protein(args = args, data = data)	

		return self.select_name(names = ['N', 'HN', 'HA', 'CA', 'C', 'O'], data = result)

	def select_water(self, args: list, data: np.recarray) -> np.recarray:
		''' Select water atoms.

		Parameters
		----------
		args : list
			List of arguments.
		data : np.recarray
			Data to be selected.

		Returns
		-------
		np.recarray
			Selected data.
		'''
		logging.debug(f"Selecting water atoms. Names: ['TIP3', 'HOH']")
		return  self.select_resnames(resnames=["TIP3", "HOH"], data=data)

	def select_lipids(self, args: list, data: np.recarray) -> np.recarray:
		''' Select lipid atoms.

		Parameters
		----------
		args : list
			List of arguments.
		data : np.recarray
			Data to be selected.

		Returns
		-------
		np.recarray
			Selected data.
		'''
		logging.debug(f"Selecting lipid atoms. Names: ['POPC', 'POPE', 'DLPE', 'DOPC', 'DMPC']")		
		return self.select_resnames(resnames = ["POPC", "POPE", "DLPE", "DOPC", "DMPC", "DMPC"], data = data)

	def select_nucleic(self, args: list, data: np.recarray) -> np.recarray:
		''' Select nucleic atoms.

		Parameters
		----------
		args : list
			List of arguments.
		data : np.recarray
			Data to be selected.

		Returns
		-------
		np.recarray
			Selected data.
		'''
		logging.debug(f"Selecting nucleic atoms. Names: ['DA', 'DT', 'DC', 'DG', 'DU', 'CYT', 'GUA', 'ADE', 'THY']")
		return self.select_resnames(resnames = ["DA", "DT", "DC", "DG", "DU", "CYT", "GUA", "ADE", "THY"], data = data)

	def __get_sels(self, selection: str) -> list:
		''' Get selections from a string.

		Parameters
		----------
		selection : str
			Selection string.

		Returns
		-------
		list
			List of selections.
		'''

		sels = []
		sp = selection.split()
		i = 0
		n = len(sp)-1

		while i <= n:
			k = sp[i].strip()

			if k in list(self.functions.keys()):
				d = {k: []}
				i += 1

				while i <= n and sp[i].strip() not in ("and", "or"):
					k = sp[i].strip()
					d[list(d.keys())[0]].append(k.strip())
					i += 1
				sels.append(d)

			elif k == "and" or k == "or":
				sels.append(k)
				i += 1

			elif k[0] == "(" and k.strip()[1:] in list(self.functions.keys()):
				k = sp[i]
				k = k[1:]
				if k in list(self.functions.keys()):
					d = {k: []}
					i += 1
					while i <= n:
						k = sp[i].strip()
						if k[-1] == ")":
							k = k[:-1]
							d[list(d.keys())[0]].append(k)
							i += 1
							break
						else:
							d[list(d.keys())[0]].append(k)
							i += 1
					sels.append(d)
			else:
				i += 1

		return sels

	def __call_fun_in_sel(self, sels: list, data: np.recarray) -> np.recarray:
		''' Call functions in selections.

		Parameters
		----------
		sels : list
			List of selections.
		data : np.recarray
			Data to be selected.

		Returns
		-------
		np.recarray
			Selected data.
		'''
		logging.debug(f"Calling functions in selections: {sels}")
		mol_data = data

		def and_fun(sel: dict, sel_data: np.recarray) -> np.recarray:
			''' Select data using the "and" operator.
			
			Parameters
			----------
			sel : dict
				Selection dictionary.
			sel_data : np.recarray
				Data to be selected.

			Returns
			-------
			np.recarray
				Selected data.
			'''
			logging.debug(f"Selecting with 'and' operator: {sel}")
			key = [*sel][0]
			sel_data = self.functions[key](sel[key], sel_data)
			return sel_data

		def or_fun(sel: dict, sel_data: np.recarray) -> np.recarray:
			''' Select data using the "or" operator.

			Parameters
			----------
			sel : dict
				Selection dictionary.
			sel_data : np.recarray
				Data to be selected.

			Returns
			-------
			np.recarray
				Selected data.
			'''
			logging.debug(f"Selecting with 'or' operator: {sel}")
			key = [*sel][0]

			sel_data = np.concatenate((sel_data, self.functions[key](sel[key], self.mol.data)))
			sel_data = np.sort(sel_data, order = 'inid').view(np.recarray)

			return sel_data

		fun = and_fun

		for sel in sels:
			if sel == "and":
				fun = and_fun
			elif sel == "or":
				fun = or_fun
			else:
				mol_data = fun(sel, mol_data)
		return mol_data

	def __call_function(self) -> np.recarray:
		''' Call the selection function.

		Returns
		-------
		np.recarray
			Selected data.
		'''
		logging.debug(f"Calling selection function for selection: {self.selection}")
		sels = self.__get_sels(selection=self.selection)
		mol_data = self.mol.data
		
		mol_data = self.__call_fun_in_sel(sels=sels, data=mol_data)

		self.sel_data = mol_data

		new_pdb = PDB()
		logging.debug(f"Creating new PDB object for selected data.")
		new_pdb.make(data = mol_data) # Depois atualziar a versão do PDB

		self.result = new_pdb

		return self.sel_data
		
	def __functions(self) -> dict:
		''' Create a dictionary of selection functions.

		Returns
		-------
		dict
			Dictionary of selection functions.
		'''

		functions = {	# Funções para seleção de resíduos
                    	"resid": self.select_resids,
						"resname": self.select_resnames,
						"residue": self.select_resinid,
						# Funções para seleção de Atomos
						"atom": self.select_atoms,
						"atomid": self.select_atoms_ids,
						"name": self.select_name,
						# Funções de seleção por coordenadas 
						"within": self.select_within,
						# Funções de seleção xxx
						"same": self.select_same,
						"not": self.select_not,
						# Funções para seleção de Cadeias
						"chain": self.select_chain,
						"segname": self.select_segname,
						# Funções para seleção por colunas
						"beta": self.select_beta,
						"occupancy": self.select_occupancy,
						# Seleções de biomoleculas
						"protein": self.select_protein,
						"backbone": self.select_backbone,
						"water": self.select_water,
						"lipids": self.select_lipids,
						"nucleic": self.select_nucleic,
						"all": self.all,
                    }
		self.functions = functions
		return functions
