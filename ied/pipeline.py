# Imports
###############################################################################
import os
import logging
import pandas as pd
from ied.ied import IED
from ied.QM.aux import AUX
from ied.MOL.PDB import PDB
from ied.QM.molden import Molden
from ied.QM.orca_out import OrcaOut
from ied.QM.overlap_matrix import OverlapMatrix
###############################################################################

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
class Pipeline: 
	def __init__(self, verbose:bool=False) -> None:
		if verbose:
			logging.getLogger().setLevel(logging.DEBUG)
			logging.debug("Verbose mode enabled inside Pipeline.")

	def __read_df(self, path:str) -> pd.DataFrame:
		"""
		Reads a DataFrame from a file.
		Parameters:
		-----------
		path: str
			Path to the file containing the DataFrame (e.g., .npy or .json).
		Returns:
		--------
			pd.DataFrame: The DataFrame read from the file.
		Raises:
		--------
			ValueError: If the file format is not supported.
			FileNotFoundError: If the file does not exist.
		"""
		ied = IED()
		if path.split(".")[-1] == "npy":
			logging.info("Reading .npy file...")
			df = ied.read_npy_to_df(path)

		elif path.split(".")[-1] == "json":
			logging.info("Reading .json file...")
			df = ied.read_json_to_df(path)
		else:
			raise ValueError("Invalid file format. Please provide a .npy or .json file.")
		
		return df


	def __summary(self) -> dict:
		"""
		Summary of the methods available in the Pipeline class.
		Returns:
		--------
			dict: A dictionary with method names as keys and method objects as values.
		"""
		summary = {}
		methods = [getattr(self, method) for method in dir(self) if callable(getattr(self, method)) and not method.startswith("__")]
		for method in methods:
			summary[method.__name__] = method
		return summary


	def __call_pipeline(self, pipeline:dict) -> None:
		"""
		Calls the methods in the pipeline dictionary.
		Parameters:
		-----------
			pipeline: dict
				A dictionary where keys are method names and values are dictionaries of parameters to pass to the methods
		Returns:
		--------
			None
		"""
		pipes = pipeline.keys()
		for pipe in pipes:

			self.__summary()[pipe](**pipeline[pipe])


	def config_file(self, file:str) -> None:
		"""
		Reads a configuration file and executes the pipeline.
		Parameters:
		-----------
			file: str
				Path to the configuration file containing the pipeline.
		Returns:
		--------
			None
		"""
		with open(file, "r") as f:
			pipeline = eval(f.read())
			self.__call_pipeline(pipeline=pipeline)
	

	def matrix(self, pdb:str, qm_sof:str, qm:str, smatrix:str="", out= "IED_out", out_format:str ="npy", gpu:bool=False, gpu_id:int=0, write:bool=True) -> tuple:
		"""
		Calculates the intermolecular electron density matrix from a PDB file and a quantum mechanics output file.
		Parameters:
		-----------
			pdb: str
				Path to the PDB file containing the structures.
			qm_sof: str
				Quantum mechanics software used (e.g., 'mopac', 'molden', 'orca', 'xTB').
			qm: str
				Path to the quantum mechanics output file (e.g., MOPAC output, Molden file, ORCA output).
			smatrix: str
				Path to the overlap matrix file (optional, default is an empty string).
			out: str
				Path to the output directory where the results will be saved (default is 'IED_out_').
			out_format: str
				Format of the output file (default is 'npy').
			procs: int
				Number of processes to use for parallel computation (default is 1).
			gpu: bool
				Whether to use GPU acceleration (default is False).
		Returns:
		--------
			None
		"""
		ied = IED()
		print(f"Path to QM output file: {qm}")

		s = OverlapMatrix()
		if smatrix != "":
			print(f"Path to S matrix file: {smatrix}")
			s.read_from_multiwfn(path=smatrix)
		
		pdb = PDB(path=pdb)
		molden = None; aux = None; orc = None

		if qm_sof == "mopac":
			print("QM software: MOPAC")
			if qm == "":
				raise ValueError("Path to MOPAC output file is required for MOPAC calculations.")

			if not os.path.exists(qm):
				raise FileNotFoundError(f"MOPAC output file not found at {qm}")
			
			if not os.path.exists(qm):
				raise FileNotFoundError(f"MOPAC output file not found at {qm}")		

			aux = AUX(); aux.read_file(path=qm, no_density_matrix=False, no_atom_charges=False)

		elif qm_sof == "molden" or qm_sof == "xTB":
			print("QM software: Molden or xTB")
			if qm == "":
				raise ValueError("Path to Molden file is required for Molden calculations.")	

			if not os.path.exists(qm):
				raise FileNotFoundError(f"Molden file not found at {qm}")

			molden = Molden(); molden.read_file(qm)
		

		elif qm_sof.strip() == "orca":
			print("QM software: ORCA")
			if qm == "":
				raise ValueError("Path to ORCA output file is required for ORCA calculations.")

			if not os.path.exists(qm):
				raise FileNotFoundError(f"ORCA output file not found at {qm}")

			orc = OrcaOut(); orc.read_file(qm)

		results = ied.matrix(
					pdb=pdb,
					aux=aux,
					molden=molden,
					orca_out=orc,
					Smatrix=s,
					path_out=out,
					out_format=out_format,
					gpu=gpu,
					gpu_id=gpu_id,
					write=write
					)
		
		return results


	def map_3D(self, pdb:str, ied:str, pdbout:str="./IED_3D_map.pdb", sel:str="all", intramol:bool=False, pep_bond:bool=False, nucleic_bond:bool=False, intrachain:bool=True, norm:str="No") -> PDB:
		a = IED()
		__pdb = PDB(pdb)
		df = self.__read_df(ied)

		results = a.map_3D(
					pdb=__pdb,
					df=df,
					sel=sel, 
					path_pdbout=pdbout,
					intramol=intramol,
					pep_bond=pep_bond,
					nucleic_bond=nucleic_bond,
					norm=norm,
					intrachain=intrachain,
					write=True
				)

		return results
	

	def two_sel(self,sel_a:str, sel_b:str,
                    pdb:str, qm_sof:str, qm:str, 
                    smatrix:str="", gpu:bool=False, gpu_id:int=0) -> tuple:
		"""Calculates the intermolecular electron density between two selections in a PDB file.

		Paremeters:
		-----------
			sel_a: str
				Selection for the first molecule (e.g., "chain A" for chain A).
			sel_b: str
				Selection for the second molecule (e.g., "chain B" for chain B).
			pdb: str
				Path to the PDB file containing the structures.
			qm_sof: str
				Quantum mechanics software used (e.g., "mopac", "molden", "orca").
			qm: str
				Path to the quantum mechanics output file (e.g., MOPAC output, Molden file, ORCA output).
			smatrix: str
				Path to the overlap matrix file (optional, default is an empty string).
		Returns:
		--------
			None
		"""
		
		if sel_a == "" or sel_b == "":
			raise ValueError("Both sel_a and sel_b must be provided.")
		if qm_sof not in ["mopac", "molden", "orca", "xTB"]:
			raise ValueError("qm_sof must be one of 'mopac', 'molden', 'orca', or 'xTB'.")
		if sel_a == sel_b:
			logging.warning(">> sel_a and sel_b must be different selections.")

		s = OverlapMatrix()
		if smatrix != "":
			if not os.path.exists(smatrix):
				raise FileNotFoundError(f"S matrix file not found at {smatrix}")
			logging.info(f"Reading S matrix file: {smatrix}")
			s.read_from_multiwfn(path=smatrix)

		ied = IED()
		logging.info(f"Reading PDB file: {pdb}")
		__pdb = PDB(path=pdb)
		molden = None; aux = None; orc = None

		if qm_sof.strip() == "mopac":
			logging.info("QM software: MOPAC")
			if qm == "":
				raise ValueError("Path to MOPAC output file is required for MOPAC calculations.")

			if not os.path.exists(qm):
				raise FileNotFoundError(f"MOPAC output file not found at {qm}")
			
			if not os.path.exists(qm):
				raise FileNotFoundError(f"MOPAC output file not found at {qm}")		

			logging.info("Reading MOPAC output file...")
			aux = AUX()
			aux.read_file(path=qm, no_density_matrix=False, no_atom_charges=False)

		elif qm_sof.strip() == "molden" or qm_sof == "xTB":
			logging.info("QM software: Molden or xTB")
			if qm == "":
				raise ValueError("Path to Molden file is required for Molden calculations.")	

			if not os.path.exists(qm):
				raise FileNotFoundError(f"Molden file not found at {qm}")

			logging.info("Reading Molden file...")
			molden = Molden()
			molden.read_file(qm)
		
		elif qm_sof.strip() == "orca":
			logging.info("QM software: ORCA")
			if qm == "":
				raise ValueError("Path to ORCA output file is required for ORCA calculations.")

			if not os.path.exists(qm):
				raise FileNotFoundError(f"ORCA output file not found at {qm}")

			logging.info("Reading ORCA output file...")
			orc = OrcaOut()
			orc.read_file(qm)

		logging.info(f"Calculating intermolecular electron density between selections '{sel_a}' and '{sel_b}' in PDB '{pdb}' using {qm_sof} software.\nRunning...")
		result = ied.two_sel(sel_a=sel_a, sel_b=sel_b, pdb=__pdb, aux=aux, molden=molden, 
										   		orca_out=orc, Smatrix=s, gpu=gpu, gpu_id=gpu_id)
		
		return result


	def plot_heatmap(self, ied:str="", norm:str="No", savefig:bool=True, plot:bool=False, annotate_res:bool=False, 
                          pdb:str="", marks:list=[], figsize="19 16", per_residue:bool=False, intramol:bool=True,
                          pep_bond:bool=True, nucleic_bond:bool=True, intrachain:bool=True, figname:str= "", cmap:str="Blues"
						  ):
		"""Plots a heatmap of the intermolecular electron density from a file.

		Paremeters:
		-----------
			ied: str
				Path to the file containing the intermolecular electron density data (e.g., .npy or .json).
			norm: str
				Normalization method to apply to the data (default is "No").
			savefig: str
				Whether to save the figure (default is "True").
			plot: str
				Whether to plot the heatmap (default is "False").
			annotate_res: str
				Whether to annotate residues in the heatmap (default is "False").
			pdb: str
				Path to the PDB file for visualization (default is an empty string).
			marks: list
				List of residue names to mark in the heatmap (default is an empty list).
			figsize: str
				Figure size for the heatmap (default is "19 16").
			per_residue: str
				Whether to plot per-residue data (default is "False").
			intramol: str
				Whether to consider intramolecular interactions (default is "True").
			pep_bond: str
				Whether to consider peptide bonds (default is "True").
			figname: str
				Name of the figure file to save (default is an empty string).
			cmap: str
				Colormap for the heatmap (default is "Blues").
			extra: dict
				Additional parameters for the heatmap (default is an empty dictionary).
		Returns:
		--------
			None
		"""
		__ied = IED()
		__ied.path_ied = ied
		__pdb = PDB(pdb)
		df = self.__read_df(ied)
		
		__ied.plot_heatmap(df=df, 
						norm=norm,
						plot=plot, 
						savefig=savefig, 
						per_residue=per_residue,
						pdb=__pdb,
						figsize=(int(figsize.split()[0]), int(figsize.split()[1])),
						annotate_res=annotate_res,
						pep_bond=pep_bond,
						intramol=intramol,
						marks=list(marks),
						figname=figname,
						nucleic_bond=nucleic_bond,
						intrachain=intrachain,
						cmap=cmap
						)


	def plot_heatmap_ref(self, ied:str, reference:str, pdb:str, reference_b:str="all",
							ncols:int= 4, norm:str="No", savefig:bool=True, plot:bool=False, 
                          	figsize="19 16", cutoff:float=.5, intramol:bool=True, multi_file:bool=False,
                          	pep_bond:bool="True", figname:str= "",nucleic_bond:bool=True, intrachain:bool=True, cmap:str="Blues", extra:dict = {}
						  ):
		"""Plots heatmaps of the intermolecular electron density from a file.

		Paremeters:
		-----------
			ied: str
				Path to the file containing the intermolecular electron density data (e.g., .npy or .json).
			reference: str
				Reference selection for the heatmap (default is an empty string).
			reference_b: str
				Reference selection for the second heatmap (default is "all").
			ncols: int
				Number of columns for the heatmap grid (default is 4).
			norm: str
				Normalization method to apply to the data (default is "No").
			savefig: bool
				Whether to save the figure (default is "True").
			plot: bool
				Whether to plot the heatmap (default is "False").
			pdb: str
				Path to the PDB file for visualization (default is an empty string).
			figsize: str
				Figure size for the heatmap (default is "19 16").
			cutoff: float
				Cutoff value for the heatmap (default is 0.5).		
			intramol: bool
				Whether to consider intramolecular interactions (default is "True").
			nucleic_bond: bool
				Whether to consider nucleic bond interactions (default is "True").
			intrachain: bool
				Whether to consider intrachain interactions (default is "True").
			multi_file: str
				Whether to plot multiple files (default is "False").
			pep_bond: bool
				Whether to consider peptide bonds (default is "True").
			figname: str
				Name of the figure file to save (default is an empty string).
			cmap: str
				Colormap for the heatmap (default is "Blues").
			extra: dict
				Additional parameters for the heatmap (default is an empty dictionary).
		Returns:
		--------
			None"""
		
		__ied = IED()
		__pdb = PDB(pdb)
		df = self.__read_df(ied)
		
		__ied.plot_heatmap_ref(df=df, 
						norm=norm,
						ncols=ncols,
						ref=reference,
						ref_b=reference_b,
						plot=plot, 
						savefig=savefig, 
						pdb=__pdb,
						figsize=(int(figsize.split()[0]), int(figsize.split()[1])),
						pep_bond=pep_bond,
						intramol=intramol,
						figname=figname,
						nucleic_bond=nucleic_bond,
						intrachain=intrachain,
						cutoff=cutoff,
						multi_file=multi_file,
						cmap=cmap,
						)
			

	def radial_distribution(self, sel:str, pdb:str, ied:str="", out:str="./IED_radial_distribution.json", selb:str = "all",
								dis_max:float=20.0, intramol:str="False", pep_bond:str="False", line:bool=False, savefig:bool=True, plot:bool=False, interachain:bool=False
								):
		""" Calculates the radial distribution of intermolecular electron density around a selection of atoms in a PDB file.

		Paremeters:
		-----------
			sel: str
				Selection of atoms for which to calculate the radial distribution (e.g., "chain A").
			pdb: str
				Path to the PDB file containing the structures.
			ied: str
				Path to the file containing the intermolecular electron density data (e.g., .npy or .json).
			out: str
				Path to the output file where the radial distribution will be saved, .json and .png (default is "./IED_radial_distribution.json").
			dis_max: float
				Maximum distance for the radial distribution (default is 20.0).
			intramol: str
				Whether to consider intramolecular interactions (default is "False").
			pep_bond: str
				Whether to consider peptide bonds (default is "False").	
		Returns:
		--------
			None
		"""
		__ied = IED()
		__pdb = PDB(pdb)
		df = self.__read_df(ied)
		
		__ied.radial_distribution(df=df,
							  		sel=sel,
									selb=selb,
									path_out=out,
									dis_max=dis_max,
									intramol=self.__str_to_bool(intramol),
									pep_bond=self.__str_to_bool(pep_bond),
									pdb=__pdb,
									line=line,
									intrachain=interachain,
									savefig=savefig,
									plot=plot
									)


	def __str_to_bool(self, s):
		if s.lower() in ['True', 'true', '1', 't', 'y', 'yes']:
			return True
		elif s.lower() in ['false', '0', 'f', 'n', 'no']:
			return False
		else:
			raise ValueError(f"Valor inválido para booleano: '{s}'")
		