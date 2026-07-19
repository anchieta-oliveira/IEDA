# Imports
###############################################################################
import re
import numpy as np
from ieda.QM.MO import MO
from ieda.QM.GTO import GTO
from ieda.core import read_file
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
This project is licensed under the MIT License.
'''

# Classes
###############################################################################
class OrcaOut: 
	def __init__(self, 
				title:str = "", 
				name:str = ""
				):
		self.title:str = title
		self.name:str = name
		self.program_version:str = ""
		self.input_file:str = ""
		self.atom_el:list = []
		self.ao_atomindex:np.array = np.array([])
		self.atom_symtype:list = []
		self.num_electrons:int = 0
		self.dipole_debye:float = ""
		self.number_scf_cycles:int = 0
		self.overlap_matrix:np.array = np.array([], dtype=np.float64)
		self.density_matrix:np.array = np.array([], dtype=np.float64)
		self.cartesian_coordinates_angstrom:np.array = np.empty((0, 3))
		self.cartesian_coordinates_au:np.array = np.empty((0, 3))
		self.atom_no:np.array  = np.array([])
		self.orbital_no:np.array  = np.array([])
		self.labels:np.array  = np.array([])
		self.atom_mass:np.array  = np.array([], dtype=np.float64)
		self.atom_frag:np.array  = np.array([], dtype=np.float64)
		self.za:np.array  = np.array([], dtype=np.float64)
		self.basis_set_information:str = ""
		self.basis_set_in_input_format:list = [] # GTOs
		self.GTOs:list = []
		self.auxiliaryj_basis_set_information:list = []
		self.shark_integral_package:str = ""
		self.dft_grid_generation:str = ""
		self.cosx_grid_generation:str = ""
		self.scf_iterations:str = ""
		self.scf_n_cycles:int = 0
		self.orbital_energies_ev:np.array = np.array([], dtype=np.float64)
		self.orbital_energies_eh:np.array = np.array([], dtype=np.float64)
		self.molecular_orbital_occupancies:np.array = np.array([])
		self.molecular_orbitals:np.array = np.array([])
		self.MOs:list = []
		self.mulliken_atomic_charges:np.array = np.array([], dtype=np.float64)
		self.sum_mulliken_atomic_charges:float = .0
		self.final_single_point_energy:float = .0
		self.cartesian_gradient:np.array = np.empty((0, 3))
		self.norm_cartesian_gradient:float = .0
		self.rms_gradient:float = .0
		self.max_gradient:float = .0
		self.timings:str = ""
		self.total_dipole_moment:np.array = np.array([], dtype=np.float64)
		self.dipole_electronic_contribution:np.array = np.array([], dtype=np.float64)
		self.dipole_nuclear_contribution:np.array = np.array([], dtype=np.float64)
		self.dipole_magnitude_au:float = .0
		self.dipole_magnitude_debye:float = .0
		self.total_run_time:str = ""
	

	def read_file(self, path:str):
		self.name = path.split("/")[-1].split(".")[-2]
		"""with open(path, "r") as file:
			lines = file.readlines()"""
		
		text = read_file(path)
		lines = text.splitlines(True)
		
		input_file = False
		cartesian_coordinates_angstrom = False
		cartesian_coordinates_au = False
		basis_set_information = False
		basis_set_in_input_format = False
		shark_integral_package = False
		cosx_grid_generation = False
		scf_iterations = False
		orbital_energies = False
		overlap_matrix = False
		dft_grid_generation = False
		molecular_orbitals = False
		mulliken_atomic_charges = False
		cartesian_gradient = False
		timings = False
		dipole_moment = False
		cc_method = False


		for line in lines:
			if "Program Version" in line:
				self.program_version = line.split()[2]

			# INPUT FILE
			elif "INPUT FILE" in line:
				input_file = True

			# INPUT FILE
			elif input_file:
				if not "END OF INPUT" in line:
					self.input_file += line
				else:
					if ("DLPNO" in self.input_file) or ("CCSD" in self.input_file) or ("dlpno" in self.input_file) or ("ccsd" in self.input_file):
						cc_method = True
					input_file = False

			# COORDINATES (ANGSTROEM)
			elif "COORDINATES (ANGSTROEM)" in line:
				cartesian_coordinates_angstrom = True

			elif cartesian_coordinates_angstrom:
				if not "----" in line and not "\n" == line:
					self.cartesian_coordinates_angstrom = np.vstack((self.cartesian_coordinates_angstrom, [float(line.split()[1]), float(line.split()[2]), float(line.split()[3])]))
					self.atom_el = np.append(self.atom_el, [line.split()[0]])
				elif "\n" == line:
					cartesian_coordinates_angstrom = False

			# CARTESIAN COORDINATES (A.U.)
			elif "CARTESIAN COORDINATES (A.U.)" in line:
				cartesian_coordinates_au = True

			elif cartesian_coordinates_au:
				if not "----" in line and not "\n" == line and not "NO" in line:
					self.atom_no = np.append(self.atom_no, [int(line.split()[0])])
					self.labels = np.append(self.labels, [line.split()[1]])
					self.za = np.append(self.za, [line.split()[2]])
					self.atom_frag = np.append(self.atom_frag, [line.split()[3]])
					self.atom_mass = np.append(self.atom_mass, [float(line.split()[4])])

					self.cartesian_coordinates_au = np.vstack((self.cartesian_coordinates_au, [float(line.split()[5]), float(line.split()[6]), float(line.split()[7])]))
				elif "\n" == line:
					cartesian_coordinates_au = False

			#BASIS SET INFORMATION
			elif "BASIS SET INFORMATION\n" == line:
				basis_set_information = True
			
			elif basis_set_information:
				if not "BASIS SET" in line:
					self.basis_set_information += line
				else:
					basis_set_information = False

			# BASIS SET IN INPUT FORMAT
			elif "BASIS SET IN INPUT FORMAT\n" == line:
				basis_set_in_input_format = True
			
			elif basis_set_in_input_format:
				gto = GTO()
				if "NewGTO" in line:
					gto.clear()					
					gto.title = line.split()[-1]
				elif not "end" in line and len(line.split()) == 2:
					gto.shell_label.append(line.split()[0])
					gto.number_of_primitives.append(int(line.split()[1]))
				elif not "end" in line and len(line.split()) == 3:
					gto.exponent_primitive.append(float(line.split()[1]))
					gto.contraction_coefficient.append(float(line.split()[2]))
				elif "end" in line:
					tmp_exp = []
					tmp_con = []
					old_exp = gto.exponent_primitive
					old_con = gto.contraction_coefficient
					gto.exponent_primitive = []
					gto.contraction_coefficient = []
					i = 0
					for n in gto.number_of_primitives:
						for j in range(n):
							tmp_exp.append(old_exp[i])
							tmp_con.append(old_con[i])
							i +=1 
						gto.exponent_primitive.append(tmp_exp) 
						gto.contraction_coefficient.append(tmp_con)
						tmp_exp = []; tmp_con = []	
					self.GTOs.append(gto.copy())
					gto.clear()
					del gto
				elif len(line.split()) > 3 and not "#" in line:
					basis_set_in_input_format = False

			# SHARK INTEGRAL PACKAGE
			elif "SHARK INTEGRAL PACKAGE" in line:
				shark_integral_package = True
			elif shark_integral_package:
				if not "Total time" in line:
					self.shark_integral_package += line
				else:
					self.shark_integral_package += line
					match = re.search(r'Number of basis functions\s+\.*\s+(\d+)', self.shark_integral_package)
					self.number_basis = int(match.group(1))
					match = re.search(r'Number of shells\s+\.*\s+(\d+)', self.shark_integral_package)
					self.num_electrons = int(match.group(1))
					shark_integral_package = False
					# Atualizar o numero de Basis se for CC
					if cc_method:
						self.number_basis = int(next((element for element in reversed(lines) if 'Basis functions' in element), None).split()[2])
			
			# OVERLAP MATRIX
			elif "OVERLAP MATRIX" in line:
				overlap_matrix = True
				self.overlap_matrix = np.zeros((self.number_basis, self.number_basis))
				lin = -2
				n_col = 0
				b = 0
			elif overlap_matrix:
				if (not "---------------" in line) and (not "Time for model" in line) and (not "model" in line) and (not "Atom" in line) and (not "SCF ITERATIONS" in line):
					lin += 1
					if lin > 0:
						line_split = line.split()[1:]
						for ic, v in enumerate(line_split):
							n_col = 1 + (6*b) + ic 
							self.overlap_matrix[lin-1, n_col-1] = float(v)
					if lin == self.number_basis:
						b += 1 
						lin = -1
				else:
					overlap_matrix = False

			# DFT GRID GENERATION
			elif "DFT GRID GENERATION" in line:
				dft_grid_generation = True
			
			elif dft_grid_generation:
				if (not "Time" in line) and ("model" in line):
					self.dft_grid_generation += line
				else:
					self.dft_grid_generation += line
					dft_grid_generation = False
			
			# COSX GRID GENERATION
			elif "COSX GRID GENERATION" in line:
				cosx_grid_generation = True
			elif cosx_grid_generation:
				if not "SCF ITERATIONS" in line:
					self.cosx_grid_generation += line
				else:
					cosx_grid_generation = False

			# SCF ITERATIONS
			elif "SCF ITERATIONS" in line:
				scf_iterations = True
			elif scf_iterations:
				if (not "ORBITAL ENERGIES" in line) and (not "**** DENSITY"):
					self.scf_iterations += line
				else:
					scf_iterations = False

			# ORBITAL ENERGIES
			elif "ORBITAL ENERGIES" in line:
				orbital_energies = True
			elif ("NO" in line and orbital_energies) or line == "\n" or "-----" in line:
				pass
			elif orbital_energies:
				if (not len(line.split()) < 4) and (not "----" in line) and (not line == "\n"):
					self.orbital_energies_ev = np.append(self.orbital_energies_ev, [float(line.split()[3])])
					self.orbital_energies_eh = np.append(self.orbital_energies_eh, [float(line.split()[2])])
					self.molecular_orbital_occupancies = np.append(self.molecular_orbital_occupancies, [float(line.split()[1])])
					self.orbital_no = np.append(self.orbital_no, [float(line.split()[0])])
				else:
					orbital_energies = False
					if "MOLECULAR ORBITALS" in line:
						molecular_orbitals = True
						self.molecular_orbitals = np.zeros((self.number_basis, self.number_basis), dtype=np.float64)
						lin = -4
						b = 0
						self.MOs = [MO() for _ in range(self.number_basis)]
						
			# MOLECULAR ORBITALS
			elif molecular_orbitals:
				if not "*****" in line:
					lin += 1
					line_split = line.split()
					for i, v in enumerate(line_split):
						if lin == -2:
							mo_id = (6*b) + i 
							self.MOs[mo_id].energy = float(v)
						elif lin == -1:
							mo_id = (6*b) + i 
							self.MOs[mo_id].occupation = int(float(v))

						elif lin >= 0:
							if i == 0:
								if len(self.ao_atomindex) < self.number_basis:
									self.ao_atomindex = np.append(self.ao_atomindex, [int(''.join([char for char in v if char.isdigit()]))])
							elif i == 1:
								if len(self.atom_symtype) < self.number_basis:
									self.atom_symtype.append(v)
								
							elif i >= 2 and lin < self.number_basis:
								mo_id = (6*b) + i -2
								self.MOs[mo_id].coefficients.append(float(v))
								self.MOs[mo_id].id = mo_id
								self.MOs[mo_id].ao_number.append(lin)
					if lin == self.number_basis:
						b += 1 
						lin = -3
	
				else:
					for orb in self.MOs:
						orb.symtype = self.atom_symtype
						orb.ao_atomindex = self.ao_atomindex
					molecular_orbitals = False
			
			# MULLIKEN ATOMIC CHARGES
			elif "MULLIKEN ATOMIC CHARGES" in line:
				mulliken_atomic_charges = True

			elif mulliken_atomic_charges:
				if not "Sum of atomic" in line:
					self.mulliken_atomic_charges = np.append(self.mulliken_atomic_charges, float(line.split()[-1]))
				else:
					mulliken_atomic_charges = False
					self.sum_mulliken_atomic_charges = float(line.split()[-1])
			
			# FINAL SINGLE POINT ENERGY
			elif "FINAL SINGLE POINT ENERGY" in line:
				self.final_single_point_energy = float(line.split()[-1])

			# CARTESIAN GRADIENT
			elif "CARTESIAN GRADIENT" in line:
				cartesian_gradient = True
			
			elif cartesian_gradient:
				if ":" in line and len(line.split()) == 6:
					self.cartesian_gradient = np.vstack((self.cartesian_gradient, [float(line.split()[-3]), float(line.split()[-2]), float(line.split()[-1])]))
				elif "Norm of the cartesian gradient" in line:
					self.norm_cartesian_gradient = float(line.split()[-1])
				elif "RMS gradient" in line:
					self.rms_gradient = float(line.split()[-1])
				elif "MAX gradient" in line:
					self.max_gradient = float(line.split()[-1])
					cartesian_gradient = False
			
			# TIMINGS
			elif "TIMINGS" in line:
				timings = True
			elif timings:
				if not "*******" in line:
					self.timings += line
				else:
					timings = False
				
			# DIPOLE MOMENT
			elif "DIPOLE MOMENT" in line:
				dipole_moment = True
			elif dipole_moment:
				if "Electronic contribution" in line:
					self.dipole_electronic_contribution = np.array([float(line.split()[-3]), float(line.split()[-2]), float(line.split()[-1])])
				elif "Nuclear contribution" in line:
					self.dipole_nuclear_contribution = np.array([float(line.split()[-3]), float(line.split()[-2]), float(line.split()[-1])])
				elif "Total Dipole Moment" in line:
					self.total_dipole_moment = np.array([float(line.split()[-3]), float(line.split()[-2]), float(line.split()[-1])])
				elif "Magnitude (a.u.)" in line:
					self.dipole_magnitude_au = float(line.split()[-1])
				elif "Magnitude (Debye)" in line:
					self.dipole_magnitude_debye = float(line.split()[-1])
					dipole_moment = False

			# TOTAL RUN TIME	
			elif "TOTAL RUN TIME:" in line:
				self.total_run_time = line
				lines.clear()
				




					




