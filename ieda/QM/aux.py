# Imports
###############################################################################
import numpy as np
from ieda.QM.MO import MO
from ieda.QM.GTO import GTO
from ieda.MOL.atom import Atom
from ieda.core import read_file
from ieda.QM.density_matrix import DensityMatrix
from ieda.QM.overlap_matrix import OverlapMatrix
##############################################################################

"""
Format: http://openmopac.net/manual/auxiliary.html
"""

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
class AUX: 
	def __init__(self, 
			    title:str = "", 
				name:str = ""
				):
		
		self.title:str = title
		self.name:str = name
		self.mopac_version:str = ""
		self.date:str = ""
		self.method:str = ""
		self.keywords:str = ""
		self.comments:str = ""
		self.atom_el:list = []
		self.atom_core:list = []
		self.atom_x_angstroms:list = []
		self.ao_atomindex:list = []
		self.atom_symtype:list = []
		self.ao_zeta:list = []
		self.atom_pqn:list = []
		self.num_electrons:int = 0
		self.empirical_formula:str = ""
		self.heat_of_formation:float = .0
		self.gradient_norm_kcal_mol_angstrom:float = .0
		self.energy_electronic:float = .0
		self.energy_nuclear:float = .0
		self.point_group:str = ""
		self.dipole_debye:float = ""
		self.area_square_angstroms:float = .0
		self.volume_cubic_angstroms:float = .0
		self.diel_ener_ev:float = .0
		self.spin_component:float = .0
		self.total_spin:float = .0
		self.number_scf_cycles:int = 0
		self.atom_x_opt_angstroms:list = []
		self.atom_charges:list = []
		self.overlap_matrix:list = []
		self.eigenvectores:list = []
		self.total_density_matrix:list = []
		self.mo_symmetry_labels:list = []
		self.eigenvalues:list = []
		self.molecular_orbital_occupancies:list = []
		self.cpu_time_sec = .0
		self.cpu_time_seconds:float = .0
		self.molecular_weight_amu:float = .0
		self.gradients_kcal_mol_angstrom:list = []
		self.set_of_mos:str = ""
		self.lmo_vectors:list = []
		self.density_matrix: list = []
		self.lmo_energy_levels:list = []

	

	def read_file(self, path:str,
			   			no_density_matrix = True,
						no_atom_x_angstroms = True,
						no_atom_el = True,
						no_atom_core = True,
						no_ao_atomindex = True,
						no_atom_symtype = True,
						no_ao_zeta = True,
						no_atom_pqn = True,
						no_atom_x_opt_angstroms = True,
						no_atom_charges = True,
						no_overlap_matrix = True,
						no_eigenvectores = True,
						no_total_density_matrix = True,
						no_mo_symmetry_labels = True,
						no_eigenvalues = True,
						no_molecular_orbital_occupancies = True,
						no_gradients_kcal_mol_angstrom = True,
						no_lmo_energy_levels = True,
						no_lmo_vectors = True
						):
		self.name = path.split("/")[-1].split(".")[-2]
		text = read_file(path)
		lines = text.splitlines(True)
			
		start_input = False
		geometry_opt = False
		final_scf = False
		
		at_el_tmp = False
		at_core_tmp = False
		at_angs_tmp = False
		ao_at_id_tmp = False
		at_sym_typ_tmp = False
		ao_zeta_tmp = False
		at_pqn_tmp = False
		at_x_opt_tmp = False
		at_chg_tmp = False
		overlap_matrix_tmp = False
		gradients_tmp = False
		density_matrix_tmp = False
		lmo_energy_levels_tmp = False
		mos_occup_tmp = False
		lmo_vectors_tmp = False

		for  line in lines:
			line_split = line.split()
			if " START OF MOPAC FILE" in line:
				start_input = True
			if start_input:
				if "MOPAC_VERSION" in line:
					self.mopac_version = line.split("=")[1].strip("\n")
				elif "DATE" in line:
					self.date = line.split("=")[1].strip("\n")
				elif "METHOD" in line:
					self.method = line.split("=")[1].strip("\n")
				elif "TITLE" in line:
					self.title = line.split("=")[1].strip("\n")
				elif "KEYWORDS" in line:
					self.keywords = line.split("KEYWORDS")[1].strip("\n")
				elif "COMMENTS" in line:
					self.comments = line.split("=")[1].strip("\n")
				elif "ATOM_EL" in line:
					at_el_tmp = True
				elif at_el_tmp and no_atom_el:
					if not "=" in line:
						self.atom_el.extend(line_split)
					else:
						at_el_tmp = False
				if "ATOM_CORE" in line and no_atom_core:
					at_core_tmp = True
				elif at_core_tmp:
					if not "=" in line:
						self.atom_core.extend(line_split)
					else:
						at_core_tmp = False
				if "ATOM_X:ANGSTROMS" in line and no_atom_x_angstroms:
					at_angs_tmp = True
				elif at_angs_tmp:
					if not "=" in line:
						self.atom_x_angstroms.append((float(line_split[0]), float(line_split[1]), float(line_split[2])))
					else:
						at_angs_tmp = False
				if "AO_ATOMINDEX" in line and no_ao_atomindex:
					ao_at_id_tmp = True
				elif ao_at_id_tmp:
					if not "=" in line:
						self.ao_atomindex.extend(list(map(lambda id: int(id), line_split)))
					else:
						ao_at_id_tmp = False
				
				if "ATOM_SYMTYPE" in line and no_atom_symtype:
					at_sym_typ_tmp = True
				elif at_sym_typ_tmp:
					if not "=" in line:
						self.atom_symtype.extend(line_split)
					else:
						at_sym_typ_tmp = False
				if "AO_ZETA" in line and no_ao_zeta:
					ao_zeta_tmp = True
				elif ao_zeta_tmp:
					if not "=" in line:
						self.ao_zeta.extend(list(map(lambda zeta: float(zeta), line_split)))
					else:
						ao_zeta_tmp = False
				
				if "ATOM_PQN" in line and no_atom_pqn:
					at_pqn_tmp = True
				elif at_pqn_tmp:
					if not "=" in line:
						self.atom_pqn.extend(list(map(lambda pqn: int(pqn), line_split)))
					else:
						at_pqn_tmp = False
				elif "NUM_ELECTRONS" in line:
					self.num_electrons = line.split("=")[1]
				elif "EMPIRICAL_FORMULA" in line:
					self.empirical_formula = line.split("EMPIRICAL_FORMULA")[1].strip("\n")
				
			if "Geometry optimization" in line:
				geometry_opt = True
				start_input = False
				
			if geometry_opt:
				pass
			
			if "Geometry optimization" in line:
				geometry_opt = False
				final_scf = True
			if final_scf:
				if "HEAT_OF_FORMATION" in line:
					self.heat_of_formation = line.split("=")[1].strip("\n")
				elif "GRADIENT_NORM" in line:
					self.gradient_norm_kcal_mol_angstrom = line.split("=")[1].strip("\n")
				elif "POINT_GROUP" in line:
					self.point_group = line.split("=")[1].strip("\n")
				elif "AREA:SQUARE ANGSTROMS" in line:
					self.area_square_angstroms = line.split("=")[1].strip("\n")
				elif "VOLUME:CUBIC ANGSTROMS" in line:
					self.volume_cubic_angstroms = line.split("=")[1].strip("\n")
				elif "DIEL_ENER:EV" in line:
					self.diel_ener_ev = line.split("=")[1].strip("\n")
				elif "SPIN_COMPONENT" in line:
					self.spin_component = line.split("=")[1].strip("\n")
				elif "TOTAL_SPIN" in line:
					self.total_spin = line.split("=")[1].strip("\n")
				elif "NUMBER_SCF_CYCLES" in line:
					self.number_scf_cycles = int(line.split("=")[1])
				elif "CPU_TIME:SEC" in line:
					self.cpu_time_sec = line.split("=")[1].strip("\n")
				elif "MOLECULAR_WEIGHT:AMU" in line:
					self.molecular_weight_amu = line.split("=")[1].strip("\n")
				elif "ATOM_X_OPT:ANGSTROMS" in line and no_atom_x_opt_angstroms:
					at_x_opt_tmp = True
				elif at_x_opt_tmp:
					if not "=" in line:
						self.atom_x_opt_angstroms.append((float(line_split[0]), float(line_split[1]), float(line_split[2])))
					else:
						at_x_opt_tmp = False
				if "ATOM_CHARGES" in line and no_atom_charges:
					at_chg_tmp = True
				elif at_chg_tmp:
					if not "=" in line:
						self.atom_charges.extend(list(map(lambda chg: float(chg), line_split)))
					else:
						at_chg_tmp = False
				if "GRADIENTS:KCAL/MOL/ANGSTROM" in line and no_gradients_kcal_mol_angstrom:
					gradients_tmp = True
				elif gradients_tmp:
					if not "=" in line:
						self.gradients_kcal_mol_angstrom.extend(list(map(lambda grad: float(grad), line_split)))
					else:
						gradients_tmp = False
					
				if "OVERLAP_MATRIX" in line and no_overlap_matrix:
					overlap_matrix_tmp = True
				elif overlap_matrix_tmp and not "#" in line:
					if not "=" in line:
						self.overlap_matrix.extend(list(map(lambda s: float(s), line_split)))
					else:
						overlap_matrix_tmp = False
				if "SET_OF_MOS" in line:
					self.set_of_mos = line.split("=")[1].strip("\n")
				
				if ("LMO_VECTORS" in line or "EIGENVECTORS" in line) and no_lmo_vectors:
					lmo_vectors_tmp = True
				elif lmo_vectors_tmp:
					if not "=" in line:
						self.lmo_vectors.extend(list(map(lambda vec: float(vec), line_split)))
					else:
						lmo_vectors_tmp = False
				
				if "DENSITY_MATRIX" in line and no_density_matrix:
					density_matrix_tmp = True
				elif density_matrix_tmp and not "#" in line:
					if not "=" in line:
						self.density_matrix.extend(list(map(lambda dm: float(dm), line_split)))
					else:
						density_matrix_tmp = False
				if ("LMO_ENERGY_LEVELS" in line or "EIGENVALUES" in line) and no_lmo_energy_levels:
					lmo_energy_levels_tmp = True
				elif lmo_energy_levels_tmp:
					if not "=" in line:
						self.lmo_energy_levels.extend(list(map(lambda eng: float(eng), line_split)))
					else:
						lmo_energy_levels_tmp = False
				
				if "MOLECULAR_ORBITAL_OCCUPANCIES" in line and no_molecular_orbital_occupancies:
					mos_occup_tmp = True
				elif mos_occup_tmp:
					if not "=" in line:
						self.molecular_orbital_occupancies.extend(list(map(lambda mo_oc: float(mo_oc), line_split)))
					else:
						mos_occup_tmp = False
				elif "CPU_TIME:SECONDS" in line:
					self.cpu_time_seconds = float(line.split("=")[1])
				
				if "END OF MOPAC FILE" in line:
					break

			
		
	def parse_Smatrix(self, S:OverlapMatrix):
		n = len(self.ao_atomindex) # Tamanho da matriz
		self.overlap_matrix = []
		for i in range(n):
			for j in range(i+1):
				self.overlap_matrix.append(S.matrix[i][j])


	def parse_Pmatrix(self, P:DensityMatrix):
		n = len(self.ao_atomindex) # Tamanho da matriz
		self.density_matrix = []
		for i in range(n):
			for j in range(i+1):
				self.density_matrix.append(P.matrix[i][j])

		

	def write(self, path:str):
		ini = """
 ####################################
 #                                  #
 #       Start of Input data        #
 #                                  #
 ####################################\n"""
		with open(path, "w") as file:
			file.write(" START OF MOPAC FILE")
			file.write(ini)
			file.write(f" MOPAC_VERSION={self.mopac_version}\n")
			file.write(f" DATE={self.date}\n")
			file.write(f" METHOD={self.method}\n")
			file.write(f" TITLE={self.title}\n")
			file.write(f" KEYWORDS={self.keywords}\n")
			file.write(f" COMMENTS={self.comments}\n")
			file.write(f" ATOM_EL=[{len(self.atom_el)}]=\n")
			n = 0
			for el in self.atom_el:
				file.write(f"  {el}")
				n += 1
				if n == 40:
					file.write(f"\n")
					n = 0
			if n !=0:
				file.write("\n")
				n = 0
			file.write(f" ATOM_CORE[0{len(self.atom_core)}]=\n")
			for core in self.atom_core:
				file.write(f"  {core}")
				n += 1
				if n == 40:
					file.write(f"\n")
					n = 0
			if n !=0:
				file.write("\n")
			file.write(f" ATOM_X:ANGSTROMS[0{len(self.atom_x_angstroms)*3}]=\n")
			for coord in self.atom_x_angstroms:
				file.write(f"{coord[0]:9.4f}{coord[1]:10.4f}{coord[2]:10.4f}\n")

			file.write(f" AO_ATOMINDEX[0{len(self.ao_atomindex)}]=\n")
			n = 0
			for at_id in self.ao_atomindex:
				file.write(f"{at_id:^5}")
				n += 1
				if n == 24:
					file.write(f"\n")
					n = 0
			if n !=0:
				file.write("\n")
			file.write(f" ATOM_SYMTYPE[0{len(self.atom_symtype)}]=\n ")
			n = 0
			for symtype in self.atom_symtype:
				file.write(f"{symtype:^3}")
				n += 1
				if n == 40:
					file.write(f"\n ")
					n = 0
			if n !=0:
				file.write("\n")
				n= 0 
			file.write(f" AO_ZETA[0{len(self.ao_zeta)}]=\n")
			for zeta in self.ao_zeta:
				file.write(f"{zeta:8.4f}")
				n += 1
				if n == 10:
					file.write(f"\n")
					n = 0

			if n !=0:
				file.write("\n")
				n = 0
			file.write(f" ATOM_PQN[0{len(self.atom_pqn)}]=\n ")
			for pqn in self.atom_pqn:
				file.write(f"{pqn:^2}")
				n += 1
				if n == 40:
					file.write(f"\n ")
					n = 0
			if n !=0:
				file.write("\n")
				n = 0
			file.write(f"NUM_ELECTRONS={self.num_electrons}\n")
			file.write(f" EMPIRICAL_FORMULA={self.empirical_formula}")
			geo_ini = """ 
 ####################################
 #                                  #
 #      Geometry optimization       #
 #                                  #
 ####################################"""
			file.write(geo_ini)
			scf_ini = """
 ####################################
 #                                  #
 #        Final SCF results         #
 #                                  #
 ####################################\n"""
			file.write(scf_ini)
			file.write(f" HEAT_OF_FORMATION:KCAL/MOL={self.heat_of_formation}\n")
			file.write(f" GRADIENT_NORM:KCAL/MOL/ANGSTROM={self.gradient_norm_kcal_mol_angstrom}\n")
			file.write(f" POINT_GROUP={self.point_group}\n")
			file.write(f" AREA:SQUARE ANGSTROMS={self.area_square_angstroms}\n")
			file.write(f" VOLUME:CUBIC ANGSTROMS={self.volume_cubic_angstroms}\n")
			file.write(f" DIEL_ENER:EV={self.diel_ener_ev}\n")
			file.write(f" SPIN_COMPONENT={self.spin_component}\n")
			file.write(f" TOTAL_SPIN={self.total_spin}\n")
			file.write(f" NUMBER_SCF_CYCLES={self.number_scf_cycles}\n")
			file.write(f" CPU_TIME:SEC={self.cpu_time_sec}\n")
			file.write(f" MOLECULAR_WEIGHT:AMU={self.molecular_weight_amu}\n")
			file.write(f" ATOM_X_OPT:ANGSTROMS[0{len(self.atom_x_opt_angstroms)*3}]=\n")
			for coord in self.atom_x_opt_angstroms:
				file.write(f"{coord[0]:9.4f}{coord[1]:10.4f}{coord[2]:10.4f}\n")
			
			file.write(f" ATOM_CHARGES[0{len(self.atom_charges)}]=\n")
			n = 0
			for chg in self.atom_charges:
				if chg > 0:
					file.write(f"{chg:+9.5f}")
				else:
					file.write(f"{chg:9.5f}")
				n += 1
				if n == 10:
					file.write(f"\n")
					n = 0
			if 0 < len(self.gradients_kcal_mol_angstrom):
				if n !=0:
					file.write("\n")
					n = 0
				file.write(f" GRADIENTS:KCAL/MOL/ANGSTROM[0{len(self.gradients_kcal_mol_angstrom)}]=\n")
				for grad in self.gradients_kcal_mol_angstrom:
					file.write(f"{grad:10.4f}")
					n += 1
					if n == 10:
						file.write(f"\n")
						n = 0

			if n !=0:
				file.write("\n")
				n = 0
			file.write(f" OVERLAP_MATRIX[0{len(self.overlap_matrix)}]=\n #  Lower half triangle only\n")
			for s in self.overlap_matrix:
				file.write(f"{s:9.4f}")
				n += 1
				if n == 10:
					file.write(f"\n")
					n = 0

			if n !=0:
				file.write("\n")
				n = 0
			file.write(f" SET_OF_MOS={self.set_of_mos}\n")
			file.write(f" LMO_VECTORS[0{len(self.lmo_vectors)}]=\n")
			n2 = 0
			for lmo in self.lmo_vectors:
				file.write(f"{lmo:9.4f}")
				n += 1
				n2 +=1
				if n == 10:
					file.write("\n")
					n = 0
				if n2 == len(self.atom_symtype):
					file.write("\n")
					n2 = 0
					n = 0 
			
			if n !=0:
				file.write("\n")
				n = 0
			file.write(f" DENSITY_MATRIX[0{len(self.density_matrix)}]=\n #  Lower half triangle only.\n")
			for density in self.density_matrix:
				file.write(f"{density:9.4f}")
				n += 1 
				if n == 10:
					file.write("\n")
					n = 0
			
			if n !=0:
				file.write("\n")
				n = 0
			file.write(f" LMO_ENERGY_LEVELS[0{len(self.lmo_energy_levels)}]=\n")
			for energy in self.lmo_energy_levels:
				file.write(f"{energy:8.3f}")
				n += 1 
				if n == 10:
					file.write("\n")
					n = 0

			if n !=0:
				file.write("\n")
				n = 0
			file.write(f" MOLECULAR_ORBITAL_OCCUPANCIES[0{len(self.molecular_orbital_occupancies)}]=\n")
			for mo_occ in self.molecular_orbital_occupancies:
				file.write(f"{mo_occ:7.4f}")
				n += 1 
				if n == 10:
					file.write("\n")
					n = 0

			if n !=0:
				file.write("\n")
				n = 0	
			file.write(f" CPU_TIME:SECONDS[1]={self.cpu_time_seconds:12.2f}\n")
			file.write(" END OF MOPAC FILE")
	
	def get_MOs_objs(self) -> list:
		i = 1
		MOs = []
		num_orb = len(self.atom_symtype)
		for energy, occ in zip(self.lmo_energy_levels, self.molecular_orbital_occupancies):
			mo = MO(id=i, energy=energy, occupation=occ)
			mo.coefficients = list(map(lambda c: c, self.lmo_vectors[(i-1)*num_orb:i*num_orb]))
			mo.spin = "Alpha"
			MOs.append(mo)
			i += 1

		return MOs


	def get_S_objs(self) -> OverlapMatrix:
		n = len(self.atom_symtype) # Tamanho da matriz
		matrix = [[0] * n for _ in range(n)]
		index = 0
		for i in range(n):
			for j in range(i + 1):
				matrix[i][j] = self.overlap_matrix[index]
				matrix[j][i] = self.overlap_matrix[index]
				index += 1
		s = OverlapMatrix(matrix=matrix)
		return s
	

	def get_P_objs(self) -> DensityMatrix:
		n = len(self.atom_symtype) # Tamanho da matriz
		matrix = np.array([[0] * n for _ in range(n)])
		index = 0
		for i in range(n):
			for j in range(i + 1):
				matrix[i][j] = self.density_matrix[index]
				matrix[j][i] = self.density_matrix[index]
				index += 1
		s = DensityMatrix(matrix=matrix)
		return s


	def get_atoms(self) -> list:
		atoms = []
		id = 0
		for coord, el, an in zip(self.atom_x_opt_angstroms, self.atom_el, self.atom_core):
			a = Atom(
							id=id,
							index=id+1,
							element=el,
							name=el,
							coordinates=coord,
						)
			a.element.guess_atomic_number()
			
			atoms.append(a)
			
			id += 1
		return atoms

	def to_molden(self):
		from ieda.QM.molden import Molden
		molden = Molden()
		# Atoms 
		id = 1
		for coord, el in zip(self.atom_x_opt_angstroms, self.atom_el):
			molden.atoms.append(
						Atom(
							id=id,
							index=id,
							element=el,
							name=el,
							coordinates=coord,
						))
			id += 1
		# GTOs
		at_seq = 1
		gto = GTO()
		gto.atom_number = self.ao_atomindex[0]
		for ao_index, at_symtyp, zeta in zip(self.ao_atomindex, self.atom_symtype, self.ao_zeta):
			if at_seq != ao_index:
				molden.GTOs.append(gto.copy())
				gto.delete()
				gto = GTO()
				gto.clear()
				gto.atom_number = ao_index
				at_seq += 1

			if at_seq == ao_index:
				if "S" == at_symtyp:
					gto.shell_label.append("s")
					gto.exponent_primitive.append([zeta*zeta/2])
					gto.contraction_coefficient.append([1])
					gto.number_of_primitives.append(1)
				elif "PX" == at_symtyp:
					gto.shell_label.append("p")
					gto.exponent_primitive.append([zeta*zeta/4])
					gto.contraction_coefficient.append([0.333])
					gto.number_of_primitives.append(3)
				elif "PY" == at_symtyp or "PZ" == at_symtyp:
					gto.exponent_primitive[1].append(zeta*zeta/4)
					gto.contraction_coefficient[1].append(0.333)
					
				elif at_symtyp == "X2": # D type = X2 XZ Z2 YZ XY
					gto.shell_label.append("d")
					gto.exponent_primitive.append([zeta])
					gto.contraction_coefficient.append([1])
					gto.number_of_primitives.append(1)
				elif at_symtyp in ["XZ", "Z2", "YZ", "XY"]:
					pass
			
		molden.GTOs.append(gto.copy())


		# MOs
		i = 1
		num_orb = len(self.atom_symtype)
		for energy, occ in zip(self.lmo_energy_levels, self.molecular_orbital_occupancies):
			mo = MO(id=i, energy=energy, occupation=occ)
			mo.coefficients = list(map(lambda c: c, self.lmo_vectors[(i-1)*num_orb:i*num_orb]))
			molden.MOs.append(mo)
			mo.spin = "Alpha"
			i += 1 
		return molden

			
