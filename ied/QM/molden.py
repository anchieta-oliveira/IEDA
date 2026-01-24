# Imports
###############################################################################
import numpy as np
from ied.QM.MO import MO
from ied.QM.aux import AUX
from ied.QM.GTO import GTO
from ied.MOL.atom import Atom
from ied.core import read_file
###############################################################################

"""
Molden (https://www.theochem.ru.nl/molden/molden_format.html)
[Molden Format]
[Title]
 Molden file ...

[Atoms] AU
element_name number atomic_number x y z
...
[GTO]
atom_sequence_number1 0
shell_label number_of_primitives 1.00
exponent_primitive_1 contraction_coefficient_1 (contraction_coefficient_1)
empty line
  1 0
s   3 1.0 
       18.7311370000         0.2149354183
        2.8253937000         0.3645711272
        0.6401217000         0.4150514277
s   1 1.0 
        0.1612778000         0.1813806839

  2 0
[MO]
 Sym= symmetry_label_1
 Ene= mo_energy_1
 Spin= (Alpha|Beta)
 Occup= mo_occupation_number_1
 ao_number_1 mo_coefficient_1
 ...
 ao_number_n mo_coefficient_n
...
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
This project is licensed under Creative Commons license (CC-BY-4.0) (Ver qual)
'''

# Classes
###############################################################################
class Molden: 
	# Gaussian Type Orbitals (GTO)
	def __init__(self):
		self.title:str = ""
		self.name_file:str = ""
		self.atoms:list = []
		self.GTOs:list = []
		self.MOs:list = []
		self.coordinates_unit: str = "Angs"

		
	def set_name(self, name:str) -> None:
		self.name = name
	
	def read_file(self, path:str) -> None:
		self.name = path.split("/")[-1].split(".")[-2]
		"""with open(path, "r") as file:
			molden_text = file.read()"""
		
		molden_text = read_file(path)

		for block in molden_text.split("["):
			# Ler atomos
			if "Atoms" in block.split("]")[0]:
				self.coordinates_unit = block.split("\n")[0].split()[1]
				for line in block.split("\n")[1:-1]:
					line_split = line.split()
					self.atoms.append(
										Atom(
											id=line_split[1],
											index=line_split[1],
											name=line_split[0],
											element=line_split[0],
											atomic_number=line_split[2],
											coordinates=(
														float(line_split[3]), 
														float(line_split[4]),
														float(line_split[5]))
											))
				
			# Ler GTO
			if "GTO" in block.split("]")[0]:
				lines = block[1:].split("\n")
				atom_id = int(lines[1].split()[0])
				gto = GTO(atom_number=atom_id)
				exp_list_tmp = []
				con_coef_list_tmp = []
				shell_list = []
				number_of_primitives_list = []

				for line in lines[2:-1]:
					line_split = line.split()
					if line_split == []:
											
						gto.shell_label = shell_list
						gto.number_of_primitives = number_of_primitives_list

						gto.exponent_primitive.append(exp_list_tmp)
						gto.contraction_coefficient.append(con_coef_list_tmp)
						self.GTOs.append(gto.copy())
						gto.clear()
						atom_id += 1 
						gto = GTO(atom_number=atom_id)

						shell_list = []
						number_of_primitives_list = []
						
					elif line_split[0] == str(atom_id):
						exp_list_tmp = []
						con_coef_list_tmp = []
						
					elif line_split[0] in ['s','p','d','f','sp','g']:

						shell_list.append(line_split[0])
						number_of_primitives_list.append(int(line_split[1]))
						if exp_list_tmp != [] and con_coef_list_tmp != []:
							gto.exponent_primitive.append(exp_list_tmp)
							gto.contraction_coefficient.append(con_coef_list_tmp)
							exp_list_tmp = []
							con_coef_list_tmp = []

					elif len(line_split) ==  2 and line_split[0] != str(atom_id):
						exp_list_tmp.append(float(line_split[0]))
						con_coef_list_tmp.append(float(line_split[1]))
	
			
			
			# Ler MOs
			if "MO" in block.split("]")[0]:
				for id, mo_text in enumerate(block.split("Sym=")[1:]):
					mo = MO(title=mo_text.split('\n')[0])
					mo.id = id +1 # Molden começa a contar de 1 
					for line in mo_text.split("\n")[1:-1]:
						line_split = line.split()
						if "Ene=" in line:
							mo.energy = float(line_split[1])
						elif "Spin=" in line:
							mo.spin = line_split[1]
						elif "Occup=" in line:
							mo.occupation = float(line_split[1])
						
						else:
							mo.ao_number.append(int(line_split[0]))
							mo.coefficients.append(float(line_split[1]))
					self.MOs.append(mo)
					del mo

	def __atomindex_atom_symtype(self) -> list:
		id_atom_mo = []
		atom_symtype = []
		for got in self.GTOs:
			for shells in got.shell_label:
				for shell in shells:
					if shell == "s":
						id_atom_mo.append(got.atom_number)
						atom_symtype.append("S")
					elif shell == "p": # PX PY PZ
						id_atom_mo.append(got.atom_number), atom_symtype.append("PX")
						id_atom_mo.append(got.atom_number), atom_symtype.append("PY")
						id_atom_mo.append(got.atom_number), atom_symtype.append("PZ")
					elif shell == "d":  # X2 XZ Z2 YZ XY; XX YY ZZ XY XZ YZ
						id_atom_mo.append(got.atom_number), atom_symtype.append("XX")
						id_atom_mo.append(got.atom_number), atom_symtype.append("YY")
						id_atom_mo.append(got.atom_number), atom_symtype.append("ZZ")
						id_atom_mo.append(got.atom_number), atom_symtype.append("XY")
						id_atom_mo.append(got.atom_number), atom_symtype.append("XZ")
						id_atom_mo.append(got.atom_number), atom_symtype.append("YZ")

					elif shell == "f": # Flata colocar 
						id_atom_mo.append(got.atom_number), id_atom_mo.append(got.atom_number) 
						id_atom_mo.append(got.atom_number), id_atom_mo.append(got.atom_number)
						id_atom_mo.append(got.atom_number), id_atom_mo.append(got.atom_number)
						id_atom_mo.append(got.atom_number)	
		return id_atom_mo, atom_symtype


	def __soma_gaussianas(self, coeficientes, alphas, x):
		resultado = 0
		for coef, alpha in zip(coeficientes, alphas):
			resultado += coef * np.exp(-alpha * (x**0 + 0**2 + 0**2))
		return resultado

	def __zeta(self) -> list:
		zetas = []

		for gto in self.GTOs:
			for coef, exp, shell in zip(gto.exponent_primitive, gto.contraction_coefficient, gto.shell_label):
				if shell == "s":
					zetas.append(self.__soma_gaussianas(coef, exp, 0))
				elif shell == "p":
					z = self.__soma_gaussianas(coef, exp, 0)
					zetas.append(z), zetas.append(z), zetas.append(z)
				elif shell == "d":
					z = self.__soma_gaussianas(coef, exp, 0)
					zetas.append(z), zetas.append(z), zetas.append(z)
					zetas.append(z), zetas.append(z), zetas.append(z) 	
		return zetas
	
	def get_ao_atomindex(self):
		ao_atomindex, atom_symtype = self.__atomindex_atom_symtype()
		return ao_atomindex
	
	def get_atom_symtype(self):
		ao_atomindex, atom_symtype = self.__atomindex_atom_symtype()
		return atom_symtype

	def to_aux(self, Smatrix:list=[]) -> AUX:
		aux = AUX()
		aux.atom_el = list(map(lambda at: at.element.symbol, self.atoms))
		aux.atom_x_angstroms = list(map(lambda at: (at.coordinates.x, at.coordinates.y, at.coordinates.z), self.atoms))
		ao_atomindex, atom_symtype = self.__atomindex_atom_symtype()
		aux.ao_atomindex = ao_atomindex
		aux.atom_symtype = atom_symtype
		aux.atom_x_opt_angstroms = list(map(lambda at: (at.coordinates.x, at.coordinates.y, at.coordinates.z), self.atoms))
		aux.lmo_vectors = [ele for li in list(map(lambda lmo: list(map(lambda c: c, lmo.coefficients)), self.MOs)) for ele in li]
		aux.lmo_energy_levels = list(map(lambda mo: mo.energy, self.MOs))
		aux.ao_zeta = self.__zeta()
		aux.molecular_orbital_occupancies = list(map(lambda mo: mo.occupation, self.MOs))
		aux.num_electrons = sum(aux.molecular_orbital_occupancies)
		aux.atom_pqn = list(map(lambda pqn: 1, aux.atom_symtype))
		if len(Smatrix) > 0:
			n = len(aux.ao_atomindex) # Tamanho da matriz
			aux.overlap_matrix = []
			for i in range(n):
				for j in range(i+1):
					aux.overlap_matrix.append(Smatrix[i][j])
		return aux

	def write(self, path:str) -> None:
		with open(path, "w") as file:
			file.write("[Molden Format]\n")
			file.write("[Title]\n")
			file.write(f"[Atoms] {self.coordinates_unit}\n")
			for atom in self.atoms:
				file.write(f"{atom.name}{atom.id:>7}{atom.element.atomic_number:>6} {atom.coordinates.x:21.14E} {atom.coordinates.y:21.14E} {atom.coordinates.z:21.14E}\n")
			file.write("[GTO]\n")

			for gto in self.GTOs:
				file.write(f"{gto.atom_number:>12} 0\n")
				for shell, n_pri, exp_list, coef_list in zip(gto.shell_label, gto.number_of_primitives, gto.exponent_primitive, gto.contraction_coefficient):
					file.write(f"{shell:>2}{n_pri:>12}   1.00000000000000\n")
					for exp, coef in zip(exp_list, coef_list):
						file.write(f"{exp:>19}   {coef:>25.15E}\n")
				file.write("\n")
			file.write("[MO]\n")
			for mo in self.MOs:
				file.write(f"Sym= {mo.id:>6}\n")
				file.write(f"Ene= {mo.energy:>23.16E}\n")
				file.write(f"Spin= {mo.spin}\n")
				file.write(f"Occup=     {mo.occupation:>1.8f}\n")
				for n, c in enumerate(mo.coefficients):
					file.write(f"{n+1:>12} {c:>23.16E}\n")
