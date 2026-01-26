# Imports
###############################################################################

"""
atom_sequence_number1 0
shell_label number_of_primitives 1.00
exponent_primitive_1 contraction_coefficient_1 (contraction_coefficient_1)
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
class GTO: 
	#A ideia é fazer classes para AO, MO, GTO, Atoms (já existe) e conectar todas pelo index do atomo
	def __init__(self, title="", 
			  	atom_number:int=0, 
				shell_label:list=[], 
				number_of_primitives:list=[], 
				exponent_primitive:list=[], 
				contraction_coefficient:list=[]):
		
		self.title:str = ""
		self.atom_number:int = atom_number
		self.shell_label:list = shell_label
		self.number_of_primitives:list = number_of_primitives
		self.exponent_primitive:list = exponent_primitive 
		self.contraction_coefficient:list = contraction_coefficient


	def copy(self):
        # Cria uma cópia da instância atual
		return GTO(
            title=self.title,
            atom_number=self.atom_number,
            shell_label=self.shell_label.copy(),
            number_of_primitives=self.number_of_primitives.copy(),
            exponent_primitive=self.exponent_primitive.copy(),
            contraction_coefficient=self.contraction_coefficient.copy()
        )
	
	def clear(self):
        # Limpa os atributos da instância atual
		self.title = ""
		self.atom_number = 0
		self.shell_label.clear()
		self.number_of_primitives.clear()
		self.exponent_primitive.clear()
		self.contraction_coefficient.clear()

	def delete(self):
        # Deleta a instância atual
		del self

	