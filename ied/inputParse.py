# Imports
###############################################################################
import re
import sys
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
class InputParser: 
	def __init__(self):
		pass

	def comand_line(self) -> dict:
		command_line = command_line = str(sys.argv[1:]).replace(",", "").replace("[", "").replace("]", "").replace("'", "")
		return self.make_input(text_arg=command_line)


	def file_args(self, file:str) -> dict:
		with open(file,"r") as file:
			text_arg = file.read()
			text_args_list = [line.split('#')[0] if '#' in line else line for line in text_arg.split('\n')]
			text_args = '\n'.join(text_args_list)
			text_args = text_args.replace("\n", "").replace("\"", "")
			return self.make_input(text_arg=text_args)


	def make_input(self, text_arg:str) -> dict:
		result = {}
		pipes =  text_arg.split("--")[1:]
		for pipe in pipes:
			dic_a = {}
			dic_b = {}
			dic_c = {}
			dic_d = {}			
			key_a = pipe.split("-")[0]
			sub_a =  re.split(r'(?<=\s)-\s', pipe)[1:]
			# Intera os args do Pipe
			for arg_a in sub_a:
				# Separa os args pelo :
				arg_aa = arg_a.split(":")
				#  Verifica se tem "++", se tiver é um argumentos, caso contrario existem subs argumentos
				if not "+++" in arg_a:
					dic_a[arg_aa[0].strip()]=arg_aa[1].strip()
				else:
					dic_b = {}
					# Entra no novo sub nivel. ++
					key_b = arg_a.split("+++")[0] # Chave do dict B
					for arg_b in arg_a.split("+++")[1:]:
						# Separa os args pelo :
						arg_bb = arg_b.split(":")
						if not "++" in arg_b:
							dic_b[arg_bb[0].strip()]=arg_bb[1].strip()
						else:
							dic_c = {}
							key_c = arg_b.split("++")[0] # Chave dict C
							for arg_c in arg_b.split("++")[1:]:
								# Separa os args pelo :
								arg_cc = arg_c.split(":")
								if not "+" in arg_c:
									dic_c[arg_cc[0].strip()]=arg_cc[1].strip()
								else:
									dic_d = {}
									# Novo nivel 
									key_d = arg_c.split("+")[0]
									for arg_d in arg_c.split("+")[1:]:
										arg_dd = arg_d.split(":")
										if not ">>>>" in arg_d:
											dic_d[arg_dd[0].strip()]=arg_dd[1].strip()
										else:
											# Novo nivel 
											pass
									dic_c[key_d.strip()]=dic_d
									pass
							dic_b[key_c.strip()]=dic_c
					# Adiciona dic_b no dict_a 
					dic_a[key_b.strip()]=dic_b
			# Adiciona pipe ao dic com seus args 	
			result[key_a.strip()]=dic_a
		return result
			
		    

		    
		    
	
