# Imports
###############################################################################
import numpy as np
##############################################################################

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
class DensityMatrix:
  def __init__(self, title="", lower_half_triangle:list = [], liner:list = [], matrix:list = []):
    self.title:str = ""
    self.lower_half_triangle = lower_half_triangle
    self.liner = liner
    self.matrix = matrix


  def make_density(self, MOs:list, nelec:int):
    nocc = nelec // 2
    c_mos = np.array([mo.coefficients for mo in MOs[:nocc]])
    nbas, _ = c_mos.shape
    dens = np.zeros((nbas, nbas))
    for mu in range(nbas):
        for nu in range(nbas):
            for occ in range(nocc):
                dens[mu, nu] += c_mos[occ, mu] * c_mos[occ, nu]
        dens[mu, nu] *= 2
    return dens


  def read_from_multiwfn(self, path:str) -> list:
    with open(path, "r") as file:
      lines = file.readlines()
    matrix = []
    matrix_lin = []
    num_col = 0
    for line in lines:
      line_split = line.split()
      if "*" in line:
         pass
      elif "\n" == line:
        break
      elif line[:6] == "      ":
        num_col = int(len(line_split))
        init_col = int(line_split[0])
        for n in range(num_col):
          matrix.append([])
        
      else:
        for id_col, v_col in enumerate(line_split[1:]):
          matrix[init_col+id_col-1].append(float(v_col))
          matrix_lin.append(float(v_col))
    self.liner = matrix_lin
    self.lower_half_triangle = matrix    
    return matrix
           
           
      

      


