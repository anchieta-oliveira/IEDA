# Imports
###############################################################################
import numpy as np
from ied.core import read_file
from scipy.special import gamma, gammainc
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
class OverlapMatrix:
  def __init__(self, title="", lower_half_triangle:list = [], liner:list = [], matrix:list = []):
    self.title:str = ""
    self.lower_half_triangle = lower_half_triangle
    self.liner = liner
    self.matrix = matrix
  

    def normalization_factor(alpha, l, m, n):
        """
        Calcula o fator de normalização para uma gaussiana primitiva.
        """
        prefactor = (2 * alpha / np.pi)**(3/4)
        normalization = (4 * alpha)**((l + m + n) / 2)
        l_factorial = np.math.factorial(l)
        m_factorial = np.math.factorial(m)
        n_factorial = np.math.factorial(n)
        return prefactor * (normalization / np.sqrt(l_factorial * m_factorial * n_factorial))


    def boys_function(n, x):
        """
        Calcula a função de Boys.
        """
        if x < 1e-8:
            return 1 / (2 * n + 1)
        else:
            return 0.5 * (x**(-n - 0.5)) * gamma(n + 0.5) * gammainc(n + 0.5, x)


    def gaussian_overlap(alpha1, A, l1, m1, n1, alpha2, B, l2, m2, n2):
        """
        Calcula a integral de sobreposição entre duas funções gaussianas primitivas com partes angulares.
        """
        p = alpha1 + alpha2
        P = (alpha1 * A + alpha2 * B) / p
        AB2 = np.dot(A - B, A - B)
        K = np.exp(-alpha1 * alpha2 / p * AB2)
        
        Sx = boys_function(l1 + l2, alpha1 * alpha2 * AB2 / p)
        Sy = boys_function(m1 + m2, alpha1 * alpha2 * AB2 / p)
        Sz = boys_function(n1 + n2, alpha1 * alpha2 * AB2 / p)
        
        prefactor = (np.pi / p)**(3/2)
        return prefactor * K * Sx * Sy * Sz


    def generate_angular_momentum(n_max):
        """
        Gera uma lista de tuplas (l, m, ml) para todos os orbitais até n_max.
        """
        angular_momentum = []
        for n in range(1, n_max + 1):
            for l in range(n):
                for m in range(-l, l + 1):
                    angular_momentum.append((l, m, m))
        return angular_momentum
    
    
    def make_overlap_matrix(centers:np.array, exponents:np.array, coefficients:np.array, angular_momentum:np.array):
        """
        Calcula a matriz de sobreposição S para um conjunto de orbitais atômicos com componentes angulares.
        """
        n = len(centers)
        S = np.zeros((n, n))
        centers = np.array(centers)
        angular_momentum = np.array(angular_momentum)

        for i in range(n):
            for j in range(i, n):  # Aproveitar a simetria da matriz de sobreposição
                S_ij = 0
                for k in range(len(exponents[i])):
                    for l in range(len(exponents[j])):
                        alpha_i = exponents[i][k]
                        alpha_j = exponents[j][l]
                        A = centers[i]
                        B = centers[j]
                        c_i = coefficients[i][k]
                        c_j = coefficients[j][l]
                        l1, m1, n1, l1, m1 = angular_momentum[i]
                        l2, m2, n2 = angular_momentum[j]
                        N_i = normalization_factor(alpha_i, l1, m1, n1)
                        N_j = normalization_factor(alpha_j, l2, m2, n2)
                        S_ij += c_i * c_j * N_i * N_j * gaussian_overlap(alpha_i, A, l1, m1, n1, alpha_j, B, l2, m2, n2)
                S[i, j] = S_ij
                S[j, i] = S_ij  # Aproveitar a simetria
        return S


    def molden_to_S(self) -> np.array:
       pass


  def read_from_multiwfn(self, path:str) -> list:
    """with open(path, "r") as file:
      lines = file.readlines()"""
        
    text = read_file(path)
    lines = text.splitlines(True)

    matrix = []
    matrix_lin = []
    num_col = 0
    for line in lines:
      line_split = line.split()
      if "*" in line:
         pass
      elif " \n" == line:
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

    n = len(matrix)
    # Inicialize uma matriz completa com zeros
    self.matrix  = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n - i):
            self.matrix[i][i + j] = self.lower_half_triangle[i][j]
            self.matrix[i + j][i] = self.lower_half_triangle[i][j]

    return self.matrix 
           
           
      

      


