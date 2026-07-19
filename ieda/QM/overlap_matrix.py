# Imports
###############################################################################
import numpy as np
from ieda.core import read_file
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
This project is licensed under the MIT License.
'''

# Classes
###############################################################################
class OverlapMatrix:
  def __init__(self, title="", lower_half_triangle:list = [], liner:list = [], matrix:list = []):
    self.title:str = ""
    self.lower_half_triangle = lower_half_triangle
    self.liner = liner
    self.matrix = matrix

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

# Gaussian Overlap Utilities
###############################################################################
def normalization_factor(alpha, l, m, n):
    N = (2 * alpha / np.pi)**(3/4)
    N *= (4 * alpha)**((l + m + n) / 2)
    N /= np.sqrt(np.math.factorial(l) * np.math.factorial(m) * np.math.factorial(n))
    return N


def boys_function(n, x):
    if x < 1e-8:
        return 1 / (2 * n + 1)
    return 0.5 * (x**(-n - 0.5)) * gamma(n + 0.5) * gammainc(n + 0.5, x)


def gaussian_overlap(alpha1, A, l1, m1, n1, alpha2, B, l2, m2, n2):
    p = alpha1 + alpha2
    AB2 = np.dot(A - B, A - B)
    K = np.exp(-alpha1 * alpha2 / p * AB2)
    Sx = boys_function(l1 + l2, alpha1 * alpha2 * AB2 / p)
    Sy = boys_function(m1 + m2, alpha1 * alpha2 * AB2 / p)
    Sz = boys_function(n1 + n2, alpha1 * alpha2 * AB2 / p)
    return (np.pi / p)**(3/2) * K * Sx * Sy * Sz


def generate_angular_momentum(n_max):
    result = []
    for n in range(1, n_max + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                result.append((l, m, m))
    return result


def make_overlap_matrix(centers, exponents, coefficients, angular_momentum):
    n = len(centers)
    S = np.zeros((n, n))
    centers = np.array(centers)
    angular_momentum = np.array(angular_momentum)
    for i in range(n):
        for j in range(i, n):
            S_ij = 0.0
            for k in range(len(exponents[i])):
                for l in range(len(exponents[j])):
                    alpha_i = exponents[i][k]
                    alpha_j = exponents[j][l]
                    A = centers[i]
                    B = centers[j]
                    c_i = coefficients[i][k]
                    c_j = coefficients[j][l]
                    li, mi, nli = angular_momentum[i]
                    lj, mj, nlj = angular_momentum[j]
                    Ni = normalization_factor(alpha_i, li, mi, nli)
                    Nj = normalization_factor(alpha_j, lj, mj, nlj)
                    S_ij += c_i * c_j * Ni * Nj * gaussian_overlap(
                        alpha_i, A, li, mi, nli,
                        alpha_j, B, lj, mj, nlj
                    )
            S[i, j] = S_ij
            S[j, i] = S_ij
    return S
    

