# Imports
###############################################################################
import numba
import logging
import zipfile
import numpy as np
from numba import cuda, float32, float64
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

# Functions
###############################################################################
@numba.njit(cache=True, parallel=True)
def intermolecular_eletron_density_two_selection(ao_ids_a:np.array, ao_ids_b:np.array, mo_coefficients:np.array, s_matrix:np.array):
    """
    Intermolecular Electron Density Analysis (IEDA) between two selections of atomic orbitals (AOs).
    This function computes the intermolecular electron density contributions using both Mulliken and OB methods.

    Parameters
    ----------
    ao_ids_a : np.array
        Array of indices for the first selection of atomic orbitals (AOs).
    ao_ids_b : np.array
        Array of indices for the second selection of atomic orbitals (AOs).
    mo_coefficients : np.array
        2D array where each row corresponds to the coefficients of a molecular orbital (MO).
    s_matrix : np.array
        Overlap matrix between atomic orbitals (AOs).
    Returns
    -------
    tuple
        A tuple containing two float values:
        - The first value is the Mulliken contribution to the intermolecular electron density.
        - The second value is the OB contribution to the intermolecular electron density.
    Notes
    -----
    - This function is optimized using Numba for parallel execution.
    - The function assumes that the input arrays are correctly formatted and compatible in dimensions.
    """
    result_comp_mulliken = .0

    for i in numba.prange(mo_coefficients.shape[0]):
        mo = mo_coefficients[i]
        tmp_comp_mulliken = .0
        for ci in ao_ids_a:
            for cj in ao_ids_b:
                tmp_comp_mulliken += (mo[ci] * mo[cj]) * s_matrix[ci, cj]
        result_comp_mulliken += 2*tmp_comp_mulliken 

    return result_comp_mulliken


@numba.njit(cache=True, parallel=True)
def matrix_intermolecular_eletron_density_numba(ats: np.array, ao_atomindex: np.array, mo_coefficients: np.array, s_matrix: np.array):
    """
    Compute the intermolecular electron density matrix using Mulliken and OB methods.
    This function calculates the intermolecular electron density contributions for each pair of atoms
    based on their associated atomic orbitals (AOs) and molecular orbital (MO) coefficients using Numba for optimization.

    Parameters
    ----------
    ats : np.array
        Array of atom indices.
    ao_atomindex : np.array    
        Array mapping each atomic orbital (AO) to its corresponding atom index.
    mo_coefficients : np.array
        2D array where each row corresponds to the coefficients of a molecular orbital (MO).
    s_matrix : np.array
        Overlap matrix between atomic orbitals (AOs).
    Returns
    -------
    tuple
        A tuple containing two 2D numpy arrays:
        - The first array is the Mulliken intermolecular electron density matrix.
        - The second array is the OB intermolecular electron density matrix.
    Notes
    -----
    - This function is optimized using Numba for parallel execution.
    - The function assumes that the input arrays are correctly formatted and compatible in dimensions.
    - The resulting matrices are symmetric, and only the upper triangular part is computed directly.
    """
    num_atoms = ats.size
    matrix_ied_m = np.zeros((num_atoms, num_atoms))

    for iat_a in numba.prange(num_atoms):
        at_a = ats[iat_a]
        #ao_ids_a = np.array([i for i, ao_id in enumerate(ao_atomindex) if ao_id == at_a])
        ao_ids_a = np.where(ao_atomindex == at_a)[0]
        

        for iat_b in range(iat_a, num_atoms):
            at_b = ats[iat_b]
            #ao_ids_b = np.array([i for i, ao_id in enumerate(ao_atomindex) if ao_id == at_b])
            ao_ids_b = np.where(ao_atomindex == at_b)[0]

            result_comp_mulliken = 0.0
            
            for i, mo in enumerate(mo_coefficients):
                tmp_comp_mulliken = 0.0

                for ci in ao_ids_a:
                    for cj in ao_ids_b:
                        """if ci == cj:
                            continue"""
                        tmp_comp_mulliken += mo[ci] * mo[cj] * s_matrix[ci, cj]
 

                result_comp_mulliken += 2*tmp_comp_mulliken
            
            # Symetric assignment
            matrix_ied_m[iat_a, iat_b] = result_comp_mulliken
            matrix_ied_m[iat_b, iat_a] = result_comp_mulliken


    return matrix_ied_m



@numba.njit(cache=True,  parallel=True, fastmath=True)
def two_selection_type_ao(ao_ids_a:np.array, ao_ids_b:np.array, mo_coefficients:np.array, s_matrix:np.array) -> float:
    """
    Calculate intermolecular electron density contributions for two selections of atomic orbitals (AOs).
    This function computes the product of molecular orbital coefficients and overlap matrix elements
    for all combinations of AOs from two selections, returning a filtered array of significant contributions.

    Parameters
    ----------
    ao_ids_a : np.array
        Array of indices for the first selection of atomic orbitals (AOs).
    ao_ids_b : np.array
        Array of indices for the second selection of atomic orbitals (AOs).
    mo_coefficients : np.array
        2D array where each row corresponds to the coefficients of a molecular orbital (MO).
    s_matrix : np.array
        Overlap matrix between atomic orbitals (AOs).
    Returns
    -------
    np.array
        A 2D numpy array where each row contains:
        - MO index
        - Product of MO coefficients and overlap matrix element
        - AO index from the first selection
        - AO index from the second selection
    Notes
    -----
    - This function is optimized using Numba for parallel execution.
    - Only contributions with a product greater than 1e-4 are included in the output.
    """

    resul = np.zeros((ao_ids_a.size * ao_ids_b.size * mo_coefficients.shape[0], 4))
    count = 0
    for i in numba.prange(mo_coefficients.shape[0]):
        mo = mo_coefficients[i]
        for ci in ao_ids_a:
            for cj in ao_ids_b:
                po = 2*(mo[ci] * mo[cj] * s_matrix[ci, cj])
                resul[count] = np.array([i, po, ci, cj])
                count += 1

    resul = resul[resul[:, 1] > 1e-4]
    return resul


# IEDA GPU  implementation
@cuda.jit(fastmath=True, cache=True)
def __kernel_ieda_two_sel(ao_ids_a, ao_ids_b, mo_coeffs, s_matrix, out_mulliken, out_bader):
    """
    GPU kernel para cálculo da densidade eletrônica intermolecular entre duas seleções AO.
    Cada bloco reduz sobre os MOs e produz um par (Mulliken, OB) como resultado.

    Parâmetros
    ----------
    ao_ids_a, ao_ids_b : 1D int32
        Índices das orbitais atômicas dos dois conjuntos.
    mo_coeffs : 2D float32
        Coeficientes moleculares (n_mos, n_ao).
    s_matrix : 2D float32
        Matriz de sobreposição AO.
    out_mulliken, out_bader : 1D float32 (length 1)
        Vetores de saída para resultados reduzidos (soma total).
    """
    n_mos = mo_coeffs.shape[0]
    n_ao_a = ao_ids_a.size
    n_ao_b = ao_ids_b.size

    tx = cuda.threadIdx.x
    bx = cuda.blockIdx.x
    bdim = cuda.blockDim.x
    gdim = cuda.gridDim.x

    partial_m = 0.0
    partial_b = 0.0

    for mo_idx in range(bx * bdim + tx, n_mos, bdim * gdim):
        tmp_m = 0.0
        tmp_b = 0.0
        for ai in range(n_ao_a):
            ci = ao_ids_a[ai]
            mo_ci = mo_coeffs[mo_idx, ci]
            mo_ci2 = mo_ci * mo_ci
            for bj in range(n_ao_b):
                cj = ao_ids_b[bj]
                mo_cj = mo_coeffs[mo_idx, cj]
                tmp_m += mo_ci * mo_cj * s_matrix[ci, cj]
                tmp_b += mo_ci2 * (mo_cj * mo_cj)
        partial_m += 2*tmp_m
        partial_b += tmp_b

    # Reduce within block
    sm_m = cuda.shared.array(256, dtype=float32)
    sm_b = cuda.shared.array(256, dtype=float32)
    sm_m[tx] = partial_m
    sm_b[tx] = partial_b
    cuda.syncthreads()

    stride = bdim // 2
    while stride > 0:
        if tx < stride:
            sm_m[tx] += sm_m[tx + stride]
            sm_b[tx] += sm_b[tx + stride]
        cuda.syncthreads()
        stride //= 2

    # Thread 0 writes block result to global memory, atomic add operation
    if tx == 0:
        cuda.atomic.add(out_mulliken, 0, sm_m[0])
        cuda.atomic.add(out_bader, 0, sm_b[0])


def intermolecular_eletron_density_two_selection_gpu(ao_ids_a:np.array, ao_ids_b:np.array, mo_coefficients:np.array, s_matrix:np.array, gpu_id) -> tuple:
    """
    Calculate intermolecular electron density contributions for two selections of atomic orbitals (AOs) using GPU acceleration.
    This function computes the intermolecular electron density contributions using both Mulliken and OB methods,
    leveraging CUDA for parallel computation.
    Parameters
    ----------
    ao_ids_a : np.array
        Array of indices for the first selection of atomic orbitals (AOs).
    ao_ids_b : np.array
        Array of indices for the second selection of atomic orbitals (AOs).
    mo_coefficients : np.array
        2D array where each row corresponds to the coefficients of a molecular orbital (MO).
    s_matrix : np.array
        Overlap matrix between atomic orbitals (AOs).
    gpu_id : int
        ID of the GPU to be used for computation.
    Returns
    -------
    tuple
        A tuple containing two float values:
        - The first value is the Mulliken contribution to the intermolecular electron density.
        - The second value is the OB contribution to the intermolecular electron density.
    Notes
    -----
    - This function utilizes Numba's CUDA JIT compiler for GPU acceleration.
    - The function assumes that the input arrays are correctly formatted and compatible in dimensions.
    """
    
    logging.info(f"Using GPU ID: {gpu_id} for IEDA calculations.")
    cuda.select_device(gpu_id)

    ao_ids_a = np.array(ao_ids_a, dtype=np.int32)
    ao_ids_b = np.array(ao_ids_b, dtype=np.int32)
    mo_coeffs = np.asarray(mo_coefficients, dtype=np.float32, order='C')
    s_matrix = np.asarray(s_matrix, dtype=np.float32, order='C')

    d_ao_a = cuda.to_device(ao_ids_a)
    d_ao_b = cuda.to_device(ao_ids_b)
    d_mo = cuda.to_device(mo_coeffs)
    d_s = cuda.to_device(s_matrix)
    d_out_m = cuda.to_device(np.zeros(1, dtype=np.float32))
    d_out_b = cuda.to_device(np.zeros(1, dtype=np.float32))


    threadsperblock = 256
    blockspergrid = min(1024, (mo_coeffs.shape[0] + threadsperblock - 1) // threadsperblock)

    __kernel_ieda_two_sel[blockspergrid, threadsperblock](
        d_ao_a, d_ao_b, d_mo, d_s, d_out_m, d_out_b
    )
    cuda.synchronize() 

    return float(d_out_m.copy_to_host()[0]), float(d_out_b.copy_to_host()[0])


def __build_ao_index_flat(ao_atomindex, ats):
    """
    Build a flattened array of atomic orbital (AO) indices for each atom, along with start indices and counts.
    This function creates a flattened representation of AO indices associated with each atom, which is useful for
    efficient access in GPU computations.

    Parameters
    ----------
    ao_atomindex : np.array
        Array mapping each atomic orbital (AO) to its corresponding atom index. 
    ats : np.array
        Array of atom indices for which the AO indices are to be extracted.
    Returns
    -------
    tuple
        A tuple containing three numpy arrays:
        - ao_index_flat: A flattened array of AO indices for all atoms in 'ats'.
        - ao_start: An array where each element indicates the starting index in 'ao_index_flat' for each atom.
        - ao_counts: An array where each element indicates the number of AOs associated with each atom.
    """
    n_atoms = ats.size

    lists = []
    for at in ats:
        aos = np.where(ao_atomindex == at)[0].astype(np.int32)
        lists.append(aos)

    ao_counts = np.array([len(l) for l in lists], dtype=np.int32)
    ao_start = np.zeros(n_atoms, dtype=np.int32)

    total = 0
    for i,c in enumerate(ao_counts):
        ao_start[i] = total
        total += c
    ao_index_flat = np.empty(total, dtype=np.int32)

    for i,l in enumerate(lists):
        st = ao_start[i]
        ao_index_flat[st:st+ao_counts[i]] = l

    return ao_index_flat, ao_start, ao_counts


@cuda.jit(fastmath=True, cache=True)
def __kernel_ieda(n_atoms, n_mos, ao_index_flat, ao_start, ao_count, mo_coeffs, s_matrix, matrix_m_out, matrix_b_out):
    """
    CUDA kernel to compute the intermolecular electron density matrix using Mulliken and OB methods.
    This kernel calculates the intermolecular electron density contributions for each pair of atoms
    based on their associated atomic orbitals (AOs) and molecular orbital (MO) coefficients.
    
    Parameters
    ----------
    n_atoms : int
        Number of atoms.
    n_mos : int
        Number of molecular orbitals (MOs).
    ao_index_flat : np.array
        Flattened array of AO indices for all atoms.
    ao_start : np.array
        Array where each element indicates the starting index in 'ao_index_flat' for each atom.
    ao_count : np.array
        Array where each element indicates the number of AOs associated with each atom.
    mo_coeffs : np.array
        2D array where each row corresponds to the coefficients of a molecular orbital (MO).
    s_matrix : np.array
        Overlap matrix between atomic orbitals (AOs).
    matrix_m_out : np.array
        Output array for the Mulliken intermolecular electron density matrix (flattened).
    matrix_b_out : np.array
        Output array for the Bader intermolecular electron density matrix (flattened).
    Notes
    -----
    - This kernel is designed to be executed on a GPU using Numba's CUDA JIT compiler.
    - Each thread computes contributions for a specific pair of atoms (i, j) in the upper triangular matrix.
    - Atomic operations are used to safely accumulate results in the output matrices.
    """

    i = cuda.blockIdx.x * cuda.blockDim.x + cuda.threadIdx.x  # linear thread id
    j = cuda.blockIdx.y * cuda.blockDim.y + cuda.threadIdx.y
    if i >= n_atoms or j >= n_atoms:
        return

    # compute only upper triangular (including diagonal)
    if j < i:
        return

    # indices in flattened AO arrays
    start_i = ao_start[i]
    count_i = ao_count[i]
    start_j = ao_start[j]
    count_j = ao_count[j]

    # local accumulators per thread
    res_m = 0.0
    res_b = 0.0

    # loop over MOs
    for im in range(n_mos):
        # tmp sums per MO
        tmp_m = 0.0
        tmp_b = 0.0

        # iterate AO pairs (ci in atom i; cj in atom j)
        for ai_idx in range(count_i):
            ci = ao_index_flat[start_i + ai_idx]  # AO index global
            mo_ci = mo_coeffs[im, ci]
            mo_ci2 = mo_ci * mo_ci
            for aj_idx in range(count_j):
                cj = ao_index_flat[start_j + aj_idx]
                mo_cj = mo_coeffs[im, cj]
                # Mulliken contribution:
                tmp_m += mo_ci * mo_cj * s_matrix[ci, cj]
                # OB contribution:
                tmp_b += mo_ci2 * (mo_cj * mo_cj)

        res_m += 2*tmp_m
        res_b += tmp_b 

    # flatten index
    idx_upper = i * n_atoms + j
    
    cuda.atomic.add(matrix_m_out, idx_upper, res_m)
    cuda.atomic.add(matrix_b_out, idx_upper, res_b)
    # symmetric partner when i != j
    if i != j:
        idx_sym = j * n_atoms + i
        cuda.atomic.add(matrix_m_out, idx_sym, res_m)
        cuda.atomic.add(matrix_b_out, idx_sym, res_b)


def matrix_intermolecular_eletron_density_numba_gpu(ats: np.array, ao_atomindex: np.array, mo_coefficients: np.array, s_matrix: np.array, dtype=np.float32, gpu_id: int = 0):
    """
    Compute the intermolecular electron density matrix using Mulliken and OB methods on a GPU.
    This function calculates the intermolecular electron density contributions for each pair of atoms
    based on their associated atomic orbitals (AOs) and molecular orbital (MO) coefficients using CUDA for acceleration.

    Parameters
    ----------
    ats : np.array 
        Array of atom indices.
    ao_atomindex : np.array
        Array mapping each atomic orbital (AO) to its corresponding atom index.
    mo_coefficients : np.array
        2D array where each row corresponds to the coefficients of a molecular orbital (MO).
    s_matrix : np.array
        Overlap matrix between atomic orbitals (AOs).
    dtype : data-type, optional
        Data type for computations (default is np.float32). Can be set to np.float64 for higher precision.
    Returns
    -------
    tuple
        A tuple containing two 2D numpy arrays:
        - The first array is the Mulliken intermolecular electron density matrix.
        - The second array is the OB intermolecular electron density matrix.
    Notes
    -----
    - This function utilizes Numba's CUDA JIT compiler for GPU acceleration.
    - The function assumes that the input arrays are correctly formatted and compatible in dimensions.
    """
    logging.info(f"Using GPU ID: {gpu_id} for IEDA calculations.")
    cuda.select_device(gpu_id)

    ao_index_flat, ao_start, ao_count = __build_ao_index_flat(ao_atomindex, ats)

    # cast to float32 for GPU
    mo = mo_coefficients.astype(dtype)
    s = s_matrix.astype(dtype)
    
    n_atoms = ats.size
    n_mos = mo.shape[0]

    # device arrays
    d_ao_flat = cuda.to_device(ao_index_flat)
    d_ao_start = cuda.to_device(ao_start)
    d_ao_count = cuda.to_device(ao_count)
    d_mo = cuda.to_device(mo)
    d_s = cuda.to_device(s)
    
    # result flattened arrays
    res_size = n_atoms * n_atoms
    d_mat_m = cuda.device_array(res_size, dtype=dtype)
    d_mat_b = cuda.device_array(res_size, dtype=dtype)

    # Grid/block sizing (tune)
    threadsperblock = (8, 8)
    blockspergrid_x = (n_atoms + threadsperblock[0] - 1) // threadsperblock[0]
    blockspergrid_y = (n_atoms + threadsperblock[1] - 1) // threadsperblock[1]
    blockspergrid = (blockspergrid_x, blockspergrid_y)

    __kernel_ieda[blockspergrid, threadsperblock](
        n_atoms, n_mos, d_ao_flat, d_ao_start, d_ao_count,
        d_mo, d_s, d_mat_m, d_mat_b
    )

    mat_m = d_mat_m.copy_to_host().reshape((n_atoms, n_atoms))
    mat_b = d_mat_b.copy_to_host().reshape((n_atoms, n_atoms))
    cuda.synchronize() 

    del d_ao_flat, d_ao_start, d_ao_count, d_mo, d_s, d_mat_m, d_mat_b
    return mat_m, mat_b


# Read file (txt or zip)
def read_file(path) -> str:
    """
    Reads the content of a file, supporting both plain text files and ZIP archives containing a single text file.
    If the file is a ZIP archive, it extracts and reads the first text file found within the archive.
    
    Parameters
    ----------
    path : str
        The path to the file to be read. Can be a plain text file or a ZIP archive.
    Returns
    -------
    str
        The content of the file as a string.
    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    ValueError
        If the ZIP archive is empty or if the path does not point to a valid file. 
    """
    if path.endswith(".zip"):
        with zipfile.ZipFile(path, "r") as zf:
            names = zf.namelist()
            if not names:
                raise ValueError(f"Null '{path}'.")
            
            inner_name = names[0]

            with zf.open(inner_name, "r") as file:
                text = file.read().decode("utf-8")
                file.close()
                return text
    else:
        with open(path, "r", encoding="utf-8") as file:
                text = file.read()
                file.close()
                return text
            
    raise ValueError("Choose lines or text to read the file.")