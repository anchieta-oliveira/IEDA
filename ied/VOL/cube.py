# Description
###############################################################################
""" Module to store and manipulate data from cubes.

This module contains the Cube class, which is used to store and manipulate
data from cubes. The class has methods to read, write, and show cubes. The data
is stored in a structured array, which can be accessed using the get_text method.

Example
-------
    >>> cube = Cube()
    >>> cube.read("file.cube")
    >>> print(cube.get_text())
    >>> cube.show("vmd")
    >>> cube.write("file.cube")

{COMMENT1 (str)} 
{COMMENT2 (str)} 
{NATOMS (int)} {ORIGIN (3x float)} {NVAL (int)} 
{XAXIS (int) (3x float)} 
{YAXIS (int) (3x float) } 
{ZAXIS (int) (3x flutuante)} 
{GEOM (int) (float) (3x flutuante)} 
      . 
      . 
{DSET_IDS (#x int)} 
      . 
      . 
{DADOS (#x scinot)} 
      . 
      .

Ref.: https://paulbourke.net/dataformats/cube/

How to Convert Bohr Radius to Angstrom
1 b, a.u. = 0.529177249 A
1 A = 1.8897259886 b, a.u.

Example: convert 15 b, a.u. to A:
15 b, a.u. = 15 × 0.529177249 A = 7.937658735 A
https://www.unitconverters.net/length/bohr-radius-to-angstrom.htm

"""

# Imports
###############################################################################
import copy
import logging
import subprocess
import tempfile

import numpy as np

from ied.MOL.PDB import PDB
from ied.MOL.atom import Atom

# Set up logging
logging.basicConfig(level = logging.INFO, format = '%(asctime)s - %(levelname)s - %(message)s')

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
class Cube:
    """ Class to store and manipulate data from cubes.

    Attributes
    ----------
    id : int
        Cube ID.
    name : str
        Cube name.
    header : str
        Cube header.
    header2 : str
        Cube header2.
    atoms : list
        List of Atom objects.
    natoms : int
        Number of atoms.
    dset_ids_bool : bool
        DSet IDs boolean.
    dset_ids : np.array
        DSet IDs.
    voxels : list
        List of voxels.
    xorg : float
        X origin.
    yorg : float
        Y origin.
    zorg : float
        Z origin.
    xdel : float
        X delta.
    ydel : float
        Y delta.
    zdel : float
        Z delta.
    values : np.array
        Array of values.
    coordinates : np.array
        Array of coordinates.
    xn : int
        X number.
    yn : int
        Y number.
    zn : int
        Z number.
    """
    
    def __init__(self,
                 id: int = 0,
	    		 name: str = "",
				):
        ''' Constructor for Cube class.

        Parameters
        ----------
        id : int
            Cube ID.
        name : str
            Cube name.
        '''

        self.id: int = id
        self.name: str = name
        self.header: str = ""
        self.header2: str = ""
        self.atoms: list = []
        self.natoms: int = 0
        self.dset_ids_bool: bool = False
        self.dset_ids: np.array = np.array([])
        self.voxels: list = []
        self.xorg: float = .0
        self.yorg: float = .0
        self.zorg: float = .0
        self.xdel: float = .0
        self.ydel: float = .0
        self.zdel: float = .0
        self.values: np.array = np.array([], dtype = np.float64)
        self.coordinates: np.array = np.array([], dtype = np.float64)
        self.xn: int = 0,
        self.yn: int = 0,
        self.zn: int = 0
    
    def read(self, path: str) -> None:
        ''' Read the cube file.

        Parameters
        ----------
        path : str
            Path to the cube file.
        '''

        self.name = path.split("/")[-1].split(".")[-2]

        with open(path, "r") as file:
            cube_lines = file.readlines()

        self.header = cube_lines[0]
        self.header2 = cube_lines[1]

        # Atom information
        self.natoms =  abs(int(cube_lines[2].split()[0]))
        self.xorg = float(cube_lines[2].split()[1]); self.yorg = float(cube_lines[2].split()[2]);  self.zorg = float(cube_lines[2].split()[3])
        self.xn =  int(cube_lines[3].split()[0]); self.yn =  int(cube_lines[4].split()[0]); self.zn =  int(cube_lines[5].split()[0])
        self.xdel = float(cube_lines[3].split()[1]); self.ydel =  float(cube_lines[4].split()[2]); self.zdel =  float(cube_lines[5].split()[3])

        # Debug information
        logging.debug(f"{self.xorg} {self.yorg} {self.zorg}")
        logging.debug(f"{self.xdel} {self.ydel} {self.zdel}")
        logging.debug(f"{self.xn} {self.yn} {self.zn}")

        # Iterate over the atoms
        for i, latom in enumerate(cube_lines[6:self.natoms+6]):
            lsp = latom.split()
            self.atoms.append(
                                Atom(
                                    id=i,
                                    index=i+1,
                                    atomic_number=int(lsp[0]),
                                    coordinates=(float(lsp[2]), float(lsp[3]), float(lsp[4]))
                                )
                            )
            
        # Check if there are DSet IDs
        if  int(cube_lines[2].split()[0]) < 0: 
            self.dset_ids_bool = True
            n_dset_ids = int(cube_lines[self.natoms+6].split()[0])
            self.dset_ids =  np.array([0] * n_dset_ids, dtype=np.float64)
            lsp = cube_lines[self.natoms+6].split()
            for i, dset in enumerate(lsp[1:]):
                self.dset_ids[i] = float(dset)
            start_data = self.natoms + 7
        else:
            start_data = self.natoms + 6
             
        float_iter = (float(x) for l in cube_lines[start_data:] for x in l.split())
        self.values = np.fromiter(float_iter, dtype=np.float64)
        self.coordinates = self.get_coordinates()

    def write(self, path: str = "./") -> None:
        ''' Write the cube file.

        Parameters
        ----------
        path : str
            Path to save the cube file.
        '''

        with open(path, 'w') as cube_file:
            cube_file.write(self.get_text())

        cube_file.close()
    
    def get_text(self) -> str:
        ''' Get the cube file text.

        Returns
        -------
        str
            Cube file text.
        '''

        cube_file = ""
        cube_file += self.header
        cube_file += self.header2

        if len(self.header2) == 0:
            cube_file += f"Totally {self.xn*self.yn*self.zn} grid points\n"
        if self.dset_ids_bool:
            cube_file += f"{-self.natoms:>5}{self.xorg:>12.6f}{self.yorg:>12.6f}{self.zorg:>12.6f}\n"
        else:
            cube_file += f"{self.natoms:>5}{self.xorg:>12.6f}{self.yorg:>12.6f}{self.zorg:>12.6f}\n"

        cube_file += f"{self.xn:>5}{self.xdel:>12.6f}{0.000000:>12.6f}{0.000000:>12.6f}\n"
        cube_file += f"{self.yn:>5}{0.000000:>12.6f}{self.ydel:>12.6f}{0.000000:>12.6f}\n"
        cube_file += f"{self.zn:>5}{0.000000:>12.6f}{0.000000:>12.6f}{self.zdel:>12.6f}\n"

        # Write atoms
        for at in self.atoms:
            cube_file += f"{at.atomic_number:>5}{at.atomic_number:>12.6f}{at.coordinates.x:>12.6f}{at.coordinates.y:>12.6f}{at.coordinates.z:>12.6f}\n"

        # Write DSet IDs
        if self.dset_ids_bool:
            cube_file += f"{len(self.dset_ids):>5}"
            for dset in self.dset_ids:
                cube_file += f"{dset:>5}"
            cube_file +="\n"

        # Write Data
        i = 0
        for v in self.values:
            cube_file += f"{v:>14.5E}"
            i += 1
            if i ==6:
                cube_file += "\n"
                i = 0

        return cube_file
        
    def get_coordinates(self) -> np.array:
        ''' Get the coordinates.

        Returns
        -------
        np.array
            Array of coordinates.
        '''
        
        # Create one-dimensional arrays of X, Y, and Z coordinates
        x_coords = np.linspace(self.xorg, self.xorg + (self.xn - 1) * self.xdel, self.xn, dtype=np.float64)
        y_coords = np.linspace(self.yorg, self.yorg + (self.yn - 1) * self.ydel, self.yn, dtype=np.float64)
        z_coords = np.linspace(self.zorg, self.zorg + (self.zn - 1) * self.zdel, self.zn, dtype=np.float64)

        # Create all combinations of X, Y, and Z coordinates using meshgrid
        X, Y, Z = np.meshgrid(x_coords, y_coords, z_coords, indexing='ij')

        # Stack the X, Y, and Z coordinates into a single array
        self.coordinates = np.vstack((X.flatten(), Y.flatten(), Z.flatten())).T

        return self.coordinates

    def guess_origem(self) -> None:
        ''' Guess the origin. '''

        # Update the origin (xorg, yorg, zorg)
        self.xorg, self.yorg, self.zorg = np.min(self.coordinates, axis=0)

    def __guess_xn_yn_zn(self) -> None:
        ''' Guess the number of voxels. '''

        xmax, ymax, zmax = np.max(self.coordinates, axis=0) 

        # Update xn, yn, and zn
        self.xn = int((xmax - self.xorg) / self.xdel) + 1
        self.yn = int((ymax - self.yorg) / self.ydel) + 1
        self.zn = int((zmax - self.zorg) / self.zdel) + 1

    def trim(self, xmin: float = None, xmax: float = None, ymin: float = None, ymax: float = None, zmin: float = None, zmax: float = None) -> None:
        ''' Trim the cube.

        Parameters
        ----------
        xmin : float
            Minimum x coordinate.
        xmax : float
            Maximum x coordinate.
        ymin : float
            Minimum y coordinate.
        ymax : float
            Maximum y coordinate.
        zmin : float
            Minimum z coordinate.
        zmax : float
            Maximum z coordinate.
        '''

        # Filtra os valores dentro dos novos limites
        mask = (
            (self.coordinates[:, 0] >= xmin) & (self.coordinates[:, 0] <= xmax) &
            (self.coordinates[:, 1] >= ymin) & (self.coordinates[:, 1] <= ymax) &
            (self.coordinates[:, 2] >= zmin) & (self.coordinates[:, 2] <= zmax)
        )

        self.coordinates = self.coordinates[mask]
        self.values = self.values[mask]

        if len(self.values) >= 0:
            self.xorg = xmin
            self.yorg = ymin
            self.zorg = zmin
            self.guess_origem()
            self.__guess_xn_yn_zn()
        else:
            logging.error("Nenhum voxel foi incluido no corte! Verifique as coordenadas mínimas!")
    
    def sum(self, cube: object) -> None:
        ''' Sum the values of two cubes.

        Parameters
        ----------
        cube : object
            Cube object.
        '''

        self.values = self.values + cube.values 

    def sub(self, cube: object) -> None:
        ''' Subtract the values of two cubes.

        Parameters
        ----------
        cube : object
            Cube object.
        '''

        self.values = self.values - cube.values

    def select(self, selection: str, mol: PDB = PDB()) -> object:
        ''' Select a region of the cube.

        Parameters
        ----------
        selection : str
            Selection string.
        mol : PDB
            PDB object.

        Returns
        -------
        object
            Cube object.
        '''

        def gets_sels() -> list:
            ''' Get the selections.
            
            Returns
            -------
            list
                List of selections.
            '''
            
            keywords = selection.split()
            sels =[]
            i = 0

            while i < len(keywords)-1:
                args_sel = []

                if keywords[i] in functions.keys():
                    o = i

                    while not keywords[i+1] in functions.keys():
                        args_sel.append(keywords[i+1])
                        i += 1

                        if i ==  len(keywords)-1:
                            break

                    sels.append({keywords[o]:args_sel})

                i +=1

            return sels

        def value(args: list, sel_values: np.array) -> np.array:
            ''' Select values based on a condition.
            
            Parameters
            ----------
            args : list
                List of arguments.
            sel_values : np.array
                Array of values.

            Returns
            -------
            np.array
                Array of values.

            Raises
            ------
            ValueError
                If the operation is not valid.
            '''

            operation = args[0]
            isov = float(args[1])

            if ">" in operation:
                mask = sel_values > isov
                sel_values[~mask] = 0 
            elif "<" in operation:
                mask = sel_values < isov
                sel_values[~mask] = 0 
            if "=" in operation:
                mask = sel_values = isov
                sel_values[~mask] = 0 
            else:
                raise ValueError("Invalid operation. Choose '>', '<' or '='.")
            
            return sel_values

        def isovalue(args:list, sel_values:np.array) -> np.array:
            return value(args=[">=", args[0]], sel_values=sel_values)


        def within(args:list, sel_values:np.array) -> np.array:
            # within 5 of mol; within 5 of center x y z; 
            dist = float(args[0])
            ref = args[2]

            if ref == "mol":
                # Extração das coordenadas dos átomos
                atom_coords = np.array([atom.coordinates.get_array() for atom in mol.atoms])
                # Calcula a distância ao quadrado para evitar a raiz quadrada
                dist_squared = dist * dist
                # Calcula a matriz de distâncias ao quadrado entre todos os voxels e todos os átomos
                # Utiliza broadcasting para expandir as dimensões e calcular todas as distâncias de uma vez
                diff = atom_coords[:, np.newaxis, :] - self.coordinates[np.newaxis, :, :]
                distances_squared = np.sum(diff**2, axis=2, dtype=np.float64)
                # Verifica se cada voxel está dentro da distância especificada de qualquer átomo
                within_distance = np.any(distances_squared <= dist_squared, axis=0)
                # Muda para zero os que estão fora da condição 
                sel_values[~within_distance] = 0

            
            elif ref == "center":
                x = float(args[3]); y = float(args[4]); z = float(args[5])
                dist_squared = dist * dist  # Evita cálculo de raiz quadrada repetidamente
                center = np.array([[x, y, z]])
                # Calcula a matriz de distâncias ao quadrado entre todos os voxels e todos os átomos
                # Utiliza broadcasting para expandir as dimensões e calcular todas as distâncias de uma vez
                diff = center[:, np.newaxis, :] - self.coordinates[np.newaxis, :, :]
                distances_squared = np.sum(diff**2, axis=2, dtype=np.float64)
                # Verifica se cada voxel está dentro da distância especificada de qualquer átomo
                within_distance = np.any(distances_squared <= dist_squared, axis=0)
                # Muda para zero os que estão fora da condição 
                sel_values[~within_distance] = 0

            return sel_values

        # Sumario de seleções     
        functions = {	# Funções para seleção de resíduos
                    	"isovalue": isovalue,
                        "value": value,
                        "within": within,
                    }
        
        # Chama as funções 
        def call_function() -> list:
            sels = gets_sels()
            print(f"Sel. DX.: {sels}")
            if len(sels) != 0:
                sel_values = np.copy(self.values)
            else:
                sel_values = []
            for i, key in enumerate(sels):
                key = [*key][0]
                if "or" in sels[i-1][key]:
                    if i == 0:
                        n = self.xn*self.yn*self.zn
                        self.values = np.array([0] * n, dtype=np.float64)
                    sel_values = sel_values + functions[key](sels[i][key], self.values)
                elif "and" in sels[i][key]: 
                    sel_values = functions[key](sels[i][key], sel_values)
                else:
                    sel_values = functions[key](sels[i][key], sel_values)

            result = sel_values
            return result
                    
        sel_values = call_function()
        cube = self.copy()
        cube.values = sel_values    
        return cube

    def make(self, 
             values: np.array = np.array, 
             xdel: float =.0, 
             ydel: float =.0, 
             zdel: float = .0, 
             xorg: float=.0, 
             yorg: float = .0, 
             zorg: float = .0, 
             xn: float = .0, 
             yn:float =.0,
             zn:float =.0
            ) -> None:
        ''' Make a cube.

        Parameters
        ----------
        values : np.array
            Array of values.
        xdel : float
            X delta.
        ydel : float
            Y delta.
        zdel : float
            Z delta.
        xorg : float
            X origin.
        yorg : float
            Y origin.
        zorg : float
            Z origin.
        xn : float
            X number.
        yn : float
            Y number.
        zn : float
            Z number.
        '''

        self.xorg = xorg
        self.yorg = yorg
        self.zorg = zorg

        # xn, yn, zn
        if xn == .0 or yn == .0 or zn == .0:
            self.__guess_xn_yn_zn()
        else:
            self.xn = xn
            self.yn = yn
            self.zn = zn

        # Delta (xdel, ydel, zdel) 
        if xdel == .0 or ydel == .0 or zdel ==.0:
            pass
        else:
            self.xdel = xdel
            self.ydel = ydel
            self.zdel = zdel

        if values == np.array:
            n = self.xn * self.yn * self.zn
            self.values = np.zeros(n, dtype=np.float64)
        else:
            self.values = self.values

        self.get_coordinates()

    def fit(self, margin: float = 1e-6) -> None:
        ''' Fit the cube.

        Parameters
        ----------
        margin : float
            Margin.
        '''

        # Seleciona valores e coordenadas maiores que zero
        mask = self.values != 0
        nw_coords = self.coordinates[mask]

        # Atualiza a origem (xorg, yorg, zorg)
        xorg_sel, yorg_sel, zorg_sel = np.min(nw_coords, axis=0)
        
        xorg = xorg_sel - (self.xdel * margin)
        yorg = yorg_sel - (self.ydel * margin)
        zorg = zorg_sel - (self.zdel * margin)

        if xorg >= self.xorg:
            self.xorg = xorg_sel
        if yorg >= self.yorg:
            self.yorg = yorg_sel
        if zorg >= self.zorg:
            self.zorg = zorg_sel
        
        xmax_old, ymax_old, zmax_old = np.max(self.coordinates, axis=0)
        xmax_sel, ymax_sel, zmax_sel = np.max(nw_coords, axis=0)

        xmax = xmax_sel + (self.xdel * margin)
        ymax = ymax_sel + (self.ydel * margin)
        zmax = zmax_sel + (self.zdel * margin)
        if xmax_old <= xmax:
            xmax = xmax_sel
        if ymax_old <= ymax:
            ymax = ymax_sel
        if zmax_old <= zmax:
            zmax = zmax_sel
        
        # Atualiza xn, yn e zn
        xn_n = int((xmax - self.xorg) / self.xdel) + 1
        yn_n = int((ymax - self.yorg) / self.ydel) + 1
        zn_n = int((zmax - self.zorg) / self.zdel) + 1
        if xn_n <= self.xn:
            self.xn = xn_n
        if yn_n <= self.yn:
            self.yn = yn_n
        if zn_n <= self.zn:
            self.zn = zn_n

        # Filtra os valores dentro dos novos limites
        new_mask = (
            (self.coordinates[:, 0] >= self.xorg) & (self.coordinates[:, 0] <= xmax) &
            (self.coordinates[:, 1] >= self.yorg) & (self.coordinates[:, 1] <= ymax) &
            (self.coordinates[:, 2] >= self.zorg) & (self.coordinates[:, 2] <= zmax)
        )
        
        filtered_coords = self.coordinates[new_mask]
        filtered_values = self.values[new_mask]
        
        # Verifica se os valores filtrados estão dentro dos limites esperados
        #assert np.all(filtered_coords[:, 0] >= self.xorg) and np.all(filtered_coords[:, 0] <= xmax)
        #assert np.all(filtered_coords[:, 1] >= self.yorg) and np.all(filtered_coords[:, 1] <= ymax)
        #assert np.all(filtered_coords[:, 2] >= self.zorg) and np.all(filtered_coords[:, 2] <= zmax)

        self.values = filtered_values
        self.coordinates = filtered_coords

    def show(self, software:str = "vmd") -> None:
        ''' Show the cube.

        Parameters
        ----------
        software : str, optional
            Visualization software. Default is "vmd". Options are "vmd" and "pymol".
        
        Raises
        ------
        ValueError
            If the software is not available.
        '''

        # Check if the software is available
        if software not in ["vmd", "pymol"]:
            raise ValueError(f"The software {software} is not available. Please choose 'vmd' or 'pymol'.")
        
        text_pdb = self.get_text()

        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".cube") as temp_file:
            # Escrever conteúdo no arquivo temporário
            temp_file.write(text_pdb)
            temp_file.flush()

            # Obter o nome do arquivo temporário
            temp_file_name = temp_file.name 

            # Imprimir o nome do arquivo temporário (opcional, apenas para fins de visualização)
            if software == "vmd":
                cmd = f"vmd -cube {temp_file_name}"
            elif software == "pymol":
                cmd = f"pymol {temp_file_name}"
            else:
                cmd = ""


            if cmd:
                subprocess.run(cmd, shell=True)

    def copy(self) -> object:
        ''' Copy the cube.

        Returns
        -------
        object
            Cube object.
        '''

        return copy.deepcopy(self)

# Functions
###############################################################################
## Private ##

## Public ##

# Aliases
###############################################################################
