# Description
###############################################################################
""" Module to store and manipulate data from DX files.

This module contains the DX class, which is used to store and manipulate
data from DX files. The class has methods to read, write, and visualize DX files.
The data is stored in a structured array, which can be accessed using the 
get_text method.

Example
-------
    >>> dx = DX()
    >>> dx.read("file.dx")
    >>> print(dx.get_text())
    >>> dx.show("vmd")
    >>> dx.write("file.dx")

"""

# Imports
###############################################################################
import copy
import logging
import subprocess
import tempfile

import numpy as np

from multiprocessing import Pool
from scipy import interpolate
from ieda.MOL.PDB import PDB

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
This project is licensed under the MIT License.
'''


# Classes
###############################################################################
class DX:
    """ Class to store and manipulate data from DX files.

    Attributes
    ----------
    id : int
        DX ID.
    name : str
        DX name.
    voxels : list
        List of voxels.
    xorg : float
        Origin in the x-axis.
    yorg : float
        Origin in the y-axis.
    zorg : float
        Origin in the z-axis.
    xdel : float
        Delta in the x-axis.
    ydel : float
        Delta in the y-axis.
    zdel : float
        Delta in the z-axis.
    values : np.array
        Array of values.
    coordinates : np.array
        Array of coordinates.
    xn : int
        Number of voxels in the x-axis.
    yn : int
        Number of voxels in the y-axis.
    zn : int
        Number of voxels in the z-axis.
    gridpositions : list
        List of grid positions.
    gridconnections : list
        List of grid connections.
    """

    def __init__(self,
                 id: int = 0,
	    		 name: str = "",
				):
        ''' Constructor for the DX class.

        Parameters
        ----------
        id : int
            DX ID.
        name : str
            DX name.
        '''
        
        self.id: int = id
        self.name: str = name
        self.voxels: list = []
        self.xorg: float = .0
        self.yorg: float = .0
        self.zorg: float = .0
        self.xdel: float = .0
        self.ydel: float = .0
        self.zdel: float = .0
        self.values: np.array = np.array([], dtype = np.float64)
        self.coordinates: np.array = np.array([], dtype = np.float64)
        self.xn: int = 0
        self.yn: int = 0
        self.zn: int = 0
        self.gridpositions: list = []
        self.gridconnections: list = []
    
    def read(self, path: str):
        ''' Read a DX file.

        Parameters
        ----------
        path : str
            Path to the DX file.
        '''

        self.name = path.split("/")[-1].split(".")[-2]

        with open(path, "r") as file:
            dx_lines = file.readlines()
            
        i_delta = 0

        for line in dx_lines: # TODO: Checar se os ifs são ifs ou elifs (otimização)
            if "gridpositions" in line:
                self.gridpositions = [int(line.split()[-3]), int(line.split()[-2]), int(line.split()[-1])]
                self.xn = int(line.split()[-3]); self.yn =  int(line.split()[-2]); self.zn = int(line.split()[-1])
                n = self.xn*self.yn*self.zn
                self.values = np.array([0] * n, dtype=np.float64)

            if "gridconnections" in line:
                self.gridconnections = [int(line.split()[-3]), int(line.split()[-2]), int(line.split()[-1])]

            if "origin" in line:
                self.xorg = float(line.split()[1])
                self.yorg = float(line.split()[2])
                self.zorg = float(line.split()[3])

            if "delta" in line:
                line = line.split()
                if float(line[1]) != 0 and i_delta == 0:
                    self.xdel = float(line[1])
                    i_delta += 1

                if float(line[2]) != 0 and i_delta == 1:
                    self.ydel = float(line[2])
                    i_delta += 1 

                if float(line[3]) != 0 and i_delta == 2:
                    self.zdel = float(line[3])
                    i_delta += 1
            
            if line[0].isdigit() :
                break
        
        data = False
        id = 0
        for line in dx_lines:
            if "class field" in line:
                data = False

            if data:
                line = line.split()
                if len(line) >= 1:
                    self.values[id] = float(line[0])
                    
                    id += 1
                if len(line) >= 2:
                    self.values[id] = float(line[1])
                    id += 1
                if len(line) == 3:
                    self.values[id] = float(line[2])
                    id += 1
        
            if "data follows" in line:
                data = True

        self.get_coordinates()
        return self.voxels
    
    def get_coordinates(self) -> np.array:
        ''' Get the coordinates of the voxels.
        
        Returns
        -------
        np.array
            Array of coordinates.
        '''

        # Criar arrays unidimensionais de coordenadas X, Y e Z
        x_coords = np.linspace(self.xorg, self.xorg + (self.xn - 1) * self.xdel, self.xn, dtype=np.float64)
        y_coords = np.linspace(self.yorg, self.yorg + (self.yn - 1) * self.ydel, self.yn, dtype=np.float64)
        z_coords = np.linspace(self.zorg, self.zorg + (self.zn - 1) * self.zdel, self.zn, dtype=np.float64)

        # Criar todas as combinações de coordenadas X, Y e Z usando meshgrid
        X, Y, Z = np.meshgrid(x_coords, y_coords, z_coords, indexing='ij')

        # Empilhar as coordenadas X, Y e Z em um único array
        self.coordinates = np.vstack((X.flatten(), Y.flatten(), Z.flatten())).T

        return self.coordinates
                
    def update_info(self) -> None:
        ''' Recalculates the origin and size of the current Volume based on it's actual volumetric data. Called after modifications such as trimming or moving. '''
        # Update origin

        # Update voxel positions

        # Update size
        pass

    def write(self, filename: str) -> bool:
        ''' Write a DX file.

        Parameters
        ----------
        filename : str
            Name of the DX file.

        Returns
        -------
        bool
            True if the file was written successfully, False otherwise.
        '''

        try:
            with open(filename, 'w') as dx_file:
                # Escreve o cabeçalho do arquivo DX
                dx_file.write("#DX potential maps\n")
                dx_file.write("#Format of the file is:\n")
                dx_file.write(f"object 1 class gridpositions counts {self.xn} {self.yn} {self.zn}\n")
                dx_file.write(f"origin {self.xorg} {self.yorg} {self.zorg}\n")
                dx_file.write(f"delta {self.xdel} 0 0\n")
                dx_file.write(f"delta 0 {self.ydel} 0\n")
                dx_file.write(f"delta 0 0 {self.zdel}\n")
                dx_file.write(f"object 2 class gridconnections counts {self.xn} {self.yn} {self.zn}\n")
                dx_file.write(f"object 3 class array type double rank 0 items [{self.xn * self.yn * self.zn}] data follows\n")
                # Escreve os valores da grade no formato DX
                i = 0
                for v in self.values:
                    dx_file.write(f"{v}")
                    i += 1
                    if i == 3:
                        dx_file.write("\n")
                        i = 0
                    else:
                        dx_file.write(" ")
                if i !=3:
                    dx_file.write("\n")
                dx_file.write(f"object \"{self.name}\" class field")
            return True
        except Exception as e:
            logging.error(f"Error writing data to file: {e}")
            return False    

    def get_text(self) -> str:
        ''' Return the text of the DX file.

        Returns
        -------
        str
            Text of the DX file.
        '''

        dx_file = ""
        dx_file += "#DX potential maps\n"
        dx_file += "#Format of the file is:\n"
        dx_file += f"object 1 class gridpositions counts {self.xn} {self.yn} {self.zn}\n"
        dx_file += f"origin {self.xorg} {self.yorg} {self.zorg}\n"
        dx_file += f"delta {self.xdel} 0 0\n"
        dx_file += f"delta 0 {self.ydel} 0\n"
        dx_file += f"delta 0 0 {self.zdel}\n"
        dx_file += f"object 2 class gridconnections counts {self.xn} {self.yn} {self.zn}\n"
        dx_file += f"object 3 class array type double rank 0 items [{self.xn * self.yn * self.zn}] data follows\n"

        # Write the grid values in DX format
        i = 0

        for v in self.values:
            dx_file += f"{v}"
            i += 1
            if i == 3:
                dx_file += "\n"
                i = 0
            else:
                dx_file += " "

        if i !=3:
            dx_file += "\n"

        dx_file += f"object \"{self.name}\" class field"

        return dx_file

    def show(self, software: str = "vmd") -> None:
        ''' Show the DX file using a visualization software.

        Parameters
        ----------
        software : str, optional
            Visualization software to use. Default is "vmd". Options are "vmd" and "pymol".

        Raises
        ------
        ValueError
            If the software is not available.
        '''

        # Check if the software is available
        if software not in ["vmd", "pymol"]:
            raise ValueError(f"Software '{software}' not available. Use 'vmd' or 'pymol'.")
        
        text = self.get_text()

        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".dx") as temp_file:
            # Write the text to a temporary file
            temp_file.write(text)
            temp_file.flush()

            # Get the name of the temporary file
            temp_file_name = temp_file.name 

            # Print the name of the temporary file (optional, only for visualization purposes)
            if software == "vmd":
                cmd = f"vmd -dx {temp_file_name}"
            elif software == "pymol":
                cmd = f"pymol {temp_file_name}"
            else:
                cmd = ""
            
            if cmd:
                logging.debug(f"Running command: {cmd}")
                subprocess.run(cmd, shell = True)
    
    def __guess_origin(self):
        ''' Guess the origin of the DX file. '''

        # Atualiza a origem (xorg, yorg, zorg)
        self.xorg, self.yorg, self.zorg = np.min(self.coordinates, axis=0)

    def __guess_xn_yn_zn(self):
        ''' Guess the number of voxels in the DX file. '''

        xmax, ymax, zmax = np.max(self.coordinates, axis=0) 

        # Atualiza xn, yn e zn
        self.xn = int((xmax - self.xorg) / self.xdel) + 1
        self.yn = int((ymax - self.yorg) / self.ydel) + 1
        self.zn = int((zmax - self.zorg) / self.zdel) + 1

    def trim(self,
             xmin: float = None,
             xmax: float = None,
             ymin: float = None,
             ymax: float = None,
             zmin: float = None,
             zmax: float = None
            ) -> None:
        ''' Trim the DX file.
        
        Parameters
        ----------
        xmin : float, optional
            Minimum value in the x-axis.
        xmax : float, optional
            Maximum value in the x-axis.
        ymin : float, optional
            Minimum value in the y-axis.
        ymax : float, optional
            Maximum value in the y-axis.
        zmin : float, optional
            Minimum value in the z-axis.
        zmax : float, optional
            Maximum value in the z-axis.
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
            self.__guess_origin()
            self.__guess_xn_yn_zn()
        else:
            logging.error("Nenhum voxel foi incluido no corte! Verifique as coordenadas mínimas!")
    
    def sum(self, dx: object) -> None:
        ''' Sum two DX files. '''

        self.values = self.values + dx.values 

    def sub(self, dx: object) -> None:
        ''' Subtract two DX files. '''

        self.values = self.values - dx.values
        #self.select(selection="value > 0")

    def corrcoef(self, dx: object) -> float:
        ''' Calculate the correlation coefficient between two DX files.

        Parameters
        ----------
        dx : object
            DX file to compare.

        Returns
        -------
        float
            Correlation coefficient between the two DX files.
        '''

        return np.corrcoef(self.values, dx.values)[0, 1]

    def interpolate_grid(self, source: object, target: object, method: str = 'linear') -> np.array:
        ''' Interpolate the values of the source grid to the target grid.
        
        Parameters
        ----------
        source : object
            Source DX file.
        target : object
            Target DX file.
        method : str, optional
            Interpolation method. Default is 'linear'. Options are 'nearest', 'linear', 'cubic', etc.

        Returns
        -------
        np.array
            Interpolated values in the target grid.
        '''

        # Defina as coordenadas do grid de origem e destino
        target_coords = target.coordinates
        src_coords = source.coordinates

        # Crie o grid de valores de origem
        source_values = np.reshape(source.values, (source.xn, source.yn, source.zn))

        # Interpole os valores para o grid de destino
        interpolated_values = interpolate.interpn(
            (np.unique(src_coords[:, 0]), np.unique(src_coords[:, 1]), np.unique(src_coords[:, 2])),
            source_values,
            target_coords,
            method = method,
            bounds_error = False,
            fill_value = 0
        )

        return interpolated_values

    def add(self, dx: object) -> object:
        ''' Add two DX grids, adjusting the grid sizes as necessary. 
        
        Parameters
        ----------
        dx : object
            Grid to add.

        Returns
        -------
        object
            A new DX object with the sum of the two grids.
        '''

        # Identificar a grade com maior extensão
        if np.prod(self.gridpositions) >= np.prod(dx.gridpositions):
            target_dx, src_dx = self, dx
        else:
            target_dx, src_dx = dx, self

        # Interpolar o grid menor para o tamanho da grade maior
        interpolated_values = self.interpolate_grid(src_dx, target_dx)

        # Somar os valores
        summed_values = target_dx.values + interpolated_values

        # Criar um novo DX com os valores somados
        summed_dx = target_dx.copy()
        summed_dx.values = summed_values
        return summed_dx

    def select(self, selection: str, mol: PDB = PDB()) -> object:
        ''' Select voxels based on a selection string.

        Parameters
        ----------
        selection : str
            Selection string.
        mol : PDB, optional
            PDB object. Default is PDB().

        Returns
        -------
        object
            A new DX object with the selected voxels.
        '''

        def gets_sels() -> list:
            ''' Get the selections from the selection string. 
            
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


        def value(args:list, sel_values:np.array) -> np.array:
            operation = args[0]
            isov = float(args[1])
            if ">" in operation:
                mask = sel_values > isov
                sel_values[~mask] = 0 
            elif "<" in operation:
                mask = sel_values < isov
                sel_values[~mask] = 0 
            elif "=" in operation:
                mask = sel_values == isov
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
                atom_coords = np.array([atom.coordinates.get_array() for atom in mol.atoms], dtype=np.float64)
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
                center = np.array([[x, y, z]], dtype=np.float64)
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
            #print(f"Sel. DX.: {sels}")
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
        dx = self.copy()
        dx.values = sel_values    
        return dx

    def fit(self, margin: float = 1e-6) -> None:
        ''' Fit the grid to the data.

        Parameters
        ----------
        margin : float, optional
            Margin to add to the grid. Default is 1e-6.
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

    def make(self, 
             values: np.array = np.array, 
             xdel: float =.0, 
             ydel: float =.0, 
             zdel: float = .0, 
             xorg: float = .0, yorg: float = .0,
             zorg: float = .0,
             xn: float = .0,
             yn: float = .0, 
             zn: float = .0
            ) -> None:
        ''' Create a DX object.

        Parameters
        ----------
        values : np.array, optional
            Array of values. Default is np.array.
        xdel : float, optional
            Delta in the x-axis. Default is .0.
        ydel : float, optional
            Delta in the y-axis. Default is .0.
        zdel : float, optional
            Delta in the z-axis. Default is .0.
        xorg : float, optional
            Origin in the x-axis. Default is .0.
        yorg : float, optional
            Origin in the y-axis. Default is .0.
        zorg : float, optional
            Origin in the z-axis. Default is .0.
        xn : float, optional
            Number of voxels in the x-axis. Default is .0.
        yn : float, optional
            Number of voxels in the y-axis. Default is .0.
        zn : float, optional
            Number of voxels in the z-axis. Default is .0.
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
        
    def copy(self) -> object:
        ''' Create a deep copy of the instance.

        Returns
        -------
        object
            Deep copy of the instance.
        '''
        
        # Realiza uma cópia profunda de todos os atributos da instância
        new_copy = DX()
        new_copy.id = self.id
        new_copy.name = self.name
        new_copy.voxels = copy.deepcopy(self.voxels)
        new_copy.xorg = self.xorg
        new_copy.yorg = self.yorg
        new_copy.zorg = self.zorg
        new_copy.xdel = self.xdel
        new_copy.ydel = self.ydel
        new_copy.zdel = self.zdel
        new_copy.values = np.copy(self.values)
        new_copy.coordinates = np.copy(self.coordinates)
        new_copy.xn = self.xn
        new_copy.yn = self.yn
        new_copy.zn = self.zn
        new_copy.gridpositions = copy.deepcopy(self.gridpositions)
        new_copy.gridconnections = copy.deepcopy(self.gridconnections)
        return new_copy

# Functions
###############################################################################
## Private ##

## Public ##

# Aliases
###############################################################################
