# Imports
###############################################################################
import json
import logging
import numpy as np
import pandas as pd
import seaborn as sns
import ieda.core as core
from ieda.MOL.PDB import PDB
import matplotlib.pyplot as plt
from collections import defaultdict
from ieda.MOL.selection import Selection
from ieda.QM.overlap_matrix import OverlapMatrix
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
class IEDA:
    def __init__(self) -> None:
        self.df = pd.DataFrame()
        self.path_ied = ""
        

    def __get_index_homo_lumo(self, MOs):
        for mo in MOs:
            if mo.occupation == 0:
                logging.info(f"HOMO: {mo.id-1} LUMO: {mo.id}")
                id_list_homo = mo.id-2 # dois pq é a posição da lista 
                id_list_lumo = mo.id-1 # para ajustar a posição da lista
                break
        return id_list_homo, id_list_lumo


    def __arrays_to_nested_dict(self, matrix, ats):
        result = defaultdict(dict)
        for i, at_a in enumerate(ats):
            key_a = int(at_a)
            for j, at_b in enumerate(ats):
                key_b = int(at_b)
                result[key_a][key_b] = float(matrix[i, j])
        return result
    

    def two_sel(self, sel_a:str, sel_b:str,
                pdb:PDB=PDB, aux=None, 
                molden=None, orca_out=None, 
                Smatrix:OverlapMatrix=OverlapMatrix,
                gpu:bool=False, 
                gpu_id:int = 0
                ) -> tuple:
        
        if aux != None:
            MOs = aux.get_MOs_objs()
            id_list_homo, id_list_lumo = self.__get_index_homo_lumo(MOs=MOs)
            mo_coefficients = np.array([mo.coefficients for mo in MOs[:id_list_homo+1]], dtype=np.float64)
            ao_atomindex = aux.ao_atomindex 
            s_matrix = aux.get_S_objs().matrix
            del aux
            
        elif molden != None:
            MOs = molden.MOs
            id_list_homo, id_list_lumo =  self.__get_index_homo_lumo(MOs=MOs)
            mo_coefficients = np.array([mo.coefficients for mo in MOs[:id_list_homo+1]], dtype=np.float64)
            ao_atomindex = molden.get_ao_atomindex() 
            
            logging.info(f"Nº AO index: {len(ao_atomindex)}")
            logging.info(f"Nº MOs: {len(MOs[0].coefficients)}")
            
            del molden
            if len(Smatrix.matrix) > 0:
                s_matrix = Smatrix.matrix
                logging.info(f"Matrix S: {len(s_matrix)}")
            else:
                n = len(MOs[0].coefficients)
                s_matrix = [[0] * n for _ in range(n)]
            del Smatrix

        elif orca_out != None:
            MOs = orca_out.MOs
            id_list_homo, id_list_lumo = self.__get_index_homo_lumo(MOs=MOs)
            id_list_homo += 2; id_list_lumo += 1 # Por que a contagem e MOs é do zero
            mo_coefficients = np.array([mo.coefficients for mo in MOs[:id_list_homo]], dtype=np.float64)
            ao_atomindex = orca_out.ao_atomindex + 1 # No orca começa a contagem de 0 e nos demais por 1; O codigo foi montado para trablhar iniciando de 1. 
            s_matrix = orca_out.overlap_matrix

            del orca_out
                
        del MOs
        s_matrix = np.array(s_matrix).astype(np.float64)

        ats_a = Selection(selection=sel_a, mol=pdb).result.data.inid +1
        logging.debug(f"Index atoms sel_a '{sel_a}': {', '.join(map(str, ats_a))}")

        ats_b = Selection(selection=sel_b, mol=pdb).result.data.inid +1 

        logging.debug(f"Index atoms sel_b '{sel_b}': {', '.join(map(str, ats_b))}")
        

        ao_ids_a = np.array([i for i, ao_id in enumerate(ao_atomindex) if ao_id in ats_a])
        ao_ids_b = np.array([i for i, ao_id in enumerate(ao_atomindex) if ao_id in ats_b])

        if gpu:
            m, b = core.intermolecular_eletron_density_two_selection_gpu(ao_ids_a=ao_ids_a, ao_ids_b=ao_ids_b, mo_coefficients=mo_coefficients, s_matrix=s_matrix, gpu_id=gpu_id)
        else:
            m, b = core.intermolecular_eletron_density_two_selection(ao_ids_a=ao_ids_a, ao_ids_b=ao_ids_b, mo_coefficients=mo_coefficients, s_matrix=s_matrix)
        logging.info(f"IED by Mulliken: {m}\nIED by OB: {b}")
        return (m, b)
    

    def matrix(self, pdb:PDB=PDB, aux=None, molden=None, orca_out=None, Smatrix:OverlapMatrix=OverlapMatrix, path_out:str=f"", out_format = "", gpu:bool=False, gpu_id:int = 0, write:bool=False) -> tuple:
        """
        Calculate the intermolecular electron density (IED) matrix between atoms in a PDB structure.
        Parameters
        ----------
        pdb : PDB object
            PDB object containing the molecular structure.
        aux : AUX object, optional 
            AUX object containing molecular orbital information. Default is None.
        molden : Molden object, optional
            Molden object containing molecular orbital information. Default is None.
        orca_out : OrcaOutput object, optional
            OrcaOutput object containing molecular orbital information. Default is None.
        Smatrix : OverlapMatrix object, optional
            OverlapMatrix object containing overlap matrix information. Default is an empty OverlapMatrix.
        path_out : str, optional
            Path to save the output files. Default is an empty string.
        out_format : str, optional
            Output file format, either 'npy' or 'json'. Default is an empty string.
        gpu : bool, optional
            Whether to use GPU acceleration. Default is False.
        gpu_id : int, optional
            GPU device ID to use if GPU acceleration is enabled. Default is 0.
        write : bool, optional
            Whether to write the output files. Default is False.
        Returns
        -------
        tuple
            A tuple containing the Mulliken and Bader IED matrices.
        Raises  
        ------
        ValueError
            If no molecular orbital data source is provided.
        Notes
        -----
        This function calculates the intermolecular electron density (IED) matrix using molecular orbital coefficients
        and overlap matrix data from various sources (AUX, Molden, OrcaOutput). It supports optional GPU acceleration
        for the calculations and can output the results in either NumPy binary or JSON format.
        """
        if path_out == "":    
            path_out = pdb.name

        if aux != None:
            MOs = aux.get_MOs_objs()
            id_list_homo, id_list_lumo = self.__get_index_homo_lumo(MOs=MOs)
            mo_coefficients = np.array([mo.coefficients for mo in MOs[:id_list_homo+1]], dtype=np.float64)
            ao_atomindex = np.array(aux.ao_atomindex)
            s_matrix = aux.get_S_objs().matrix
            del aux
            
        elif molden != None:
            MOs = molden.MOs
            id_list_homo, id_list_lumo =  self.__get_index_homo_lumo(MOs=MOs)
            mo_coefficients = np.array([mo.coefficients for mo in MOs[:id_list_homo+1]], dtype=np.float64)
            ao_atomindex = np.array(molden.get_ao_atomindex())
            
            logging.info(f"Nº AO inex: {len(ao_atomindex)}")
            logging.info(f"Nº MOs: {len(MOs[0].coefficients)}")
            
            del molden
            if len(Smatrix.matrix) > 0:
                s_matrix = Smatrix.matrix
                logging.info(f"Matrix S: {len(s_matrix)}")

            else:
                n = len(MOs[0].coefficients)
                s_matrix = [[0] * n for _ in range(n)]
            del Smatrix

        elif orca_out != None:
            MOs = orca_out.MOs
            id_list_homo, id_list_lumo = self.__get_index_homo_lumo(MOs=MOs)
            id_list_homo += 2; id_list_lumo += 1 # Por que a contagem e MOs é do zero
            mo_coefficients = np.array([mo.coefficients for mo in MOs[:id_list_homo]], dtype=np.float64)
            ao_atomindex = orca_out.ao_atomindex + 1 # No orca começa a contagem de 0 e nos demais por 1; O codigo foi montado para trablhar iniciando de 1. 
            s_matrix = orca_out.overlap_matrix

            del orca_out
                
        del MOs
        #list_ats = [at.id+1 for at in pdb.atoms]
        list_ats = pdb.data.inid + 1
        s_matrix = np.array(s_matrix).astype(np.float64)
        
        if gpu:
            mulliken_matrix, bader_matrix = core.matrix_intermolecular_eletron_density_numba_gpu(ats=list_ats, ao_atomindex=ao_atomindex, mo_coefficients=mo_coefficients, s_matrix=s_matrix, gpu_id=gpu_id)
        
        else:
            mulliken_matrix, bader_matrix = core.matrix_intermolecular_eletron_density_numba(ats=np.array(list_ats), ao_atomindex=ao_atomindex, mo_coefficients=mo_coefficients, s_matrix=s_matrix)

        if write:
            if out_format == "npy": 
                logging.info(f"{path_out}_matrix_IED_mulliken.npy...")
                np.save(f"{path_out}_matrix_IED_mulliken.npy", mulliken_matrix)

                logging.info(f"{path_out}_matrix_IED_OB.npy...")
                np.save(f"{path_out}_matrix_IED_OB.npy", bader_matrix)

            elif out_format == "json":
                resultado_m = self.__arrays_to_nested_dict(mulliken_matrix, list_ats)
                resultado_b = self.__arrays_to_nested_dict(bader_matrix, list_ats)

                with open(f"{path_out}_matrix_IED_mulliken.json", 'w') as arquivo_m:
                    logging.info(f"Writing .json Mulliken")
                    json.dump(resultado_m, arquivo_m, indent=2)
                
                with open(f"{path_out}_matrix_IED_OB.json", 'w') as arquivo_b:
                    logging.info(f"Writing .json OB")
                    json.dump(resultado_b, arquivo_b, indent=2)

        return mulliken_matrix, bader_matrix
                

    def __annotate_atom_ranges(self, pdb:PDB, ax=None, marks:list=[]):
        """
        Annotate a group of consecutive yticklabels with a group name.

        Arguments:
        ----------
        pdb : PDB object
            PDB object containing the information about atoms and residues.
        ax : matplotlib.axes object (default None)
            The axis instance to annotate.
        """

        if ax is None:
            ax = plt.gca()

        # Obter os resíduos e seus intervalos de átomos
        residues = []
        last_res = (pdb.data.resseq[0], pdb.data.resnames[0], pdb.data.chainids[0])
        atom_group = [pdb.data.serials[0]]

        for ser, resn, chain, at in zip(pdb.data.resseq[1:], pdb.data.resnames[1:], pdb.data.chainids[1:], pdb.data.serials[1:]):
            current_res = (ser, resn, chain)
            if current_res == last_res:
                atom_group.append(at)
            else:
                residues.append((last_res, atom_group))
                last_res = current_res
                atom_group = [at]

        residues.append((last_res, atom_group))

    # Índices de átomos no eixo Y
        indices_at = pdb.data.serials.tolist()
        ymax = max(indices_at)

        # --- Anotações ---
        for (residue_info, atom_indices) in residues:
            resseq, resname, chain = residue_info
            first_atom_index = atom_indices[0]
            last_atom_index = atom_indices[-1]
            label = f"{resname} {resseq} {chain}"

            # Barras verticais ao lado esquerdo
            ax.plot([-1, -1], [first_atom_index, last_atom_index],
                    color='gray', linewidth=2, clip_on=False)

            # Nome ao lado esquerdo
            ax.annotate(label,
                        xy=(-0.5, (first_atom_index + last_atom_index) / 2),
                        xytext=(-50, 0),
                        textcoords='offset points',
                        ha='right', va='center',
                        fontsize=16, color='black', clip_on=False)

            # Linha horizontal acima dos átomos
            ax.plot([first_atom_index, last_atom_index],
                    [ymax + 2, ymax + 2],
                    color='gray', linewidth=2, clip_on=False)

            # Nome rotacionado acima
            ax.annotate(label,
                        xy=((first_atom_index + last_atom_index) / 2, ymax - 0.5),
                        xytext=(0, -70),
                        textcoords='offset points',
                        ha='center', va='center',
                        fontsize=16, rotation=-90,
                        color='black', clip_on=False)

        # Ajustar limites do gráfico
        ax.set_xlim(left=-1)

        # --- Marcação opcional de resíduos específicos ---
        if len(marks) >= 1:
            for start_idx, end_idx in marks:
                if start_idx - 1 < len(residues) and end_idx - 1 < len(residues):
                    start_atoms = residues[start_idx - 1][1]
                    end_atoms = residues[end_idx - 1][1]
                    xy = (float(end_atoms[0]), float(start_atoms[0]))
                    height = float(len(start_atoms))
                    width = float(len(end_atoms))

                    ax.add_patch(Rectangle( # type: ignore
                        xy=xy, height=height, width=width,
                        fill=False, edgecolor='red', lw=2, alpha=.9
                    ))

        plt.draw()


    def __annotate_residue(self, pdb:PDB, ax=None, marks:list=[]):
        """
        Annotate a group of consecutive yticklabels with a group name.

        Arguments:
        ----------
        pdb : PDB object
            PDB object containing the information about atoms and residues.
        ax : matplotlib.axes object (default None)
            The axis instance to annotate.
        """

        if ax is None:
            ax = plt.gca()
            
        if len(marks) >=1:
            from matplotlib.patches import Rectangle        
            for mark in marks:
                xy = (float(mark[1]), float(mark[0]))
                
                height = 1
                width = 1
                # Add a rectangle patch to highlight the residue
                ax.add_patch(Rectangle(xy=xy, height=height, width=width, fill=False, edgecolor='red', lw=1, alpha=.9))
                
        plt.draw()


    def __transform_data_by_residue(self, df: pd.DataFrame, pdb):
        atom_to_residue = {
            at: (ser, resn, chain)
            for at, ser, resn, chain in zip(pdb.data.inid, pdb.data.resseq, pdb.data.resnames, pdb.data.chainids)
        }
        
        residue_data = defaultdict(lambda: defaultdict(float))

        for atom_index, atom_interactions in df.items():
            residue_1 = atom_to_residue.get(int(atom_index))
            if residue_1:
                residue_key_1 = f"{residue_1[1]} {residue_1[0]} {residue_1[2]}"
                for interacting_atom_index, value in atom_interactions.items():
                    residue_2 = atom_to_residue.get(int(interacting_atom_index))
                    if residue_2:
                        residue_key_2 = f"{residue_2[1]} {residue_2[0]} {residue_2[2]}"
                        residue_data[residue_key_1][residue_key_2] += value

        residue_df = pd.DataFrame(residue_data).T.fillna(0.0)
        return residue_df


    def __replace_none_with_zero(self, obj):
        if isinstance(obj, dict):
            return {k: self.__replace_none_with_zero(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self.__replace_none_with_zero(elem) for elem in obj]
        elif obj is None:
            return 0
        else:
            return obj


    def plot_heatmap(self, df:pd.DataFrame, norm:str="No", savefig:bool=True, plot:bool=False, annotate_res:bool=False, 
                          pdb:PDB=PDB, marks:list=[], figsize=(39, 36), per_residue:bool=False, intramol:bool=True, cmap:str="Blues",
                          pep_bond:bool=True, nucleic_bond:bool=True, intrachain:bool=True, figname:str=""
                          ):
        """
        Plot heatmap from a DataFrame.
        
        Parameters
        ----------
        df : pd.DataFrame
            DataFrame containing the data to plot.
        norm : str, optional
            Normalization type. Default is "No".
        savefig : bool, optional
            Whether to save the figure. Default is True.
        plot : bool, optional
            Whether to display the plot. Default is False.
        annotate_res : bool, optional
            Whether to annotate residues. Default is False.
        pdb : PDB object, optional
            PDB object for residue information. Default is an empty PDB.
        marks : list, optional
            List of residue ranges to mark. Default is an empty list.
        figsize : tuple, optional
            Figure size. Default is (39, 36).
        per_residue : bool, optional
            Whether to aggregate data per residue. Default is False.
        intramol : bool, optional
            Whether to include intramolecular interactions. Default is True.
        pep_bond : bool, optional
            Whether to include peptide bond interactions. Default is True.
        nucleic_bond : bool, optional  
            Whether to include nucleic bond interactions. Default is True.
        intrachain : bool, optional
            Whether to include intrachain interactions. Default is True.
        figname : str, optional
            Filename to save the figure. Default is an empty string.
        cmap : str, optional
            Colormap for the heatmap. Default is "Blues".
        Returns
        -------
            None
        Notes
        -----
        This function plots a heatmap from the provided DataFrame, with options for normalization,
        residue annotation, and saving the figure.
        """

        df.columns = df.columns.astype(int)
        df.index =  df.index.astype(int)

        if not intramol:
            df = self.intramol_remove(df=df, pdb=pdb)

        if not pep_bond:
            df = self.remove_pep_bond(df=df, pdb=pdb)
        
        if not nucleic_bond:
            logging.debug("Removing nucleic bond.")
            df = self.remove_nucleic_bond(df=df, pdb=pdb)
        
        if not intrachain:
            logging.debug("Removing intrachain.")
            df = self.remove_intrachain(df=df, pdb=pdb)        
        
        if per_residue:
            df = self.__transform_data_by_residue(df=df, pdb=pdb)
            #df = df.applymap(lambda x: 0 if x >= 0 else x)
            #df = df * -1
            #df = ((df - np.min(df)) / (np.max(df) - np.min(df)))

        if norm != "No":
            df = self.normalize_dataframe(df=df, type_norm=norm)
        
        sns.set_theme()
        sns.set_theme("poster")
        f, ax = plt.subplots(figsize=figsize)
        sns.color_palette("Blues", as_cmap=True)
        sns.heatmap(df, annot=False, linewidths=.5, ax=ax, cmap=cmap)

        if annotate_res and not per_residue:
            self.__annotate_atom_ranges(pdb=pdb, ax=ax, marks=marks)
        
        elif annotate_res and per_residue:
            self.__annotate_residue(pdb=pdb, ax=ax, marks=marks)

        plt.tight_layout()
        if savefig:
            if figname == "":
                plt.savefig(f"{self.path_ied.split('.')[-2]}_heatmap.png")
            else:
                plt.savefig(figname)

        if plot:
            plt.show()
        
        plt.close()


    def plot_multiple_heatmaps(self, dfs:list, residue_names:list, ncols:int=3, figsize:tuple=(15, 5),  dpi:int=200, plot:bool=False, savefig:bool=True, figname:str="",cmap:str="Blues", **kwargs):
        """
        Plot multiple heatmaps in a grid layout.
        Parameters
        ----------
        dfs : list of pd.DataFrame
            List of DataFrames to plot.
        residue_names : list of str
            List of residue names corresponding to each DataFrame.
        ncols : int, optional
            Number of columns in the grid. Default is 3.
        figsize : tuple, optional
            Figure size. Default is (15, 5).
        dpi : int, optional
            Dots per inch for the figure. Default is 200.
        plot : bool, optional
            Whether to display the plot. Default is False.
        savefig : bool, optional
            Whether to save the figure. Default is True.
        figname : str, optional
            Filename to save the figure. Default is an empty string.
        cmap : str, optional
            Colormap for the heatmaps. Default is "Blues".
        **kwargs : additional keyword arguments
            Additional arguments to pass to plt.subplots().
        Returns
        -------
            None
        Notes
        -----
        This function plots multiple heatmaps in a grid layout, allowing for customization of the layout and
        appearance of the plots.
        """
        num_plots = len(dfs)
        nrows = (num_plots - 1) // ncols + 1
        
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize, sharey='row', dpi=dpi, **kwargs)
        
        for i, (df, residue_name) in enumerate(zip(dfs, residue_names)):
            row = i // ncols
            col = i % ncols
            ax = axes[row, col]
            sns.color_palette("Blues", as_cmap=True)
            sns.heatmap(df, annot=False, linewidths=.5, ax=ax, cmap=cmap)
            ax.set_title(residue_name)
            ax.set_ylabel('')

        #plt.tight_layout()
        cbar_ax = fig.add_axes([0.94, 0.15, 0.02, 0.7])  # Posição da barra de cores
        fig.colorbar(axes[0, 0].collections[0], cax=cbar_ax)
        plt.tight_layout(rect=[0, 0, 0.9, 1])
        
        if savefig:
            if figname == "":
                plt.savefig(f"{self.path_ied.split('.')[-2]}_heatmap.png")
            else:
                plt.savefig(figname)
        if plot:
            plt.show()
            
        plt.close()


    def plot_heatmap_ref(self, df:pd.DataFrame, ref:str, ref_b:str="all", norm:str="No", savefig:bool=True, 
                              plot:bool=False, pdb:PDB=PDB, figsize=(39, 36), ncols=3, cutoff:float=.0, intramol:bool=True, 
                              pep_bond:bool=False, nucleic_bond:bool=True, intrachain:bool=True, figname:str="", multi_file:bool=False, cmap:str="Blues"
                              ):
        """
        Plot heatmaps for residues based on a reference selection.
        Parameters
        ----------
        df : pd.DataFrame
            DataFrame containing the data to plot.
        ref : str
            Reference selection for rows.
        ref_b : str, optional
            Reference selection for columns. Default is "all".
        norm : str, optional
            Normalization type. Default is "No".
        savefig : bool, optional
            Whether to save the figure. Default is True.
        plot : bool, optional
            Whether to display the plot. Default is False.
        pdb : PDB object, optional
            PDB object for residue information. Default is an empty PDB.
        figsize : tuple, optional
            Figure size. Default is (39, 36).
        ncols : int, optional
            Number of columns in the grid. Default is 3.
        cutoff : float, optional
            Cutoff value for residue inclusion. Default is 0.0.
        intramol : bool, optional
            Whether to include intramolecular interactions. Default is True.
        pep_bond : bool, optional
            Whether to include peptide bond interactions. Default is False.
        nucleic_bond : bool, optional
            Whether to include nucleic bond interactions. Default is True.
        intrachain : bool, optional
            Whether to include intrachain interactions. Default is True.
        figname : str, optional
            Filename to save the figure. Default is an empty string.
        multi_file : bool, optional
            Whether to save multiple files. Default is False.
        cmap : str, optional
            Colormap for the heatmaps. Default is "Blues".
        Returns
        -------
            None
        Notes
        -----
        This function plots heatmaps for residues based on a reference selection, with options for normalization,
        residue inclusion based on a cutoff, and saving the figures in multiple files.
        """
        if norm != "No":
            df = self.normalize_dataframe(df=df,type_norm=norm)

        if not intramol:
            logging.debug("Removing intramolecular interactions.")
            df = self.intramol_remove(df=df, pdb=pdb)

        if not pep_bond:
            logging.debug("Removing peptide bond.")
            df = self.remove_pep_bond(df=df, pdb=pdb)

        if not nucleic_bond:
            logging.debug("Removing nucleic bond.")
            df = self.remove_nucleic_bond(df=df, pdb=pdb)
        
        if not intrachain:
            logging.debug("Removing intrachain.")
            df = self.remove_intrachain(df=df, pdb=pdb)     


        sel = Selection(selection=ref, mol=pdb)
        selb = Selection(selection=ref_b, mol=pdb)
        at_ids = sel.result.data.inid
        df = df.iloc[at_ids, :] 
        dfs = []
        res_names = []
        
        df.index = [f"{serial}:{name}" for serial, name in zip(sel.result.data.serials, sel.result.data.names)]

        for res in np.unique(selb.result.data.resinid):
            at_ids_res = np.where(selb.result.data.resinid == res)[0]
            df_temp = df.iloc[:,at_ids_res]
            df_temp.columns = [f"{serial}:{name}" for serial, name in zip(selb.result.data.serials[at_ids_res], selb.result.data.names[at_ids_res])]
            if df_temp.values.sum() > cutoff:
                dfs.append(df_temp.copy())
                res_names.append(f"{selb.result.data.resnames[at_ids_res][0]} {selb.result.data.resseq[at_ids_res][0]} {selb.result.data.chainids[at_ids_res][0]}")

        logging.info(f"{len(res_names)} residues IED...")
        if multi_file:
            for df, name in zip(dfs, res_names):
                sns.set_theme("poster")
                f, ax = plt.subplots(figsize=figsize)
                sns.heatmap(df, annot=False, linewidths=.5, ax=ax, cmap=cmap, cbar=True)
                ax.set_title(name)
                ax.set_ylabel(ref)
                if figname == "":
                    plt.savefig(f"{'_'.join(name.split())}_heatmap.png")
                else:
                    plt.savefig(figname+"_"+name+"_heatmap.png")
        else: 
            self.plot_multiple_heatmaps(dfs, res_names, savefig=savefig, plot=plot, ncols=ncols, figsize=figsize, figname=figname)


    def read_json_to_df(self, path:str) -> pd.DataFrame:
        """
        Parses a JSON file and converts it to a pandas DataFrame.
        Parameters
        ----------
        path : str
            Path to the JSON file.
        Returns
        -------
        pd.DataFrame
            DataFrame containing the data from the JSON file.
        Notes
        -----
        The JSON file should contain a nested dictionary structure where the outer keys represent row indices
        and the inner keys represent column indices. The values are the data points for the DataFrame.
        The DataFrame will have integer indices and columns, corresponding to the atom indices.
        The DataFrame will be stored in the `self.df` attribute.
        The path to the JSON file will be stored in the `self.path_ied` attribute.
        Example
        -------
        >>> ied = IED()
        >>> df = ied.read_json_to_df("path/to/matrix.json")  
        >>> print(df.shape)
        (100, 100)
        >>> print(df.head())
           0         1         2         3         4
        0  0.0       0.1       0.2       0.3       0.4
        1  0.1       0.0       0.5       0.6       0.7
        2  0.2       0.5       0.0       0.8       0.9
        3  0.3       0.6       0.8       0.0       1.0
        4  0.4       0.7       0.9       1.0       0.0  
        """
        self.path_ied = path
        with open(path, 'r') as arquivo:
            resultado = json.load(arquivo)

        resultado = self.__replace_none_with_zero(resultado)
        df = pd.DataFrame(resultado)
        self.json = resultado
        self.df = df
        return df
    

    def read_npy_to_df(self, path: str) -> pd.DataFrame:
        """
        Reads a numpy file and converts it to a pandas DataFrame.  
        Parameters
        ----------
        path : str
            Path to the numpy file.
        Returns
        -------
        pd.DataFrame
            DataFrame containing the data from the numpy file.
        Raises
        ------
        ValueError
            If the loaded matrix is not 2D.
        Notes
        -----
        The numpy file should contain a 2D matrix. If the matrix is not 2D, a ValueError will be raised.
        The DataFrame will have the same shape as the loaded matrix.
        The DataFrame will have integer indices and columns, corresponding to the atom indices.
        The DataFrame will be stored in the `self.df` attribute.
        The path to the numpy file will be stored in the `self.path_npy` attribute.
        Example
        -------
        >>> ied = IED()
        >>> df = ied.read_npy_to_df("path/to/matrix.npy")  
        >>> print(df.shape)
        (100, 100)
        >>> print(df.head())
           0         1         2         3         4
        0  0.0       0.1       0.2       0.3       0.4
        1  0.1       0.0       0.5       0.6       0.7
        2  0.2       0.5       0.0       0.8       0.9
        3  0.3       0.6       0.8       0.0       1.0
        4  0.4       0.7       0.9       1.0       0.0
        """

        self.path_ied = path
        matriz = np.load(path) 

        if matriz.ndim != 2:
            logging.error(f"The loaded matrix has shape {matriz.shape}, which is not 2D.")
            raise ValueError(f"The matrix must be 2D, but it has shape {matriz.shape}")

        df = pd.DataFrame(matriz)
        self.df = df
        return df
    

    def normalize_dataframe(self, df: pd.DataFrame, type_norm: str = "MinMaxLog") -> pd.DataFrame:
        """
        Normalize the values in a DataFrame using the specified normalization method.
        Parameters
        ----------
        df : pd.DataFrame
            DataFrame to be normalized.
        type_norm : str, optional
            Type of normalization to apply. Options are "MinMaxLog" and "MinMax". Default is "MinMaxLog".
        Returns
        -------
        pd.DataFrame
            Normalized DataFrame.
        Notes
        -----
        The "MinMaxLog" normalization applies a logarithmic transformation followed by min-max scaling.
        The "MinMax" normalization applies min-max scaling directly to the data.
        If the minimum and maximum values are equal, the DataFrame is filled with zeros to avoid division by zero.
        Example
        -------
        >>> ied = IED()
        >>> df = pd.DataFrame([[1, 2], [3, 4]])
        >>> normalized_df = ied.normalize_dataframe(df, type_norm="MinMaxLog")
        """
        df = df.copy()

        if type_norm == "MinMaxLog":
            # log do DF inteiro
            x = np.log(df.values + 1e-10)

            min_val = x.min()
            max_val = x.max()

            if min_val == max_val:
                df[:] = 0.0
            else:
                df[:] = (x - min_val) / (max_val - min_val)

            return df

        elif type_norm == "MinMax":
            x = df.values

            min_val = x.min()
            max_val = x.max()

            if min_val == max_val:
                df[:] = 0.0
            else:
                df[:] = (x - min_val) / (max_val - min_val)

            return df

        else:
            return df

    
    def intramol_remove(self, df:pd.DataFrame, pdb:PDB) -> pd.DataFrame:
        """ Remove intramolecular interactions from the dataframe.
        This function sets the values of the dataframe to 0 for all pairs of atoms that belong to the same residue.

		Parameters
		----------
            df : pd.DataFrame
                DataFrame containing the interaction values.
            
            pdb :PDB 
                PDB object containing the structure information.
        Returns
		-------
            pd.DataFrame
                DataFrame with intramolecular interactions removed.
        """
        df_values = df.values.copy()  
        resinid = pdb.data.resinid
        inid = pdb.data.inid

        for res in np.unique(resinid):
            idx = inid[resinid == res]

            df_values[np.ix_(idx, idx)] = 0


        df.iloc[:, :] = df_values

        return df
    

    def remove_pep_bond(self, df:pd.DataFrame, pdb:PDB) -> pd.DataFrame:
        """Remove peptide bonds from the dataframe.
        This function sets the values of the dataframe to 0 for all pairs of atoms that belong to consecutive residues in the protein chain.

        Parameters
        ----------
            df : pd.DataFrame
                DataFrame containing the interaction values.
            pdb : PDB
                PDB object containing the structure information.
        Returns
        -------
            pd.DataFrame
                DataFrame with peptide bonds removed.
        """
        pro = Selection(selection="protein", mol=pdb).result
        
        if pro.natoms == 0:
            return df
        
        df_values = df.values.copy() 

        resseqs = np.unique(pro.data.resseq)
        resinid = pro.data.resseq
        inid = pro.data.inid

        for i, res in enumerate(resseqs[:-1]):
            res_next = resseqs[i + 1]
            if res_next == res + 1:

                idx_res = inid[resinid == res]
                idx_next = inid[resinid == res_next]

                df_values[np.ix_(idx_res, idx_next)] = 0
                df_values[np.ix_(idx_next, idx_res)] = 0

        # atualiza o DataFrame original
        df.iloc[:, :] = df_values
        return df


    def remove_nucleic_bond(self, df: pd.DataFrame, pdb: PDB) -> pd.DataFrame:
        """
        Remove bonds between consecutive nucleotides (phosphodiester bonds)
        from the interaction DataFrame.
        Parameters
        ----------
            df : pd.DataFrame
                DataFrame containing the interaction values.
            pdb : PDB
                PDB object containing the structure information.
        Returns
        -------
            pd.DataFrame
                DataFrame with nucleic bonds removed.
        """
        nuc = Selection(selection="nucleic", mol=pdb).result

        if nuc.natoms == 0:
            return df
        
        df_values = df.values.copy() 
        resseqs = np.unique(nuc.data.resseq)
        resinid = nuc.data.resseq
        inid = nuc.data.inid

        for i, res in enumerate(resseqs[:-1]):
            res_next = resseqs[i + 1]
            
            if res_next == res + 1:
                idx_res = inid[resinid == res]
                idx_next = inid[resinid == res_next]

        
                df_values[np.ix_(idx_res, idx_next)] = 0
                df_values[np.ix_(idx_next, idx_res)] = 0

        
        df.iloc[:, :] = df_values
        return df
    
    def remove_intrachain(self, df:pd.DataFrame, pdb:PDB) -> pd.DataFrame:
        """ Remove intramolecular interactions from the dataframe.
        This function sets the values of the dataframe to 0 for all pairs of atoms that belong to the same chain.
        
        Parameters
        ----------
            df : pd.DataFrame
                DataFrame containing the interaction values.
            pdb :PDB
                PDB object containing the structure information.
        Returns
        -------
            pd.DataFrame
                DataFrame with intramolecular interactions removed.
        """
        chains = np.unique(pdb.data.chainids)
        df_values = df.values.copy()

        for chain in chains:
            sel_chain = Selection(selection=f"chain {chain}", mol=pdb).result
            ids_chain = sel_chain.data.inid

            # usa indexação matricial para zerar tudo de uma vez
            df_values[np.ix_(ids_chain, ids_chain)] = 0

        # reatribui os valores modificados
        df.iloc[:, :] = df_values
        return df
    

    def map_3D(self, pdb:PDB, df:pd.DataFrame, path_pdbout:str="./3D_map.pdb", sel:str="all", intramol:bool=False, pep_bond:bool=False, nucleic_bond:bool=False, intrachain:bool=True, norm:str="No", write:bool=False) -> PDB:
        """
        Map interaction data onto a PDB structure as B-factors.
        Parameters
        ----------
            pdb : PDB
                PDB object containing the structure information.
            df : pd.DataFrame
                DataFrame containing the interaction values.
            path_pdbout : str, optional
                Path to save the output PDB file. Default is "./3D_map.pdb".
            sel : str, optional
                Selection string to specify which atoms to map. Default is "all".
            intramol : bool, optional
                Whether to include intramolecular interactions. Default is False.
            pep_bond : bool, optional
                Whether to include peptide bond interactions. Default is False.
            nucleic_bond : bool, optional
                Whether to include nucleic bond interactions. Default is False.
            intrachain : bool, optional
                Whether to include intrachain interactions. Default is True.
            norm : str, optional
                Normalization type. Default is "No".
            write : bool, optional
                Whether to write the output PDB file. Default is False.
        Returns
        -------
            PDB
                PDB object with mapped B-factors.
        Notes
        -----
        This function maps interaction data from a DataFrame onto a PDB structure by setting the B-factors
        of the atoms based on the interaction values. It allows for various filtering options and normalization.
        """
        pdb_map = pdb.__copy__()
        
        pdb_map.data.bfactors = np.zeros(pdb.natoms)

        df.columns = df.columns.astype(int)
        df.index =  df.index.astype(int)

        if not intramol:
            logging.debug("Removing intramol.")
            df = self.intramol_remove(df=df, pdb=pdb)

        if not pep_bond:
            logging.debug("Removing pepitid bond.")
            df = self.remove_pep_bond(df=df, pdb=pdb)

        if not nucleic_bond:
            logging.debug("Removing nucleic bond.")
            df = self.remove_nucleic_bond(df=df, pdb=pdb)
        
        if not intrachain:
            logging.debug("Removing intrachain.")
            df = self.remove_intrachain(df=df, pdb=pdb)


        if sel == "all":
            row_sums = df.sum(axis=1)

            for i, _ in enumerate(pdb_map.data.bfactors):
                pdb_map.data.bfactors[i] = row_sums.iloc[i] - df.iloc[i, i]

            if norm != "No":
                pdb_map.data.bfactors = self.normalize_dataframe(df=pd.DataFrame(pdb_map.data.bfactors), type_norm=norm)[0].to_numpy()

        elif sel != "all":
            sel_pdb = Selection(selection=sel, mol=pdb).result
            ids_sel = sel_pdb.data.inid

            row_sums = df.sum(axis=1)
           
            for i in pdb_map.data.inid:
                if i in ids_sel:
                    pdb_map.data.bfactors[i] = row_sums.iloc[i] - df.iloc[i, i]

            if norm != "No":
                pdb_map.data.bfactors = self.normalize_dataframe(df=pd.DataFrame(pdb_map.data.bfactors), type_norm=norm)[0].to_numpy()

            logging.info(f"Writing PDB file {path_pdbout}.")        
        
        if write:
            pdb_map.write(path=path_pdbout)
            #pdb_map.show()
        return pdb_map


    def radial_distribution(self, sel:str, pdb:PDB, df:pd.DataFrame, path_out:str="./IED_radial_distribution.json", selb:str = "all", 
                               dis_max:float=20.0, intramol:bool=True, pep_bond:bool=True, line:bool=True, plot:bool=False, savefig:bool=True, intrachain:bool=False
                                ):
        """Radial distribution of IED around a selection of atoms.

        Parameters
        ----------
            sel : str
                Selection string to select the atoms of interest.
            pdb : PDB
                PDB object containing the structure information.
            df : pd.DataFrame
                DataFrame containing the IED values between atom pairs.
            path_out : str, optional
                Path to save the output JSON file (default is "./IED_radial_distribution.json").
            dis_max : float, optional
                Maximum distance to consider for the radial distribution (default is 20.0).
            remove_intramol : bool, optional
                Whether to remove intramolecular interactions (default is True).
            remove_pep_bond : bool, optional
                Whether to remove peptide bond interactions (default is True).
        Returns
        -------
            None
        Notes
        ----- 
            The output JSON file will contain a list of dictionaries, each with the keys distance, IED, and atoms (atomid:atom name:resname:).
            A scatter plot of the radial distribution will be saved as a PNG file with the same name as the JSON file.
        """
        
        sel_ids = Selection(selection=sel, mol=pdb).result.data.inid
        selb_ids = Selection(selection=selb, mol=pdb).result.data.inid
        dis_mat = pdb.get_distance_matrix()
        
        if not intramol: 
            df = self.intramol_remove(df=df, pdb=pdb)
            logging.debug("Removing intramol.")

        if not pep_bond:
            df = self.remove_pep_bond(df=df, pdb=pdb)
            logging.debug("Removing pep bond.")
        
        if not intrachain:
            df = self.remove_intrachain(df=df, pdb=pdb)
            logging.info("Removing intrachain.")
            

        distances = []
        ieds = []
        ied_atoms = []

        for i in sel_ids:
            for j in selb_ids:
                if i == j:
                    continue
                dist = dis_mat[i, j]
                if dist <= dis_max and  df.iloc[i, j] > 10e-5:
                    distances.append(dist)
                    ieds.append(df.iloc[i, j])
                    ied_atoms.append((f"{pdb.data.serials[i]}:{pdb.data.names[i]}:{pdb.data.resnames[i]} - {pdb.data.serials[j]}:{pdb.data.names[j]}:{pdb.data.resnames[j]}"))

        
        data = [{"distance": float(d), "ied": float(i), 'atoms': at} for d, i, at in zip(distances, ieds, ied_atoms)]
        data = sorted(data, key=lambda x: x['distance'])

        with open(path_out, "w") as f:
            json.dump(data, f, indent=2)

        # Plotar
        sns.set(style="whitegrid")
        plt.figure(figsize=(8, 6))

        # Scatter dos pontos reais
        sns.set_palette(["#DD8452", "#4C72B0" ])
        if line:
            sns.histplot(x=distances, weights=ieds, bins=50, stat="density", kde=True, alpha=0.3)

        sns.scatterplot(x=distances, y=ieds, s=10)
        #sns.kdeplot(x=distances, color="red", linewidth=2, label="KDE")



        plt.xlabel("Distance (Å)")
        plt.ylabel("Int")
        plt.title("Radial Distribution of IED")
        plt.legend()
        plt.tight_layout()

        if savefig:
            plt.savefig(path_out.replace(".json", ".png"), dpi=300)

        if plot:
            plt.show()

        plt.close()
        logging.info(f"Radial distribution saved to {path_out} and plot saved as {path_out.replace('.json', '.png')}")


