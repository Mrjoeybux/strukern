from abc import ABC, abstractmethod
from difflib import get_close_matches
from rdkit import Chem
import glob
import os
import json
import numpy as np
import ctypes as ct
import multiprocessing as mp
datasets = ["Ahneman", "Musk", "3d", "3.5d", "images", "Arabidopsis", "vertebrates", "benchmark_graph"]
DATAPATH = "../../data/"

class DataLoader:

    def __init__(self, dataset_name: str):
        self.path = DATAPATH
        self.datasets = {"3d" : self.load_3d,
                         "3.5d": self.load_3point5d,
                         "cifar": self.load_cifar,
                         "mnist": self.load_mnist,
                         "3DQSAR": self.load_3DQSAR,
                         "aids": self.load_AIDS,
                         "mutagenicity": self.load_mutagenicity,
                         "protein": self.load_mutagenicity,
                         "musk": self.load_musk
                         }
        if dataset_name in self.datasets:
            self.load_func = self.datasets[dataset_name]
        else:
            msg = "{0} is not in the list of available datasets!\n".format(dataset_name)
            matches = get_close_matches(dataset_name, datasets)
            if not len(matches) == 0:
                msg += "Did you mean?:\n"
                for match in matches:
                    msg += str(match) + "\n"
            raise TypeError(msg)

    def load_3d(self):
        path = os.path.join(self.path, "catalysts/gxl/")

    def load_3point5d(self):
        pass
        path = os.path.join(self.path, "catalysts/gxl/")
        allfiles = []
        with open(path + "ee.json", "r") as f:
            label_dict = json.load(f)
        labels = []
        max_num_conformations = 32
        for i in [1, 2, 3, 4, 5, 6, 7, 8]:
            for letter in ["a", "b", "c", "d", "e"]:
                files = [os.path.split(file)[1] for file in glob.glob(path + "3.5d/" + str(i) + letter + "*.gxl")]
                ### Pad with empty strings for ease of use in creating the shared array
                if len(files) < max_num_conformations:
                    files += [""]*(max_num_conformations - len(files))
                allfiles.append(files)
                labels.append(label_dict[str(i)][letter])
        return np.asarray(allfiles), np.array(labels)

    def load_cifar(self):
        pass

    def load_mnist(self):
        pass

    def load_3DQSAR(self):
        pass

    def load_AIDS(self):
        pass

    def load_mutagenicity(self):
        pass

    def load_Protein(self):
        pass

    def load_musk(self):
        pass

    def load(self):
        return self.load_func()


class Mol2GXL(ABC):

    def __init__(self):
        self.mol = None

    @abstractmethod
    def generate_atom_descriptors(self, atom_index: int) -> dict:
        pass

    @abstractmethod
    def generate_bond_descriptors(self, bond_index1: int, bond_index2: int) -> dict:
        pass

    def _write_atom(self, atom_index):
        descriptors = self.generate_atom_descriptors(atom_index)
        stream = "<node id=\"_" + str(atom_index) + "\">"
        for name in descriptors.keys():
            stream += "<attr name=\"" + name + "\"><TYPENAME>" + str(descriptors[name]) + "</TYPENAME></attr>"
        stream += "</node>"
        return stream

    def _write_bond(self, bond_index):
        descriptors = self.generate_bond_descriptors(bond_index)
        bond_index1 = self.mol.GetBondWithIdx(bond_index).GetBeginAtomIdx()
        bond_index2 = self.mol.GetBondWithIdx(bond_index).GetEndAtomIdx()
        stream = "<edge from=\"_" + str(bond_index1) + "\" to=\"_" + str(bond_index2) + "\">"
        for name in descriptors.keys():
            stream += "<attr name=\"" + name + "\"><TYPENAME>" + str(descriptors[name]) + "</TYPENAME></attr>"
        stream += "</edge>"
        return stream

    def _set_mol(self, mol):
        self.mol = mol

    def generate_initial(self, ID: str):
        stream = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        stream += "<!DOCTYPE gxl SYSTEM \"http://www.gupro.de/GXL/gxl-1.0.dtd\">\n"
        stream += "<gxl xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n"
        stream += "<graph id=\"" + ID + "\" edgeids=\"false\" edgemode=\"undirected\">"
        return stream

    def write(self, mol: Chem.rdchem.Mol, path: str, ID: str) -> None:
        self._set_mol(mol)
        stream = self.generate_initial(ID)
        for atom_index in range(self.mol.GetNumAtoms()):
            stream += self._write_atom(atom_index)
        for bond_index in range(self.mol.GetNumBonds()):
            stream += self._write_bond(bond_index)
        stream += "</graph></gxl>"
        f = open(path, "w")
        f.write(stream)
        f.close()


class Catalyst2GXL(Mol2GXL):

    def generate_atom_descriptors(self, atom_index: int) -> dict:
        pt = self.mol.GetConformer().GetAtomPosition(atom_index)
        atom = self.mol.GetAtomWithIdx(atom_index)
        return {"AtomicNum": atom.GetAtomicNum(),
                       "x": pt.x,
                       "y": pt.y,
                       "z": pt.z,
                       "charge": atom.GetProp("_TriposPartialCharge")}

    def generate_bond_descriptors(self, bond_index: int) -> dict:
        return {"BondType": self.mol.GetBondWithIdx(bond_index).GetBondType()}


def threeD_to_gxl():
    filepath = "../../data/catalysts/mol/3d/"
    savepath = "../../data/catalysts/gxl/3d"
    files = glob.glob(filepath + "*.mol2")
    catalyst_writer = Catalyst2GXL()
    for filename in files:
        ID = os.path.split(filename)[1][:-5]
        savename = os.path.join(savepath, ID + ".gxl")
        mol = Chem.rdmolfiles.MolFromMol2File(filename, removeHs=False, sanitize=False)
        catalyst_writer.write(mol, savename, ID)

def threepointfiveD_to_gxl():
    filepath = "../../data/catalysts/mol/3.5d/"
    savepath = "../../data/catalysts/gxl/3.5d"
    files = glob.glob(filepath + "*.mol2")
    catalyst_writer = Catalyst2GXL()
    for filename in files:
        ID = os.path.split(filename)[1][:-5]
        savename = os.path.join(savepath, ID + ".gxl")
        mol = Chem.rdmolfiles.MolFromMol2File(filename, removeHs=False, sanitize=False)
        catalyst_writer.write(mol, savename, ID)

def to_shared_memory(source_arr, array_type=ct.c_double, with_type=False):
    source_arr = source_arr.reshape(-1)
    if array_type is ct.c_char:
        str_array = source_arr.tostring()
        mp_array = mp.RawArray(array_type, len(str_array))
        arr_type = source_arr.dtype
    elif array_type is ct.c_int:
        mp_array = mp.RawArray(array_type, source_arr.shape[0])
        arr_type = np.dtype(array_type)
    else:
        mp_array = mp.RawArray(array_type, source_arr.shape[0])
        arr_type = np.dtype(array_type)
    np_array = np.frombuffer(mp_array, dtype=arr_type)
    np_array[:] = source_arr
    if with_type:
        return mp_array, arr_type
    else:
        return mp_array

def data(sm_X, X_type, sz, args=None):
    data = np.frombuffer(sm_X, dtype=X_type).reshape(sz, -1)
    print(data.shape)
    return data[args]

if __name__ == "__main__":
    #threeD_to_gxl()
    #threepointfiveD_to_gxl()
    dl = DataLoader("3.5d")
    x, y = dl.load()
    sz, dim = x.shape
    sm_x, arr_type = to_shared_memory(x, ct.c_char, with_type=True)
    new_x = data(sm_x, arr_type, sz, [0, 1])
    print(new_x)