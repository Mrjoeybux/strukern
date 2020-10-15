from abc import ABC, abstractmethod
import sys

sys.path.append("..")
from python_stuff.kernel_functions import GEDKernel, JPEGImageKernel, LocalityImproved, ZlibStringKernel, \
    MultiInstanceKernel
from difflib import get_close_matches
from rdkit import Chem
import glob
import os
import json
import numpy as np
import ctypes as ct
import multiprocessing as mp
import pickle
import mnist
from lxml import etree
import configparser
import argparse
from sklearn.model_selection import train_test_split, ParameterGrid, ParameterSampler
import json
import itertools
import logging
__logger = logging.getLogger('experiments')
__DATAPATH = "../../data/"
__LOGPATH = "../../results/logs/"
__MODELPATH = "../../results/models/"
_SEED = 42
_TEST_SIZE = 0.3

data_dict = {"3d": __DATAPATH + "catalysts/gxl/3d/",
             "3.5d": __DATAPATH + "catalysts/gxl/3.5d/",
             "HIV": __DATAPATH + "hiv/gxl/",
             "musk1": __DATAPATH + "musk/",
             "musk2": __DATAPATH + "musk/",
             "mutagenicity": __DATAPATH + "mutagenicity/data/",
             "protein": __DATAPATH + "protein/data/",
             "mnist": "N/A"}

kernels = ["WL", "GED", "MI", "JPEG", "OTHERIMAGEKERNEL", "LI", "ZLIB"]
kernel_dict = {"GED": GEDKernel,
               "MI": MultiInstanceKernel,
               "JPEG": JPEGImageKernel,
               "LI": LocalityImproved,
               "ZLIB": ZlibStringKernel,
               }


class DataLoader:

    def __init__(self, dataset_name: str, dataset_path: str = None):
        self.dataset_name = dataset_name
        self.test_size = _TEST_SIZE
        self.seed = _SEED
        self.datasets = {"3d": self.load_3d,
                         "3.5d": self.load_3point5d,
                         "cifar": self.load_cifar,
                         "mnist": self.load_mnist,
                         "3dqsar": self.load_3DQSAR,
                         "aids": self.load_AIDS,
                         "mutagenicity": self.load_mutagenicity,
                         "protein": self.load_Protein,
                         "musk1": self.load_musk1,
                         "musk2": self.load_musk2,
                         "HIV": self.load_HIV,
                         "Ahneman": self.load_ahneman}
        if dataset_name in self.datasets:
            self.load_func = self.datasets[dataset_name]
        else:
            msg = "{0} is not in the list of available datasets!\n".format(dataset_name)
            matches = get_close_matches(dataset_name, list(data_dict.keys()))
            if not len(matches) == 0:
                msg += "Did you mean?:\n"
                for match in matches:
                    msg += str(match) + "\n"
            raise TypeError(msg)

    def load_3d(self, path):
        allfiles = []
        with open(path.split("3d/")[0] + "ee.json", "r") as f:
            label_dict = json.load(f)
        labels = []
        for i in [1, 2, 3, 4, 5, 6, 7, 8]:
            for letter in ["a", "b", "c", "d", "e"]:
                ### Pad with empty strings for ease of use in creating the shared array
                allfiles.append(path + str(i) + letter + ".gxl")
                labels.append(label_dict[str(i)][letter])
        return np.asarray(allfiles), np.array(labels), None, None, None, None

    def load_3point5d(self, path):
        allfiles = []
        with open(path.split("3.5d/")[0] + "ee.json", "r") as f:
            label_dict = json.load(f)
        labels = []
        max_num_conformations = 32
        for i in [1, 2, 3, 4, 5, 6, 7, 8]:
            for letter in ["a", "b", "c", "d", "e"]:
                files = [os.path.split(file)[1] for file in glob.glob(path + str(i) + letter + "*.gxl")]
                ### Pad with empty strings for ease of use in creating the shared array
                if len(files) < max_num_conformations:
                    files += [""] * (max_num_conformations - len(files))
                allfiles.append(files)
                labels.append(label_dict[str(i)][letter])
        return np.asarray(allfiles), np.array(labels), None, None, None, None

    def load_cifar(self, path):
        # TODO
        pass

    def load_mnist(self, path):
        return mnist.train_images(), mnist.train_labels(), mnist.test_images(), mnist.test_labels(), None, None

    def load_3DQSAR(self, path):
        # TODO
        pass

    def load_AIDS(self, path):
        # TODO
        pass

    def load_mutagenicity(self, path):
        X_all = []
        y_all = []
        classnames = {"nonmutagen": 0, "mutagen": 1}
        for file in ["train.cxl", "test.cxl", "valid.cxl"]:
            X, y = self.process_gxl_collection_file(path + file, classnames)
            X_all.append(X), y_all.append(y)
        return X_all[0], y_all[0], X_all[1], y_all[1], X_all[2], y_all[2]

    def load_Protein(self, path):
        X_all = []
        y_all = []
        classnames = {"1": 1, "2": 2, "3": 3, "4": 4, "5": 5, "6": 6}
        for file in ["train.cxl", "test.cxl", "valid.cxl"]:
            X, y = self.process_gxl_collection_file(path + file, classnames)
            X_all.append(X), y_all.append(y)
        return X_all[0], y_all[0], X_all[1], y_all[1], X_all[2], y_all[2]

    def load_musk1(self, path):
        fname = path + "clean1.data"
        X = []
        y = []
        current_conformations = []
        with open(fname, "r") as f:
            for i, line in enumerate(f):
                values = line.split(",")
                if i == 0:
                    prev_molecules = values[0]
                if values[0] != prev_molecules:
                    if values[0][:4] == "MUSK":
                        y.append(1)
                    else:
                        y.append(0)
                    X.append(current_conformations)
                    current_conformations = []
                    prev_molecules = values[0]
                descriptors = [int(x) for x in values[2:-1]]
                descriptors.append(int(values[-1].split(".\n")[0]))
                current_conformations.append(np.array(descriptors))
        X.append(current_conformations)
        if values[0][:4] == "MUSK":
            y.append(1)
        else:
            y.append(0)
        # return X, np.array(y), None, None, None, None
        return np.array(X), np.array(y), None, None, None, None

    def load_musk2(self, path):
        fname = path + "clean2.data"
        X = []
        y = []
        current_conformations = []
        with open(fname, "r") as f:
            for i, line in enumerate(f):
                values = line.split(",")
                if i == 0:
                    prev_molecules = values[0]
                if values[0] != prev_molecules:
                    if values[0][:4] == "MUSK":
                        y.append(1)
                    else:
                        y.append(0)
                    X.append(current_conformations)
                    current_conformations = []
                    prev_molecules = values[0]
                descriptors = [int(x) for x in values[2:-1]]
                descriptors.append(int(values[-1].split(".\n")[0]))
                current_conformations.append(np.array(descriptors))
        X.append(current_conformations)
        if values[0][:4] == "MUSK":
            y.append(1)
        else:
            y.append(0)
        # return X, np.array(y), None, None, None, None
        return np.array(X), np.array(y), None, None, None, None

    def load_ahneman(self):
        # TODO
        pass

    def load_HIV(self, path):
        labelfile = path + "labels.dat"
        y = []
        X = []
        with open(labelfile, "r") as f:
            for i, line in enumerate(f):
                y.append(int(line))
                X.append("mol_" + str(i) + ".gxl")
        return np.array(X), np.array(y), None, None, None, None

    def load(self):
        ## Returns train_X, train_y, test_X, test_y, validate_X, validate_y
        ## if (test_X and test_y) == None then there is no predefined test set
        ## if (validate_X and validate_y) == None then there is no predefined validation set
        ## self.split creates a validation and/or test set
        # return self.load_func(data_dict[self.dataset_name])
        return self.split(self.load_func(data_dict[self.dataset_name]))

    def process_gxl_collection_file(self, filepath, classnames):
        X = []
        y = []
        with open(filepath) as f:
            xml = f.read()
        root = etree.fromstring(xml)
        enzymes = root.getchildren()[0]
        for child in enzymes.getchildren():
            info = str(etree.tostring(child)).split("<print file=\"")[1]
            X.append(info[:info.find(".gxl") + 4])
            fileclassinfo = info.split("class=\"")[-1][:]
            y.append(classnames[fileclassinfo[:fileclassinfo.find("\"")]])
        return np.array(X), np.array(y)

    def split(self, data_tuple):
        if (data_tuple[4] is None) and (data_tuple[5] is None):
            if (data_tuple[2] is None) and (data_tuple[3] is None):
                return self._no_predefined_split(data_tuple)
            else:
                return self._predefined_train_and_test_split(data_tuple)
        else:
            return data_tuple

    def _no_predefined_split(self, data_tuple):
        X_train_val, X_test, y_train_val, y_test = train_test_split(data_tuple[0],
                                                                    data_tuple[1],
                                                                    random_state=self.seed,
                                                                    test_size=self.test_size)
        X_train, X_val, y_train, y_val = train_test_split(X_train_val,
                                                          y_train_val,
                                                          random_state=self.seed,
                                                          test_size=self.test_size)
        return X_train, y_train, X_test, y_test, X_val, y_val

    def _predefined_train_and_test_split(self, data_tuple):
        X_train, X_val, y_train, y_val = train_test_split(data_tuple[0],
                                                          data_tuple[1],
                                                          random_state=self.__SEED,
                                                          __TEST_SIZE=self.__TEST_SIZE)
        return X_train, y_train, data_tuple[2], data_tuple[3], X_val, y_val


class KernelLoader:
    def __init__(self, config_key: str, kernel_key: str, dataset_key: str):
        self.config_key = config_key
        self.kernel_key = kernel_key
        self.dataset_key = dataset_key

    def load(self):
        kwargs = dict(config._sections[self.config_key])
        print(kwargs)
        if ("GED" in self.kernel_key) or ("GED" in self.config_key):
            kwargs["graph_dir"] = data_dict[self.dataset_key]
        return kernel_dict[self.kernel_key](**kwargs), self._load_param_range()

    def _load_param_range(self):
        if self.kernel_key == "MI":
            kernel_key = self._parse_MI_config_key()
        else:
            kernel_key = self.kernel_key
        vals = config.get('DISCRETE_PARAM_VALS', kernel_key)
        if vals == 'None':
            return None
        return json.loads(vals)

    def _parse_MI_config_key(self):
        if "GED" in self.config_key:
            return "GED"
        if "RBF" in self.config_key:
            return "RBF"


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


class HIV2GXL(Mol2GXL):
    def generate_atom_descriptors(self, atom_index: int) -> dict:
        # pt = self.mol.GetConformer().GetAtomPosition(atom_index)
        atom = self.mol.GetAtomWithIdx(atom_index)
        return {"AtomicNum": atom.GetAtomicNum()}

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


def hiv_to_gxl():
    filepath = "../../data/hiv/csv/"
    savepath = "../../data/hiv/gxl/"
    hiv_writer = HIV2GXL()
    count = 0
    with open(filepath + "HIV.csv", "r") as f, open(savepath + "labels.dat", "w") as g:
        for i, line in enumerate(f):
            if (i == 0) or (i % 2 == 1):
                continue
            smiles = line.split(",")[0]
            label = line.split(",")[2].split("\n")[0]
            mol = Chem.MolFromSmiles(smiles)  # , sanitize=False)
            hiv_writer.write(mol, savepath + "mol_" + str(count) + ".gxl", str(count))
            g.write(label + "\n")
            count += 1


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


def main():
    global config
    config = configparser.ConfigParser()
    config.read("../build/config.ini")
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--dataset',
                        choices=list(data_dict.keys()), required=True)
    parser.add_argument('--kernel', choices=kernels, required=True)

    parser.add_argument('--config_key', choices=[key for key in dict(config._sections) if key != "DISCRETE_PARAM_VALS"],
                        required=True)
    parser.add_argument('--random_search', nargs='?', default=0, type=int)
    parser.add_argument('--num_inner_cv_folds', nargs='?', default=5, type=int)
    parser.add_argument('--mode', nargs='?', default='regression', type=str, choices=['regression', 'classification'])
    parser.add_argument('--num_random_restarts', nargs='?', default=10, type=int)
    args = vars(parser.parse_args(sys.argv[1:]))
    print(args['random_search'])
    data = DataLoader(args["dataset"]).load()
    kernel_f, discrete_param_range = KernelLoader(args["config_key"], args["kernel"], args["dataset"]).load()
    # print(kernel_f)
    # print(type(discrete_param_range), discrete_param_range)
    import strukern
    print(args)
    if args['random_search'] == 0:
        param_combinations = discrete_param_search(discrete_param_range)
    else:
        param_combinations = discrete_param_search(discrete_param_range, random=True, num_samples=args['random_search'])
    print(create_log_and_model_filenames(args['dataset'], args['kernel'], args['config_key']))
    # params = {"bandwidth": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]}
    # X = data[0]
    # K = kernel_f.compute_kernel_matrix(X, params)


def discrete_param_search(discrete_parameter_range, random=False, num_samples=20):
    if discrete_parameter_range is None:
        return None
    if random:
        return list(ParameterSampler(discrete_parameter_range, n_iter=num_samples))
    else:
        return list(ParameterGrid(discrete_parameter_range))


def create_log_and_model_filenames(dataset, kernel, kernel_config):
    if kernel_config == 'None':
        return __LOGPATH + dataset + "_" + kernel + ".log", \
               __MODELPATH + dataset + "_" + kernel + ".dat"
    if kernel == "MI":
        return __LOGPATH + dataset + "_" + kernel_config + ".log", \
               __MODELPATH + dataset + "_" + kernel_config + ".dat"
    return __LOGPATH + dataset + "_" + kernel + "_" + kernel_config + ".log", \
           __MODELPATH + dataset + "_" + kernel + "_" + kernel_config + ".dat"


def run_experiment(data_tuple, kernel_function, parameter_combinations, log_fname, model_fname, num_inner_cv_folds, mode, num_random_restarts):
    pass


if __name__ == "__main__":
    main()
    """pass
    #threeD_to_gxl()
    #threepointfiveD_to_gxl()
    dl = DataLoader(sys.argv[1])
    data = dl.load()
    for i in range(6):
        shape = data[i].shape if data[i] is not None else "N/A"
        print("data {}, shape {}".format(i, shape))"""

    # X = dl.load()
    # print(X)
    # hiv_to_gxl()
    # print(X.shape)

    """x, y = dl.load()
    sz, dim = x.shape
    sm_x, arr_type = to_shared_memory(x, ct.c_char, with_type=True)
    new_x = data(sm_x, arr_type, sz, [0, 1])
    print(new_x)"""
