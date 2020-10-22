import sys
from src import strukern
import numpy as np
import sys
import time
#import imageio
import matplotlib.pyplot as plt
import rdkit.Chem
def test_eigen():
    string = strukern.stringkernels.TestCompressionKernel()
    #string = strukern.DiracStringKernel()
    a = ["hi", "hello", "bye", "hi"]
    K = string.computeKernelMatrix(a, {"compressionlevel": 1})
    """for i in range(K.shape[0]):
        K[i, i] = -1*i
        if i < K.shape[0] - 1:
            K[i, i + 1] = 1.0
            K[i + 1, i] = 1.0
    K[3, 0] = -5.0
    K[0, 3] = -5.0"""
    Knorm = string.normaliseKrein(K)
    print(K)
    print("\n")
    print(Knorm)

def test_graph():
    #help(strukern)
    adjlist = {0: [1, 2], 1: [0], 2: [0]}
    nodeattr = {0: [666], 1: [999], 2: [42]}
    edgeattr = {(0, 1): [333], (0, 2): [111]}
    g = strukern.Graphii(adjlist, nodeattr, edgeattr, Directed=False)
    print(g)
    #print(g.getNodeAttribute(i))


def test_mnist_image():
    path = "/home/mrjoeybux/coding/strukern/src/mnist_image_converted.dat"
    vals = []
    f = open(path, "r")
    for line in f:
        vals.append(int(line.split("\n")[0]))
    X = np.array(vals)
    X = np.resize(X, (28, 28))
    method = strukern.imagekernels.ImageCompressionMethod
    jp = strukern.imagekernels.JPEGCompressionKernel(method.Vertical)
    params = strukern.KernelParams()
    params.JPEGCompressionQuality = 100
    #print(params.JPEGCompressionQuality)
    print(jp.dot(X, X, params))
    n = 1000
    times = []


def test_multi_instance_ged():
    if __name__ == '__main__':
        gk = strukern.graphkernels
        ged_meth = gk.GEDMethods.BIPARTITE
        ged_ec = gk.GEDEditCosts.MABTS
        ged = gk.GEDKernel(ged_ec, ged_meth, "")
        multi_meth = strukern.MultiInstanceMethod.Sum
        ged_multi = strukern.GEDMultiInstance(ged, multi_meth)

def test_li():
    x = "hello friend"
    y = "howdy friend"
    li = strukern.stringkernels.LocalityImprovedKernel()
    params = strukern.datastructures.KernelParams()
    params.LocalityImproved = {"sub_window_length": 3,
                               "d1": 1,
                               "d2": 1}
    print(li.dot(x, y, params))


def test_ppm():
    x = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.\nUt enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. \nDuis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. \nExcepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."

    y = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.\nUt enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. \nDuis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. \nExcepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."

    z = "ipsum Loren dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.\nUt enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. \nDuis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. \nExcepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."
    import random
    base_size = 100
    times = []
    sizes = []
    size = 100
    passages = [x, z]
    X = [random.choice(passages) for i in range(size)]
    params = strukern.datastructures.KernelParams()
    params.StringCompressionLevel = 5

    distance_measure = strukern.datastructures.CompressionDistanceMeasure.NCD
    #distance_measure = strukern.datastructures.StringCompressionDistanceMeasure.DIFF

    ppm = strukern.stringkernels.PPMCompressionKernel(distance_measure)
    start = time.time()
    K = ppm.computeKernelMatrix(passages, params)
    end = time.time() - start
    print("done in {0:.3f} seconds!".format(end))
    print(K)
    vals, vecs = np.linalg.eigh(K)
    #print(vals)


def test_tthresh():
    shape = (10, 10, 10)
    x1 = np.random.rand(*shape)
    x2 = np.ones(shape)
    x3 = np.random.rand(*shape)
    measure = strukern.datastructures.CompressionDistanceMeasure.NCD
    cub = strukern.tensorkernels.TThreshTensorCompressionKernel(measure)
    params = strukern.datastructures.KernelParams()

    params.TensorCompressionLevel = 1
    params.TensorConcatDim = 2
    val = cub.computeRectangularKernelMatrix([x1, x2], [x3], params)
    #val = cub.computeKernelMatrix([x1, x2], params)
    #val = cub.dot(x1, x3, params)
    print(val)


def smiles_to_graph(smiles):
    mol = rdkit.Chem.MolFromSmiles(smiles)
    graph_mol = strukern.datastructures.Graphvv()
    for atom in mol.GetAtoms():
        desc = []
        desc.append(atom.GetAtomicNum())
        desc.append(atom.GetIsAromatic())
        desc.append(atom.GetHybridization())
        desc.append(atom.GetChiralTag())
        desc.append(atom.GetNumRadicalElectrons())
        desc.append(atom.GetFormalCharge())
        graph_mol.add_node(atom.GetIdx(), desc)
    for bond in mol.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        desc = []
        desc.append(bond.GetBondType())
        desc.append(bond.GetIsAromatic())
        graph_mol.add_edge(begin, end, desc)
    return graph_mol

from scipy.optimize import linear_sum_assignment as munkres

class OptimalMoleculeAssignmentKernel():
    def __init__(self, base_kernel):
        self.name = "oma_" + base_kernel.name
        self.kf = base_kernel

    def cost_matrix(self, x1, x2, params):
        # x1, x2 - rdkit Molecules
        # this method finds the optimal assignment of atoms in x1 to atoms in x2
        n, m = x1.GetNumAtoms(), x2.GetNumAtoms()
        cost_mat = np.zeros((n, m))
        for i in range(n):
            for j in range(m):
                cost_mat[i, j] = self.kf.dot(x1, x2, i, j, params)
        return cost_mat

    def dot(self, x1, x2, params):
        n, m = x1.GetNumAtoms(), x2.GetNumAtoms()
        if n <= m:
            cost_mat = self.cost_matrix(x1, x2, params)
        else:
            cost_mat = self.cost_matrix(x2, x1, params)
        print(cost_mat)
        idxs = munkres(-1*cost_mat)
        return cost_mat[idxs].sum()


class GaussianKernel:

    def dot(self, x1, x2, params):
        sq_diff = np.linalg.norm(x1 - x2)**2
        sigma = params["sigma"]
        return np.exp((-1.0*sq_diff)/(2*sigma**2))


class SubStructFeatureKernel():
    def __init__(self):
        super(SubStructFeatureKernel, self).__init__()
        self.kf = GaussianKernel()
        self.name = "substruct_feature"

    def dot(self, mol1, mol2, mol1_atomidx, mol2_atomidx, params=None):
        self.mol1, self.mol2 = mol1, mol2
        atom1 = mol1.GetAtomWithIdx(mol1_atomidx)
        atom2 = mol2.GetAtomWithIdx(mol2_atomidx)
        return self._neighbour_kernel(atom1, atom2, params)

    def _atom_kernel(self, atom1, atom2, params):
        if atom1.GetAtomicNum() != atom2.GetAtomicNum():
            return 0.0
        atom1_features = self._generate_atom_features(atom1)
        atom2_features = self._generate_atom_features(atom2)
        return self.kf.dot(atom1_features, atom2_features, params)

    def _neighbour_and_bond_generator(self, mol, atom):
        for neigh in atom.GetNeighbors():
            yield neigh, mol.GetBondBetweenAtoms(atom.GetIdx(), neigh.GetIdx())

    def _bond_kernel(self, bond1, bond2, params):
        bond1_features = self._generate_bond_features(bond1)
        bond2_features = self._generate_bond_features(bond2)
        if bond1_features[0] == bond2_features[0]:
            if bond1_features[1] == bond2_features[1]:
                return self.kf.dot(bond1_features, bond2_features, params)
        return 0.0

    def _generate_atom_features(self, atom):
        desc = []
        desc.append(atom.GetAtomicNum())
        # desc.append(atom.GetImplicitValence())
        # desc.append(atom.GetExplicitValence())
        desc.append(atom.GetIsAromatic())
        desc.append(atom.GetHybridization())
        desc.append(atom.GetChiralTag())
        desc.append(atom.GetNumRadicalElectrons())
        desc.append(atom.GetFormalCharge())
        return np.array(desc)

    def _generate_bond_features(self, bond):
        desc = []
        desc.append(bond.GetBondType())
        desc.append(bond.GetIsAromatic())
        # desc.append(bond.IsInRing())
        return np.array(desc)

    def _R0(self, atom1, atom2, params):
        cost_mat = []
        for neigh1, bond1 in self._neighbour_and_bond_generator(self.mol1, atom1):
            row = []
            for neigh2, bond2 in self._neighbour_and_bond_generator(self.mol2, atom2):
                row.append(self._atom_kernel(neigh1, neigh2, params)*self._bond_kernel(bond1, bond2, params))
            cost_mat.append(row)
        cost_mat = np.array(cost_mat)
        idxs = munkres(-1*cost_mat)
        return cost_mat[idxs].sum()/max(cost_mat.shape[0], cost_mat.shape[1])

    def _neighbour_kernel(self, atom1, atom2, params):
        radius = params["molecular_radius"] if "molecular_radius" in params else 2
        decay = params["decay"] if "decay" in params else 1
        Rl = 0
        for i in range(1, radius + 1):
            Rl += self._Rl(atom1, atom2, i, params)*(decay**i)
        return self._atom_kernel(atom1, atom2, params) + self._R0(atom1, atom2, params) + Rl

    def _Rl(self, atom1, atom2, i, params):
        if i == 1:
            r1 = 0
            count = 0
            for neigh1 in atom1.GetNeighbors():
                for neigh2 in atom2.GetNeighbors():
                    r1 += self._R0(neigh1, neigh2, params)
                    count += 1
            return r1/count
        rl_1 = 0
        count = 0
        for neigh1 in atom1.GetNeighbors():
            for neigh2 in atom2.GetNeighbors():
                rl_1 += self._Rl(neigh1, neigh2, i - 1, params)
                count += 1
        return rl_1/count

def test_graph():
    """x0_attr = [1, 2, 3]
    x0_index = 0
    edge_attr = [0, 1, 5]
    edge_attr1 = [0, 1, 6]
    x1_attr = [5, 2, 4]
    x1_index = 1

    A = strukern.datastructures.Graphvv()
    for i in range(3):
        A.add_node(i, [0, i])
    A.add_edge(0, 1, [0, 0, 1])
    A.add_edge(1, 2, [0, 0, 2])
    B = strukern.datastructures.Graphvv()
    for i in range(4):
        B.add_node(i, [0, i])
    B.add_edge(0, 1, [0, 0, 1])
    B.add_edge(1, 2, [0, 0, 2])
    B.add_edge(1, 3, [0, 0, 3])

    for i in range(3):
        print(A.get_neighbours(i))
    print()
    for i in range(4):
        print(B.get_neighbours(i))
    """

    benzene = "C1=CC=CC=C1"
    benzene_graph = smiles_to_graph(benzene)
    cyclohexane = "N[C@@H](C)C(=O)O"
    cyclohexane_graph = smiles_to_graph(cyclohexane)


    oma = strukern.graphkernels.OptimalMolecularAssignmentKernel()
    params = strukern.datastructures.KernelParams()
    params.MolecularRadius = 2
    params.AssignmentDecay = 0.1
    params.RBFSigma = 10
    print(oma.dot(benzene_graph, cyclohexane_graph, params))
    params = {"sigma": 10, "molecular_radius": 2, "decay": 0.1}
    oma_py = OptimalMoleculeAssignmentKernel(SubStructFeatureKernel())
    print(oma_py.dot(rdkit.Chem.MolFromSmiles(benzene), rdkit.Chem.MolFromSmiles(cyclohexane), params))


if __name__ == "__main__":
    inp = int(sys.argv[1])
    if inp == 0:
        test_eigen()
    elif inp == 1:
        test_graph()
    elif inp == 2:
        test_mnist_image()
    #elif inp == 3:
    #   test_large_image()
    elif inp == 4:
        test_multi_instance_ged()
    #elif inp == 5:
    #   test_colour_image()
    elif inp == 6:
        test_li()
    elif inp == 7:
        test_ppm()
    elif inp == 8:
        test_tthresh()
    elif inp == 9:
        test_graph()