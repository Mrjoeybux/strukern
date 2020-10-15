import numpy as np
import sys
from abc import abstractmethod

sys.path.append("../build/")
import strukern


class KernelFunction:
    def __init__(self):
        self.ID = None

    def compute_kernel_matrix(self, X, params=None):
        n = len(X)
        kmat = np.zeros((n, n))
        for i, xi in enumerate(X):
            for j, xj in enumerate(reversed(X), start=1):
                kmat[i, -j] = self.dot(xi, xj, params)
                if i + j != n:
                    kmat[-j, i] = kmat[i, -j]
                else:
                    break
        return kmat

    @abstractmethod
    def dot(self, x1, x2, params=None):
        pass

    def set_unique_ID(self, ID):
        self.ID = ID


class StrukernKernel(KernelFunction):
    def __init__(self, **kwargs):
        super(StrukernKernel, self).__init__()
        self.strukern_params = strukern.datastructures.KernelParams()
        self.strukern_kernel = self.init_strukern_kernel(**kwargs)

    @abstractmethod
    def init_strukern_kernel(self, **kwargs):
        pass

    @abstractmethod
    def bandwidth_to_params(self, bandwidth: dict) -> None:
        pass

    def dot(self, x1, x2, params: dict = None):
        self.bandwidth_to_params(params)
        return self.strukern_kernel.dot(x1, x2, self.strukern_params)

    def compute_kernel_matrix(self, X, params: strukern.datastructures.KernelParams = None):
        self.bandwidth_to_params(params)
        return self.strukern_kernel.computeKernelMatrix(X, self.strukern_params)


EDIT_COST_CONSTRUCTOR = {"CG": strukern.datastructures.GEDEditCosts.CatalystGauss,
                         "CE": strukern.datastructures.GEDEditCosts.CatalystEuclid,
                         "Cons": strukern.datastructures.GEDEditCosts.Constant}

GED_METHOD_CONSTRUCTOR = {"Bi": strukern.datastructures.GEDMethods.Bipartite}


class GEDKernel(StrukernKernel):
    def init_strukern_kernel(self, **kwargs):
        if "graph_dir" in kwargs:
            dir = kwargs["graph_dir"]
        else:
            raise ValueError("Please specify the directory which contains the .gxl files!")
        self.edit_costs = EDIT_COST_CONSTRUCTOR[
            kwargs["edit_costs"]] if "edit_costs" in kwargs else strukern.datastructures.GEDEditCosts.Constant
        method = GED_METHOD_CONSTRUCTOR[
            kwargs["ged_method"]] if "ged_method" in kwargs else strukern.datastructures.GEDMethods.Bipartite
        return strukern.graphkernels.GEDKernel(self.edit_costs, method, dir)

    def bandwidth_to_params(self, bandwidth: dict) -> None:
        if self.edit_costs in [strukern.datastructures.GEDEditCosts.CatalystGauss,
                               strukern.datastructures.GEDEditCosts.CatalystEuclid]:
            params = bandwidth["bandwidth"]
            edit_cost_params = {"base_node_ins": params[0],
                                "base_node_del": params[1],
                                "base_node_rel": params[2],
                                "base_edge_ins": params[3],
                                "base_edge_del": params[4],
                                "base_edge_rel": params[5]}
            self.strukern_params.MolecularEditCosts = edit_cost_params


JPEG_METHOD_CONSTRUCTOR = {"V": strukern.datastructures.ImageCompressionMethod.Vertical,
                           "H": strukern.datastructures.ImageCompressionMethod.Horizontal,
                           "B": strukern.datastructures.ImageCompressionMethod.Both}

JPEG_IMAGETYPE_CONSTRUCTOR = {"BW": strukern.datastructures.ImageType.BW,
                              "RGB": strukern.datastructures.ImageType.RGB}


class JPEGImageKernel(StrukernKernel):
    def init_strukern_kernel(self, **kwargs):
        method = JPEG_METHOD_CONSTRUCTOR[kwargs["image_compression_method"]] if "image_compression_method" in kwargs \
            else strukern.datastructures.ImageCompressionMethod.Vertical
        image_type = JPEG_IMAGETYPE_CONSTRUCTOR[kwargs["image_type"]] if "image_type" in kwargs \
            else strukern.datastructures.ImageType.BW
        return strukern.imagekernels.JPEGCompressionKernel(method, image_type)

    def bandwidth_to_params(self, bandwidth: dict) -> None:
        self.strukern_params.JPEGCompressionQuality = bandwidth["bandwidth"]


class LocalityImproved(StrukernKernel):
    def init_strukern_kernel(self, **kwargs):
        return strukern.stringkernels.LocalityImprovedKernel()

    def bandwidth_to_params(self, bandwidth: dict) -> None:
        li_params = bandwidth["bandwidth"]
        self.strukern_params.LocalityImproved = li_params


class ZlibStringKernel(StrukernKernel):
    def init_strukern_kernel(self, **kwargs):
        return strukern.stringkernels.ZlibCompressionKernel()

    def bandwidth_to_params(self, bandwidth: dict) -> None:
        zlib_compression_level = bandwidth["bandwidth"]
        self.strukern_params.ZlibCompressionLevel = zlib_compression_level


MULTI_INSTANCE_METHOD_CONSTRUCTOR = {"Max": strukern.datastructures.MultiInstanceMethod.Max,
                                     "Min": strukern.datastructures.MultiInstanceMethod.Min,
                                     "OA": strukern.datastructures.MultiInstanceMethod.OA,
                                     "Sum": strukern.datastructures.MultiInstanceMethod.Sum}


class RBFKernel(StrukernKernel):
    def init_strukern_kernel(self, **kwargs):
        return strukern.vectorkernels.RBFKernel()

    def bandwidth_to_params(self, bandwidth: dict) -> None:
        rbf_sigma = bandwidth["bandwidth"]
        self.strukern_params.RBFSigma = rbf_sigma


class MultiInstanceKernel(StrukernKernel):
    def __init__(self, **kwargs):
        self.multi_instance_bound_kernels = {"GED": strukern.graphkernels.GEDMultiInstance,
                                             "RBF": strukern.vectorkernels.VectorMultiInstance}
        self.constructors = {"GED": GEDKernel,
                             "RBF": RBFKernel}
        self.base_kernel = None
        super(MultiInstanceKernel, self).__init__(**kwargs)

    def init_strukern_kernel(self, **kwargs):
        if "base_kernel" not in kwargs:
            raise ValueError("Please specify the base_kernel argument!")
        base_kernel = kwargs["base_kernel"]
        if base_kernel not in self.constructors.keys():
            raise ValueError("Please specify a valid strukern base kernel for use with the multi-instance kernel!")
        method = MULTI_INSTANCE_METHOD_CONSTRUCTOR[kwargs["multi_instance_method"]] if "multi_instance_method" in kwargs \
            else strukern.datastructures.MultiInstanceMethod.Min
        self.base_kernel = self.constructors[base_kernel](**kwargs)
        return self.multi_instance_bound_kernels[base_kernel](self.base_kernel.strukern_kernel, method)

    def bandwidth_to_params(self, bandwidth: dict) -> None:
        self.base_kernel.bandwidth_to_params(bandwidth)
        self.strukern_params = self.base_kernel.strukern_params


if __name__ == "__main__":
    meth = strukern.datastructures.ImageCompressionMethod.Vertical
    image_type = strukern.datastructures.ImageType.RGB
    jpeg = JPEGImageKernel(compression_method=meth, image_type=image_type)
    multi_rbf = MultiInstanceKernel(base_kernel="RBF")
    """import os, glob
    multi_ged = MultiInstanceKernel(base_kernel="GED",
                                    graph_dir="../../data/catalysts/gxl/3.5d/",
                                    edit_costs=strukern.datastructures.GEDEditCosts.CatalystGauss)
    allfiles = []
    params = strukern.datastructures.KernelParams()
    param_dict = {"base_node_ins": 1.0,
              "base_node_del": 1.0,
              "base_node_rel": 1.0,
              "base_edge_ins": 1.0,
              "base_edge_del": 1.0,
              "base_edge_rel": 1.0}
    params.CatalystEUCLID = param_dict
    for i in [1, 2, 3, 4, 5, 6, 7, 8]:
        for letter in ["a", "b", "c", "d", "e"]:
            files = [os.path.split(file)[1] for file in glob.glob("../../data/catalysts/gxl/3.5d/" + str(i) + letter + "*.gxl")]
            allfiles.append(files)
    print(multi_ged.strukern_kernel.dot(allfiles[0], allfiles[1], params))"""
