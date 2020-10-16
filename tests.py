import sys
from src import strukern
import numpy as np
import sys
import time
#import imageio
import matplotlib.pyplot as plt
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