import strukern
import numpy as np
import sys
import time
import imageio
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
    """for i in range(n):
        start = time.time()
        jp.dot(X, X, {"compressionlevel" : 1})
        times.append(time.time() - start)
    print("Mean Execution time: {0:.3g} +/- {1:.3g}".format(np.mean(times), np.std(times)))
    print("Total Execution time: {0:.3g}".format(np.sum(times)))"""

def test_large_image():
    path = "/home/mrjoeybux/coding/strukern/src/build/random_black_and_white.png"
    X = imageio.imread(path)
    X = X[:, :, 0]
    method = strukern.imagekernels.ImageCompressionMethod
    jp = strukern.imagekernels.JPEGCompressionKernel(method.Vertical)
    n = 10
    times = []
    for i in range(n):
        start = time.time()
        jp.dot(X, X, {"compressionlevel" : 1})
        times.append(time.time() - start)
    print("Mean Execution time: {0:.3g} +/- {1:.3g}".format(np.mean(times), np.std(times)))
    print("Total Execution time: {0:.3g}".format(np.sum(times)))


def test_colour_image():
    path = "/home/mrjoeybux/coding/strukern/src/build/random_black_and_white.png"
    X = imageio.imread(path)
    method = strukern.datastructures.ImageCompressionMethod.Both
    params = strukern.datastructures.KernelParams()
    params.JPEGCompressionQuality = 1
    kernel = strukern.imagekernels.JPEGCompressionKernel(method)
    print(kernel.dot(X, X, params))

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

if __name__ == "__main__":
    inp = int(sys.argv[1])
    if inp == 0:
        test_eigen()
    elif inp == 1:
        test_graph()
    elif inp == 2:
        test_mnist_image()
    elif inp == 3:
        test_large_image()
    elif inp == 4:
        test_multi_instance_ged()
    elif inp == 5:
        test_colour_image()
    elif inp == 6:
        test_li()