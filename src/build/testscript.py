import strukern
import numpy as np
import sys

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

def test_image():
    path = "/home/mrjoeybux/coding/strukern/src/mnist_image_converted.dat"
    vals = []
    f = open(path, "r")
    for line in f:
        vals.append(int(line.split("\n")[0]))
    X = np.array(vals)
    X = np.resize(X, (28, 28))
    jp = strukern.imagekernels.vJPEGCompressionKernel()
    assert jp.compress(X, 1) == 610.0

if __name__ == "__main__":
    inp = int(sys.argv[1])
    if inp == 0:
        test_eigen()
    elif inp == 1:
        test_graph()
    elif inp == 2:
        test_image()