path = "/home/mrjoeybux/coding/strukern/src/mnist_image.dat"
f = open(path, "r")
vals = []
for line in f:
    split = line.split("\n")[0].split(" ")
    nums = [x for x in split if x != ""]
    vals += nums
f.close()
f = open("/home/mrjoeybux/coding/strukern/src/mnist_image_converted.dat", "w")
for x in vals:
    f.write(x + "\n")
f.close()