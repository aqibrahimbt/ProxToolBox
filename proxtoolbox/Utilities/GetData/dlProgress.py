import sys

def dlProgress(counter, blocksize, size):
    p = counter*blocksize*100.0/size
    sys.stdout.write("\rProgress: %d%%" % p)
    sys.stdout.flush()

