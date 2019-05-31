# Author : Ali Snedden
# Date   : May, 2019
# License: MIT
import random
import numpy as np
import time
import sys

def write_matrix(Matrix = None, FileName = None): 
    """
    ARGS:
        Matrix   : Matrix to write to filename
        FileName : 
    DESCRIPTION:
    RETURN:
    DEBUG:
    NOTES: 
    FUTURE:
    """
    f = open(FileName, "w+")
    for i in range(Matrix.shape[0]):
        for j in range(Matrix.shape[1]):
            f.write("{} ".format(Matrix[i,j]))
        f.write("\n")
    f.close()


def print_help(arg):
    """
    ARGS:
    DESCRIPTION:
    RETURN:
    DEBUG:
    NOTES: 
    FUTURE:
    """
    sys.stdout.write(
        "\nUSAGE : ./matrix_generator.py Ax Ay Bx By OutDir\n\n"
        "     Ax : dimension of A in x\n"
        "     Ay : dimension of A in y\n"
        "     Bx : dimension of B in x\n"
        "     By : dimension of B in y\n\n"
        "     OutDir : Directory to output data\n\n"
        )
    sys.exit(arg)
    


def main():
    """
    ARGS:
        Matrix   : Matrix to write to filename
        FileName : 
    DESCRIPTION:
    RETURN:
    DEBUG:
    NOTES: 
    FUTURE:
    """
    startTime = time.time()
    random.seed(42)
    nArg = len(sys.argv)
    ######### Get Command Line Options ##########
    if(nArg == 2 and (sys.argv[1][0:3] == "--h" or sys.argv[1][0:2] == "-h")):
        print_help(0)
    elif(nArg != 6):
        print_help(1)
    Ax = int(sys.argv[1])
    Ay = int(sys.argv[2])
    Bx = int(sys.argv[3])
    By = int(sys.argv[4])
    outDir = sys.argv[5]
    A=np.zeros([Ax,Ay])
    B=np.zeros([Bx,By])
    
    for i in range(A.shape[0]):
      for j in range(A.shape[1]):
        A[i,j] = random.randint(0,10)
    
    for i in range(B.shape[0]):
      for j in range(B.shape[1]):
        B[i,j] = random.randint(0,10)
    
    print("Numpy multiplication time: ")
    npStartTime = time.time()
    AB = np.dot(A,B)
    print("Run time : {:.4f} s".format((time.time() - npStartTime)))
    # 
    # array([[  7.,   0.,   3., ...,   5.,  10.,   9.],
    #        [  0.,   7.,   7., ...,   0.,  10.,   9.],
    #        [ 10.,  10.,   9., ...,   0.,   0.,   9.],
    #        ..., 
    #        [  3.,   3.,   8., ...,   8.,   5.,   2.],
    #        [  7.,   9.,   3., ...,   1.,   1.,  10.],
    #        [  5.,   8.,  10., ...,   9.,   2.,   0.]])
    
    # array([[ 10.,   9.,   7., ...,   8.,   9.,   9.],
    #        [  0.,   1.,   2., ...,   8.,  10.,   1.],
    #        [  1.,   9.,   5., ...,   0.,   8.,   2.],
    #        ..., 
    #        [  7.,   4.,   7., ...,   4.,   7.,   3.],
    #        [  3.,   7.,   9., ...,   6.,   3.,   1.],
    #        [  5.,   6.,   2., ...,   1.,   7.,   5.]])
    # 
    
    write_matrix(Matrix=A, FileName="{}/A.txt".format(outDir))
    write_matrix(Matrix=B, FileName="{}/B.txt".format(outDir))
    write_matrix(Matrix=AB, FileName="{}/AB.txt".format(outDir))
    
    print("Ended : %s"%(time.strftime("%D:%H:%M:%S")))
    print("Tota Time To Run : {:.4f} h".format((time.time() - startTime)/3600.0))

if __name__ == "__main__":
    main()

