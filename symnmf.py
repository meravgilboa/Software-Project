import sys
import pandas as pd
import numpy as np
import math
import symnmfmodule as sn

np.random.seed(1234)

def initialize_H(n, k, W):
    m = np.mean(W)
    H = np.random.uniform(0, 2*math.sqrt(m/k), size=(n,k))
    return H.tolist()

#Python Interface Of Our Code
def main():
    if (len(sys.argv) != 4):
        print("An Error Has Occured")
        sys.exit(1)
    try:
        #Retreive inputs from CMD 
        k = int(sys.argv[1])
        goal = sys.argv[2]
        file_name = sys.argv[3]
        #Get the points from the file
        data = pd.read_csv(file_name, header=None)
        points = [point.tolist() for i, point in data.iterrows()]
        n = len(points)
        d = len(points[0])
        #Check that the input is valid
        if (k <= 0 or k >=len(points)):
            print("An Error Has Occured")
            sys.exit(1)
        #Run the algorithm subject to goal value
        if (goal == "symnmf"): #Perform full symnnmf and output H
            W = sn.norm(points, n, d)
            H = initialize_H(n, k, W)
            res = sn.symnmf(W, H, n, k)
        elif (goal == "sym"): # Caclculate and ouput the symilarity matrix
            res = sn.sym(points, n, d)
        elif (goal == "ddg"): #Calculate and output the diagonal degree matrix
            res = sn.ddg(points, n, d)
        elif (goal == "norm"): #calculate and output the normalized similarity matrix
            res = sn.norm(points, n, d)
        else:
            print("An Error Has Occured")
            sys.exit(1)
    except Exception as e:
        print("An Error Has Occured")
        sys.exit(1)
    #Print the result
    for row in res:
        print(",".join(str("{:.4f}".format(round(x, 4))) for x in row))       
    
if __name__ == '__main__':
    main()