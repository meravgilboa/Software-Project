import sys
import pandas as pd
import numpy as np
import kmeans as km
import math
import symnmfmodule as sn
import symnmf as snpy
import sklearn.metrics as skm

def main():
    if len(sys.argv) != 3:
        print ("An Error has occurred")
        sys.exit(1)
    try:
        #Retreive inputs from CMD 
        k = int(sys.argv[1])
        file_name = sys.argv[2]
        #Get the points from the file
        data = pd.read_csv(file_name, header=None)
        points = [point.tolist() for i, point in data.iterrows()]
        n = len(points)
        d = len(points[0])
        #Check that the input is valid
        if (k <= 0 or k >= len(points)):
            print("An Error Has Occured")
            sys.exit(1)
        #Run the algorithm subject to goal value
        W = sn.norm(points, n, d)
        H  = snpy.initialize_H(n, k, W)
        res = sn.symnmf(W, H, n, k) #H is now the solution
        clusters_list = sn.derive(res, n, k)
        kmeans_list = km.assign_clusters(k, points)
        nmf_sol = skm.silhouette_score(points, clusters_list)
        kmeans_sol = skm.silhouette_score(points, kmeans_list)

    except Exception as e:
        print ("An Error has occurred:", e)
        return
    

    print("nmf: {:.4f}".format(nmf_sol))
    print("kmeans: {:.4f}".format(kmeans_sol))

if __name__ == '__main__':
    main()