import math
import random
import sys


class Point:
    def __init__(self, arr):
        self.arr = arr

    # calculate distance between two points
    def distance(self, other):
        sum = 0
        for i in range(len(self.arr)):
            sum += (self.arr[i] - other.arr[i])**2
        return math.sqrt(sum)

#read data from file, I'm assuming the input is correct they are all the in the same dimension
def read_data(filename):
    points = []
    f = open(filename)
    for line in f:
        line = line.strip('\n')
        line = line.split(',')
        values = [float(i) for i in line]
        point = Point(values)
        points.append(point)
    f.close()
    return points

def assign_points(points, centroids):
    clusters = [[] for i in range(len(centroids))]
    for point in points:
        min_dist = sys.maxsize
        distances = [point.distance(centroid) for centroid in centroids]
        cluster_index = distances.index(min(distances))
        clusters[cluster_index].append(point)
    return clusters

def update_centroids(clusters):
    centroids = []
    for cluster in clusters:
        sum = [0 for i in range(len(cluster[0].arr))]
        for point in cluster:
            for i in range(len(point.arr)):
                sum[i] += point.arr[i]
        centroid = Point([x/len(cluster) for x in sum])
        centroids.append(centroid)
    return centroids

def diff(centroids, new_centroids, threshold):
    for i in range(len(centroids)):
        if centroids[i].distance(new_centroids[i]) > threshold:
            return False
    return True

def kmeans(K, points, iter = 300):
    epsilon = 0.0001
    centroids = points[:K]
    iteration = 0
    while True:
        #for each iteration calculate distance and assign points to clusters
        clusters = assign_points(points, centroids)
        new_centroids = update_centroids(clusters)

        if diff(centroids,new_centroids,epsilon) or (iteration == iter):
            break

        #update centroids
        iteration += 1
        centroids = new_centroids
    return centroids

def assign_clusters(k, points_lst):
    points = []
    for p in points_lst:
        point = Point([float(i) for i in p])
        points.append(point)
    centroids = kmeans(k, points)
    assigned = []
    for point in points:
        distances = [point.distance(centroid) for centroid in centroids]
        cluster_index = distances.index(min(distances))
        assigned.append(cluster_index)
    return assigned

def main(k, iter, filename):
    #assert iuput
    if (k <1 or not isinstance(k,int)):
        print("Invalid number of clusters!")
        return
    if (iter <1 or not isinstance(iter,int) or iter>1000):
        print("Invalid maximum iteration!")
        return
    
    #read data from file
    points = read_data(filename)

    #initialize k centroids to be the first k points
    centroids = points[:k]

    iteration = 0
    while True:
        #for each iteration calculate distance and assign points to clusters
        clusters = assign_points(points, centroids)
        new_centroids = update_centroids(clusters)

        if diff(centroids,new_centroids,0.001 or iteration == iter):
            break

        #update centroids
        iteration += 1
        centroids = new_centroids
            
    for centroid in centroids:
        formatted_coords = [format(coord, '.4f') for coord in centroid.arr]
        print(','.join(formatted_coords))

    return

if __name__ == '__main__':
    args = sys.argv
    if len(args) == 3:
        k = int(args[1])
        filename = args[2]
        main(k, 200, filename)
    elif len(args) == 4:
        k = int(args[1])
        iter = int(args[2])
        filename = args[3]
        main(k, iter, filename)
    else:
        print("An Erorr Has Occured!")
        exit(1)