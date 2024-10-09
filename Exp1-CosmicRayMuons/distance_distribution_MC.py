'''
We perform a Monte Carlo simulation of the following problem:

There are two rectangular horizontal boards, each having 43 x 55 cm,
and one of the boards is placed heigh h directly above the other.

Muons come in from above. We detect a muon if it passes through both boards.
The angular distribution of the muons is uniform (so it's proportional to cos^2(theta)).

We will perform a reweighted Monte Carlo simulation where we randomly select a point on the top board,
then randomly generate a direction for the muon, and check if the muon passes through the bottom board.
if it does, we record the distance between the point on the top board and the point on the bottom board.
'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import scipy.special
import math
import random
from tqdm import tqdm

L = 43
W = 55

def generate_point_on_top_board():
    '''
    Generates a random point on the top board
    '''
    x = random.uniform(0, L)
    y = random.uniform(0, W)
    return x, y

def generate_direction():
    '''
    Generates a random direction for the muon. 
    
    generate a point on the unit sphere, 
    and then accept it with probability cos(theta)^2
    '''
    # generate a random point on the unit sphere
    x = np.random.normal(0, 1)
    y = np.random.normal(0, 1)
    z = np.random.normal(0, 1)
    norm = np.sqrt(x**2 + y**2 + z**2)
    x /= norm
    y /= norm
    z /= norm
    # accept the point with probability cos(theta)^2 = z^2
    theta = np.arccos(z)
    if np.random.uniform(0, 1) < z**2:
        return theta, np.arctan2(y, x)
    else:
        return generate_direction()

def MC_sample(d):
    '''
    generate a point on top board and an angle,
    then check if the muon passes through the bottom board
    and return the distance between the two points if it does
    Parameters:
        d: the distance between the two boards
    '''
    x, y = generate_point_on_top_board()
    theta, phi = generate_direction()
    # check if the muon passes through the bottom board
    # the equation of the line is given by:
    # new point is x + cos(phi)*tan(theta)*d, y + sin(phi)*tan(theta)*d
    x_new = x + np.cos(phi)*np.tan(theta)*d
    y_new = y + np.sin(phi)*np.tan(theta)*d
    if 0 <= x_new <= L and 0 <= y_new <= W:
        return np.sqrt((x_new - x)**2 + (y_new - y)**2 + d**2)
    else:
        return None
    
def MC_simulation(d, N, verbose=False):
    '''
    Perform a Monte Carlo simulation of the problem
    Parameters:
        d: the distance between the two boards
        N: the number of samples
    Returns:
        the distance distribution
    '''
    distances = []
    while len(distances) < N:
        if verbose and np.random.uniform(0, 1) < 0.00001:
            print(len(distances))
        distance = MC_sample(d)
        if distance is not None:
            distances.append(distance)
    return distances


def get_mean_dist(d, N, verbose=False):
    data = MC_simulation(d, N, verbose)
    return np.mean(data), np.std(data)/np.sqrt(len(data))

if __name__ == '__main__':
    data = MC_simulation(d=100, N=10000)
    # print the mean distance
    print(np.mean(data))
    print(np.std(data)/np.sqrt(len(data)))


    # plt.hist(data, bins=100, density=True)
    # plt.show()