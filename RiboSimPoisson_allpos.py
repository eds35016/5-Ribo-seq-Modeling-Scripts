'''
Created on Aug 29, 2022

@author: daily
'''

#Requires numpy and matplotlib packages
import sys
import random
import numpy as np
from numpy.random import default_rng
import matplotlib.pyplot as plt

if len(sys.argv) != 4:
    print('Usage: RiboSimPoisson_allpos.py read_count transcript_length initiation_rate_avg')
    exit()

#Set parameters
read_count = int(sys.argv[1])
#Transcript length counted in amino acids
transcript_length = int(sys.argv[2])
initiation_rate_avg = int(sys.argv[3])

#Create arrays
read_array = np.zeros(shape=(read_count, transcript_length))
ribo_dist = []
pois_dist = []

rng = default_rng()

#Fill read array with transcript arrays consisting of ribosome positions (0 = No ribosome in that position, 1 = Ribosome present)
#Final ribosome placements for each read determined through shifting a random number of times after the first passthrough
for j in range(len(read_array)):
    rand_int = random.randint(1, transcript_length)
    initiation_step = -1
    for i in range(len(read_array[j]) + rand_int):
        #Calculate new initiation rate, variance uses poisson distribution
        while initiation_step < 1:
            initiation_step = rng.poisson(initiation_rate_avg, 1)[0]
            pois_dist.append(initiation_step)
        read_array[j] = np.roll(read_array[j], 1)
        read_array[j][:1] = 0
        if (i/initiation_step).is_integer():
            read_array[j][:1] = 1
            initiation_step = -1
    if (j/100).is_integer():
        print("Creating reads: ", j, "/", len(read_array))
print("Creating reads:", len(read_array), "/", len(read_array))

#Plot distribution of ribosome initiation rates
#pois_array = np.asarray(pois_dist)
#plt.hist(pois_dist, bins=len(np.unique(pois_array)))
#plt.show()

#Find the ribosomes on each read and add their positions within their respective reads to an array
for j in range(len(read_array)):
    for i in range(len(read_array[j])):
        if read_array[j][i] != 0:
            ribo_dist.append(i + 1)
    if (j/100).is_integer():
        print("Parsing reads: ", j, "/", len(read_array))
print("Parsing reads:", len(read_array), "/", len(read_array))

ribo_array = np.asarray(ribo_dist)
with np.printoptions(threshold=np.inf):
    print()
    print("All Ribosome Positions:")
    print()
    print(ribo_array)

#Prediction model

plt.hist(ribo_array, bins=len(np.unique(ribo_array)), range=[1, transcript_length])
plt.show()

print()
print("Done.")
