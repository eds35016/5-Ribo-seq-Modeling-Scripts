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
    print('Usage: RiboSimPoisson.py read_count transcript_length initiation_rate_avg')
    exit()

#Set parameters
read_count = int(sys.argv[1])
#Transcript length counted in amino acids
transcript_length = int(sys.argv[2])
initiation_rate_avg = int(sys.argv[3])

#Create arrays
read_array = np.zeros(shape=(read_count, transcript_length))
first_ribo_array = np.zeros(read_count)
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

#Trim reads to first ribosome and add first ribosome position of each read to first_ribo_array (2 = trimmed position)
for j in range(len(read_array)):
    for i in range(len(read_array[j])):
        if read_array[j][i] == 0:
            read_array[j][i] = 2
            first_ribo_array[j] += 1
        else:
            first_ribo_array[j] += 1
            break
    if (j/100).is_integer():
        print("Trimming reads: ", j, "/", len(read_array))
print("Trimming reads:", len(read_array), "/", len(read_array))

with np.printoptions(threshold=np.inf):
    print()
    print('First Ribosome Positions:')
    print()
    print(first_ribo_array)

sum_reads = 0

for i in range(len(first_ribo_array)):
    sum_reads += first_ribo_array[i]

avg_FRP = sum_reads/len(first_ribo_array)
sd_FRP = np.std(first_ribo_array)

print()
print("Mean FRP: ", avg_FRP)
print("FRP SD: ", sd_FRP)

#Prediction model

plt.hist(first_ribo_array, bins=len(np.unique(first_ribo_array)), range=[1, transcript_length])
plt.show()

print()
print("Done.")
