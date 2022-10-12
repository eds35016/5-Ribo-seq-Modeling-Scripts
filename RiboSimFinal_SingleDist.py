'''
Created on Oct 3, 2022

@author: daily
'''
#Requires numpy, matplotlib, and sklearn packages
import sys
import random
import numpy as np
from numpy.random import default_rng
import matplotlib.pyplot as plt
from sklearn.linear_model import PoissonRegressor

if len(sys.argv) != 4:
    print('Usage: RiboSimFinal_SingleDist.py read_count transcript_length initiation_rate_avg')
    exit()

#Set parameters
read_count = int(sys.argv[1])
#Transcript length counted in amino acids
transcript_length = int(sys.argv[2])
initiation_rate_avg = int(sys.argv[3])

rng = default_rng()

read_array = np.zeros(shape=(read_count, transcript_length))
first_ribo_array = np.zeros(read_count)

#Fill read array with transcript arrays consisting of ribosome positions (0 = No ribosome in that position, 1 = Ribosome present)
#Final ribosome placements for each read determined through shifting a random number of times after the first passthrough
for j in range(len(read_array)):
    rand_int = random.randint(1, transcript_length)
    initiation_step = -1
    for i in range(len(read_array[j]) + rand_int):
        #Calculate new initiation rate, variance uses poisson distribution
        while initiation_step < 1:
            initiation_step = rng.poisson(initiation_rate_avg, 1)[0]
        read_array[j] = np.roll(read_array[j], 1)
        read_array[j][:1] = 0
        if (i/initiation_step).is_integer():
            read_array[j][:1] = 1
            initiation_step = -1
    if (j/100).is_integer():
        print("Creating reads: ", j, "/", len(read_array))
print("Creating reads:", len(read_array), "/", len(read_array))

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
    print()
    
count_array = np.unique(first_ribo_array, return_counts=True)

#Generate poisson regression model for the dataset
clf = PoissonRegressor()
mod_array = count_array[0].reshape(-1,1)
clf.fit(mod_array, count_array[1])

yvals = clf.predict(mod_array)
coef = clf.coef_
intercept = clf.intercept_

#Plot model
plt.plot(count_array[0], count_array[1], '.', count_array[0], yvals, '-')
title_string = "Model 1: IR Avg of " + str(initiation_rate_avg) + " bp"
plt.title(title_string)
plt.ylim(bottom=0)
plt.show()

print(coef)
print(intercept)
print()

print("Done.")
