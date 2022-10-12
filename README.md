# 5-Ribo-seq-Modeling-Scripts
Scripts to generate theoretical 5' RiboSeq data and theoretical RiboSeq data in silico

Requirements:

Python 3+
NumPy Package
MatPlotLib Package
Sklearn Package

RiboSim.py: Original 5' Ribo-seq simulated data generator using a normal distribution for ribosome initiation rate variance. Generates a histogram of first ribosome positions.

  Params: 
  
      read_count: number of transcript reads to simulate
      
      transcript_length: length of the original transcript
      
      initiation_rate_avg: rate at which new ribosomes are added to the transcript on average (ie. new ribosome added to the transcript after the previous ribosome has moved 150 codons from the start on average)
      
      initiation_rate_sd: Standard deviation of the normal distribution defining the variance in initiation rate
      
RiboSimPoisson.py: Same as RiboSim.py, but uses a poisson distribution to generate the ribosome initiation rate variance. This means that the standard deviation is equivalent to the square root of the initiation rate average.

  Params: 
  
      read_count: number of transcript reads to simulate
      
      transcript_length: length of the original transcript
      
      initiation_rate_avg: rate at which new ribosomes are added to the transcript on average

RiboSimPoisson_allpos.py: Same as RiboSimPoisson.py, but outputs a histogram of all ribosome positions on each transcript, not just the first one.

  Params: 
  
      read_count: number of transcript reads to simulate
      
      transcript_length: length of the original transcript
      
      initiation_rate_avg: rate at which new ribosomes are added to the transcript on average

RiboSimFinal.py: Generates three separate sets of simulated data based on the method from RiboSimPoisson.py. The initiation_rate_avg of the three datasets will differ by the amount set for model_deviation. Then outputs individual dotplots and a poisson regression line for each individual dataset and one with all models combined for comparison.

  Params: 

      read_count: number of transcript reads to simulate
      
      transcript_length: length of the original transcript
      
      initiation_rate_avg: rate at which new ribosomes are added to the transcript on average for the middle model (model 2)
      
      model_deviation: the amount of deviation in initiation rate between each dataset

RiboSimFinal_SingleDist.py: Same as RiboSimFinal.py, except it generates only one dataset instead of three

  Params: 

      read_count: number of transcript reads to simulate
      
      transcript_length: length of the original transcript
      
      initiation_rate_avg: rate at which new ribosomes are added to the transcript on average for the middle model
