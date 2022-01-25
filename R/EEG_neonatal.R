# file EEG_neonatal.R
# copyright (C) 2021-2022 Artur Gramacki and Jaroslaw Gramacki
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


# These two libraries must by first installed.
# Install from https://bioconductor.org
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rhdf5")

library(edf)
library(rhdf5)

source("EEG_neonatal_FUNS.R")

# We assume the following directory structure is created under the 'dir':
#   |---annotations
#   |---edf
#   |---Python
#   |---R
#   |---working
#     |---acc_loss
#     |---best_models
#     |---hists
#     |---inputs
#     |---logs
#     |---results    
#     |---ROC

# Set dir to indicating the appropriate directory for your local computer. 
# If you are working with a directory structure as shown above, do not change this variable. 
dir = "../"

# Symbols of human experts
we <- c( "A", "B", "C")

# Annotations file names, as downloaded from https://zenodo.org/record/4940267
ann.f <- c("annotations_2017_A_fixed.csv",
          "annotations_2017_B.csv",
          "annotations_2017_C.csv")

# infant IDs which have seizures 
s.IDs <- c(1,4,5,7,9,11,13,14,15,16,17,19,20,21,22,25,31,34,36,38,39,40,41,44,47,50,51,52,62,63,66,67,69,71,73,75,76,77,78,79)
# infant IDs which are seizure free
ns.IDs <- c(3,10,18,27,28,29,30,32,35,37,42,45,48,49,53,55,57,58,59,60,70,72)

timestamp()

for (i in 1:3) { 
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 1, chunks = 1, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 1, chunks = 2, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 1, chunks = 5, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 1, chunks = 10, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 1, chunks = 20, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 2, chunks = 1, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 2, chunks = 2, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 2, chunks = 5, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 2, chunks = 10, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 2, chunks = 20, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 5, chunks = 1, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 5, chunks = 2, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 5, chunks = 5, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 5, chunks = 10, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 5, chunks = 20, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)

  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 10, chunks = 1, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 10, chunks = 2, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 10, chunks = 5, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 10, chunks = 10, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 10, chunks = 20, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)

  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 20, chunks = 1, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 20, chunks = 2, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 20, chunks = 5, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 20, chunks = 10, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 20, chunks = 20, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)

  # We set chunks = 10000 and this way we are sure that the maximum possible set of chunks will be selected. 
  # Our dataset simply doesn't have seizures as long as 10,000 seconds.  
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 1, chunks = 10000, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 2, chunks = 10000, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 5, chunks = 10000, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 10, chunks = 10000, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
  out <-  generate_samples(which.expert = we[i], annotations_file = ann.f[i], seizure.IDs = s.IDs, non.seizure.IDs = ns.IDs, window = 20, chunks = 10000, down.sampling.factor = 4, preprocessing = FALSE, dir = dir, random = FALSE, write.txt.files = FALSE, write.hdf5.files = TRUE)
}

timestamp()

# For testing purposes we generate one file in txt format. 
# However, in Python we use data saved in HDF5 format which is a binary format which is not directly readable.
out <-  generate_samples(
  which.expert = "A",
  annotations_file = "annotations_2017_A_fixed.csv",
  seizure.IDs = s.IDs,
  non.seizure.IDs = ns.IDs,
  window = 1,
  chunks = 1,
  down.sampling.factor = 4,
  preprocessing = FALSE,
  dir = dir,
  random = FALSE,
  write.txt.files = TRUE,
  write.hdf5.files = TRUE
)

# For testing purposes 
ann <-
  read.csv(
    paste(dir, "annotations/", ann.f[1], sep = ""),
    sep = ",",
    header = TRUE,
    stringsAsFactors = F,
    check.names = FALSE,
    encoding = 'UTF-8'
  )
head(ann)

for (i in 1:79) {
  out <- select_seizure_chunks(ann, 256, i)
  if (!is.na(out$total.seizures))
    cat(out$seizures$seizure_duration, "", sep = ";")
}

