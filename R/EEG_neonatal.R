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
#       |---test_files
#   |---working
#     |---acc_loss
#     |---best_models
#     |---hists
#     |---inputs
#     |---logs
#     |---results    
#     |---ROC
#     |---waveforms

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

# Write signal lengths (in min, sec) to a txt file
read_signal_lengths(1:79, dir)

# The 3 below `generate_eeg_waveforms` function calls will generate 1379 plots representing
# all the seizures annotated by 3 experts, see Table 7. 
# The lengths of the waveforms were arbitrarily set at 10 seconds. However, the
# user can generate waveforms with different lengths controled by the `window` parameter.

# infant IDs where 1 or 2 experts annotated seizures, expert A
s12A.IDs <- c(2,8,23,33,54,68)
# infant IDs where 1 or 2 experts annotated seizures, expert B
s12B.IDs <- c(8,24,64,68,74)
# infant IDs where 1 or 2 experts annotated seizures, expert C
s12C.IDs <- c(6,12,23,26,33,43,46,54,56,61,64,65,74)

time_elapsed <- system.time({
p.ID <- c(s.IDs, s12A.IDs)
generate_eeg_waveforms(
  mode = "seizure",
  which.expert = we[1],
  annotations.file = ann.f[1],
  patient.IDs = p.ID,
  window = 10,
  down.sampling.factor = 4,
  non.seizures.chunks = 10,
  preprocessing = FALSE,
  random = FALSE,
  format = "png",
  save.to.file = TRUE,
  dir = dir
)

p.ID <- c(s.IDs, s12B.IDs)
generate_eeg_waveforms(
  mode = "seizure",
  which.expert = we[2],
  annotations.file = ann.f[2],
  patient.IDs = p.ID,
  window = 10,
  down.sampling.factor = 4,
  non.seizures.chunks = 10,
  preprocessing = FALSE,
  random = FALSE,
  format = "png",
  save.to.file = TRUE,
  dir = dir
) 

p.ID <- c(s.IDs, s12C.IDs)
generate_eeg_waveforms(
  mode = "seizure",
  which.expert = we[3],
  annotations.file = ann.f[3],
  patient.IDs = p.ID,
  window = 10,
  down.sampling.factor = 4,
  non.seizures.chunks = 10,
  preprocessing = FALSE,
  random = FALSE,
  format = "png",
  save.to.file = TRUE,
  dir = dir
) 
})
min <- time_elapsed[3] %/% 60
sec <- round(time_elapsed[3] %% 60, 0)
cat("Execution time: ", min,  "min ", sec, "sec", sep = "")


####################################################################
# Preparation of data that will be fed to the CNN network.
# This data was not used during training of the CNN network.
# This data is used to demonstrate how to use a pre-trained model to 
# classify data other than that used to build the CNN network,
# see Table 7 in the paper.
# We arbitrarily decided to set window = 10 and chunks = 20.
# To avoid a conflict with the filenames (see `working\inputs` directory), 
# the expert was marked `CC` instead of `C`.
####################################################################
out <-
  generate_samples(
    which.expert = "CC",
    annotations_file = ann.f[3],
    seizure.IDs = s12C.IDs,
    non.seizure.IDs = ns.IDs,
    window = 10,
    chunks = 20,
    down.sampling.factor = 4,
    preprocessing = FALSE,
    dir = dir,
    random = FALSE,
    write.txt.files = FALSE,
    write.hdf5.files = TRUE
  )

####################################################################
# Example 1
# 25 plots will be generated 
# (expert A annotated 25 seizures in neonate # 1)
####################################################################
generate_eeg_waveforms(
  mode = "seizure",
  which.expert = "A",
  annotations.file = "annotations_2017_A_fixed.csv",
  patient.IDs = 1,
  window = 10,
  down.sampling.factor = 4,
  non.seizures.chunks = 10,
  preprocessing = FALSE,
  random = TRUE,
  format = "png",
  save.to.file = TRUE,
  dir = dir
)

####################################################################
# Example 2
# For testing purposes we generate one file in txt format. 
# However, in our codes in Python we use data in the binary HDF5 format 
# which is not directly readable. 
# The generated data will be stored in the working\inputs directory
####################################################################
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

####################################################################
# Example 3
# For testing purposes we read the first annotations file 
# and display a short summary of annotations.
####################################################################
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

total.seizures <- 0
seizure.durations <- c()
for (i in s.IDs) {
  out <- select_seizure_chunks(ann, 256, i)
  if (!is.na(out$total.seizures)) {
    cat(out$seizures$seizure_duration, "", sep = ";")
    seizure.durations <- c(seizure.durations, out$seizures$seizure_duration)
    total.seizures <- total.seizures + out$total.seizures
  }
}
total.seizures
seizure.durations
sort(seizure.durations)
range(seizure.durations)

####################################################################
# Example 4
# Ranges for regenarating EEG waveforms in Figures 4, 5, 6
####################################################################
sec.range <- c(104, 121)
sec.range <- c(6847, 6863)

# Ranges in Figure 4
sec.range <- c(44, 49)
sec.range <- c(555, 560)
sec.range <- c(2000, 2005)
sec.range <- c(3001, 3006)
sec.range <- c(3560, 3565)

# Ranges in Figure 5
sec.range <- c(885, 889)
sec.range <- c(1222, 1226)
sec.range <- c(3000, 3004)
sec.range <- c(3335, 3339)

# Ranges in Figure 6
sec.range <- c(55, 56)
sec.range <- c(89, 90)
sec.range <- c(560, 561)
sec.range <- c(1334, 1335)
sec.range <- c(2222, 2223)
sec.range <- c(3017, 3018)
sec.range <- c(3331, 3332)
sec.range <- c(3876, 3877)
sec.range <- c(3991, 3992)
sec.range <- c(4002, 4003)

# one selected channel
# Figure 4a, first seizure
dev.off()
out <- generate_selected_waveforms(
  patient = 1,
  down.sampling.factor = 4,
  channel = 6,
  sec.range = c(104, 121),
  preprocessing = FALSE,
  dir
)
plot(out$points, type = "l", xlab = "", ylab = out$channel.names, yaxt = "n")
write.table(out, "clipboard", row.names = FALSE, col.names = FALSE)

# Figure 4a, second seizure
dev.off()
out <- generate_selected_waveforms(
  patient = 1,
  down.sampling.factor = 4,
  channel = 6,
  sec.range = c(6847, 6863),
  preprocessing = FALSE,
  dir
)
plot(out$points, type = "l", xlab = "", ylab = out$channel.names, yaxt = "n")

# Figure 4b, first seizure
dev.off()
out <- generate_selected_waveforms(
  patient = 10,
  down.sampling.factor = 4,
  channel = 6,
  sec.range = c(44, 49),
  preprocessing = FALSE,
  dir
)
plot(out$points, type = "l", xlab = "", ylab = out$channel.names, yaxt = "n")

# etc. for other plots in Figures 4, 5, 6

####################################################################
# Example 5
# Plot all waveforms for a selected patient and a selected range
####################################################################
# all channels
dev.off()
out <- generate_selected_waveforms(
  patient = 44,
  down.sampling.factor = 4,
  channel = "all",
  sec.range = c(350, 380),
  preprocessing = TRUE,
  dir = dir
)
n <- length(out$points)
par(mfrow = c(n-1, 1))
par(mai = c(0.0, 1.0, 0.0, 0.4)) # bottom, left, top, right)
for (j in 1:(n-1)) {
  if (j == 1) xaxt <- "s" else xaxt <- "n"
  plot(out$points[, n] + 1, out$points[, j], type = "l", xlab = "", ylab = "", yaxt = "n", xaxt = xaxt) 
  mtext(text = colnames(out$points[j]), las = 1, adj = 0, cex = 0.8, side = 2, outer = FALSE, line = 5)
}

generate_eeg_waveforms(
  mode = "seizure",
  which.expert = "A",
  annotations.file = "annotations_2017_A_fixed.csv",
  patient.IDs = 1,
  window = 10,
  down.sampling.factor = 4,
  non.seizures.chunks = 10,
  preprocessing = FALSE,
  random = TRUE,
  format = "png",
  save.to.file = TRUE,
  dir = dir
)


