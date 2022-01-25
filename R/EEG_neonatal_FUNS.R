# file EEG_neonatal_FUNS.R
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

select_seizure_chunks <- function(data, f, k) {
  # data - seizure annotations file
  # f - base frequency 
  # k - which patient (1:79)
  
  aa <- 0
  seizures <- data.frame()
  sec.1 <- which(data[,k] == 1)
  # only for patients with at least one seizure 
  if (length(sec.1) != 0 ) {
    # https://stackoverflow.com/questions/23095415/how-to-find-if-the-numbers-are-continuous-in-r
    # s: each list item is the starting and ending second of an epileptic seizure
    s <- unname(tapply(sec.1, cumsum(c(1, diff(sec.1)) != 1), range))
    for (i in 1:length(s)) {
      (s.duration.in.secs <- s[[i]][2] - s[[i]][1] + 1)
      #(s.samples <- (s[[i]][1] - 1) * f + 1)  
      (s.samples <- c( 
        (s[[i]][1] - 1) * f + 1, 
        (s[[i]][1] - 1)* f  + (s.duration.in.secs * f)
      ))
      seizures[i, 1] <- k
      seizures[i, 2] <- s.duration.in.secs
      seizures[i, 3] <- s[[i]][1]
      seizures[i, 4] <- s[[i]][2]
      seizures[i, 5] <- s.samples[1]
      seizures[i, 6] <- s.samples[2]
    }  
    colnames(seizures) <- c("patient", "seizure_duration", "from_sec", "to_sec", "from_sample", "to_sample")
    list(
      seizures = seizures, 
      total.seizures = length(s)
    )
  } else {
    list(seizures = NA, total.seizures = NA)
  }  
}  

generate_montage <- function(matrix) {
  # The code based on read_data_montage.m  
  # https://github.com/ktapani/Neonatal_Seizure_Detection
  str <- data.frame()
  str[1,1] = 'Fp2'; str[1,2] = 'F4';    # Fp2-F4
  str[2,1] = 'F4'; str[2,2] = 'C4';     # F4-C4
  str[3,1] = 'C4'; str[3,2] = 'P4';     # C4-P4
  str[4,1] = 'P4'; str[4,2] = 'O2';     # P4-O2
  str[5,1] = 'Fp1'; str[5,2] = 'F3';    # Fp1-F3
  str[6,1] = 'F3'; str[6,2] = 'C3';     # F3-C3
  str[7,1] = 'C3'; str[7,2] = 'P3';     # C3-P3
  str[8,1] = 'P3'; str[8,2] = 'O1';     # P3-O1
  str[9,1] = 'Fp2'; str[9,2] = 'F8';    # Fp2-F8
  str[10,1] = 'F8'; str[10,2] = 'T4';   # F8-T4
  str[11,1] = 'T4'; str[11,2] = 'T6';   # T4-T6
  str[12,1] = 'T6'; str[12,2] = 'O2';   # T6-O2
  str[13,1] = 'Fp1'; str[13,2] = 'F7';  # Fp1-F7
  str[14,1] = 'F7'; str[14,2] = 'T3';   # F7-T3
  str[15,1] = 'T3'; str[15,2] = 'T5';   # T3-T5
  str[16,1] = 'T5'; str[16,2] = 'O1';   # T5-O1
  str[17,1] = 'Fz'; str[17,2] = 'Cz';   # Fz-Cz
  str[18,1] = 'Cz';  str[18,2] ='Pz';   # Cz-Pz
  
  (label <- colnames(matrix))
  mtx.mont <- matrix(NA, nrow = nrow(matrix), ncol = 20)
  
  for (jj in 1:18) {
    ref1 = rep(0, 21)
    ref2 = rep(0, 21)
    for (ii in 1:21) {
      ref1[ii] <- as.numeric(grepl(str[jj,1], label[ii]))
      ref2[ii] <- as.numeric(grepl(str[jj,2], label[ii]))
    }
    qq1 = which(ref1 == 1)[1]
    qq2 = which(ref2 == 1)[1]
    mtx.mont[,jj] = matrix[, qq1] - matrix[, qq2] 
  }
  df.mont <- as.data.frame(mtx.mont)
  df.mont[,19] = matrix[,22]
  for (kk in 1:18) {
    colnames(df.mont)[kk] <- paste(str[kk,1], "-", str[kk,2], sep ="")
  }
  colnames(df.mont)[19] <- c("t")
  df.mont
}

generate_samples <- function(which.expert,
                             annotations_file,
                             seizure.IDs,
                             non.seizure.IDs,
                             window,
                             chunks,
                             down.sampling.factor = 4,
                             preprocessing = FALSE,
                             dir,
                             random = TRUE,
                             write.txt.files = FALSE,
                             write.hdf5.files = TRUE) {
  
  ####################################################################
  # which.expert - "A", "B" or "C" (or any other symbol)
  # annotations_file - file name where annotations are saved
  # seizure.IDs - infant IDs which have seizures. 
  #                We assume that EDF file names are eeg1.edf, eeg2.edf, eeg3.edf etc.
  # non.seizure.IDs - infant IDs which are seizure free. 
  #                   We assume that EDF file names are eeg1.edf, eeg2.edf, eeg3.edf etc.
  # window - window size in seconds
  # chunks - number of window's chunks  
  # down.sampling.factor - down sampling factor of the original edf file   
  #                        (must be divisible by "f", that is, 1,2,4,8 etc. 
  #                        In practice no more than 4, 1 - no down-sampling)
  # preprocessing - see preprocess.m in https://github.com/ktapani/Neonatal_Seizure_Detection
  # dir - root dir, see the directory structure given in the paper
  # random - if TRUE generate random non-seizure chunks from non-seizure EDF files 
  # write.txt.files - if TRUE safe the seizure and nos-seizure final files as TXT files
  #                   (we generate TXT files for illustrative purposes only. As for the content, 
  #                   they are fully compatible with HDF5 binary files)     
  # write.hdf5.files - if TRUE safe the final file in the HDF5 format
  #
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
  ####################################################################
  
  # preprocessing, see preprocess.m
  Bn <- c(0.963823074384152,  -0.649404428450631,   0.963823074384152)
  An <- c(1.000000000000000,  -0.649404428450631,   0.927646148768304)
  NumL <- c(7.011441810666725e-05,   5.609153448533380e-04,   1.963203706986683e-03,   3.926407413973366e-03,   4.908009267466706e-03,   3.926407413973366e-03,   1.963203706986683e-03,   5.609153448533380e-04,   7.011441810666725e-05)
  DenL <- c(1.000000000000000e+00,  -4.233841605746147e+00,   8.384782529806564e+00,  -9.942427619247514e+00,   7.647300224702220e+00,  -3.881875628015500e+00,   1.264260508787803e+00,  -2.407066741662927e-01,   2.045755491417291e-02)
  NumH <- c(0.990425738509653,  -3.961702954038612,   5.942554431057918,  -3.961702954038612,   0.990425738509653)
  DenH <- c(1.000000000000000,  -3.980759293223619,   5.942462766030441,  -3.942646613397998,   0.980943143502392)

  # read annotations file  
  ann <- read.csv(paste(dir, "annotations/", annotations_file, sep=""), sep=",", header = TRUE, stringsAsFactors = F, check.names=FALSE, encoding = 'UTF-8')
  
    # Read EDF params from the first file.
  # n.sigs - number of recorded channels in the EDF file (+1 because then we add "t" column) 
  # f.edf - sampling frequency of the original EDF files
  filename <- paste(dir, "edf/eeg", 1, ".edf", sep = "")
  edf <- read.edf(filename = filename, read.annotations = FALSE, header.only = FALSE)
  n.sigs <- edf[["header.global"]][["n.signals"]] + 1
  f.edf <- edf[["header.signal"]][[1]][["n.samples"]]
  (f <- f.edf /down.sampling.factor)
  
  cat("---------------------------------------------------------------------------", "\n")
  cat("annotations file name:     ", annotations_file, "\n", sep = "")
  cat("seizure patients:         ", seizure.IDs, "\n", sep = " ")
  cat("non-seizure patients:     ", non.seizure.IDs, "\n", sep = " ")    
  if (write.hdf5.files) cat("hdf5 file being generated: ", "expert_", which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.hdf5", "\n", sep = "")
  if (write.txt.files) cat("txt file being generated:  ", "expert_", which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", "\n", sep = "")
  cat("---------------------------------------------------------------------------", "\n\n")
  
  SEIZURE <- data.frame()
  NON.SEIZURE <- data.frame() 
  FINAL <- data.frame()
  
  unlink(paste(dir, "working/inputs/signal_lengths_in_mins.txt", sep = ""))
  unlink(paste(dir, "working/inputs/seizures_",     which.expert, "_", f, "Hz.txt", sep = ""))
  unlink(paste(dir, "working/inputs/seizures_",     which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", sep = ""))
  unlink(paste(dir, "working/inputs/non_seizures_", which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", sep = ""))
  
  unlink(paste(dir, "working/inputs/SEIZURE_expert_",      which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", sep = ""))
  unlink(paste(dir, "working/inputs/NON.SEIZURE_expert_",  which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", sep = ""))
  unlink(paste(dir, "working/inputs/expert_",              which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", sep = ""))

  # Counter of seizure chunks. 
  S <- 0
  
  ####################################################################
  # First, we analyze patients with at least one annotated seizure. 
  ####################################################################
  for (i in seizure.IDs) {
    filename <- paste(dir, "edf/eeg", i, ".edf", sep = "")
    cat(filename, sep = "")
    edf <- read.edf(filename = filename, read.annotations = FALSE, header.only = FALSE)
    
    # Calculate the number of samples. Take the first signal because each has the same number of samples. 
    (len <- length(edf$signal[[1]]$data))
    
    write.table(
      x = t(c(i, round(len / f.edf / 60, 2))),
      file = paste(dir, "working/inputs/signal_lengths_in_mins.txt", sep = ""),
      append = TRUE,
      col.names = FALSE,
      row.names = FALSE
    )
    
    edf.mtx <- matrix(NaN, nrow = len, ncol = n.sigs)
    for (n in 1:(n.sigs - 1)) {
      edf.mtx[, n] <- edf[["signal"]][[n]][["data"]]
    }
    head(edf.mtx)
    
    sig.names <- NaN
    # Read the names of the signals.
    for (s in 1:(n.sigs - 1)) {
      sig.names[s] <- edf$header.signal[[s]]$label
    }
    colnames(edf.mtx) <- c(sig.names, "t")
    head(edf.mtx)
    
    # Add time stamps to the last column. Take data from any channel, the same everywhere. 
    edf.mtx[, n.sigs] <- edf[["signal"]][[1]][["t"]]
    head(edf.mtx)
    
    # down-sampling
    edf.mtx <- edf.mtx[seq(1, nrow(edf.mtx), down.sampling.factor),]    
    
    # montage
    edf.mtx.m <- generate_montage(edf.mtx)
    # last column - reserved for class label
    # penultimate column - timestamps
    head(edf.mtx.m)
    cc <- ncol(edf.mtx.m)
    
    # preprocessing, see preprocess.m in https://github.com/ktapani/Neonatal_Seizure_Detection  
    if (preprocessing) {
      edf.mtx.p <- matrix(NA, nrow = nrow(edf.mtx.m), ncol = ncol(edf.mtx.m))
      colnames(edf.mtx.p) <- colnames(edf.mtx.m)
      #for (m in 1:18) {
      for (m in 1:(cc - 2)) {
        edf.mtx.p[, m] = signal::filter(Bn, An, edf.mtx.m[, m]); # 50Hz notch filter
        edf.mtx.p[, m] = signal::filter(NumL, DenL, edf.mtx.p[, m]); # Low pass IIR Butterworth, cutoff 30Hz
        edf.mtx.p[, m] = signal::filter(NumH, DenH, edf.mtx.p[, m]); # High pass IIR Butterwoth, cutoff 0.3Hz
      }
      # edf.mtx.p[, 19] <- edf.mtx.m[, 19]
      edf.mtx.p[, cc - 1] <- edf.mtx.m[, cc - 1]
    } else {
      edf.mtx.p <- edf.mtx.m
    }
    head(edf.mtx.p)
    
    out <- select_seizure_chunks(data = ann, f = f, k = i)
    out$seizures
    out$total.seizures
    cat(" (total number of seizures annotated: ", out$total.seizures, ")\n", sep = "")
    
    if (i == seizure.IDs[1]) col.names = TRUE else col.names = FALSE
    write.table(
      x = out$seizures,
      file = paste(dir, "working/inputs/seizures_", which.expert, "_", f, "Hz.txt", sep = ""),
      append = TRUE,
      col.names = col.names,
      row.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
    seizure.mtx <- matrix(NaN, nrow = 0, ncol = ncol(edf.mtx.p)) 
    
    m <- 0
    out2 <- data.frame()
    for (j in 1:out$total.seizures) {
      (r <- out$seizures[j, 2] %/% window)
      if (r > 0) {
        if (r > chunks) { # so as not to select more chunks than possible 
          fr <- chunks
        } else {
          fr <- r
        }
        for (k in 1:fr) {
          S <- S + 1
          m <- m + 1
          (from <- out$seizures[j, 5] + (k - 1) * f * window) 
          (to <- out$seizures[j ,5] + (k * window * f) - 1)  
          out2[m, 1] <- out$seizures[j, 1]
          out2[m, 2] <- out$seizures[j, 2]
          out2[m, 3] <- from
          out2[m, 4] <- to
        }  
      } else {
        # do nothing
      }
    }
    colnames(out2) <- c("patient", "seizure_duration", "from_sample",	"to_sample")
    
    if (i == seizure.IDs[1]) col.names = TRUE else col.names = FALSE
    write.table(
      x = out2,
      file = paste(dir, "working/inputs/seizures_", which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", sep = ""), 
      append = TRUE,
      col.names = col.names,
      row.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
    
    if (nrow(out2) > 0) {
      for (p in 1:nrow(out2)) {
        (from <- out2[p, 3])
        (to <- out2[p, 4])
        temp <- edf.mtx.p[from:to,]
        seizure.mtx <- rbind(seizure.mtx, temp)
      }
      head(seizure.mtx)
      nrow(seizure.mtx)
      
      seizure.mtx[, cc] <- rep(1, nrow(seizure.mtx))
      colnames(seizure.mtx)[cc] <- "seizure"
      
      SEIZURE <- rbind(SEIZURE, seizure.mtx)
    }
    
  } ### for (i in seizure.IDs)
  
  if (write.txt.files) {
    write.table(
      x = round(SEIZURE, 4),
      file = paste(dir, "working/inputs/SEIZURE_", which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", sep = ""), 
      dec = ".",
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE
    )
  }
  
  ####################################################################  
  # Second, we analyze patients with NO ONE annotated seizure. 
  # We select RANDOMLY as many chunks as we have selected from the seizured files. 
  ####################################################################
  # Set how many chunks to take from each non-seizured file, so that there are about the same number of chunks.
  (nn <- ceiling(S / length(non.seizure.IDs)))

  for (q in non.seizure.IDs) {
    filename <- paste(dir, "/edf/eeg", q, ".edf", sep = "")
    cat(filename, "\n", sep = "")
    edf <- read.edf(filename = filename, read.annotations = FALSE, header.only = FALSE)
    
    # Calculate the number of samples. Take the first signal because each has the same number of samples. 
    (len <- length(edf$signal[[1]]$data))
    
    write.table(
      x = t(c(q, round(len / f.edf / 60, 2))),
      file = paste(dir, "working/inputs/signal_lengths_in_mins.txt", sep=""),
      append = TRUE,
      col.names = FALSE,
      row.names = FALSE
    )

    edf.mtx <- matrix(NaN, nrow = len, ncol = n.sigs)
    for (n in 1:(n.sigs - 1)) {
      edf.mtx[, n] <- edf[["signal"]][[n]][["data"]]
    }
    head(edf.mtx)

    sig.names <- NaN
    # Read the names of the signals.
    for (s in 1:(n.sigs - 1)) {
      sig.names[s] <- edf$header.signal[[s]]$label
    }
    colnames(edf.mtx) <- c(sig.names, "t")
    head(edf.mtx)

    # Add time stamps to the last column. Take data from any channel, the same everywhere. 
    edf.mtx[, n.sigs] <- edf[["signal"]][[1]][["t"]]
    head(edf.mtx)

    # down-sampling 
    edf.mtx <- edf.mtx[seq(1, nrow(edf.mtx), down.sampling.factor),]

    # montage
    edf.mtx.m <- generate_montage(edf.mtx)
    head(edf.mtx.m)
    (cc <- ncol(edf.mtx.m))

    # preprocessing  
    if (preprocessing) {
      edf.mtx.p <- matrix(NA, nrow = nrow(edf.mtx.m), ncol = ncol(edf.mtx.m))
      colnames(edf.mtx.p) <- colnames(edf.mtx.m)
      #for (m in 1:18) {
      for (m in 1:(cc - 2)) {
        edf.mtx.p[, m] = signal::filter(Bn, An, edf.mtx.m[, m]); # 50Hz notch filter
        edf.mtx.p[, m] = signal::filter(NumL, DenL, edf.mtx.p[, m]); # Low pass IIR Butterworth, cutoff 30Hz
        edf.mtx.p[, m] = signal::filter(NumH, DenH, edf.mtx.p[, m]); # High pass IIR Butterwoth, cutoff 0.3Hz
      }
      #edf.mtx.p[, 19] <- edf.mtx.m[, 19]
      edf.mtx.p[, cc - 1] <- edf.mtx.m[, cc - 1]
    } else {
      edf.mtx.p <- edf.mtx.m
    }
    head(edf.mtx.p)
    
    # allocate a matrix
    non.seizure.mtx <- matrix(NaN, nrow = nn * f * window, ncol = ncol(edf.mtx.p)) 
    colnames(non.seizure.mtx) <- colnames(edf.mtx.m)
    
    (nrows.edf <- nrow(edf.mtx.p))
    mm <- matrix(NaN, nrow = nn, ncol = 4)
    for (k in 1:nn) {
      if (random) {
        r <- sample(1:(nrows.edf - f * window), 1) # randomly 
      } else {
        set.seed(42 + k)
        r <- sample(1:(nrows.edf - f * window), 1) # seed
      }  
      temp <- as.matrix(edf.mtx.p[r:(r + f * window - 1), ])
      non.seizure.mtx[((k - 1) * f * window + 1):(k * f * window),] <- temp
      
      mm[k,] = t(c(q, k, r, (r + f * window - 1)))
    } # for (k in 1:nn)
    head(non.seizure.mtx)

    NON.SEIZURE <- rbind(NON.SEIZURE, non.seizure.mtx)
    
    
    colnames(mm) <- c("patient", "sequence_number", "from_sample",	"to_sample")
    
    if (q == non.seizure.IDs[1]) col.names = TRUE else col.names = FALSE
    write.table(
      x = mm,
      file = paste(dir, "working/inputs/non_seizures_", which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", sep = ""), 
      append = TRUE,
      sep = "\t",
      col.names = col.names,
      row.names = FALSE
    )
    
    colnames(NON.SEIZURE) <- colnames(edf.mtx.m)
    head(NON.SEIZURE)
    
  } # for (ii in non.seizure.IDs)
  
  # In the last column 0 - no seizure 
  NON.SEIZURE[, cc] <- rep(0, nrow(NON.SEIZURE))
  colnames(NON.SEIZURE)[cc] <- "seizure"
  head(NON.SEIZURE)
  nrow(NON.SEIZURE)
  
  if (write.txt.files) {
    write.table(
      x = round(NON.SEIZURE, 4),
      file = paste(dir, "working/inputs/NON.SEIZURE_", which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", sep = ""), 
      dec = ".",
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE
    )
  }
  
  ####################################################################
  # Combining the partial results into the final one file 
  ####################################################################
  FINAL <- rbind(SEIZURE, NON.SEIZURE)
  FINAL.mtx <- as.matrix(FINAL)
  
  if (write.hdf5.files) {
    fname <- paste(dir, "working/inputs/expert_", which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.hdf5", sep ="")
    unlink(fname)
    h5write(obj = FINAL.mtx, file = fname, name = "FINAL.mtx", createnewfile = TRUE)
    h5closeAll()
  }  
    
  if (write.txt.files) {
    write.table(
      x = round(FINAL, 4),
      file = paste(dir, "working/inputs/expert_", which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", sep = ""),
      sep = "\t", 
      append = FALSE,
      col.names = FALSE,
      row.names = FALSE
    )
  } 
  
  list(
    SEIZURE = SEIZURE,
    NON.SEIZURE =NON.SEIZURE,
    FINAL = FINAL
  )
}



