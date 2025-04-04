# Testing -----------------------------------------------------------------
library(EGM)
options(wfdb_path = '/mmfs1/home/dseaney2/wfdb/bin')
file_path <- "/mmfs1/projects/cardio_darbar_chi/common/cohorts/ml_genetics/annotations/ludb/files"
leads <- c('i','ii','iii','avr','avl','avf','v1','v2','v3','v4','v5','v6')

signal <- list()
annotations <- list()

symbols <- c()
for (sample in 1:200) {
  for (lead in 1:12) {
    ann <- read_annotation(record = sample,
                           record_dir = file_path,
                           annotator = leads[lead])
    
    # symbols <- unique(c(symbols,ann$type))
    
    # Any ( without N/t/p or ), any N/t/p without ( or ), any ) without N/t/p or (
    open_ind <- which(ann$type == '(')
    wave_ind <- which(ann$type == 'N' | ann$type == 't' | ann$type == 'p')
    close_ind <- which(ann$type == ')')
    
    # Note any ( without wave
    if (sum(!open_ind %in% (wave_ind - 1))) {
      ind <- which(!open_ind %in% (wave_ind - 1))
      print(paste('(    without wave on sample',sample,'lead',lead))
    }
    
    # Note any ( without )
    if (sum(!open_ind %in% (close_ind - 2))) {
      ind <- which(!open_ind %in% (close_ind - 2))
      print(paste('(    without )    on sample',sample,'lead',lead))
    }
    
    
    
    # Note any wave without (
    if (sum(!wave_ind %in% (open_ind + 1))) {
      ind <- which(!wave_ind %in% (open_ind + 1))
      print(paste('wave without (    on sample',sample,'lead',lead))
    }
    
    # Note any wave without )
    if (sum(!wave_ind %in% (close_ind - 1))) {
      ind <- which(!wave_ind %in% (close_ind - 1))
      print(paste('wave without )    on sample',sample,'lead',lead))
    }
    
    
    
    # Note any ) without (
    if (sum(!close_ind %in% (open_ind + 2))) {
      ind <- which(!close_ind %in% (open_ind + 2))
      print(paste(')    without (    on sample',sample,'lead',lead))
    }
    
    # Note any ) without wave
    if (sum(!close_ind %in% (wave_ind + 1))) {
      ind <- which(!close_ind %in% (wave_ind - 1))
      print(paste(')    without wave on sample',sample,'lead',lead))
    }
    
  }
  
  if (!sample %% 10) {
    print(paste(sample))
  }
}

# Load --------------------------------------------------------------------
ecg_list <- list()
leads <- c('i','ii','iii','avr','avl','avf','v1','v2','v3','v4','v5','v6')
record_dir <- 'LUDB/data/'

for (i in 1:200) {
  ecg <- read_wfdb(record = i, record_dir = record_dir)
  ecg_list[[i]] <- ecg
  names(ecg_list)[i] <- paste0('sample_', i)
  
  ecg_list[[i]]$annotation <- list()
  
  for (lead in 1:12) {
    ann <- read_annotation(record = i, record_dir = record_dir, annotator = leads[lead])
    ecg_list[[i]]$annotation[[lead]] <- ann
    names(ecg_list[[i]]$annotation)[lead] <- paste0(leads[lead])
  }
  
  if (!i%%10) {
    print(paste(i))
  }
}

old_list <- ecg_list

# Cleaning ----------------------------------------------------------------

# Problem samples: 7, 8, 34, 95, 104, 111, 116, 198
# sinus, AF, sinus, AF, sinus, sinus, sinus, sinus

# 7:   lead 7-9
# 8:   lead 11
# 34:  lead 7-9
# 95:  lead 2-4, 6-7, 9-12
# 104: lead 2-4, 6
# 111: lead 2-3, 6, 9-12
# 116: lead 7-8, 10
# 198: lead 8
leads <- c('i','ii','iii','avr','avl','avf','v1','v2','v3','v4','v5','v6')

sample <- 7
lead <- 7
a <- read_wfdb(paste0('LUDB/data/',sample), annotator = leads[lead])
plot_func(y=a$signal$v1, color = ann_wfdb_to_complete(a$annotation))
a$annotation

# Cleanup 7 -----------------------------------------------------------------
source('C:/Users/darre/OneDrive/Documents/UICOM Research/ECG Segmentation/ML_functions.R')
setwd('C:/Users/darre/OneDrive/Documents/UICOM Research/')
leads <- c("i","ii","iii","avr","avl","avf","v1","v2","v3","v4","v5","v6")

lead <- 7
sample <- 7

test <- ecg_list[[sample]]$annotation[[lead]]
# test <- read_annotation(record = sample, record_dir = rec_dir, annotator = leads[lead])
# test
# missed_waves[missed_waves$Sample == sample & missed_waves$Lead == leads[lead],]

# matrix_sample = which(all_sinus == 7)
# plot_func(sinus_samples12[matrix_sample,,lead],sinus_annotations12[matrix_sample,,lead])

test$type[c(1,18,27,35,43,60)] <- '('
test$type[c(19)] <- 'N'
test$type[c()] <- ')'

index <- c(979,2351,2859,3360,4355)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

ecg_list[[sample]]$annotation[[lead]] <- test

# corrected <- ann_wfdb_to_complete(test)
# plot_func(sinus_samples12[matrix_sample,,lead],corrected)
# sinus_annotations12[matrix_sample,,lead] <- corrected


lead <- 8
sample <- 7
test <- ecg_list[[sample]]$annotation[[lead]]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(1,10,28,37,45,53,62)] <- '('
test$type[c(2,11,29,54)] <- 'N'
# test$type[c(38,46)] <- ')'

index <- c(2859,3359,4355)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

a <- ann_wfdb_to_complete(test)
ecg_list[[sample]]$annotation[[lead]] <- test

lead <- 9
sample <- 7
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(1,10,19,27,54)] <- '('
test$type[c(2,11,28)] <- 'N'

index <- c(1775,3860)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test


# Cleanup 8 ---------------------------------------------------------------
lead <- 11
sample <- 8

test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c()] <- '('
test$type[c(55)] <- 'N'
test$type[c(56)] <- ')'

index <- c(1106)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

ecg_list[[sample]]$annotation[[lead]] <- test

# Cleanup 34 --------------------------------------------------------------
lead <- 7
sample <- 34
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(7,12,23,28,33)] <- '('
test$type[c()] <- 'N'
test$type[c()] <- ')'

index <- c(1250,1752,2753,3253,3752)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test



lead <- 8
sample <- 34
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(1,6,11,28,33,38)] <- '('
test$type[c(39)] <- 'N'
test$type[c()] <- ')'

index <- c(750,1248,1750,3251,3749)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test


lead <- 9
sample <- 34
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(1,6,18,29)] <- '('
test$type[c(7)] <- 'N'
test$type[c()] <- ')'

index <- c(750,2251,3251)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test

# Sample 95 ---------------------------------------------------------------
sample <- 95
lead <- 2

test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(1,6,11,16)] <- '('
test$type[c()] <- 'N'
test$type[c()] <- ')'

index <- c(828,1327,1824,2329)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test


sample <- 95
lead <- 3

test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(1,6,11,16,21)] <- '('
test$type[c()] <- 'N'
test$type[c()] <- ')'

index <- c(828,1327,1824,2329,2787)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test


sample <- 95
lead <- 4
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(25)] <- '('
test$type[c()] <- 'N'
test$type[c()] <- ')'

index <- c(2811)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test


sample <- 95
lead <- 6
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(1,6,11,16,21)] <- '('
test$type[c()] <- 'N'
test$type[c()] <- ')'

index <- c(828,1327,1824,2329,2811)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test


sample <- 95
lead <- 7

test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(1,6,11,16)] <- '('
test$type[c()] <- 'N'
test$type[c()] <- ')'

index <- c(828,1327,1824,2329)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test


sample <- 95
lead <- 9

test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(1,6,11,16)] <- '('
test$type[c()] <- 'N'
test$type[c()] <- ')'

index <- c(828,1327,1824,2329)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test


sample <- 95
lead <- 10

test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(1,6,11,16)] <- '('
test$type[c()] <- 'N'
test$type[c()] <- ')'

index <- c(828,1327,1824,2329)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test


sample <- 95
lead <- 11

test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(1,6,11,16)] <- '('
test$type[c()] <- 'N'
test$type[c()] <- ')'

index <- c(828,1327,1824,2329)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test


sample <- 95
lead <- 12

test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(1,6,11,16)] <- '('
test$type[c()] <- 'N'
test$type[c()] <- ')'

index <- c(828,1327,1824,2329)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test

# Cleanup 104 -------------------------------------------------------------
lead <- 2
sample <- 104
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))


test$type[c(1,6,11,16,21,26,31,36,41)] <- '('
test$type[c()] <- 'N'
test$type[c()] <- ')'

index <- c(750,1221,1690,2159,2625,3074,3534,3993,4447)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test


lead <- 3
sample <- 104
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))


test$type[c(1,6,11,16,21,26,31,36,41)] <- '('
test$type[c()] <- 'N'
test$type[c()] <- ')'

index <- c(750,1221,1690,2159,2625,3074,3534,3993,4447)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test


lead <- 4
sample <- 104
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(1,6,11,16,21,26,31,36,41)] <- '('
test$type[c()] <- 'N'
test$type[c()] <- ')'

index <- c(750,1208,1669,2137,2602,3052,3514,3976,4430)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test


lead <- 6
sample <- 104
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(1,6,11,16,21,26,31,36,41)] <- '('
test$type[c()] <- 'N'
test$type[c()] <- ')'

index <- c(768,1201,1690,2159,2623,3074,3533,3972,4446)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test

# Cleanup 111 -------------------------------------------------------------

lead <- 2
sample <- 111
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(37)] <- '('
test$type[c(38)] <- 'N'
test$type[c()] <- ')'

# VERIFY THESE- these were previously commented 
# index <- c()
# type <- array('N',length(index))
# time <- paste0("0:0",index/500)
# zeros <- array(0,length(index))
# test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
# test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test


lead <- 3
sample <- 111
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(1,24)] <- '('
test$type[c(25)] <- 'N'
test$type[c()] <- ')'

index <- c(956)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test


lead <- 6
sample <- 111
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(25,37)] <- '('
test$type[c(26,38)] <- 'N'
test$type[c()] <- ')'

# PREVIOUSLY COMMENTED
# index <- c(956)
# type <- array('N',length(index))
# time <- paste0("0:0",index/500)
# zeros <- array(0,length(index))
# test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
# test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test


lead <- 9
sample <- 111
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))


test$type[c(43)] <- '('
test$type[c()] <- 'N'
test$type[c()] <- ')'

index <- c(4459)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

ecg_list[[sample]]$annotation[[lead]] <- test


lead <- 10
sample <- 111
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(25,42)] <- '('
test$type[c(43)] <- 'N'
test$type[c()] <- ')'

index <- c(2958)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test


lead <- 11 
sample <- 111
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(1,7,25,30,35,40)] <- '('
test$type[c(2,8)] <- 'N'
test$type[c()] <- ')'

index <- c(2958,3455,3957,4459)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test


lead <- 12 
sample <- 111
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c(7,24,36)] <- '('
test$type[c(25)] <- 'N'
test$type[c()] <- ')'

index <- c(1464,3960)
type <- array('N',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
ecg_list[[sample]]$annotation[[lead]] <- test

# Cleanup 116 -------------------------------------------------------------
lead <- 7 
sample <- 116
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

# VERFIY
# test$type[c(1,7,25,30,35,40)] <- '('
# test$type[c(2,8)] <- 'N'
# test$type[c()] <- ')'
# VERIFY:
test <- test[-c(8:12,14:18,19:23)]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

ecg_list[[sample]]$annotation[[lead]] <- test


lead <- 8
sample <- 116
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))
#   VERIFY
test <- test[-c(8:12,14:17)]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

ecg_list[[sample]]$annotation[[lead]] <- test


lead <- 10
sample <- 116
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

test$type[c()] <- '('
test$type[c()] <- 'N'
test$type[c(81)] <- ')'

index <- c(3164)
type <- array('(',length(index))
time <- paste0("0:0",index/500)
zeros <- array(0,length(index))
test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
test <- test[order(as.numeric(test$sample)),]

plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

ecg_list[[sample]]$annotation[[lead]] <- test

# Cleanup 198 -------------------------------------------------------------
lead <- 8
sample <- 198
test <- ecg_list[[sample]]$annotation[[lead]]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))


test$type[c(61)] <- '('
test$type[c(62)] <- 'p'
test$type[c()] <- ')'

# index <- c(3164)
# type <- array('(',length(index))
# time <- paste0("0:0",index/500)
# zeros <- array(0,length(index))
# test <- rbind(test,data.frame(time, index, type, zeros, zeros, zeros), use.names = FALSE)
# test <- test[order(as.numeric(test$sample)),]
plot_func(ecg_list[[sample]]$signal[[lead]], ann_wfdb_to_complete(test))

ecg_list[[sample]]$annotation[[lead]] <- test



# Verify ------------------------------------------------------------------
ludb_set <- ecg_list

for (sample in 1:200) {
  for (lead in 1:12) {
    
    ann <- ludb_set[[sample]]$annotation[[lead]]
    
    # Any ( without N/t/p or ), any N/t/p without ( or ), any ) without N/t/p or (
    open_ind <- which(ann$type == '(')
    wave_ind <- which(ann$type == 'N' |
                        ann$type == 't' | ann$type == 'p')
    close_ind <- which(ann$type == ')')
    
    # Note any ( without wave
    if (sum(!open_ind %in% (wave_ind - 1))) {
      ind <- which(!open_ind %in% (wave_ind - 1))
      print(paste('(    without wave on sample', sample, 'lead', lead))
    }
    
    # Note any ( without )
    if (sum(!open_ind %in% (close_ind - 2))) {
      ind <- which(!open_ind %in% (close_ind - 2))
      print(paste('(    without )    on sample', sample, 'lead', lead))
    }
    
    
    
    # Note any wave without (
    if (sum(!wave_ind %in% (open_ind + 1))) {
      ind <- which(!wave_ind %in% (open_ind + 1))
      print(paste('wave without (    on sample', sample, 'lead', lead))
    }
    
    # Note any wave without )
    if (sum(!wave_ind %in% (close_ind - 1))) {
      ind <- which(!wave_ind %in% (close_ind - 1))
      print(paste('wave without )    on sample', sample, 'lead', lead))
    }
    
    
    
    # Note any ) without (
    if (sum(!close_ind %in% (open_ind + 2))) {
      ind <- which(!close_ind %in% (open_ind + 2))
      print(paste(')    without (    on sample', sample, 'lead', lead))
    }
    
    # Note any ) without wave
    if (sum(!close_ind %in% (wave_ind + 1))) {
      ind <- which(!close_ind %in% (wave_ind - 1))
      print(paste(')    without wave on sample', sample, 'lead', lead))
    }
    
  }

  if (!sample %% 10) {
    print(paste(sample))
  }
}

# Write files -------------------------------------------------------------



save(ludb_set,file = 'ludb_set.RData')
load('ludb_set.RData')

# LUDB set testing --------------------------------------------------------
# starts <- rep(NA,200)
# ends <- rep(NA,200)
# for (i in 1:200){
#   table <- ludb_set[[i]]$annotation$i$sample
#   starts[i] <- table[1]
#   ends[i] <- table[length(table)]
# }

# Samples with first annotation after 2 sec (1000th index):
# 25  39  59 124 127 192
# 1035 1010 ~1100 ~1030 1366 ~1015

# Samples with last annotation before 8 sec (4000th index):
#  1     41  47    54  59   125   127   134   156   158   160
# 3995 3800 3000* 3925 3930 3850 3600*  3975  3995 3900  3700*

# sample = 160
# plot_func(ludb_set[[sample]]$signal$i,wfdb_ann2continuous2(ludb_set[[sample]]$annotation$i))


# delete? -----------------------------------------------------------------



# Download LUDB dataset under wes_ml/annotation: wget -r -N -c -np https://physionet.org/files/ludb/1.0.1/
# Read in files
#   Create list 'ludb' with $signal, $annotations, $info
#   Within both signal and annotations, create $i, ii, ...

# 200 total samples

# library(EGM)
# options(wfdb_path = '/mmfs1/home/dseaney2/wfdb/bin')
# file_path <- "/mmfs1/projects/cardio_darbar_chi/common/cohorts/ml_genetics/annotations/ludb/files"
# leads <- c('i','ii','iii','avr','avl','avf','v1','v2','v3','v4','v5','v6')
# 
# signal <- list()
# annotations <- list()
# 
# for (sample in 1:200) {
#   sig <- read_signal(record = sample,
#                      record_dir = file_path)
#   sig <- sig[,-1]
#   signal[[sample]] <- sig
#   signal[[paste0("sample_", sample)]] <- sig
# }
# 
# for (sample in 1:200) {
#   for (lead in 1:12) {
#     ann <- read_annotation(record = sample,
#                            record_dir = file_path,
#                            annotator = leads[lead])
#     
#   }
# }



# Separate into training and testing data (70/30)

# Training: 
#   LUDB preprocessing: (will need to handle for spectrogram too?)
#     Pick random 4 sec interval starting between 2 to 4 seconds. Set all other values to 0
#     Train using 


# Cubic spline for upsampling as needed