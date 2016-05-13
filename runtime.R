wd <- '~/projects/nonparasatt/'

# Run script and get runtime
runtime <- system.time(source(paste0(wd, 'run.R')))

# Make readable
runtime <- as.data.frame(as.table(runtime))
names(runtime) <- c("Type of Time", "Time (seconds)")

# Write CSV
write.csv(runtime, 
          file = paste0(wd, "submission/runtime.txt"), 
          row.names = F)
