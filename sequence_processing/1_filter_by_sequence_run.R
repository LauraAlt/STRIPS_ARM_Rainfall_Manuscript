a <- data.table::fread('sequence_data/Armstrong_Sequences_Metadata.csv', sep = ',')
data.table::setkey(a, sample_name)

seq_files <- dir('sequence_data/sequence_files')
for(sample in seq_along(a$sample_name)){
  dir.create(file.path('sequence_data/sequence_files', a$miseq_run[sample]), showWarnings = FALSE)
  files <- seq_files[grep(gsub('_','_',a$sample_name[sample]), seq_files)]
  for(file in files){
    file.rename(file.path('sequence_data/sequence_files', file), file.path('sequence_data/sequence_files', a$miseq_run[sample], file))
    
  }
}

seq_folders <- dir('sequence_data/sequence_files', pattern = 'Run', full.names = TRUE)
