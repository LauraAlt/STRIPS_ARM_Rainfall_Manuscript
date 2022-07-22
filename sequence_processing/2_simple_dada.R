library(simple.dada)

for( folder in "sequence_data/sequence_files/STRIPS_Run_1/"){
  simple_dada(folder, cores = 3)
}

for( folder in "sequence_data/sequence_files/STRIPS_Run_2/"){
  simple_dada(folder, cores = 3)
}

for( folder in "sequence_data/sequence_files/STRIPS_Run_3/"){
  simple_dada(folder, cores = 3)
}
