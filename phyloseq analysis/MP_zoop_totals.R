rm(list = ls())

metadata <- read.csv("sam_data.csv")

nrow(metadata[metadata$Organism == 'Calanoid', ])
nrow(metadata[metadata$Organism == 'Cyclopoid', ])
nrow(metadata[metadata$Organism == 'Daphnia', ])
nrow(metadata[metadata$Organism == 'Ceriodaphnia', ])
nrow(metadata[metadata$Organism == 'Bosmina', ])
nrow(metadata[metadata$Organism == 'Holopedium', ])
nrow(metadata[metadata$Sample_Type == 'Zooplankton', ])
nrow(metadata[metadata$Sample_Type == 'Water', ])

sum.run<-metadata %>% count(Time.Point, Lake, Organism)
