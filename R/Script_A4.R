##Required Packages----
install.packages("ggseqlogo")

library(Biostrings)
library(tidyverse)
library(DECIPHER)
library(ggseqlogo)
library(phangorn)
library(ape)

##Data Aqcuisition: Hyla ctyb sequences were downloaded (all results) from NCBI, after filtering for genomic DNA data only.

raw_seqs <- readDNAStringSet("data/sequence.fasta")

head(names(raw_seqs), 20)

seq_lengths <- width(raw_seqs)
summary(seq_lengths)


headers <- names(raw_seqs)

species <- headers %>% 
  str_split(" ") %>% 
  map_chr(~ paste(.x[2], .x[3]))

head(species, 20)

#convert to df

seq_df <- tibble(
  header = headers,
  species = species %>% str_squish(),
  length = width(raw_seqs)
)

head(seq_df)

#Lets look at the amount of species we have in this dataset.

table(seq_df$species)

#Lets filter the data set.
#Lets filter by length, no less than 700, and no more than 1500

seq_length_filtered <- seq_df %>% 
  filter(length >= 700, length <= 1500)

nrow(seq_length_filtered)
length(unique(seq_length_filtered$species))

#To further filter, we will keep the longest sequence for each species, which should be the most complete cytb

seq_longest_per_species <- seq_length_filtered %>% 
  group_by(species) %>% 
  arrange(desc(length)) %>% 
  slice(1) %>% 
  ungroup()

nrow(seq_longest_per_species)
distinct(seq_longest_per_species, species)

#I notice there is some junk, so I will need to keep filtering

seq_taxa_filtered <- seq_longest_per_species %>% 
  filter(!species %in% "UNVERIFIED: Hyla") %>% 
  filter(!str_detect(species, "Hyla cf\\."))


seq_final <- seq_taxa_filtered %>% 
  group_by(species) %>% 
  slice_max(order_by = length, n = 1) %>% 
  ungroup()

nrow(seq_final)

seq_final$species

#Great, we have removed <700 and >1500, we have filtered the data down to one sequence per species. This now leaves us with 16 species to look at for our analysis!

#Now, I must retrive the DNA data from the raw_seqs stringset, that corresponds to my final species list.

final_headers <- seq_final$header

final_seqs <- raw_seqs[final_headers]
final_seqs

#Now, let's align these sequences

aligned_seqs <- AlignSeqs(final_seqs)

aligned_seqs
BrowseSeqs(aligned_seqs)

##Build Neighbor Joining Tree----

#Lets create a neighbor joining tree to start. It will serve as a quick/exploratory tree, just as a reference point.

#Convert alignment to Phangorn object (must first be recognized as a matrix)

alignment_matrix <- as.matrix(aligned_seqs)

Phang_alignment <- phyDat(alignment_matrix, type = "DNA")
Phang_alignment


dist_matrix <- dist.ml(Phang_alignment)

nj_tree <- NJ(dist_matrix)
plot(nj_tree, cex = 0.4)

#This tree is completely illegible, so I will rename the tree labels, and adjust the margins.
#I recall I have seq_final as a tibble, which I can use for renaming

label_map <- setNames(seq_final$species, seq_final$header)

nj_tree_clean <- nj_tree

nj_tree_clean$tip.label <- label_map[nj_tree_clean$tip.label]

par(mar = c(1, 1, 2, 1))
#Adjusting bottom, left, top, and right margins

plot(
  nj_tree_clean,
  cex = 0.7,
  label.offset = 0.002,
)
title("Neighbor-Joining Tree (species labels only)", cex.main = 1.0)


##Build Maximum Likelihood Tree----


#Lets evaluate and choose the best substitution model
#I will use modelTest function, which tests my alignment using many different evolutionary models
model_test <- modelTest(Phang_alignment)

#Now out of this long input, I need to find the one with the lowest AIC score, this will tell me which model fits my data the best.

best_model <- model_test$Model[which.min(model_test$AIC)]
best_model

#The model with the lowest AIC score is "TPM2u+G(4)+I"

#Now, lets fit the ML model, using the best model

fit <- pml(nj_tree, Phang_alignment)

?optim.pl
#Note, I cannot easily use "best_model" as input for the next command, as it is not one of the accepted inputs, so I must express my model in terms of an acceptable input.

best_model
#I must edit the name of my model to input it in.

best_model_base_name <- sub("\\+.*", "", best_model)
best_model_base_name

#optGamma enables the G, optInv enables I. Which were both present in the best model.
#optBf and optEdge enable optimization of base frequencies, and branch lengths, respectively.

fit_opt <- optim.pml(
  fit,
  model = best_model_base_name,  
  optGamma = TRUE,
  optInv   = TRUE,
  optBf    = TRUE,
  optEdge  = TRUE
)

fit_opt
plot(fit_opt$tree, cex = 0.6)

#Again the tree is not legible, so I will perform similar edits to prepare it.

#Recall:
label_map <- setNames(seq_final$species, seq_final$header)

ml_tree_clean <- fit_opt$tree
ml_tree_clean$tip.label <- label_map[ml_tree_clean$tip.label]

par(mar = c(1, 1, 2, 1))
#Adjusting bottom, left, top, and right margins
#Redundant, but just for sanity check

plot(ml_tree_clean, cex = 0.7)
title("Maximum Likelihood Tree (species labels only)")

####if par is broken run, then try again:
#dev.off()

##Adding bootstrap values to ML tree----

##Let's add bootstrap values to the ML tree, finally producing our fully polished tree!

#For reproducibility:
set.seed(123) 

#Lets set 100 bootstrap replicates
#Note: the optNni option allows for nearest-neighbor interchange

bs_reps <- bootstrap.pml(fit_opt, bs = 100, optNni = TRUE)



#upon lengthy troubleshooting, I found that the plotBS() is using the long ugly names from my original tree.
#This was causing an error resulting in "unmatched tip labels" so I have to fix this too.
#The following is adapted from phangorn documentation, to sequentially rename bootstrapped tree labels. (What a headache!)

bs_reps_renamed <- lapply(bs_reps, function(tr) {
  tr$tip.label <- label_map[tr$tip.label]
  tr
})

#The class must also be switched to "multiphylo" so it's usable
class(bs_reps_renamed) <- "multiPhylo"


ml_tree_bs <- plotBS(ml_tree_clean, bs_reps_renamed,)

#Again, lastly I will edit the tree to make it more presentable.

plot(ml_tree_bs, cex = 0.7, label.offset = 0.002)
nodelabels(
  ml_tree_bs$node.label,
  frame = "none",
  cex = 0.7,
  adj = c(1.0, -0.2)
)
title("Maximum Likelihood Tree with Bootstrap Support")



