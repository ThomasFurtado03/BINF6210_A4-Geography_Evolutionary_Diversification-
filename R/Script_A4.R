##Required Packages----
install.packages("ggseqlogo")
install.packages("caper")

library(phytools)
library(ggplot2)
library(Biostrings)
library(tidyverse)
library(DECIPHER)
library(ggseqlogo)
library(phangorn)
library(ape)
library(caper)


#Data collection and inspection----

##Data Aqcuisition: Hyla ctyb sequences were downloaded (all results) from NCBI, after filtering for genomic DNA data only.

raw_seqs <- readDNAStringSet("data/sequence.fasta")

head(names(raw_seqs), 20)

seq_lengths <- width(raw_seqs)
summary(seq_lengths)


headers <- names(raw_seqs)
headers
species <- headers %>% 
  str_split(" ") %>% 
  map_chr(~ paste(.x[2], .x[3]))

head(species, 20)

##Data exploration and filtering----

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
#Lets filter by length, no less than 700, and no more than 1500. (This is the expected range for cytb!)

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

#Lets use an exploratory figure to evaluate the data set so far!
ggplot(seq_final, aes(x = length)) +
  geom_histogram(binwidth = 20, fill = "steelblue") +
  theme_minimal() +
  labs(
    title = "Distribution of CytB Sequence Lengths",
    x = "Length (bp)",
    y = "Count"
  )
#There is a little variances in the sizes, with a clear cluster around 900bp. Significant differences could be from extra flanking regions, or slightly different amplification fragments may have been used.
#Nonetheless, the data looks prepared for further analysis






##Alignment and visualization----

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


##Collection and addition of geographic information----

#Geographic regions were manually assigned at the species level using global distribution information from GBIF (gbif.org). A custom region table was created linking each Hyla species in the dataset to its known geographic distribution.

region_broad <- tibble(
  species = c(
    "Hyla annectans",
    "Hyla arborea",
    "Hyla carthaginiensis",
    "Hyla chinensis",
    "Hyla felixarabica",
    "Hyla hallowellii",
    "Hyla heinzsteinitzi",
    "Hyla intermedia",
    "Hyla meridionalis",
    "Hyla molleri",
    "Hyla orientalis",
    "Hyla perrini",
    "Hyla sanchiangensis",
    "Hyla sarda",
    "Hyla savignyi",
    "Hyla tsinlingensis"
  ),
  region_broad = c(
    "Asia",         # annectans
    "Europe",       # arborea
    "North Africa", # carthaginiensis
    "Asia",         # chinensis
    "Middle East",  # felixarabica
    "Asia",         # hallowellii
    "Middle East",  # heinzsteinitzi
    "Europe",       # intermedia
    "Europe",       # meridionalis
    "Europe",       # molleri
    "Middle East",  # orientalis
    "Middle East",  # perrini
    "Asia",         # sanchiangensis
    "Europe",       # sarda
    "Middle East",  # savignyi
    "Asia"          # tsinlingensis
  )
)


#Now finally, lets join this with my data so we can continue and answer geographic questions.

seq_with_regions <- seq_final %>%
  left_join(region_broad, by = "species")

table(seq_with_regions$region_broad)


#Lets turn this into a clear and simple histogram!

ggplot(seq_with_regions, aes(x = region_broad)) +
  geom_bar(fill = "steelblue") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Number of Hyla Species per Broad Geographic Region",
    x = "Region",
    y = "Species Count"
  )

#Just as a quick check, Summing all of these numbers up gives us 16, which matches what we have filtered down for!

##ML tree coloured by region----


#I will assign each region a generic colour
region_colours <- c(
  "Europe" = "#1f78b4",
  "Asia" = "#33a02c",
  "Middle East" = "#ff7f00",
  "North Africa" = "#e31a1c"
)

#Ensure consistent margins
par(mar = c(1, 1, 2, 1))

#Lets plot the bootstrapped tree with colours indicating which region they are from!
#(Run all lines in subsequent block)
plot(
  ml_tree_bs,
  cex = 0.8,
  label.offset = 0.010
)
tiplabels(
  pch = 21,
  bg = region_colours[seq_with_regions$region_broad],
  cex = 1.0,
  label.offset = 0.0100
)
legend(
  "topright",
  legend = names(region_colours),
  col = region_colours,
  pch = 19,
  bty = "n",
  xpd = NA,
  cex = 0.90,
  inset = c(0, -0.07)
)
title("Maximum Likelihood Tree of Hyla Species\nColoured by Broad Geographic Region",
      cex.main = 0.85 )

#Europe seems to be split, as some cluster and some do not. The

#Middle East is scattered, phylogenetically dispersed around the tree. 

#As for NA, there is insufficient data for comparison.

#Looking at the distribution by region, The species from Asia look to be very compact in the centre of the tree and all share a clade. I will test their relationship further.



##Performing D test----

#The D test tests a if binary trait is phylogenetically clutered, or randomly distributed across the tree.

#I will use this to answer my initial research question: 
#"Do Hyla species form distinct geographic clades corresponding to major regions?"


#I must first build a data frame for caper
trait_df <- data.frame(
  species = seq_with_regions$species,
  region = seq_with_regions$region_broad
)

#Make rownames which match the tree tips
rownames(trait_df) <- trait_df$species

#Binary Trait: from asia, or not from asia. Lets make it express the result numerically to easily use in the test.

trait_df$is_asia <- as.numeric(trait_df$region == "Asia")

#Let's check the results quickly
table(trait_df$is_asia)

#Nice, we have 11 non-asian species, and 5 asian species!

#Now I must prepare comparitive data as input for phylo.d

#But before this, I find I must "root my tree" for the analysis.
ml_tree_rooted <- midpoint(ml_tree_clean)
is.rooted(ml_tree_rooted)

ml_tree_rooted

?comparative.data
#Must enable vcv to calculate the covariance matrix
comp_data_asia <- comparative.data(
  phy = ml_tree_rooted,
  data = trait_df,
  names.col = species,
  vcv = TRUE,
  warn.dropped = TRUE
)

#Now lets run the D test!
?phylo.d

dtest_result_asia <- phylo.d(comp_data_asia, binvar = is_asia)
dtest_result_asia


#D = -2.94 (VERY negative). This means that the Asian species are strongly phlogenetically clustered.
#This means that the Asian species share a clear phylogenetic signal and likely diverged from a common Asian ancestor rather than being scattered across the tree!

#The test also reveals that it follows a Brownian clustering model quite closely (0.979)



#Now lets run the same test using the European species!
trait_df$is_europe <- as.numeric(trait_df$region == "Europe")
table(trait_df$is_europe)

comp_data_europe <- comparative.data(
  phy = ml_tree_rooted,
  data = trait_df,
  names.col = species,
  vcv = TRUE,
  warn.dropped = TRUE
)

dtest_result_europe <- phylo.d(comp_data_europe, binvar = is_europe)
dtest_result_europe


##Interpretation of results and literature consultation----


#These results are very interesting.. D is close to 1, which means that the trait is evolving randomly on the tree, and shows no real phylogenetic clustering. In other words, European species do NOT group together more than expected by chance!

#The p(Brownian) is not significant, and therefore does not show a pattern of phylogenetic conservatism.

#This means that the situation is as follows: Asia shows phylogentetic clustering, while others (Europe) are randomly distributed. Possibly indicating that geography incosistently affects diversification across Hyla.



#Asian Hyla species showed strong phylognetic clustering (D << 0), which likely reflects a combination of true biological history, and sampling structure.
#Upon consulting the literature, most Asian Hyla actually belong to a single and recent  monophyletic radiation, primarily 2 distinct clades existing in East Asia, both of which diverged only 5 million years ago. (Dufresnes et al., 2016). 

##Dufresnes, C., Litvinchuk, S. N., Borzée, A., Jang, Y., Li, J.-T., Miura, I., Perrin, N., & Stöck, M. (2016). Phylogeography reveals an ancient cryptic radiation in east-asian tree frogs (Hyla japonica group) and complex relationships between Continental and island lineages. BMC Evolutionary Biology, 16(1). https://doi.org/10.1186/s12862-016-0814-x 

#Since this is true, This means that there is a higher likelihood that closely related species will end up being sampled from Asian regions. Thus explaining the strong correlation detected by the D test. This means that these strong results from Asian species are true, as they are in fact closely related, but can ALSO be explained by sampling bias.


#Now how do we explain the results for Europe? Why are they not as tightly packed together as Asia, even though they are both their own distinct regions?


#The literature shows that Phylogenetic analysis revealed there to be THREE highly diverged lineages: H. arborea (Greece, North France, Central Europe), H. molleri (Iberian Peninsula [Spain/Portugal]) and H. orientalis (Northeastern Europe) (Stöck et al., 2012). 

#In addition to 3 highly diverged and old lineages, H. arborea exhibits molecular signatures of postglacial range expansion(Stöck et al., 2012).

#Given these pieces of information, the results seem a lot more reasonable. The answer is that Europe's geography is diverse, widespread, and consists of hybrid zones and glacial range expansions throughout the Earth's history. So it is quite reasonable that the species studied form Europe do not match as closely as the Asian ones do. And it comes down to the geographic history of the sampled regions! 

#Due to this, this investigation shows that geographic diversification patterns in Hyla differ STRONGLY among regions. As a result of specific geographic history, and sampling structure.


##Stöck, M., Dufresnes, C., Litvinchuk, S. N., Lymberakis, P., Biollay, S., Berroneau, M., Borzée, A., Ghali, K., Ogielska, M., & Perrin, N. (2012). Cryptic diversity among western palearctic tree frogs: Postglacial range expansion, range limits, and secondary contacts of three European tree frog lineages (Hyla arborea group). Molecular Phylogenetics and Evolution, 65(1), 1–9. https://doi.org/10.1016/j.ympev.2012.05.014 


