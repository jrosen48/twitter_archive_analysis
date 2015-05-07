#---------------------------------
# Setting up
#---------------------------------

require(curl)
require(ggplot2)
require(plyr)
require(dplyr)
require(SnowballC)
require(tm)
require(lsa)
require(cluster)

#---------------------------------
# Getting data
#---------------------------------

setwd("~/Dropbox/Statistics/my_archive/twitter_archive_analysis/csvs")

path = "~/Dropbox/Statistics/my_archive/twitter_archive_analysis/csvs"

out.file <- ""
file.names <- dir(path, pattern =".csv")

for(i in 1:32){ #change to load more files
  file <- read.table(file.names[i],header=TRUE, sep=";", stringsAsFactors=FALSE)
  out.file <- rbind(out.file, file)
} 

my_csv_list <- list()
my_csv_text_list <- list()
for (i in 1:length(file.names)){
  my_csv_list[[i]] <- read.table(file.names[i], sep=",", header = TRUE, fill = TRUE)
  my_csv_text_list[[i]] <- my_csv_list[[i]][8]
}

str(my_csv_text_list)

#---------------------------------
# Processing text
#---------------------------------

# StemCompletion modification
stemCompletion_mod <- function(x, dict=myCorpus.copy) {
  PlainTextDocument(stripWhitespace(paste(stemCompletion(unlist(strsplit(as.character(x), " ")),dictionary=dict, type="shortest"),sep="", collapse=" ")))
}

# Remove URL
removeURL <- function(x) gsub("http[[:alnum:]]*", "", x)

# Remove @Username
removeUsername <- function(x) gsub("@\\w+", "", x)

# RemoveHashtag
removeHashtag <- function(x) gsub("#\\w+", "", x)

custom_stopwords <- c()

my_list <- list()
for (i in 1:length (my_csv_text_list)){ #change this
  #temp_df <- hashtags_list[[i]][7]
  temp_df <-  my_csv_text_list[[i]]
  print(paste("Currently processing hashtag", "#", i))
  myCorpus <- Corpus(DataframeSource(temp_df)) # change to grouped_df for all terms
  myCorpus <- tm_map(myCorpus, content_transformer(tolower))
  myCorpus <- tm_map(myCorpus, removeUsername)
  #myCorpus <- tm_map(myCorpus, removeHashtag)
  myCorpus <- tm_map(myCorpus, removePunctuation)
  myCorpus <- tm_map(myCorpus, removeNumbers)
  myCorpus <- tm_map(myCorpus, removeWords, stopwords("english"))
  myCorpus <- tm_map(myCorpus, removeWords, custom_stopwords, mc.cores=1)
  myCorpus <- tm_map(myCorpus, removeURL, mc.cores = 1)
  myCorpus <- tm_map(myCorpus, stemDocument, mc.cores = 1)
  myCorpus.copy <- myCorpus
  #myCorpus <- tm_map(myCorpus, stemCompletion_mod, mc.cores = 1)
  myCorpus <- tm_map(myCorpus, stripWhitespace, mc.cores = 1)
  myCorpus <- tm_map(myCorpus, PlainTextDocument)
  name <- paste("doc", i, sep = "")
  my_list[[name]] <- myCorpus
}


# Create a Term Document Matrix

myTDM <- list()
for (i in 1:length(my_csv_text_list)){
  name <- paste("doc", i, sep = "")
  myTDM[[name]] <- TermDocumentMatrix(my_list[[name]], 
                                      control = list(weighting = function(x) weightTfIdf(x, normalize = TRUE)))
}

# can replace above
# myTDM[[name]] <- TermDocumentMatrix(my_list[[name]], control = list(weighting = weightTfIdf, normalize = TRUE))

myTDM_all <- myTDM[[1]]
for (i in 2:length(my_csv_text_list)){
  myTDM_all <- c(myTDM_all, myTDM[[i]])
}

myTDM_all
myTDM

# Removing sparse terms

myTDM_common <- removeSparseTerms(myTDM_all, .995) # ignore maybe
myTDM_common

# Creating clusters
d <- dist(myTDM_common, method = "euclidean")
kfit <- kmeans(d, 11)

# Combining clusters with words and number of words
y <- rowSums(as.matrix(myTDM_common))
x <- kfit$cluster

my_df <- data.frame(x, y)
my_df <- add_rownames(my_df, var = "word")
my_df <- arrange(my_df, x, desc(y))

# Generating total frequencies for all terms and for each hashtag
all_terms <- rowSums(inspect(myTDM_common[,]))
all_terms <- all_terms[sort(names(all_terms))]

myTDM_list <- list()
for (i in 1:length(myTDM)){
  myTDM_list[[i]] <- removeSparseTerms(myTDM[[i]], .99)
  myTDM_list[[i]] <- rowSums(inspect(myTDM_list[[i]][,]))
}

length(myTDM)
myTDM_list

# Use myTDM_list to create overlapping common terms

tags_vec_list <- list()
for (i in 1:length(myTDM_list)){
  tags_vec_list[[i]] <- vector(length = dim(myTDM_common)[1]) 
  names(tags_vec_list[[i]]) <- names(all_terms) 
  for (j in 1:length(all_terms)){
    if (names(all_terms[j]) %in% names(myTDM_list[[i]])) {
      tags_vec_list[[i]][j] <- all_terms[j] 
    }
    else {tags_vec_list[[i]][j] <- 0} 
  }
}

# Creating vectors for each cluster with the total number of times each term in the cluster occurs across all terms

my_df_clusters <- list()
for (i in 1:dim(table(my_df[2]))){
  my_df_clusters[[i]] <- filter(my_df, x == i)
  my_df_clusters[[i]] <- my_df_clusters[[i]][,c(1,3)]
}

my_df_clusters

my_clust_list <- list()
for (i in 1:dim(table(my_df[2]))){
  my_clust_list[[i]] <- vector(length = dim(myTDM_common)[1]) 
  names(my_clust_list[[i]]) <- names(all_terms) 
  for (j in 1:length(all_terms)){
    if (names(all_terms[j]) %in% my_df_clusters[[i]][,1]) {
      my_clust_list[[i]][j] <- all_terms[j] 
    }
    else {my_clust_list[[i]][j] <- 0} 
  }
}

my_clust_list

# Computing similarities 

my_df_tags_vec <- list()
my_cosine_list <- list()

for (i in 1:length(tags_vec_list)) {
  my_df_tags_vec[[i]] <- vector()
  for (j in 1:dim(table(my_df[2]))) {
    my_df_tags_vec[[i]] <- append(my_df_tags_vec[[i]], cosine(tags_vec_list[[i]], my_clust_list[[j]]))
    my_cosine_list[[i]] <- my_df_tags_vec[[i]]
  }
}

my_cosines <- as.data.frame(do.call(rbind, my_cosine_list))
for (i in length(my_cosines)){
  names(my_cosines[i]) <- paste("c", i)
}

my_cosines <- my_cosines[1:36,]
my_cosines
write.csv(my_cosines, "my_archive_cosines.csv")

# Interpretation

print(my_df_clusters)
print(my_cosines)

#SAMR
#TPACK
#NGSSChat
#EdTechBridge
#TeacherEdChat
#SciChat

# View(my_cosine_list[[1]])
# View(my_cosine_list[[2]])
# View(my_cosine_list[[3]])
# View(my_cosine_list[[4]])
# View(my_cosine_list[[5]])
# View(my_cosine_list[[6]])

# Dodged bar charts
# ggplot(diamonds, aes(clarity, fill=cut)) + geom_bar(position="dodge")


