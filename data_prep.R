##################
### Setting up ###
##################

rm(list=ls())

### Load packages
library(AnthroTools)  # for salience calculations

### Load data
flw1 <- read.csv("FreeList_CERC_V0.1_FIN.csv", sep = ";") # wave 1 free-list data
flw2 <- read.csv("FreeList_Wave2_BGD.csv", sep = ";") # wave 2 free-list BGD data
dfallw <- read.csv("MoralizingGods_data set_V1.0_wIDs.csv", sep = ";") # wave 1 and 2 game and demographic, incl. original IDs for merging across waves

### Subset 
fldata <- data.frame(
  Group = c(flw1$Culture, flw2$CULTURE),
  ID = c(flw1$FLID, flw2$LOCID),
  Order = c(flw1$Order, flw2$ORDER),
  BGD_ORIG = c(flw1$BGD_ORIG_NC, flw2$BGD),
  BGD_GEN = c(flw1$BGD, flw2$BGD_GEN_TB),
  BGD_SPEC1 = c(flw1$BGD_SPEC_TB_cl,flw2$BGD_SPEC_AO),
  BGD_SPEC2 = c(rep(NA, length(flw1$BGD_SPEC_TB_cl)),flw2$BGD_SPEC_CS),
  DIS_BGD = c(flw1$DIS_BGD_AB, flw2$DIS_BGD)
)

fldata$BGD_GEN <- stringr::str_to_title(fldata$BGD_GEN) # streamline so first letter of general codes are uppercase

### Check inter-coder reliability
aggregate(DIS_BGD ~ BGD_GEN + Group, data=fldata, FUN=sum)

### General codes: Morality and Virtue

table(unique(fldata$BGD_GEN)) # check that general codes are streamlined

bgd.sal <- CalculateSalience(fldata, Order = "Order", Subj = "ID",
                             CODE = "BGD_GEN", GROUPING = "Group", Rescale = FALSE, Salience = "BGD.S")

labs <- c("Group", "ID", "Order", "BGD_GEN", "BGD.S") # select relevant variables
bgd.fl <- bgd.sal[labs]
bgd.df <- bgd.fl[complete.cases(bgd.fl), ] # Only complete cases in free-list data

bgd.freq <- FreeListTable(bgd.df, CODE = "BGD_GEN", Order = "Order", Salience = "BGD.S", # Calculate frequency of FL items
                         Subj = "ID", GROUPING = "Group",
                         tableType = "Frequency")

rownames(bgd.freq) <- NULL # remove superfluous rownames

mvlabs <- c("Subject", "Morality", "Virtue", "Group")

bgd.morality <- bgd.freq[mvlabs]

### Specific codes: Game relevant responses

fldata$BGD_SPEC1 <- stringr::str_to_title(fldata$BGD_SPEC1) # streamline so first letter of codes are uppercase
fldata$BGD_SPEC2 <- stringr::str_to_title(fldata$BGD_SPEC2) # streamline so first letter of codes are uppercase

gamecode.words <- c("Dishonesty", "Dishonest", "Deception", "Truthfulness.n", "Not Truthful", "Helpful.n", "Lies",
                    "Stealing", "Theft", "Greed", "Selfish", "Not Sharing", "Egoism", "Cooperation.n", "Kindness.n",
                    "Justice.n", "Unjustice", "Exploiting Others", "Compassion.n", "Betray", "Betrayal",
                    "Loyalty.n", "Loyal.n", "Disloyal", "Disloyalty")

fldata$gamecode <- ifelse(fldata$BGD_SPEC1 %in% gamecode.words | fldata$BGD_SPEC2 %in% gamecode.words, 1,0)

unique(fldata$BGD_SPEC1[fldata$gamecode==0]) # check non-game codes
unique(fldata$BGD_SPEC2[fldata$gamecode==0]) # check non-game codes

hist(fldata$gamecode)

bgd.spec <- fldata[complete.cases(fldata$BGD_SPEC1),]

spec.n <- aggregate(gamecode ~ ID + Group, data=bgd.spec, FUN=sum)

# get total number of free-list responses in BGD_SPEC
bgd.n <- aggregate(BGD_GEN ~ ID + Group, data=fldata, function(x) {sum(!is.na(x))})

colnames(bgd.n) <- c("ID", "Group", "BGD.TOTAL")

# merge number of "Morality" and "Virtue" with total free-list responses
colnames(bgd.morality) <- c("ID", "Morality", "Virtue", "Group")

bgd.prop <- merge(bgd.morality, bgd.n, by = c("ID", "Group"))

bgd.data <- merge(bgd.prop, spec.n, by = c("ID", "Group"))

# compute proportions
bgd.data$MPROP <- with(bgd.data, Morality/BGD.TOTAL)
bgd.data$MVPROP <- with(bgd.data, (Morality+Virtue)/BGD.TOTAL)
bgd.data$SPECPROP <- with(bgd.data, gamecode/BGD.TOTAL)
str(bgd.data)
hist(bgd.data$MPROP)
hist(bgd.data$MVPROP)
hist(bgd.data$SPECPROP)

### select game and demo vars
gdvars <- c("SITE", # field site
            "Original.ID", # original participant id
            "CID",  # participant id
            "RAG.ORDER", # order of RAGs
            "RAG.DISTANT1", # coins to DISTANT (in local vs. distant RAG)
            "RAG.DISTANT2",  # coins to DISTANT (in self vs. distant RAG)
            "HONEST13", # game check (RAGs)
            "DG.ORDER", # order of DGs
            "DG.DISTANT1", # coins to DISTANT (in local vs. distant DG)
            "DG.DISTANT2",  # coins to DISTANT (in self vs. distant DG)
            "DG.HONEST", # game check (DG)
            "CHILDREN", # number of children
            "MAT1",  # material security
            "BGPUNISH", # moralistic god punitiveness 
            "BGDIE", # moralistic god punitiveness
            "BGFEEL", # moralistic god omniscience
            "BGSEE", # moralistic god omniscience
            "BGSTLIMP", # moralistic god punish stealing
            "BGLIEIMP", # moralistic god punish lying
            "BGMURDIMP" # moralistic god punish murder
            )

gddf <- dfallw[gdvars]

dff <- with(gddf, 
            data.frame(
              Group = SITE, # field site
              ID = Original.ID, # original participant id
              CID = CID, # participant id
              RAG.ORDER = ifelse(grepl('^1', RAG.ORDER), 1, 0), # if local game first, then 1
              RAG.LOCAL = RAG.DISTANT1, # RAG LOCAL vs. DISTANT
              RAG.SELF = RAG.DISTANT2, # RAG SELF vs. DISTANT
              RAG.CHECK = HONEST13, # RAG game check
              DG.ORDER = ifelse(grepl('^1', DG.ORDER), 1, 0), # if local game first, then 1
              DG.LOCAL = DG.DISTANT1, # DG LOCAL vs. DISTANT
              DG.SELF = DG.DISTANT2, # DG SELF vs. DISTANT
              DG.CHECK = DG.HONEST, # DG game check
              CHILDREN = CHILDREN, # number of children
              MAT = MAT1, # material security
              DIEPUN = rowMeans(gddf[, 14:15], na.rm = TRUE), # compute punishment score
              OMNI.BG = rowMeans(gddf[, 16:17], na.rm = TRUE), # compute omniscience score
              MINDEX = rowMeans(gddf[, 18:20], na.rm = TRUE) # compute morality item scale index
            )
          )
str(dff)

unique(dff$Group) # Check groups in demo df
unique(bgd.data$Group) # Check groups in free-list df

dff$Group <- gsub("Co.Tanna", "Coastal Tanna", dff$Group) # streamline site name
dff$Group <- gsub("In.Tanna", "Inland Tanna", dff$Group) # streamline site name
dff$Group <- gsub("Yasawa", "Yasawa Fiji", dff$Group) # streamline site name
dff$Group <- gsub("Lovu", "Lovu Fiji", dff$Group) # streamline site name
bgd.data$Group <- gsub("Peru", "Huatasani", bgd.data$Group) # streamline site name
bgd.data$Group <- gsub("Congo", "Kananga", bgd.data$Group) # streamline site name

df <- merge(bgd.data, dff, by = c("ID", "Group"), all.x = TRUE, all.y = TRUE) # merge free-list and game/demo data

unique(df$Group) # Check groups

df$NEW.ID <- 1:nrow(df) # provide new unique ID number

# stack dataset for modeling
gamecols <- c("RAG.LOCAL", "RAG.SELF","DG.LOCAL", "DG.SELF" )

ddf <- tidyr::pivot_longer(
  df,
  all_of(gamecols),
  names_to = "TYPE",
  values_to = "Y")

ddf$SELF <- ifelse(grepl("SELF", ddf$TYPE), 1, 0) # if SELF game, then 1

ddf$RAG <- ifelse(grepl("RAG", ddf$TYPE), 1, 0) # if RAG game, then 1

# select final columns and streamline names

dddf <- with(ddf, 
             data.frame(
               NEW.ID = NEW.ID,
               GROUP = Group,
               BGD.MORALITY = Morality,
               BGD.TOTAL = BGD.TOTAL,
               BGD.SPECCODE = gamecode,
               BGD.MORPROP = MPROP,
               BGD.MORVIRTPROP = MVPROP,
               BGD.SPECPROP = SPECPROP,
               RAG.ORDER = RAG.ORDER,
               RAG.CHECK = RAG.CHECK,
               DG.ORDER = DG.ORDER,
               DG.CHECK = DG.CHECK,
               CHILDREN = CHILDREN,
               MAT = MAT,
               DIEPUN = DIEPUN,
               OMNI.BG = OMNI.BG,
               MORINDEX = MINDEX,
               TYPE = TYPE,
               SELF = SELF,
               RAG = RAG,
               Y = Y
             ))

# export prepped data
write.csv(dddf, file = "data_prepped.csv", row.names = FALSE)

### End
