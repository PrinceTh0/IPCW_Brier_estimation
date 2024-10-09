
remove(list=ls()) 

load("breast.RData")
# contains data frame named "DF",
# as obtained from pre-processing using R package SEERaBomb

#path to output pre-processed ds
out.path <- "..."

library(forcats)

## Pre-processing, part 1

# consider female patients only
DF <- DF[DF$sex == 2, ]
# exclude patients aged 76 or older
DF <- DF[DF$agedx <= 75, ]
# exclude patients with unclear diagnostic confirmation
DF <- DF[DF$dxconf < 9, ]
# restrict data to patients with first malignant tumor
DF <- DF[DF$seqnum <= 1, ]
DF <- DF[DF$firstprm == 1, ]
# restrict data to patients without distant metastases
DF <- DF[DF$adjm6value == 0 & !is.na(DF$adjm6value), ]


## Define outcome variable and event indicator

survtime <- DF$surv
event <- DF$"dthclass"

## Define covariates

COVARS <- DF[,
             c(
               "marstat", # marital status
               "racrecy", # race recode
               "nhia", # Spanish origin
               "agedx", # age at diagnosis
               "yrdx", # year of diagnosis
               "siteo2", # primary site of tumor
               "lateral", # laterality
               "beho3", # behavior
               "grade",
               "dxconf", # diagnostic confirmation
               "eod10pn", # lymph nodes positive
               "eod10ne", # nodes examined
               "eod10sz", # tumor size (1973 - 2003)
               "cssize", # tumor size (2004 - 2010)
               "surgprim", # surgery of primary site
               "radiatn", # type of radiation therapy
               "radsurg", # radiation sequence
               "numprims", # number of primaries
               "firstprm", # first malignant primary
               "erstatus", # ER status
               "prstatus", # PR status
               "insrecpub", # insurance
               "adjtm6value", # ADJUSTED AJCC 6TH T Stage
               "adjnm6value", # ADJUSTED AJCC 6TH N Stage
               "adjm6value" # ADJUSTED AJCC 6TH M Stage
             )]

################################################################################

# Pre-processing, part 2

breast <- data.frame(survtime, event, COVARS)

# restrict years of diagnosis to 1998 - 2010
breast04 <- breast[breast$yrdx >= 1998,]
# exclude patients with unknown event indicator
breast04 <- breast04[breast04$event != 9,]
# exclude patients aged 17 or younger
breast04 <- breast04[breast04$agedx >= 18,]
# restrict data to patients with known number of positive lymph nodes
breast04 <- breast04[breast04$eod10pn <= 90,]
# restrict data to patients with known number of examined lymph nodes
breast04 <- breast04[breast04$eod10ne <= 90,]
# restrict data to patients with known tumor size
breast04 <- breast04[breast04$cssize < 990 | breast04$eod10sz < 990,]
breast04$tumorsize <- breast04$cssize
breast04$tumorsize[is.na(breast04$tumorsize)] <- breast04$eod10sz[is.na(breast04$tumorsize)]
breast04$tumorsize[breast04$tumorsize==888] <- NA
breast04 <- breast04[,!(names(breast04) %in% c("yrdx", "cssize", "eod10sz", "insrecpub","adjm6value"))]

breast04$survtime <- breast04$survtime + 1
names(breast04)[1] <- "time.discrete"
names(breast04)[2] <- "state"

# convert surgery of primary site to a factor (yes/no)
breast04$surgprimRecode <- breast04$surgprim
breast04$surgprimRecode[breast04$surgprim %in% c(10:90)] <- 10

# delete cases with missing values and delete columns that are not needed for analysis
breast04 <- breast04[complete.cases(breast04),]
breast04 <- breast04[,!(names(breast04) %in% c("beho3", "firstprm", "surgprim"))]

# exclude patients with undetermined cell type
breast04 <- breast04[breast04$grade < 9,]
# restrict marital status to single, married, separated, divorced, widowed
breast04 <- breast04[breast04$marstat <= 5,]
breast04 <- breast04[breast04$racrecy <= 9,]
# restrict data to tumors with known laterality
breast04 <- breast04[breast04$lateral <= 2,]
# exclude diagnostic confirmations other than positive histology
breast04 <- breast04[breast04$dxconf == 1,]
# exclude patients with unclear or inconsistent TNM stage classification
breast04 <- breast04[breast04$adjnm6value < 60,]
breast04 <- breast04[breast04$adjtm6value < 60,]
# exclude patients with unknown radiation therapy and/or sequence
breast04 <- breast04[breast04$radiatn < 8,]
breast04 <- breast04[breast04$radsurg < 9,]
# exclude patients with unknown ER or PR status
breast04 <- breast04[breast04$erstatus <= 2,]
breast04 <- breast04[breast04$prstatus <= 2,]

# delete columns that are not needed for analysis
breast04 <- breast04[,!(names(breast04) %in% c("dxconf"))]

# convert categorical variables into factor variables
for (k in 1:ncol(breast04))
  if(!(names(breast04)[k] %in% c("time.discrete", "state",
                                 "agedx", "eod10pn", "tumorsize",
                                 "eod10ne", "numprims"))) breast04[,k] <- as.factor(breast04[,k])

###additional data prepping


#variable adjtm6value, description:
#https://seer.cancer.gov/seerstat/variables/seer/ajcc-stage/6th/breast.html#t
#0: no evidence of primary tumor (but either positive lymph nodes or distant mets)
#5: in situ and Paget's disease with no invasive underlying tumor
#11: microinvasion 0.1 cm or less in greatest dimension
#12: tumor more than 0.1 cm but not more than 0.5 cm in greatest dimension
#15: tumor more than 0.5 cm but not more than 1 cm in greatest dimension
#18: tumor more than 1 cm but not more than 2 cm in greatest dimension
#20: tumor more than 2 cm but not more than 5 cm in greatest dimension
#30: tumor more than 5 cm in greatest dimension
#41: direct extension to the chest wall
#42: edema (including peau d'orange) or ulceration of the skin of the breast, or satellite skin nodules confined to the same breast
#43: both 41 & 42
#44: inflammatory carcinoma

#recode nhia to Hispanic vs. Non-Hispanic origin
breast04$nhiarec <- fct_collapse(breast04$nhia, "1" = c("1","2","3","4","5","6","7","8"))
#recode radsurg to combine intraoperative therapy categories
breast04$radsurgrec <- fct_collapse(breast04$radsurg, "5" = c("5","6"))
#recode adjtm6value to combine categories (0 and 5), (41, 42 and 43)
breast04$adjtm6valuerec <- fct_collapse(breast04$adjtm6value, "00a" = c("0","5"), "01a" = c("41","42","43"))

#drop old unused variables
breast04 <- breast04[,!(names(breast04) %in% c("nhia","radsurg","adjtm6value"))]

#save to .RData
save(breast04, list = "breast04", file = file.path(out.path,"breast04.RData"))

#clean workspace, remove breast and DF
rm(breast)
rm(DF)
gc()
