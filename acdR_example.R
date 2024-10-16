
#
### Example of using acdR
#
# Requires:
#    - RSQLite
#    - AcdFiaBench.db

library( dplyr )
library( acdR )
options(dplyr.summarise.inform = FALSE)

# constants
METRIC <- 0
IMPERIAL <- 1

###################

db.dir <- "DIRECTORY WITH YOUR FIA DB FILE"

# connect to FIA Benchmarking DB
db <- dbConnect( RSQLite::SQLite(), paste0(db.dir,"AcdFiaBench.db") )

# get list of selected plots (UnitCountyPlot is the key)
selected.plots <- dbGetQuery( db, "select * from FiaSelectedPlots") %>% filter( MeasurePeriod > 0 )

# pick an arbitrary plot to project
plot <- selected.plots[1,]

# set stand init data
stand.init <- dbGetQuery( db,
                          paste0( "select * from FVS_StandInit where UnitCountyPlot == '",
                                  trimws(plot$UnitCountyPlot), "' and TrtPeriod == '", plot$TrtPeriod, "'")) %>% arrange( Inv_Year )

# get tree list
treelist <- dbGetQuery( db, paste0( "select * from FVS_TreeInit where Stand_CN == '", stand.init$Stand_CN[1], "'"))

# calculate ba expansion factor for the plot design
ba.factor <- ifelse( stand.init$Basal_Area_Factor[1] == 0, 1.0, -stand.init$Basal_Area_Factor[1] )


# compute trees per acre using plot design info
treelist$tpa <- treelist$TREE_COUNT * ifelse( treelist$DBH < stand.init$Brk_DBH[1] & stand.init$Brk_DBH[1] != 999,
                                              stand.init$Inv_Plot_Size[1],
                                              ba.factor ) / stand.init$Num_Plots[1]

# make sure tree variables are properly set (e.g., missing set to 0)
treelist$HT <- ifelse( is.na(treelist$HT), 0.0, treelist$HT )
treelist$cr <- ifelse( !is.na(treelist$CRRATIO), treelist$CRRATIO/100.0, 0 )

# set various ACD options
region <- "ME"
age <- ifelse( is.na(stand.init$Age[1]), as.integer(10), as.integer(stand.init$Age[1]) )
use_sbw <- as.integer(0)
use_hw <- as.integer(0)
use_thin <- as.integer(0)
treelist$form <- as.integer(0)
treelist$risk <- as.integer(0)
cdef <- -1.0
use_ingrowth <- as.integer(0)
cut_point <- 0.5
MinDBH <- 3 / 2.54

periods <- as.integer(10)

grown.trees <- grow_acd( periods,
                         region,
                         age,
                         IMPERIAL,
                         stand.init$CSI[1]*3.28, # convert CSI to imperial units
                         as.numeric(stand.init$ElevFt[1]),
                         cdef,
                         use_sbw, use_hw, use_thin,
                         use_ingrowth, cut_point, MinDBH,
                         as.integer(treelist$Plot_ID),
                         as.integer(treelist$Tree_ID),
                         as.integer(treelist$SPECIES),
                         treelist$DBH,
                         treelist$HT,
                         treelist$tpa,
                         treelist$cr,
                         treelist$form,
                         treelist$risk )


