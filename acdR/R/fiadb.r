


#' Establish connection to FIADB database
#'
#' @param fiadb -- fully qualified file name of an FIADB SQLite database
#'
#' @return
#' Returns an connection to the SQLite DB.
#'
#' @examples
#' # connect to Maine FIADB database
#' # fia.connect( "c:/temp/ME.db" )
#' @export

fia_connect <- function( fiadb )
{
    db <- dbConnect( RSQLite::SQLite(), fiadb )
    return( db )
}




#' Get a STAND_CN's data from an FIADB database for growing with ACD.
#'
#' @param fiadbcon -- Connection to an FIADB database (see `fia.connect()`)
#' @param stand.condition.code -- FIA stand condition code (STAND_CN) to retrieve
#'
#' @return
#' Returns a list with two members if successful:
#' \itemize{
#'     \item stand : a data.frame with stand/plot level data from the STAND_CN
#'     \item tree : a data.frame with tree records from the STAND_CN
#' }
#'
#' If an error occurs, the function returns `NULL`.
#'
#' @examples
#' # connect to Oregon FIADB database
#' # or.db <- fia.connect( "c:/temp/OR.db" )
#' # fiadata <- get.fia( or.db, 1234567890 )
#'
#' @export

get_fia <- function( fiadb.con, stand.condition.code )
{

    stand.condition.code <- as.character( stand.condition.code )
    standinit <- dbGetQuery( fiadb.con,
                             paste( "select STAND_CN, STAND_ID, VARIANT, INV_DAY, INV_YEAR, INV_MONTH, LATITUDE, LONGITUDE, REGION, LOCATION,",
                                    "DISTRICT, ECOREGION, BASAL_AREA_FACTOR, INV_PLOT_SIZE, BRK_DBH,
                                    AGE, ASPECT, SLOPE, ELEVFT, NUM_PLOTS, MAX_sDI,",
                                    "SITE_SPECIES, SITE_INDEX, PHYSIO_REGION, STATE, COUNTY, FOREST_TYPE, FOREST_TYPE_FIA FROM",
                                    "FVS_STANDINIT_PLOT where STAND_CN = ", stand.condition.code ))

    if( nrow(standinit) == 0 )
        return( NULL )

    ba.factor <- ifelse( standinit$BASAL_AREA_FACTOR == 0, 1.0, 24.07219 )

    treeinit <- dbGetQuery( fiadb.con,
                            paste( "select STAND_CN, STAND_ID, PLOT_ID, TREE_ID,",
                                   "TREE_COUNT, SPECIES, DIAMETER, HT, CRRATIO, DEFECT_CUBIC, DEFECT_BOARD",
                                   "AGE, BH_YEARS FROM ",
                                   "FVS_TREEINIT_PLOT where STAND_CN = ", stand.condition.code )) %>%
        #filter( TREE_ID < 1000000) %>%
        mutate( AGE=ifelse(!is.na(AGE), AGE, standinit$AGE ),
                TPA=TREE_COUNT*ifelse( DIAMETER < standinit$BRK_DBH & standinit$BRK_DBH != 999,
                                       299.8611,
                                       ba.factor ),
                DEFECT_CUBIC=ifelse(!is.na(DEFECT_CUBIC),DEFECT_CUBIC,0.0))

    if( nrow(treeinit) == 0 )
        return( NULL )

    # compute growth year
    growth.year <- standinit$INV_YEAR + ifelse( standinit$INV_MONTH >= 7, 1, 0 )

    # construct list of plot condition data for use in growth
    plot.data <- list()

    plot.data[["stand"]] <- data.frame( PLT_CN=stand.condition.code,
                                        n.plots = standinit$NUM_PLOTS,
                                        site = standinit$SITE_INDEX,
                                        site.species = standinit$SITE_SPECIES,
                                        aspect = standinit$ASPECT, slope = standinit$SLOPE,
                                        elevation= standinit$ELEVFT,
                                        age = standinit$AGE,
                                        growth.year = growth.year,
                                        latitude = standinit$LATITUDE,
                                        longitude = standinit$LONGITUDE,
                                        national.forest = standinit$LOCATION,
                                        eco.region = standinit$ECOREGION ) %>%
                            mutate_at(vars(5:12), ~replace(., is.na(.), 0))

    fia.tree <- treeinit %>% mutate_at(vars(c("DIAMETER","HT","CRRATIO")), ~replace(., is.na(.), 0)) %>% mutate( CR=CRRATIO/100 )
    plot.data[["tree"]] <- fia.tree[,c("PLOT_ID","TREE_ID","AGE","SPECIES","DIAMETER","HT","TPA","CR")]

    return( plot.data )

}


#' grow_acd_fia() : grow an FIA tree list with Acadian Variant of FVS optionally using metric or imperial units.
#'
#' @param db         : SQLite DB connection
#' @param stand_cn   : FIA STAND_CN identifier
#' @param periods    : Number of periods to project the tree list
#' @param cdef       : Cumulative Defoliation Percentage (between 0 and 1) (default = 0)
#' @param use_sbw    : Use Spruce Budworm Modifiers (F = No (default)), T = Yes)
#' @param use_hw     : Use Hardwood Risk and Form Modifiers (F = No (default)), T = Yes)
#' @param use_thin   : Use Thinning Modifiers (F = No (default)), T = Yes)
#' @param use_ingrowth :  Use ingrowth functionality (F = No (default)), T = Yes)
#' @param cut_point  : Ingrowth probability threshold (default = 0.5)
#' @param MinDBH     : Minimum DBH for ingrowth trees (cm or inches) (default = 1.18in)
#' @param form       : Optional Form Code (Castle et al. 2017) 1 through 8 (default = NA)
#' @param risk       : Optional Risk Code (Castle et al. 2017) 1 through 4 (default = NA)
#' @param csi_default : Default climate site index used if site is missing (feet) (defaults to 52.5ft)
#'
#' @description
#' Extracts an FIA condition from an open FIA DB data base connection (see \link{fia_connect}). Projects the extracted
#' tree list using the Acadian Variant of FVS.
#'
#' @return
#' Returns a data.frame with the grown tree list in following variables (all values in imperial units):
#' \itemize{
#'     \item STAND_CN
#'     \item PLOT_ID
#'     \item TREE_ID
#'     \item SPECIES
#'     \item DBH
#'     \item HT
#'     \item TPA
#'     \item CR
#' }
#'
#' @examples
#' # open FIA DB connection
#' db <- fia_connect( "C:/Models/FIA/WV/SQLite_FIADB_WV.db" )
#'
#' # project the FIA condition 7340111010661 10 years
#' grown.plot <- grow_acd_fia( db, 7340111010661, 10 )
#'
#' @export

grow_acd_fia <- function( db, stand_cn, periods, cdef = 0,
                          use_sbw = F, use_hw = F, use_thin = F,
                          use_ingrowth = F, cut_point = 0.5, MinDBH = 1.18,
                          form = NA, risk = NA,
                          csi_default = 52.5 )
{
    ZERO <- as.integer(0)
    ONE <- as.integer(1)

    # get plot condition data from FIA DB
    plot.data <- acdR::get_fia( db, stand_cn )

    if( is.null(plot.data ) )
        return( NULL )

    # insure variable types are correct
    n <- nrow( plot.data$tree )
    use_sbw <- as.integer(use_sbw)
    use_hw <- as.integer(use_hw)
    use_thin <- as.integer(use_thin)
    if( !is.na(form) )
    {
        form <- as.integer(form)
    } else {
        form <- rep(ZERO,n)
    }
    if( !is.na(risk) )
    {
        risk <- as.integer(risk)
    } else {
        risk <- rep(ZERO,n)
    }


    # project the plot condition
    g <- grow_acd( as.integer(periods), "ME", as.integer(plot.data$stand$age), ONE,
                   ifelse(!is.na(plot.data$stand$site),plot.data$stand$site,csi_default),
                   plot.data$stand$elevation, cdef,
                   use_sbw, use_hw, use_thin, use_ingrowth, cut_point, MinDBH,
                   as.integer(plot.data$tree$PLOT_ID),
                   as.integer(plot.data$tree$TREE_ID),
                   as.integer(plot.data$tree$SPECIES),
                   plot.data$tree$DIAMETER, plot.data$tree$HT, plot.data$tree$TPA, plot.data$tree$CR, form, risk )

    if( nrow(g) == 0 )
        return( NULL )

    # build grown tree list data.frame
    results <- data.frame( STAND_CN = as.character(stand_cn),
                           PLOT_ID = g$plot.id,
                           TREE_ID = g$tree.id,
                           SPECIES = g$species,
                           DBH = g$dbh,
                           HT = g$ht,
                           CR = g$cr,
                           TPA = g$expf )

    return( results )
}
