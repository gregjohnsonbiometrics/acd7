// Interface between the Acadian Variant of FVS to R
//
// G.P. Johnson
// Greg Johnson Biometrics LLC
// (c) 2023, 2024
//

#include <Rcpp.h>

#include <string>
#include "stand.hpp"

//' grow_acd() : grow tree list with Acadian Variant of FVS optionally using metric or imperial units.
//'
//' @param periods    : Integer | Number of periods to project the tree list
//' @param Region     : String  | Region code: 'ME' (Maine), 'NB' (New Brunswick)
//' @param year       : Integer | Current year or age of stand
//' @param units      : Integer | 0 = metric, 1 = imperial
//' @param csi        : Numeric | Site Index (meters or feet)
//' @param elev       : Numeric | Elevation (meters or feet)
//' @param cdef       : Numeric | Cumulative Defoliation Percentage (between 0 and 1)
//' @param use_sbw    : Integer | Use Spruce Budworm Modifiers (0 = No, 1 = Yes)
//' @param use_hw     : Integer | Use Hardwood Risk and Form Modifiers (0 = No, 1 = Yes)
//' @param use_thin   : Integer | Use Thinning Modifiers (0 = No, 1 = Yes)
//' @param use_ingrowth : Integer | Use ingrowth functionality
//' @param cut_point  : Numeric | ingrowth probability threshold
//' @param MinDBH     : Numeric | Minimum DBH for ingrowth trees (cm or inches)
//' @param plot_id    : vector  | Plot ID
//' @param tree_id    : vector  | Tree ID
//' @param spp        : vector  | FIA Numeric Species Code
//' @param dbh        : vector  | diameter at breast height (cm or inches)
//' @param ht         : vector  | total height (meters or feet)
//' @param expf       : vector  | expansion factor, trees per hectare or per acre
//' @param cr         : vector  | crown ratio (between 0 and 1)
//' @param form       : vector  | Form Code (Castle et al. 2017) 1 through 8
//' @param risk       : vector  | Risk Code (Castle et al. 2017) 1 through 4
//' @md
//'
//' @description
//' Projects the tree list using the Acadian Variant of FVS with Metric Units.
//'
//' FIA species codes are defined here: https://www.fia.fs.usda.gov/program-features/urban/docs/Appendix%203a%20CORE%20Species%20List.pdf.
//'
//' @return
//' Returns a data.frame with the grown tree list in following variables:
//' \itemize{
//'     \item plot.id
//'     \item tree.id
//'     \item species
//'     \item dbh
//'     \item ht
//'     \item expf
//'     \item cr
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame grow_acd(
    Rcpp::IntegerVector  periods,
    Rcpp::StringVector   Region,
    Rcpp::IntegerVector  year,
    Rcpp::IntegerVector  units,
    Rcpp::NumericVector  csi,
    Rcpp::NumericVector  elev,
    Rcpp::NumericVector  cdef,
    Rcpp::IntegerVector  use_sbw,
    Rcpp::IntegerVector  use_hw,
    Rcpp::IntegerVector  use_thin,
    Rcpp::IntegerVector  use_ingrowth,
    Rcpp::NumericVector  cut_point,
    Rcpp::NumericVector  MinDBH,

    Rcpp::IntegerVector  plot_id,
    Rcpp::IntegerVector  tree_id,
    Rcpp::IntegerVector  spp,
    Rcpp::NumericVector  dbh,
    Rcpp::NumericVector  ht,
    Rcpp::NumericVector  expf,
    Rcpp::NumericVector  cr,
    Rcpp::IntegerVector  form,
    Rcpp::IntegerVector  risk )
{

    // perform sanity checks (all tree vectors must be the same size and greater than 0)
    if(  plot_id.size() > 0                &&
         plot_id.size() ==  tree_id.size() &&
         plot_id.size() ==  spp.size()     &&
         plot_id.size() ==  dbh.size()     &&
         plot_id.size() ==  ht.size()      &&
         plot_id.size() ==  expf.size()    &&
         plot_id.size() ==  cr.size()      &&
         plot_id.size() ==  form.size()    &&
         plot_id.size() ==  risk.size() )
    {
        std::streambuf* stderrbuf = std::cerr.rdbuf(Rcpp::Rcerr.rdbuf());

        const double ft_m =  (units[0] == 0) ? 1.0 : 0.3048;
        const double in_cm = (units[0] == 0) ? 1.0 : 2.54;
        const double ac_ha = (units[0] == 0) ? 1.0 : 2.47;

        try {
            // create STAND object
            STAND s( (std::string) Region[0], year[0], csi[0]*ft_m, elev[0]*ft_m, cdef[0],
                     use_sbw[0], use_hw[0], use_thin[0], use_ingrowth[0], cut_point[0], MinDBH[0]*in_cm );

            // create tree list
            for( R_xlen_t i = 0; i < plot_id.size(); i++ )
            {
                TREE t(  plot_id[i],  tree_id[i],  spp[i],  dbh[i]*in_cm,  ht[i]*ft_m,  expf[i]*ac_ha,  cr[i],  form[i],  risk[i] );
                s.trees.emplace_back( t );
            }

            // grow stand/plot
            s.grow( periods[0] );

            auto n = s.trees.size();

            std::vector<int> gplot_id(n, 0);
            std::vector<int> gtree_id(n, 0);
            std::vector<int> gspp(n, 0);
            std::vector<double> gdbh(n, 0.0);
            std::vector<double> ght(n, 0.0);
            std::vector<double> gcr(n, 0.0);
            std::vector<double> gexpf(n, 0.0);

            // load grown results back
            for( size_t i = 0; i < n; i++ )
            {
                auto &t = s.trees[i];
                if( t.expand_tree_id == 0 )
                {
                    gplot_id[i] = t.plot_id;
                    gtree_id[i] = t.tree_id;
                    gspp[i] = t.spp;
                    gdbh[i] = t.dbh / in_cm;
                    ght[i] = t.ht / ft_m;
                    gexpf[i] = t.tph / ac_ha;
                    gcr[i] = t.cr;
                }
            }

            std::cerr.rdbuf(stderrbuf);
                Rcpp::DataFrame grown_trees = Rcpp::DataFrame::create(
                Rcpp::Named("plot.id") = gplot_id,
                Rcpp::Named("tree.id") = gtree_id,
                Rcpp::Named("species") = gspp,
                Rcpp::Named("dbh") = gdbh,
                Rcpp::Named("ht") = ght,
                Rcpp::Named("expf") = gexpf,
                Rcpp::Named("cr") = gcr );

            return Rcpp::wrap(grown_trees);
        } catch(  std::exception &e ) {
            std::cerr << "Exception in growth()\n" << e.what() << "\n";
            std::cerr.rdbuf(stderrbuf);
            return Rcpp::wrap(0);
        }
    } else {
        return Rcpp::wrap(0);
    }

}


