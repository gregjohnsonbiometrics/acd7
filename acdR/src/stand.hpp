// Definition of stand class
//
// G.P. Johnson
// Greg Johnson Biometrics LLC
// (c) 2024
//

#ifndef STAND_HPP
#define STAND_HPP

#include "tree.hpp"
#include <vector>
#include <string>
#include <exception>

enum class INGROWTH_MODEL_TYPE {
    GNLS = 0,
    NLME = 1
};

class STAND {

    public:
        std::string Region;     // Region code: ME (Maine), NB (New Brunswick)
        int    year;            // current year of the stand
        double csi;             // climate site index
        double ccf;             // crown competition factor
        double elevation;       // elevation (meters)
        double average_dbh;     // average dbh
        double average_dbh_10;  // average dbh for trees with dbh >= 10 cm
        double average_dbh_sw;  // average softwood dbh 
        double average_dbh_hw;  // average hardwood dbh
        double average_dbh_10_sw;  // average softwood dbh for trees with dbh >= 10 cm
        double average_dbh_10_hw;  // average hardwood dbh for trees with dbh >= 10 cm
        double dbh_sd;          // standard deviation of dbh
        double dbh_10_sd;       // standard deviation of dbh for trees with dbh >= 10 cm
        double average_height_sw;  // average softwood height
        double average_height_hw;  // average hardwood height
        double average_sg;      // average specific gravity
        double average_sg_10;   // average specific gravity  for trees with dbh >= 10 cm
        double CDEF = -1;       // cumulative defoliation %
        double topht;           // top height
        double ba_sw;           // softwood basal area per hectare
        double ba_hw;           // hardwood basal area per hectare
        double ba;              // basal area per hectare
        double tph;             // trees per hectare
        double qmd;             // quadratic mean diameter
        int    n_species;       // number of species in the tree list

        double min_dbh;         // minimum dbh
        double max_dbh;         // maximum dbh
        double min_dbh_10;      // minimum dbh
        double rd;              // relative density
        double rd_10;           // relative density for trees with dbh >= 10 cm
        double sdi_10;          // stand density index for trees with dbh >= 10 cm
        double sdi;             // stand density index for all trees
        double bf_ba;           // balsam fir basal area
        double ithw_ba;         // intolerant hardwood basal area

        bool use_sbw_mod  = false;
        bool use_hw_mod   = false;
        bool use_thin_mod = false;
        bool use_ingrowth = false;
        double cut_point  = 0.5;
        double MinDBH     = 0.0;

        bool initialized = false;

        // thinning data
        double percent_ba_removed = 0.0;
        double ba_pre_thin = 0.0;
        double qmd_ratio = 0.0;
        double thin_year = -1;

        std::vector<TREE> trees;

        void compute_ba_tph_bal();
        void compute_ccfl();
        void compute_ccf();
        void predict_hcb();
        void compute_tree_statistics();
        void compute_sdi_rd();
        void compute_topht();
        void compute_n_species();

        STAND( std::string t_Region, int t_year, double t_csi, double t_elev, double t_cdef, 
               bool use_sbw_t, bool use_hw_t, bool use_thin_t, 
               int t_use_ingrowth, double t_cut_point, double t_MinDBH ); 

        // initialize stand and tree variables, impute missing data
        void initialize();

        // grow the stand
        void grow( int n_years );

        // growth functions
        void diameter_growth();
        void height_growth();
        void crown_recession();
        void survival();

        // apply growth and mortality
        void apply_growth_mortality();

        // ingrowth 
        double ingrowth( INGROWTH_MODEL_TYPE type );

    private:

        unsigned long max_tree_id = 0ULL;
        unsigned long find_max_tree_id();

        bool expand_tree_list( double threshold );
        bool unexpand_tree_list();

        std::unordered_map<int,double> ba_spp;
        std::unordered_map<int,double> ba_grp_spp;
        std::unordered_map<int,std::unordered_map<int,double>> plot_species_ba; // map is plot: spp : ba
        void build_ba_spp_map();
        void ingrowth_composition( double IPH );

        void compute_survival_prob();
        double mort_sbw();
        double mort_thin();

};

#endif

