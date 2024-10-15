// Definition of tree class
//
// G.P. Johnson
// Greg Johnson Biometrics LLC
// (c) 2024
//

#ifndef TREE_HPP
#define TREE_HPP

#include <exception>
#include <iostream>
#include "parameters.hpp"

struct FORM_CLASS {
    double STM = 0.0;
    double LSW = 0.0;
    double MST = 0.0;
    double LF = 0.0;
};

class TREE {
    public:
        unsigned long plot_id;
        unsigned long tree_id;
        unsigned long expand_tree_id = 0;

        int    spp;     // FIA species code
        double dbh;     // diamter at breast height (cm)
        double ht;      // total height (meters) 
        double tph;     // trees per hectare represented by this tree
        double cr;      // crown ratio
        int Form;       // Form Code (Castle et al. 2017) : F1 : Ideal form,                        F2 : Acceptable form,           F3 : Poor form; large branches, 
                        //                                  F4 : Unacceptable form; multiple stems, F5 : Poor form; multiple stems, F6 : Poor form; single stem lean,
                        //                                  F7 : Acceptable significant fork,       F8 : Poor form; multiple stems
        int Risk;       // Risk Code (Castle et al. 2017) : R1 : Low Risk, R2 : Moderate Risk, R3 : High Risk, R4 : Very high risk
        double ba;      // tree basal area per hectare
        double bal;     // basal area in larger trees
        double bal_hw;  // hardwood basal area in larger trees
        double bal_sw;  // softwood basal area in larger trees
        double ccfl;    // crown competition factor in larger trees
        double ccfl_hw; // crown competition factor in larger hardwood trees
        double ccfl_sw; // crown competition factor in larger softwood trees

        double mcw;     // maximum crown width
        double lcw;     // largest crown width
        double mca;     // maximum crown area
        double hcb;     // height to crown base

        TREE( unsigned long plot_id_t, unsigned long tree_id_t, int spp_t, double dbh_t, double ht_t, double tph_t, double cr_t, int form_t, int risk_t );
        ~TREE() = default;

        // utilities
        void decode_form_and_risk();

        // accessor functions
        bool is_softwood() { return species_id.softwood; };

        double get_sg() { return _attributes->sg; };
        double get_shade() { return _attributes->shade; };
        double get_drought() { return _attributes->drought; };
        double get_wd() { return _attributes->wd; };
        double get_waterlog() { return _attributes->waterlog; };
        double get_survival() { return _p_survival; };

        void set_ddbh( double ddbh ) { _ddbh = ddbh; };
        void set_dht( double dht ) { _dht = dht; };
        void set_dtph( double dtph ) { _dtph = dtph; }; 
        double get_ddbh() { return _ddbh; };
        double get_dht() { return _dht; };
        double get_dhcb() { return _dhcb; };

        void compute_attributes();
        
        // crown functions
        void compute_mcw(); // maximum crown width
        void compute_lcw(); // largest crown width
        void compute_mca(); // maximum crown area

        // form and risk probabilities
        double risk_probability();
        FORM_CLASS form_probability();

        // imputation
        //void ht_pred( double csi, double ccf, bool override_ht = false );
        void ht_pred( double ccf, int region, bool override_ht = false );
        void hcb_pred( double ccf );

        // growth functions
        void dDBH( std::string Region, double csi, double sba, double percent_ba_removed, double ba_pre_thin, double qmd_ratio, 
                   int thin_year, int year, double average_dbh_sw, double topht, double CDEF );                   
        void dHT( std::string Region, double csi, double percent_ba_removed, double ba_pre_thin, double qmd_ratio, 
                  int thin_year, int year, double average_dbh_sw, double topht, double CDEF  );                  
        void dHCB( double ccf, double percent_ba_removed, double ba_pre_thin, double qmd_ratio, int thin_year, int year );
        bool survival_prob( std::string Region, double csi, double ba, double qmd, 
                            double percent_ba_removed, double ba_pre_thin, double qmd_ratio, int thin_year, int year,
                            double average_height_hw, double average_height_sw, double CDEF,
                            bool use_sbw_mod, bool use_hw_mod, bool use_thin_mod );

        // increment and apply change estimates to tree dimensions
        void apply_growth_mortality();
        void reset();

    private:
        int species_index = -1;
        SPP_ID species_id;
        bool _formB = false;
        bool _low_risk = true; 

        double _ddbh = 0.0;
        double _dht = 0.0;
        double _dhcb = 0.0;
        double _dtph = 0.0;
        double _p_survival = 1.0;

        const SPECIES_ATTRIB *_attributes;

        // pointers to parameter estimates
        std::array<double,6> ddbh_p = {0.0,0.0,0.0,0.0,0.0,0.0};
        std::array<double,6> dht_p = {0.0,0.0,0.0,0.0,0.0,0.0};

        const CROWN_PARMS  *mcw_p = nullptr;
        const CROWN_PARMS  *lcw_p = nullptr;
        //const HTPRED_PARMS *htpred_p = nullptr;
        std::array<double,6> htpred_p = {0.0,0.0,0.0,0.0,0.0,0.0};
        const double       *hcb_p = nullptr;
        std::array<double,6> dhcb_p = {0.0,0.0,0.0,0.0,0.0,0.0};
        std::array<double,5> mort_beta = {0.0,0.0,0.0,0.0,0.0};

        // modifiers
        double dDBH_thin( double percent_ba_removed, double ba_pre_thin, double qmd_ratio, int thin_year, int year );
        double dDBH_sbw( std::string Region, double average_dbh_sw, double topht, double CDEF = -1.0 );
        double dDBH_hw_form_risk();

        double dHT_thin( double percent_ba_removed, double ba_pre_thin, double qmd_ratio, int thin_year, int year );
        double dHT_sbw( double topht, double average_dbh_sw, double CDEF = -1.0 );

        double dHCB_thin( double percent_ba_removed, double ba_pre_thin, double qmd_ratio, int thin_year, int year );

        double surv_sbw( std::string Region, double average_height_hw, double average_height_sw, double CDEF );
        double surv_hw( double ba );
        double surv_thin( double percent_ba_removed, double ba_pre_thin, double qmd_ratio, int thin_year, int year );

};

#endif
