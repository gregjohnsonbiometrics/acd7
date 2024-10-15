// Definition of tree class
//
// G.P. Johnson
// Greg Johnson Biometrics LLC
// (c) 2024
//

#include "tree.hpp"
#include "utility.hpp"

#include <math.h>
#include <exception>
#include <unordered_map>

// TREE constructor
TREE::TREE( unsigned long plot_id_t, unsigned long tree_id_t, int spp_t, double dbh_t, double ht_t, double tph_t, double cr_t, int form_t, int risk_t ) :
    plot_id(plot_id_t), tree_id(tree_id_t), spp(spp_t), dbh(dbh_t), ht(ht_t), tph(tph_t), cr(cr_t), Form(form_t), Risk(risk_t)
{
    try {
        if( species_map.find( spp_t ) != species_map.end() )
        {
            // get species index
            species_id = species_map.at( spp_t );
            species_index = species_id.spp_index;

            // if species index is -1, use mapped species
            if( species_index == -1 )
            {
                species_index = species_map.at( species_crosswalk.at( spp_t ).mapped_code ).spp_index;
            }

        } else {
            // species not in map, try crosswalk
            auto cw_spp = species_crosswalk.at( spp_t ).mapped_code;
            species_id = species_map.at( cw_spp );
            species_index = species_map.at( cw_spp ).spp_index;
        }

        _attributes = get_species_attrib( species_index );

    } catch(const std::out_of_range &e) {
        std::cerr << "Exception in TREE(): species " << spp_t << " not found.\n";
        throw;
    } catch( ... ) {
        std::cerr << "Exception in TREE(): could not create tree.\n";
        throw;
    }

    bal = 0.0;
    bal_hw = 0.0;
    bal_sw = 0.0;
    ccfl = 0.0;
    hcb = cr > 0.0 && ht > 0.0 ? (1.0 - cr) * ht : 0.0;
 
    // get species index for other hardwoods or other softwoods depending on species
    auto si_other = get_species_index( is_softwood() ? 9991 : 9990 ); 

    // Get parameters. If a parameter has not been estimated for the species (or mapped species), assign other hardwood or softwood parameters
    if( ddbh_parms.find(spp_t) != ddbh_parms.end() )
        ddbh_p = ddbh_parms.at(spp_t);
    else
        ddbh_p = is_softwood() ? ddbh_parms.at(9991) : ddbh_parms.at(9990);

    if( dht_parms.find(spp_t) != dht_parms.end() )
        dht_p = dht_parms.at(spp_t);
    else
        dht_p = is_softwood() ? dht_parms.at(9991) : dht_parms.at(9990);

    if( dhcb_parms.find(spp_t) != dhcb_parms.end() )
        dhcb_p = dhcb_parms.at(spp_t);
    else
        dhcb_p = is_softwood() ? dhcb_parms.at(9991) : dhcb_parms.at(9990);

    if( htpred_parms.find(spp_t) != htpred_parms.end() )
        htpred_p = htpred_parms.at(spp_t);
    else
        htpred_p = is_softwood() ? htpred_parms.at(9991) : htpred_parms.at(9990);

    //htpred_p = &htpred_parms[species_index];
    hcb_p = hcb_parms[species_index] != 0.0 ? &hcb_parms[species_index] : &hcb_parms[si_other];
    mcw_p = mcw_parms[species_index].a1 != 0.0 ? &mcw_parms[species_index] : &mcw_parms[si_other];
    lcw_p = lcw_parms[species_index].a1 != 0.0 ? &lcw_parms[species_index] : &lcw_parms[si_other];
    if( mort_parms.find( spp_t ) != mort_parms.end() )
        mort_beta = mort_parms.at( spp_t );
    else
        mort_beta = is_softwood() ? mort_parms.at( 9991 ) : mort_parms.at( 9990 );

    compute_attributes();
    decode_form_and_risk();
}


//////////////////////////////////////////////////////////////////////
//  TREE class methods

// Reset all growth variables
void TREE::reset()
{
    _ddbh = 0.0;
    _dht = 0.0;
    _dhcb = 0.0;
    _dtph = 0.0;
    _p_survival = 1.0;
}

// Decode alphanumeric Castle et al. 2017 Form and Risk Codes
void TREE::decode_form_and_risk()
{
    // check for valid species, form and risk codes
    if( Form >= 1 && Form <= 8 && Risk >= 1 && Risk <= 4 )
    {
        // Convert NHRI form classes (Form 'A' == 0, 'B' == 1)
        _formB = !( Form == 1 || Form == 7 || Form == 3 || Form == 4 );

        // Convert NHRI risk classes (LR == 0, HR == 1)
        _low_risk = ( Risk == 1 || Risk == 2 );
    } else {
        _formB = false;
        _low_risk = true;
    }
}

// Compute tree attributes: tree ba, mcw, lcw, mca
void TREE::compute_attributes() 
{
    ba = dbh * dbh * 0.00007854 * tph;
    compute_mcw();
    compute_lcw();
    compute_mca();
}
 
// #### Crown prediction ####

// Maximum crown width
// For species without a specific parameter estimate, default to 999
void TREE::compute_mcw()
{
    try {
        mcw = mcw_p->a1 * std::pow( dbh, mcw_p->a2 );
    } catch( std::exception &e ) {
        std::cerr << "Exception in compute_mcw() for dbh " << dbh << "\n" << e.what() << "\n";
        throw;
    }
}

   
// Largest Crown Width
// For species without a specific parameter estimate, default to 999
void TREE::compute_lcw()
{
    try {
        lcw = mcw/(lcw_p->a1 * std::pow( dbh, lcw_p->a2 ));
    } catch( std::exception &e ) {
        std::cerr << "Exception in compute_lcw() for dbh " << dbh << "\n" << e.what() << "\n";
        throw;
    }
}   

// Maximum Crown Area
// Requires mcw to be populated
void TREE::compute_mca()
{
    constexpr double pi = 3.14159265358979323846;

    try {
        mca = 100.0*((pi*(mcw*mcw/4.0))/10000.0)*tph;
    } catch( std::exception &e ) {
        std::cerr << "Exception in compute_mca() for tph " << tph << ", and mcw " << mcw << "\n" << e.what() << "\n";
        throw;
    }
}   


// Height prediction 
// Total height prediction function (updated 8/31/2012) using species as a random effect
// If override_ht is true, all height will be replace with imputed values
// void TREE::ht_pred( double csi, double ccf, bool override_ht )
// {
//     auto &c = htpred_fixed_parms;
//     try {
//         if( ht <= 0.0 || override_ht )
//             ht = 1.37 + ((c[0]+htpred_p->c0) + std::pow(csi,c[1])) * pow( 1.0-std::exp(-c[2]*dbh), c[3]+htpred_p->c3 + c[4]*std::log(ccf+1.0)+c[5]*bal);
 
//     } catch( std::exception &e ) {
//         std::cerr << "Exception in ht_pred() for dbh " << dbh << ", bal " << bal << ", csi " << csi << ", ccf" << ccf << "\n" << e.what() << "\n";
//         throw;
//     }
// }

// Height prediction (Johnson 06/06/2024) 
// Total height prediction function 
// region is an indicator variable 0 = ME, 1 = NB
// If override_ht is true, all height will be replace with imputed values
void TREE::ht_pred( double ccf, int region, bool override_ht )
{
    try {
        if( ht <= 0.0 || override_ht )
            ht = 1.37 + (htpred_p[0]+htpred_p[1]*region) * 
                   std::pow( 1.0-std::exp( -htpred_p[2] * dbh - htpred_p[4] *(bal+1.0)), htpred_p[3] ) * std::pow(std::log(ccf),htpred_p[5]);
    } catch( std::exception &e ) {
        std::cerr << "Exception in ht_pred() for dbh " << dbh << ", bal " << bal << ", ccf" << ccf << "\n" << e.what() << "\n";
        throw;
    }
}


// Height to crown base prediction (updated 9/11/12 using species as random effect)
// Requires ccf and bal to be computed prior to call
// Recomputes crown ratio (cr)
void TREE::hcb_pred( double ccf ) 
{     
    auto &a = hcb_fixed_parms;

    double dhr = dbh / ht;

    try{
        hcb = ht/(1.0 + std::exp((a[0] + *hcb_p) + a[1]*dbh+ a[2]*ht + a[3]*dhr + a[4]*std::log(ccf+1.0) + a[5]*(bal+1.0)));
        cr = (ht - hcb)/ht;
    } catch( std::exception &e ) {
        std::cerr << "Exception in hcb_pred() for dbh " << dbh << ", bal " << bal << ", ht " << ht << ", ccf" << ccf << "\n" << e.what() << "\n";
        throw;
    }
}

// Revised diameter increment function from Greg Johnson (02/21/2024)
void TREE::dDBH( std::string Region, double csi, double sba, double percent_ba_removed, double ba_pre_thin, double qmd_ratio, 
                 int thin_year, int year, double average_dbh_sw, double topht, double CDEF )
{   
    _ddbh = 0.0;

    try {
        double tdbh = std::max(dbh,1.0);              
        _ddbh = std::exp( ddbh_p[0] +
                        ( ddbh_p[1] * std::log(tdbh + 1.0)) +
                        ( ddbh_p[2] * tdbh) + 
                        ( ddbh_p[3] * std::log(cr)) +
                        ( ddbh_p[4] * bal/std::log(tdbh + 1.0)) +
                        ( ddbh_p[5] * std::log(csi)));
    
        double thin_modifier = dDBH_thin( percent_ba_removed, ba_pre_thin, qmd_ratio, thin_year, year );
        double sbw_modifier = dDBH_sbw( Region, average_dbh_sw, topht, CDEF );
        double hw_risk_modifier = dDBH_hw_form_risk();

        _ddbh *= thin_modifier * sbw_modifier * hw_risk_modifier;

    } catch( std::exception &e ) {
            std::cerr << "Exception in dDBH() for dbh " << dbh << ", bal " << bal << ", cr " << cr << ", csi " << csi << "\n" << e.what() << "\n";
            throw;
    }
}

// Diameter thinning modifier function (3/15/16 based on results by Christian Kuehne)
double TREE::dDBH_thin( double percent_ba_removed, 
                        double ba_pre_thin, 
                        double qmd_ratio, 
                        int    thin_year, 
                        int    year )
{
    double dbh_thin_modifier = 1.0;

    // there must be a valid thin_year (non-negative, occuring in the past)
    if( thin_year >= 0 && thin_year <= year && percent_ba_removed > 0.0 && qmd_ratio > 0.0 && ba_pre_thin > 0.0 )
    {
        int time_since_thinning = year - thin_year; 

        try {
            if( spp == 12 )
            {
                // balsam fir (FIA spp == 12)
                double y0 = -0.2566;
                double y1 = -22.7609;
                double y2 =  0.7745;
                double y3 =  1.0511;
                dbh_thin_modifier = 1.0 + ( std::exp(y0+y1/((100.0*percent_ba_removed*qmd_ratio)+0.01)) * 
                                                std::pow(y2,time_since_thinning) * 
                                                std::pow(time_since_thinning,y3) );
            } else if( spp == 97 ) {

                // red spruce 
                double y0 =  -0.5010;
                double y1 = -20.1147;
                double y2 =   0.8067;
                double y3 =   1.1905;
                dbh_thin_modifier = 1.0 + ( std::exp(y0+y1/((100.0*percent_ba_removed*qmd_ratio)+0.01)) * 
                                                std::pow(y2,time_since_thinning) * 
                                                std::pow(time_since_thinning,y3) );
            }

            // constrain thinning modifer to: 0.75 <= modifier <= 1.25
            dbh_thin_modifier = std::max( 0.75, std::min( dbh_thin_modifier, 1.25 ) );

        } catch( std::exception &e ) {
            std::cerr << "Exception in dDBH_thin() for percent_ba_removed " << percent_ba_removed <<
                         ", qmd_ratio " << qmd_ratio << ", time_since_thinning " << time_since_thinning << "\n" << e.what() << "\n";
            throw;
        }

    } else {
        dbh_thin_modifier = 1.0;   
    }

    return dbh_thin_modifier;
}  

// Diameter increment modifier for SBW (Cen et al. 2016)
// This is expecting average dbh sw for 10+ cm trees
double TREE::dDBH_sbw( std::string Region, double average_dbh_sw, double topht, double CDEF )
{
    double dbh_sbw_modifier = 1.0;

    if( (spp == 12 || spp == 94 || spp == 95 || spp == 97) && CDEF >= 0.0 )
    {
        double b1 = 0.0;
        double b2 = 0.0;
        double b3 = 0.0;
        double b4 = 0.0;
        double b5 = 0.0;
        double b6 = 0.0;
        double b7 = 0.0;

        if( Region == "ME" )
        {
            b2 =  0.0019;
            b3 = -0.0327;
            b4 = -0.0412;
            b5 =  0.3950;

            if( spp == 12)
            {
                // balsam fir
                b1 =  0.1187;
                b6 = -1.2813;
                b7 = -0.0016;
            } else if( spp == 97 || spp == 95 ) { 
                // red spruce or black spruce
                b1 =  0.0675;
                b6 = -0.9477;
                b7 = -0.0006;
            } else if( spp == 94 ) {
                // white spruce
                b1 =  0.0321;
                b6 = -0.3715;
                b7 = -0.0183; 
            }
        }              
        else if( Region== "NB" )
        {                                                                     
            b2 = -0.0190;
            b3 = -0.0277;
            b4 = -0.0027;
            b5 =  0.0000;
            if( spp == 12)
            {
                // balsam fir
                b1 =  0.0701;
                b6 = -0.8200;
                b7 = -0.0018;
            } else if( spp == 97 || spp == 95 ) { 
                // red spruce or black spruce
                b1 =  0.0320;
                b6 = -0.6861;
                b7 = -0.0012;
            } else if( spp == 94 ) {
                // white spruce
                b1 =  0.0487;
                b6 = -0.7839;
                b7 = -0.0006; 
            }
        }

        try {
            double dDBHa = b1*dbh * std::exp( b2*bal_hw + b3*bal_sw + b4*topht+b5*cr + b6*(dbh/average_dbh_sw) );

            double dDBHb = b1*dbh * std::exp( b2*bal_hw + b3*bal_sw + b4*topht + b5*cr + b6*(dbh/average_dbh_sw) + b7*CDEF );

            dbh_sbw_modifier = dDBHb / dDBHa;
        } catch( std::exception &e ) {
            std::cerr << "Exception in dDBH_sbw() for dbh " << dbh << ", bal_sw " << bal_sw << ", bal_hw " << bal_hw << ", topht " << topht << 
                         ", average_dbh_sw " << average_dbh_sw << ", CDEF " << CDEF << "\n" << e.what() << "\n";
            throw;
        }
    }

    return dbh_sbw_modifier;
}


// Diameter increment HW form and risk modifier
double TREE::dDBH_hw_form_risk()
{  
     double dbh_hw_form_risk_modifier = 1.0;

    // check for valid species, form and risk codes
    if( Form >= 1 && Form <= 8 && Risk >= 1 && Risk <= 4 &&
        /*       RO             YB            RM            PB            QA */
        ( spp == 833 || spp == 371 || spp == 316 || spp == 375 || spp == 746 ) ) // removed SM from the test spp == 318 ||
    {
        double b0 =  -2.9487;
        double b1 =  -0.1090;
        double b2 =   1.2111;
        double b3 =  -0.0430;
        double b4 =   0.0000;
        double b5 =   0.0000;
        double b6a =  0.2176;
        double b6b = -0.0250*_formB + 0.2176*_low_risk; // Form == B and Risk = Low

        switch( spp ) {
            case 746: b4 = -0.1059; // QA
                      break;
            case 316: b4 = -0.6377; // RM
                      break;
            case 833: b4 = -0.3453; // RO
                      break;
            case 371: b4 = -0.2494; // YB
                      break;  
            case 375: b4 =  0.0000; // PB
                      break;                                          
        };

        switch( spp ) {
            case 746: b5 = 0.0476; // QA
                      break;
            case 316: b5 = 0.0477; // RM
                      break;
            case 833: b5 = 0.0511; // RO
                      break;
            case 371: b5 = 0.0251; // YB
                      break; 
            case 375: b5 = 0.0000; // PB
                      break;                                                                  
        };

        try {
            double a = b0 + b1*dbh + b2*std::log(dbh) + b3*bal + b4 + b5*dbh;

            double dDBH_a = std::exp( a + b6a );
            double dDBH_b = std::exp( a + b6b );

            dbh_hw_form_risk_modifier = dDBH_b / dDBH_a;
        } catch( std::exception &e ) {
            std::cerr << "Exception in dDBH_hw_form_risk() for dbh " << dbh << ", bal " << bal << "\n" << e.what() << "\n";
            throw;
        }
    } else {
        // invalid species, form or risk codes
       
    }

    return dbh_hw_form_risk_modifier;
}  

// Estimate Mortality
// G. Johnson alternative individual tree mortality/survival estimator 06/19/2024 
// Currently estimated for 18 species and all conifers, and all hardwoods
bool TREE::survival_prob( std::string Region, double csi, double ba, double qmd,
                          double percent_ba_removed, double ba_pre_thin, double qmd_ratio, int thin_year, int year,
                          double average_height_hw, double average_height_sw, double CDEF,
                          bool use_sbw_mod, bool use_hw_mod, bool use_thin_mod )
{  
    try {
        _p_survival = 1.0 - std::exp( -std::exp( -mort_beta[0] + 
                                                  mort_beta[1]*(std::pow(dbh,mort_beta[2])/(bal+1.0) )));

        // compute modifiers
        double sbw_modifier = 1.0;
        double hw_modifier = 1.0;
        double thin_modifier = 1.0;
        if( use_sbw_mod )
            sbw_modifier = surv_sbw( Region, average_height_hw, average_height_sw, CDEF );

        if( use_hw_mod )
            hw_modifier = surv_hw( ba );

        if( use_thin_mod )
            thin_modifier = surv_thin( percent_ba_removed, ba_pre_thin, qmd_ratio, thin_year, year );

        _p_survival *= sbw_modifier * (1.0/thin_modifier) * hw_modifier;
 
        return true;

    } catch( const std::exception &e ) {
        std::cerr << "Exception in survival_prob() for dbh " << dbh << ", bal " << bal << ", shade " << get_shade() << "\n" << e.what() << "\n";
        throw;
    }

    return false;
}

// SBW survival modifier: tree record (Cen et al. 2016).
// Estimate modifier and apply to probability of survival
//
// Cen Chen, Aaron Weiskittel, Mohammad Bataineh, and David A. MacLean. 2017. Even low levels of spruce budworm defoliation affect mortality and ingrowth but net growth is more driven by competition. Canadian Journal of Forest Research. 47(11): 1546-1556. https://doi.org/10.1139/cjfr-2017-0012
// C Chen, A Weiskittel, M Bataineh, DA MacLean. 2017. Evaluating the influence of varying levels of spruce budworm defoliation on annualized individual tree growth and mortality in Maine, USA and New Brunswick, Canada
// Forest Ecology and Management 396:184-194. https://doi.org/10.1016/j.foreco.2017.03.026 . 
double TREE::surv_sbw( std::string Region,
                       double average_height_hw,
                       double average_height_sw,
                       double CDEF )
{
    double surv_sbw_modifier = 1.0;

    if( (spp == 12 || spp == 94 || spp == 95 || spp == 97) && CDEF >= 0.0 )
    {
        double b1 = 0.0;
        double b2 = 0.0;
        double b3 = 0.0;
        double b4 = 0.0;
        double b5 = 0.0;
        double b6 = 0.0;
        double b7 = 0.0;
        double b8 = 0.0;

        if( Region == "ME" )
        {
            b1 = -6.5208;
            b2 = -0.4866;
            b4 =  0.0316;
            b6 = -0.0175;
            b7 =  0.0274;

            if( spp == 12)
            {
                // balsam fir
                b3 = -0.0355; 
                b5 =  1.5087;
                b8 =  0.0040;
            } else if( spp == 97 || spp == 95 ) { 
                // red spruce or black spruce
                b3 = -0.1231;
                b5 =  1.5087;
                b8 =  0.0056;
            } else if( spp == 94 ) {
                // white spruce
                b3 = -0.1755; 
                b5 =  1.5087;
                b8 =  0.0207;
            }
        }
        else if( Region == "NB" )
        {
            b1 = -6.8310;
            b2 =  0.0000;
            b4 =  0.2025;
            b6 =  0.0000;
            b7 =  0.0000;
            if( spp == 12)
            {
                // balsam fir
                b3 = -0.2285;
                b5 =  2.1703;
                b8 =  0.0029;
            } else if( spp == 97 || spp == 95 ) { 
                // red spruce or black spruce
                b3 = -0.2285;
                b5 =  2.0809;
                b8 =  0.0101;
            } else if( spp == 94 ) {
                // white spruce
                b3 = -0.2285;
                b5 =  1.5802;
                b8 =  0.0021;
            }
        }

        try {
            double x = b1 + b2*cr + b3*dbh + b4*average_height_sw + b5*(ht/average_height_sw) + b6*bal_sw + b7*bal_hw;
            double mort_a = (1.0 - std::exp(-std::exp( x )));
            double mort_b = (1.0 - std::exp(-std::exp( x + b8*CDEF )));
            surv_sbw_modifier = (mort_a > 0.0) ? (1.0 - mort_b)/(1.0 - mort_a) : 1.0;
        } catch( std::exception &e ) {
            std::cerr << "Exception in surv_sbw() for dbh " << dbh << ", ht " << ht << ", cr " << cr << ", average_height_sw" << average_height_sw <<
                         ", average_height_hw" << average_height_hw << ", bal_sw " << bal_sw << ", bal_hw " << bal_hw << ", CDEF " << CDEF << 
                         "\n" << e.what() << "\n";
            throw;
        }
    }

    // apply modifier (constrain to (0,1))
    return (surv_sbw_modifier <= 1.0) ? surv_sbw_modifier : 1.0;
}


// HW survival modifier: tree record from Castle et al. (2017). Adjust survival estimate for HW form and risk
double TREE::surv_hw( double ba )
{
    double surv_hw_modifier = 1.0;

    enum class NHRIFORM {
        STM,
        SWP,
        MST,
        OTHER
    };

    // check for valid species and form codes
    if( (Form >= 1 && Form <= 8 ) &&
          //     'RO',         'YB',         'RM',         'PB',          'QA'
        ( spp == 833 || spp == 371 || spp == 316 || spp == 375 || spp == 746 ) ) // removed SM from the test spp == 318 || 
    {
        double b0 = 15.1991;
        double b1 = -0.1509;
        double b2 = -0.1232;
        double b3 = -1.4053;
        double b4 =  3.3082;
        double b5 =  0.0000;
        double b6 =  0.0000;

        // Convert NHRI form classes
        switch( Form ) {
            case 1: b5 = 3.3082; // STM
                    break;
            case 2: b5 = 2.2518; // SWP
                    break;
            case 5: b5 = 0.0000; // MST;
                    break;
            case 8: b5 = 0.0000; // OTHER;
                    break;                    
        }

        // assign random effects by species
        if( spp == 746 )
        {   
            // QA
            b4 = -2.7907;
            b6 =  0.0791;
        } else if( spp == 316 ) {
            // RM
            b4 = -3.9809;
            b6 =  0.8343;
        } else if( spp == 833 ) {
            // RO
            b4 = -0.7937;
            b6 =  0.8944;
        } else if( spp == 371 ) {
            // YB
            b4 = 5.2531;
            b6 = 0.1528;
        }

        try {
            double x = b0 + b1*dbh + b2*bal + b3*std::sqrt(ba) + b4 + b6*dbh;

            double mort_a = std::exp( x ) / ( 1.0 + std::exp( x ) );
            double mort_b = std::exp( x + b5 ) / ( 1.0 + std::exp( x + b5 ) );

            surv_hw_modifier = (mort_a != 0.0 ) ? mort_b / mort_a : 1.0;

        } catch( std::exception &e ) {
            std::cerr << "Exception in surv_hw() for dbh " << dbh << ", bal " << bal << ", ba " << ba << "\n" << e.what() << "\n";
            throw;
        }
    }

    // apply modifier (constrain to (0,1))
    return (surv_hw_modifier <= 1.0) ?  surv_hw_modifier : 1.0;
}


// Thinning survival modifier: tree record. Adjust survival estimate for thinning.
double TREE::surv_thin( double percent_ba_removed, 
                      double ba_pre_thin, 
                      double qmd_ratio, 
                      int thin_year, 
                      int year )
{

    double surv_thin_modifier = 1.0;

    if( thin_year >= 0 && thin_year <= year && percent_ba_removed > 0.0 && qmd_ratio > 0.0 && ba_pre_thin > 0.0 )
    {
        int time_since_thinning = year - thin_year; 
      
        try {
            if( spp == 12 ) 
            {
                // balsam fir 
                double y0 = 1.7414;
                double y1 = 7.0805;
                double y2 = 0.6677; 
                double y3 = 0.8474;
                surv_thin_modifier = 1.0 + (std::exp( y0 + (y1/(((100.0*percent_ba_removed+ba_pre_thin)*qmd_ratio)+0.01))) *
                                    std::pow(y2,time_since_thinning) * std::pow(time_since_thinning,y3) );
            } else if( spp == 97 ) {
                // red spruce
                double y0=   10.5057;
                double y1= -650.8260;
                double y2=    0.6948; 
                double y3=    0.6429;
                surv_thin_modifier = 1.0 + (std::exp( y0 + (y1/(((100.0*percent_ba_removed)+ba_pre_thin)+0.01))) * 
                                    std::pow(y2,time_since_thinning) * std::pow(time_since_thinning,y3) );
            }
        } catch( std::exception &e ) {
            std::cerr << "Exception in surv_thin() for percent_ba_removed " << percent_ba_removed <<
                         ", qmd_ratio " << qmd_ratio << ", time_since_thinning " << time_since_thinning << "\n" << e.what() << "\n";
            throw;
        }
    }

    // apply modifier (constrain modifier to (0,1))
    double x = 1.0 / surv_thin_modifier;
    return (x <= 1.0) ? x : 1.0;
}

// apply growth to current tree and update crown dimensions
// this function resets the tree growth variables
void TREE::apply_growth_mortality()
{
    try {  
        // increment tree dimensions
        dbh += _ddbh;
        ht += _dht; 
        hcb += _dhcb;
        if( hcb > ht ) hcb = ht;
        cr = (ht - hcb)/ht;
        tph -= (_dtph <= tph) ? _dtph : tph;

        ba = dbh * dbh * 0.00007854 * tph;

        // update crown attributes
        compute_mcw();
        compute_lcw();
        compute_mca();
        
    } catch( ... ) {
        throw;
    }

    // clean up
    reset();
}

// Alternative dHT equation from Greg Johnson (02/26/2024)
void TREE::dHT( std::string Region, double csi, double percent_ba_removed, double ba_pre_thin, double qmd_ratio, 
                int thin_year, int year, double average_dbh_sw, double topht, double CDEF )
{
    // model form by Greg Johnson
    _dht = 0.0;
    try {
        try {
            _dht = dht_p[0] * dht_p[1] * dht_p[2] * std::pow( cr, dht_p[5] ) * std::pow( csi/30.0, dht_p[5] ) * 
                    std::exp( -dht_p[1] * ht - dht_p[4] * (ccfl/100.0) ) * std::pow((1.0 - std::exp(-dht_p[1]*ht)),dht_p[2]-1.0);
        } catch( std::exception &e ) {
            std::cerr << "Exception in dHT() for ht " << ht << ", cr " << cr << ",ccfl " << ccfl << ", csi " << csi << "\n" << e.what() << "\n";
            throw;
        }

        double thin_modifier = dHT_thin( percent_ba_removed, ba_pre_thin, qmd_ratio, thin_year, year );
        double sbw_modifier = dHT_sbw( topht, average_dbh_sw, CDEF );

        _dht *= thin_modifier * sbw_modifier;

    } catch( ... ) {
        throw;
    }    
}


// Thinning height modifier (Kuehne et al. 2016)
// thin_year (the year thinning occurred) is negative for no thinning
// Operative for balsam fir and red spruce
double TREE::dHT_thin( double percent_ba_removed, 
                       double ba_pre_thin, 
                       double qmd_ratio, 
                       int    thin_year, 
                       int    year )
{
    double height_thin_modifier = 1.0;

    // there must be a valid thin_year (non-negative, occuring in the past, and less than 5 years ago)
    if( thin_year >= 0 && thin_year <= year && (year - thin_year) < 5 )
    {
        int time_since_thinning = year - thin_year; 

        try {
            if( spp == 12 )
            {
                // balsam fir (FIA spp == 12)
                double y0 = -1.8443;
                double y1 =  5.2969;
                double y2 =  1.0532;
                double y3 =  0.0000;
                height_thin_modifier = 1.0 - ( std::exp(y0+y1/((100.0*percent_ba_removed)+0.01)) * 
                                                std::pow(y2,time_since_thinning) * 
                                                std::pow(time_since_thinning,y3) );
            } else if( spp == 97 ) {

                // red spruce 
                double y0 = -1.8426;
                double y1 =  6.2781;
                double y2 =  1.1596;
                double y3 =  0.0000;
                height_thin_modifier = 1.0 -( std::exp(y0+y1/((100.0*percent_ba_removed)+0.01)) *
                                            std::pow(y2,time_since_thinning) * 
                                            std::pow(time_since_thinning,y3) );
            }

            // constrain thinning modifer to: 0.75 <= modifier <= 1.25
            height_thin_modifier = std::max( 0.75, std::min( height_thin_modifier, 1.25 ) );
        } catch( std::exception &e ) {
            std::cerr << "Exception in dHT_thin() for percent_ba_removed" << percent_ba_removed <<
                         ", qmd_ratio " << qmd_ratio << ", time_since_thinning " << time_since_thinning << "\n" << e.what() << "\n";

            throw;
        }

    } else {
        height_thin_modifier = 1.0;   
    }

    return height_thin_modifier;
}


// SBW Height modifier (Cen et al. 2016)
// CDEF of -1 is missing or not supplied, no modifier will be computed
// This is expecting average dbh sw for 10+ cm trees
double TREE::dHT_sbw( double topht,
                      double average_dbh_sw,
                      double CDEF )
{
    double height_sbw_modifier = 1.0;

    if( (spp == 12 || spp == 94 || spp == 95 || spp == 97) && CDEF >= 0.0 )
    {
        double b1 =  0.0000;
        double b2 = -0.0011;
        double b3 =  0.0316;
        double b4 =  2.4512;
        double b5 =  0.0000;
        double b6 =  0.0000;

        if( spp == 12 )
        {
            // balsam fir
            b1 =  0.0013; 
            b5 =  0.3676; 
            b6 = -0.0017; 
        } else if( spp == 97 || spp == 95 ) { 
            // red spruce or black spruce
            b1 =  0.0009; 
            b5 =  0.2881; 
            b6 = -0.0014;
        } else if( spp == 94 ) {
            // white spruce
            b1 = 0.0005;
            b5 = 0.6800; 
            b6 = 0.0001;
        }

        try {
            double dHTa = b1 * dbh * std::exp( b2*dbh*dbh + b3*topht + b4*cr + b5 *(dbh/average_dbh_sw) );

            double dHTb = b1 * dbh * std::exp( b2*dbh*dbh + b3*topht + b4*cr + b5*(dbh/average_dbh_sw) + b6*CDEF);

            height_sbw_modifier = dHTb / dHTa;
        } catch( std::exception &e ) {
             std::cerr << "Exception in dHT_sbw() for dbh " << dbh << ", topht " << topht << ", cr " << cr <<
                          ", average_dbh_sw " << average_dbh_sw << ", CDEF " << CDEF << "\n" << e.what() << "\n";
        }
    }

    return height_sbw_modifier;
}

// ### Crown recession ####
// dHCB equation using late-stage hcb ratio: G. Johnson 07/2024
void TREE::dHCB( double ccf, 
                 double percent_ba_removed, 
                 double ba_pre_thin, 
                 double qmd_ratio, 
                 int    thin_year, 
                 int    year )
{
    try {
        try {
            //_dhcb = dhcb_p[0] * ((ht - hcb) + std::pow(_dht, dhcb_p[1] ));
            _dhcb = dhcb_p[0] * std::pow( hcb/dhcb_p[5], dhcb_p[2] ) * 
                    ( (ht - hcb) + std::pow( _dht, dhcb_p[1] ) ) *
                    std::pow( (1.0 - std::exp( -dhcb_p[3]*(ccf+1.0) )), dhcb_p[4] );
                    
            double thin_modifier = dHCB_thin( percent_ba_removed, ba_pre_thin, qmd_ratio, thin_year, year );

            _dhcb *= thin_modifier;
        } catch( std::exception &e ) {
             std::cerr << "Exception in dHCB() for hcb " << hcb << ", dht " << _dht << ", ccf " << ccf <<
                          ", shade " << _attributes->shade << "\n" << e.what() << "\n";
        }
    } catch( ... ) {
        throw;
    }
}

// Crown recession thinning modifier
// thin_year (the year thinning occurred) is negative for no thinning
// Operative for balsam fir and red spruce
double TREE::dHCB_thin( double percent_ba_removed, 
                        double ba_pre_thin, 
                        double qmd_ratio, 
                        int    thin_year, 
                        int    year )
{
    double hcb_thin_modifier = 1.0;

    // there must be a valid thin_year (non-negative and occuring in the past) for balsam fir and red spruce
    if( thin_year >= 0 && thin_year <= year && (spp == 12 || spp == 97) )
    {
        int time_since_thinning = year - thin_year; 

        try {
            if( spp == 12 )
            {
                // balsam fir (FIA spp == 12)
                double y0 =  -0.4208;
                double y1 = -17.0998;
                double y2 =   0.7986;
                double y3 =   0.0521;
                hcb_thin_modifier = 1.0 - (std::exp(y0+y1/((100.0*percent_ba_removed*qmd_ratio)+0.01)) *
                                        std::pow(y2,time_since_thinning) * 
                                        std::pow(time_since_thinning,y3));
            } else if( spp == 97 ) {

                // red spruce 
                double y0 =  -1.0778;
                double y1 = -14.7694;
                double y2 =   0.7758;
                double y3 =   1.1164;
                hcb_thin_modifier = 1.0 - (std::exp(y0+y1/((100.0*percent_ba_removed*qmd_ratio)+0.01)) *
                                        std::pow(y2,time_since_thinning) * 
                                        std::pow(time_since_thinning,y3));
            }

            // constrain hcb modifer to: modifier <= 1.00
            hcb_thin_modifier = std::min( std::abs(hcb_thin_modifier), 1.0 );
        } catch( std::exception &e ) {
            std::cerr << "Exception in dHCB_thin() for percent_ba_removed" << percent_ba_removed <<
                         ", qmd_ratio " << qmd_ratio << ", time_since_thinning " << time_since_thinning << "\n" << e.what() << "\n";
            throw;
        }
    } else {
        // nothing to do here 
    }

    return hcb_thin_modifier;
}


// Risk Classification from Castle et al. (2017; CJFR 47: 1457-1467)
// Returns the probability of a tree being high risk
double TREE::risk_probability()
{
    double b0 = -0.6886;
    double b1 = -0.0001;
    double b2 =  0.0000;
    double b3 =  0.0000;

    double HR = 0.0;

    //          RM            RO            SM            YB
    if( spp == 316 || spp == 833 || spp == 318 || spp == 371 )
    {
        // Red maple (Acer rubrum L.) was used as the reference level with which other species were compared.
        if( spp == 833 )
        {
            // red oak
            b2 = -0.0184;
            b3 = -0.0393;
        } else if( spp == 318  ) {
            // sugar maple
            b2 = -0.1513;
            b3 = -0.0164;        
        } else if( spp == 371 ) {
            // yellow birch
            b2 = -0.9851;
            b3 =  0.0196;
        }

        try {
            HR = std::exp(b0 + b1*dbh + b2 + b3*dbh ) / (1.0 + std::exp(b0 + b1*dbh + b2 + b3*dbh ));
        } catch( std::exception &e ) {
            std::cerr << "Exception in risk_probability() for dbh " << dbh << "\n" << e.what() << "\n";
            throw;
        }
    }

    return HR;
}

// Form Classification from Castle et al. (2017; CJFR 47: 1457-1467)
// Returns the probability of single straight stem (STM), 
// extensive sweep and lean (LSW), multiple stems (MST),
// significant fork on first 5 m (LF)
FORM_CLASS TREE::form_probability()
{
    FORM_CLASS t;

    double b0_stm = -0.9491;
    double b1_stm =  0.0174;
    double b2_stm =  0.0000;

    double b0_lsw = -1.1143;
    double b1_lsw = -0.0322;
    double b2_lsw =  0.0000;

    double b0_mst = -0.4110;
    double b2_mst =  0.0000;

    double b0_lf = -4.0677;
    double b1_lf =  0.0322;
    double b2_lf =  0.0000;

    // Red maple (Acer rubrum L.) was used as the reference level with which other species were compared.
    if( spp == 833 )
    {
        // red oak
        b2_stm = -0.2826;
        b2_lsw =  0.7910;
        b2_mst = -0.5009;
        b2_lf  =  0.1139;
    } else if( spp == 318 ) {
        // sugar maple
        b2_stm =  0.7541;
        b2_lsw = -0.2325;
        b2_mst = -1.1347;
        b2_lf  =  0.6278;
    } else if( spp == 371 ) {
        // yellow birch
        b2_stm = -0.0208;
        b2_lsw =  0.2980;
        b2_mst = -0.7557;
        b2_lf  =  1.0681;
    }

    //          RM            RO            SM            YB
    if( spp == 316 || spp == 833 || spp == 318 || spp == 371 )
    {
        try {
            t.STM = std::exp( b0_stm + b1_stm*dbh + b2_stm ) / (1.0 + std::exp( b0_stm + b1_stm*dbh + b2_stm ));
            t.LSW = std::exp( b0_lsw + b1_lsw*dbh + b2_lsw ) / (1.0 + std::exp( b0_lsw + b1_lsw*dbh + b2_lsw ));
            t.MST = std::exp( b0_mst              + b2_mst ) / (1.0 + std::exp( b0_mst              + b2_mst ));
            t.LF  = std::exp( b0_lf  + b1_lf*dbh  + b2_lf  ) / (1.0 + std::exp( b0_lf  + b1_lf*dbh  + b2_lf ));

            double xx = 1.0 / (t.STM + t.LSW + t.MST + t.LF);
            t.STM *= xx;
            t.LSW *= xx;
            t.MST *= xx;
            t.LF  *= xx;
        } catch( std::exception &e ) {
            std::cerr << "Exception in form_probability() for dbh " << dbh << "\n" << e.what() << "\n";
            throw;
        }
    }

    return t;
}

