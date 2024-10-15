// Parameter data structures
//
// G.P. Johnson
// Greg Johnson Biometrics LLC
// (c) 2024
//

#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <unordered_map>
#include <array>
#include <string>

struct SPP_CROSSWALK {
    std::string fvs_code;
    int mapped_code;
};

struct SPP_ID { 
    int spp_index;              // index into parameter arrays
    std::string spp_code;       // alpha species code
    bool softwood;              // softwood/hardwood flag
    std::string common_name;    // common species name
}; 

extern const std::unordered_map<int, SPP_CROSSWALK> species_crosswalk;

extern const std::unordered_map<int,SPP_ID> species_map;

const std::string get_common_name( int fia_species );
const int get_species_index( int fia_species );
 
struct CROWN_PARMS {
    double a1;
    double a2;
};

struct HTPRED_PARMS {
    double c0;
    double c3;
};

struct SPECIES_ATTRIB {
    double sg;
    double wd;
    double shade;
    double drought;
    double waterlog;
};

SPECIES_ATTRIB const *get_species_attrib( int species_index );

extern const std::array<SPECIES_ATTRIB,71> species_attrib;
extern const std::array<CROWN_PARMS,71> mcw_parms;
extern const std::array<CROWN_PARMS,71> lcw_parms;
extern const std::array<double,71> hcb_parms;
extern const std::array<double,6> hcb_fixed_parms;

extern const std::unordered_map<int, std::array<double,6>> htpred_parms;

extern const std::unordered_map<int,std::array<double,6>> ddbh_parms;

extern const std::unordered_map<int,std::array<double,6>> dht_parms;

extern const std::unordered_map<int,std::array<double,6>> dhcb_parms;

extern const std::unordered_map<int,std::array<double,5>> mort_parms;
#endif
