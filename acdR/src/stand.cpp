// Implementation of stand class
//
// G.P. Johnson
// Greg Johnson Biometrics LLC
// (c) 2024
//
#include "stand.hpp"

#include <vector>
#include <set>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
#include <random>       // std::default_random_engine, std::uniform_real_distribution

#include <math.h>
#include <exception>

// sorting machinery for generated indices to decreasing sorted vector by dbh
template <typename T>
std::vector<size_t> sort_indices(const std::vector<T> &v, bool use_dbh ) 
{
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  std::stable_sort(idx.begin(), idx.end(),
       [&v, use_dbh](size_t i1, size_t i2) { bool t = use_dbh ? v[i1].dbh >= v[i2].dbh : v[i1].ht >= v[i2].ht;                                    
                                    return t; });

  return idx;
}

STAND::STAND( std::string t_Region, int t_year, double t_csi, double t_elev, double t_cdef, 
        bool use_sbw_t, bool use_hw_t, bool use_thin_t, 
        int t_use_ingrowth, double t_cut_point, double t_MinDBH ) : 
        Region(t_Region), year(t_year), csi(t_csi), elevation(t_elev), CDEF(t_cdef), 
        use_sbw_mod(use_sbw_t), use_hw_mod(use_hw_t), use_thin_mod(use_thin_t),
        use_ingrowth( t_use_ingrowth ), cut_point( t_cut_point), MinDBH( t_MinDBH)
{
    if( Region != "ME" && Region != "NB")
    {
        std::cerr << "Invalid Region (" << Region << ")\n";
        throw std::invalid_argument( "Invalid Region");
    }

    if( csi <= 0.0 )
    {
        std::cerr << "Invalid CSI (" << csi << ")\n";
        throw std::invalid_argument( "Invalid CSI" );
    }
};


// compute:
//      total basal area, 
//      total trees per hectare, 
//      basal area in larger trees (bal_sw and bal_hw) 
// for the tree list
void STAND::compute_ba_tph_bal()
{
    double _bal_sw = 0.0;
    double _bal = 0.0;
    double last_sw_dbh = 9999.0;
    double last_dbh = 9999.0;
    double pending_sw = 0.0;
    double pending = 0.0;

    ba = 0.0;
    ba_sw = 0.0;
    ba_hw = 0.0;
    bf_ba = 0.0;
    ithw_ba = 0.0;
    tph = 0.0;

    // traverse the vector computing and storing bal
    for( auto i : sort_indices(trees, true) )
    {
        ba += trees[i].ba;
        tph += trees[i].tph;

        // balsam fir basal area
        if( trees[i].spp == 12 ) bf_ba += trees[i].ba;

        // intolerant hardwood basal area
        if( !trees[i].is_softwood() && trees[i].get_shade() < 2.0 ) ithw_ba += trees[i].ba;

        if( trees[i].dbh < last_dbh ) 
        {
            trees[i].bal = _bal;
            pending = _bal;
            _bal += trees[i].ba;
            last_dbh = trees[i].dbh;
        } else if( trees[i].dbh == last_dbh ) {
            trees[i].bal = pending;
            _bal += trees[i].ba;
        } else {
            // error, throw exception (need to define exception)
            throw;
        }

        if( trees[i].is_softwood() )
        {
            ba_sw += trees[i].ba;

            if( trees[i].dbh < last_sw_dbh )
            {
                trees[i].bal_sw = _bal_sw;
                pending_sw = _bal_sw;
                _bal_sw += trees[i].ba;
                last_sw_dbh = trees[i].dbh;
            } else if( trees[i].dbh == last_sw_dbh ) {
                trees[i].bal_sw = pending_sw;
                _bal_sw += trees[i].ba;
            } else {
                // error, throw exception (need to define exception)
                throw;
            }
        } else {
            ba_hw += trees[i].ba;
            trees[i].bal_sw = _bal_sw;
        }

        trees[i].bal_hw = trees[i].bal - trees[i].bal_sw;
    }

    qmd = ( tph > 0.0 ) ? std::sqrt( ba / tph / 0.00007854 ) : 0.0;
}


// compute crown competition factor in larger trees (ccfl_sw and ccfl_hw) for the tree list
// requires mcw and lcw to be populated
void STAND::compute_ccfl()
{
    double _ccfl_sw = 0.0;
    double _ccfl = 0.0;
    double last_sw_dbh = 9999.0;
    double last_dbh = 9999.0;
    double pending_sw = 0.0;
    double pending = 0.0;

    // traverse the vector computing and storing ccfl
    for( auto i : sort_indices(trees, true) )
    {
        if( trees[i].dbh < last_dbh ) 
        {
            trees[i].ccfl = _ccfl;
            pending = _ccfl;
            _ccfl += trees[i].mca;
            last_dbh = trees[i].dbh;
        } else if( trees[i].dbh == last_dbh ) {
            trees[i].ccfl = pending;
            _ccfl += trees[i].mca;
        } else {
            // error, throw exception
            throw std::runtime_error( "Exception in compute_ccfl()" );
        }

        if( trees[i].is_softwood() )
        {
            if( trees[i].dbh < last_sw_dbh )
            {
                trees[i].ccfl_sw = _ccfl_sw;
                pending_sw = _ccfl_sw;
                _ccfl_sw += trees[i].mca;
                last_sw_dbh = trees[i].dbh;
            } else if( trees[i].dbh == last_sw_dbh ) {
                trees[i].ccfl_sw = pending_sw;
                _ccfl_sw += trees[i].mca;
            } else {
                throw std::runtime_error( "Exception in compute_ccfl()" );
            }
        } else {
            trees[i].ccfl_sw = _ccfl_sw;
        }

        trees[i].ccfl_hw = trees[i].ccfl - trees[i].ccfl_sw;
    }
}

// compute crown competition factor
void STAND::compute_ccf()
{
    try {
        ccf = 0.0;
        for( auto &t : trees )
            ccf += t.mca;
    } catch( ... ) {
         throw std::runtime_error( "Exception in compute_ccf()" );
    }            
}

// Predict height to crown base if cr is missing, otherwise computes it directly from cr
void STAND::predict_hcb()
{
    try {
        for( auto &t : trees )
            if( t.hcb == 0.0 )
            {
                if( t.cr > 0.0 )
                {
                    t.hcb = (1.0 - t.cr) * t.ht;
                }
                else {
                    t.hcb_pred( ccf ); 
                    t.cr = 1.0 - t.hcb/t.ht;
                }
            }       
    } catch( ... ) {
        throw std::runtime_error( "Exception in predict_hcb()" );
    }

}

// compute:
//      average softwood and hardwood dbh
//      average softwood and hardwood height
//      min and max dbh
//      standard deviation of dbh
void STAND::compute_tree_statistics()
{
    average_dbh = 0.0;
    average_dbh_10 = 0.0;
    average_dbh_sw = 0.0;
    average_dbh_hw = 0.0;
    average_dbh_10_sw = 0.0;
    average_dbh_10_hw = 0.0;
    average_height_sw = 0.0;
    average_height_hw = 0.0;
    average_sg = 0.0;
    average_sg_10 = 0.0;
    min_dbh = 9999.0;
    max_dbh = 0.0;
    sdi = 0.0;
    sdi_10 = 0.0;
    double t_tph_sw = 0.0;
    double t_tph_hw = 0.0;
    double t_tph_10 = 0.0;
    double dbh2 = 0.0;
    double dbh_102 = 0.0;

    try {
        for( auto &t : trees )
        {
            average_dbh += t.dbh * t.tph;
            sdi += std::pow( t.dbh / 25.4, 1.6 ) * t.tph;

            if( t.dbh >= 10.0 )
            {
                average_dbh_10 += t.dbh * t.tph;
                dbh_102 += t.dbh * t.dbh * t.tph;
                t_tph_10 += t.tph;
                sdi_10 += std::pow( t.dbh / 25.4, 1.6 ) * t.tph;
                average_sg_10 += t.get_sg() * t.tph;
            }  

            dbh2 += t.dbh * t.dbh * t.tph;
            average_sg += t.get_sg() * t.tph;

            if( t.is_softwood() )
            {
                average_dbh_sw += t.dbh * t.tph;
                if( t.dbh >= 10.0 )
                    average_dbh_10_sw += t.dbh * t.tph;
                    
                average_height_sw += t.ht * t.tph;
                t_tph_sw += t.tph;
            } else {
                average_dbh_hw += t.dbh * t.tph;
                if( t.dbh >= 10.0 )
                    average_dbh_10_hw += t.dbh * t.tph;

                average_height_hw += t.ht * t.tph;
                t_tph_hw += t.tph;
            }

            if( t.dbh > max_dbh ) max_dbh = t.dbh;
            if( t.dbh < min_dbh ) min_dbh = t.dbh;
        }

        if( tph > 0.0 )
        {
            average_dbh /= tph;
            dbh2 /= tph;
            dbh_sd = std::sqrt((dbh2 - average_dbh*average_dbh)*tph/(tph - 1.0));

            average_sg /= tph;
        }

        if( t_tph_10 > 0.0 )
        {
            average_dbh_10 /= t_tph_10;
            dbh_102 /= t_tph_10;
            dbh_10_sd = std::sqrt((dbh_102 - average_dbh_10*average_dbh_10)*t_tph_10/(t_tph_10 - 1.0));
            average_sg_10 /= t_tph_10;
        }

        if( t_tph_sw > 0.0 )
        {
            average_dbh_sw /= t_tph_sw;
            average_dbh_10_sw /= t_tph_sw;
            average_height_sw /= t_tph_sw;
        }

        if( t_tph_hw > 0.0 )
        {
            average_dbh_hw /= t_tph_hw;
            average_dbh_10_hw /= t_tph_hw;
            average_height_hw /= t_tph_hw;
        }
    } catch( std::exception &e ) {
        std::cerr << "Exception in compute_tree_statistics()\n" << e.what() << "\n";
    }
}

// compute topht for a tree list (average height ot largest 100 trees per hectare )
void STAND::compute_topht()
{
    double sum_tph = 0.0;
    double sum_ht = 0.0;

    for( auto i : sort_indices(trees, false) )
    {
        if( sum_tph + trees[i].tph <= 100.0 )
        {
            sum_ht += trees[i].ht * trees[i].tph;
            sum_tph += trees[i].tph;
        } else if( sum_tph < 100.0 ) {
            sum_ht += trees[i].ht * (100.0 - sum_tph);
            sum_tph = 100.0;
        } else {
            // nothing to do here
        }
    }

    if( sum_tph > 0.0 )
        sum_ht /= sum_tph;
    else
        sum_ht = 0.0;

    topht = sum_ht;
}

// compute the number of unique species in the tree list
void STAND::compute_n_species()
{
    std::set<int> species;

    for( auto &t : trees )
        species.emplace( t.spp );

    n_species = species.size();
}

// SDI and RD
void STAND::compute_sdi_rd()
{
    double dbh_range = min_dbh < max_dbh ? max_dbh - min_dbh : 0.0;
    double dbh_10_range = min_dbh_10 < max_dbh ? max_dbh - min_dbh_10 : 0.0;
    double SDImax = 0.0;
    double SDImax2 = 0.0;
    double SDImax_all = 0.0;
    rd = 0.0;

    try {
        double meanSG = std::max( average_sg, 0.80 ); 
        double meanSG_10 = std::max( average_sg_10, 0.80 ); 

        // for trees with dbh >= 10 cm
        // Weiskittel & Kuehne 2019
        SDImax = 475.2079 - 1.5908*(ba_hw/ba) - 236.9051*std::log(meanSG_10) + 50.3299*std::sqrt(dbh_10_range) + 
                 13.5202*(double)n_species + 0.0685*elevation - 2.8537*std::sqrt(elevation) + 222.7836*(1.0/csi);

        SDImax2 = 1347.445 - 1003.870*meanSG; // Weiskittel and Kuehne (2019)

        SDImax = SDImax > 0.0 ? SDImax2 : SDImax;
        rd_10 = sdi_10 / SDImax;

        // for all trees 
        // Weiskittel & Kuehne 2019
        SDImax_all = 475.2079 - 1.5908*(ba_hw/ba) - 236.9051* std::log(meanSG)+ 50.3299* std::sqrt(dbh_range) +
                     13.5202*(double)n_species + 0.0685*elevation - 2.8537* std::sqrt(elevation) + 222.7836*(1.0/csi); 
        SDImax_all = SDImax_all > 0.0 ? SDImax_all : SDImax2;

        rd = sdi / SDImax_all;

    } catch( std::exception &e ) {
        std::cerr << "Exception in compute_sdi_rd() for average_sg " << average_sg << ", average_sg_10 " << average_sg_10 <<
                     ", ba" << ba << ", ba_hw" << ba_hw << ", dbh_10_range " << dbh_10_range << ", n_species " << n_species <<
                     ", elevation " << elevation << ", csi " << csi << "\n" << e.what() << "\n";
        throw;
    }

}

// process tree list to expand tree records with greater than threshold tph
bool STAND::expand_tree_list( double threshold )
{
    // check to see if there is a valid tree list
    if( trees.size() == 0 )
        return false;

    // create random number generator
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-0.005,0.005);

    try {
        // process tree list
        auto n_trees = trees.size();
        for( size_t i = 0; i < n_trees; i++ )
        {
            // if tph of the currrent tree record is greater than threshold, expand it
            if( trees[i].tph > threshold )
            {
                // compute number of new tree records required
                auto n_new_trees = std::trunc(trees[i].tph / threshold) - 1;

                double cum_tph = threshold;
                size_t j;
                for( j = 0; j < n_new_trees; j++ )
                {
                    TREE new_tree = trees[i];
                    new_tree.expand_tree_id = j+1;
                    new_tree.dbh += distribution(generator);
                    new_tree.ht += new_tree.ht > 0.0 ? distribution(generator) : 0.0;
                    new_tree.tph = threshold;
                    new_tree.compute_attributes();
                    cum_tph += threshold;

                    trees.emplace_back( new_tree );
                }

                if( cum_tph < trees[i].tph )
                {
                    TREE new_tree = trees[i];
                    new_tree.expand_tree_id = ++j;
                    new_tree.dbh += distribution(generator);
                    new_tree.ht += new_tree.ht > 0.0 ? distribution(generator) : 0.0;
                    new_tree.tph = trees[i].tph - cum_tph;
                    new_tree.compute_attributes();

                    trees.emplace_back( new_tree );
                }

                trees[i].tph = threshold;
                trees[i].expand_tree_id = ++j;
                trees[i].compute_attributes();

            }
        }
    } catch( std::exception &e ) {
        std::cerr << "Exception while expanding tree list.\n" << e.what() << "\n";
        throw;
    }

    return true;
}

// collapse tree records that were spread using STAND::expand_tree_list() due to exceeding tph threshold
bool STAND::unexpand_tree_list()
{
    // search list for spread trees (expand_tree_id > 0)
    for( auto t_itr = trees.begin(); t_itr != trees.end(); t_itr++ )
    {
        if( t_itr->expand_tree_id > 0 )
        {
            // use first expanded tree as accumulator
            // weight by tph
            t_itr->dbh *= t_itr->tph; 
            t_itr->ht  *= t_itr->tph;
            t_itr->hcb *= t_itr->tph;
            t_itr->cr  *= t_itr->tph;

            // aggregate spread tree data (assuming all tree records above this point are compressed)
            auto t2_itr = t_itr;
            for( ++t2_itr; t2_itr != trees.end(); t2_itr++ )
            {
                if( t2_itr->plot_id == t_itr->plot_id && 
                    t2_itr->tree_id == t_itr->tree_id && 
                    t2_itr->expand_tree_id > 0 &&
                    t2_itr->expand_tree_id != t_itr->expand_tree_id )
                {
                    t_itr->tph += t2_itr->tph;
                    t_itr->dbh += t2_itr->dbh * t2_itr->tph;
                    t_itr->ht  += t2_itr->ht * t2_itr->tph;
                    t_itr->hcb += t2_itr->hcb * t2_itr->tph;
                    t_itr->cr  += t2_itr->cr * t2_itr->tph;
                    t2_itr->expand_tree_id = -1;
                    t2_itr->tph = 0.0;
                }
            }

            // compute averages
            if( t_itr->tph > 0.0 )
            {
                t_itr->dbh /= t_itr->tph;
                t_itr->ht  /= t_itr->tph;
                t_itr->hcb /= t_itr->tph;
                t_itr->cr  /= t_itr->tph;
                t_itr->expand_tree_id = 0;

                // recompute tree attributes
                t_itr->compute_attributes();

            }
        }
    }

    // clean up tree vector
    std::erase_if( trees, [](TREE t){return t.tph == 0.0;} );

    return true;
}

// ingrowth species group cross walk
const std::unordered_map<int,int> cross_walk = {
    {  531 /*AB*/,  531 },
    {  746 /*QA*/,  746 },
    {  318 /*SM*/,  318 },
    {  241 /*WC*/,  241 },
    {  379 /*GB*/,  379 },
    {  375 /*PB*/,  379 },
    {  371 /*YB*/,  379 },
    {   12 /*BF*/,   12 },
    {  316 /*RM*/,  316 },
    {   97 /*RS*/,   97 },
    {   95 /*BS*/,   97 },
    {   94 /*WS*/,   97 },
    {  129 /*WP*/,  129 },
    { 9990 /*OH*/, 9990 },
    { 9991 /*OS*/, 9991 },
};

// build basal area by species map for use in ingrowth distribution
// birches and spruce are grouped together
void STAND::build_ba_spp_map()
{
    // build species basal area map
    for( auto &t : trees )
    {
        if( cross_walk.find( t.spp ) != cross_walk.end() )
        {
            ba_spp[t.spp] += t.ba;
            ba_grp_spp[cross_walk.at(t.spp)] += t.ba;
            plot_species_ba[t.plot_id][t.spp] += t.ba;
        } else {
            // any species not explicitly enumerated, defaults to other softwoods/hardwoods
            ba_spp[t.is_softwood() ? 9991 : 9990] += t.ba;
            ba_grp_spp[t.is_softwood() ? 9991 : 9990] += t.ba;
            plot_species_ba[t.plot_id][t.is_softwood() ? 9991 : 9990] += t.ba;
        }
    }    
}

// allocate ingrowth to species and create new trees in the tree list 
void STAND::ingrowth_composition( double IPH )
{
    std::array<double,35> b = {
        /* b10 0*/  -2.5645, /* b11 1*/   0.0020,  /* b12 2*/  2.6624,  /* b13 3*/ -0.0010,  /* b14 4*/ -0.0127, 
        /* b20 5*/  -3.0291, /* b21 6*/   0.0027,  /* b22 7*/  2.7779,  /* b23 8*/  0.0211,  /* b24 9*/  0.0221,
        /* b30 10*/ -0.6566, /* b31 11*/  0.0123, /* b32 12*/  1.7669, /* b33 13*/ -0.0421, /* b34 14*/ -0.0283,
        /* b40 15*/ -1.2500, /* b41 16*/ -0.0132, /* b42 17*/  2.0470, /* b43 18*/ -0.0514, /* b44 19*/  0.0351,
        /* b50 20*/ -5.1074, /* b51 21*/ -0.0117, /* b52 22*/  3.8817, /* b53 23*/  0.0501, /* b54 24*/  0.0726,
        /* b60 25*/ -2.9832, /* b61 26*/ -0.0020, /* b62 27*/  2.4837, /* b63 28*/  0.0673, /* b64 29*/ -0.0167, 
        /* b70 30*/ -4.7182, /* b71 31*/  0.0070, /* b72 32*/  3.2269, /* b73 33*/  0.1000, /* b74 34*/  0.0188 }; 

    // ingrowth allocation by species group
    std::unordered_map<int,double> species_group_percent;

    double total_percent = 0.0;

    // compute ingrowth allocation by species/species group
    for( auto &[spp, sba] : ba_grp_spp )
    {
        double percent = 0.0;
        double pba = sba/ba;

        if( spp == 379 || spp == 375 || spp == 371 )   // BCH: GB 379 PB 375 YB 371 
        {             
            percent = b[0] + b[1]*ba + b[2]*pba + b[3]*csi + b[4]*MinDBH;

        } else if( spp == 12 ) // BF
        {
            percent = b[5] + b[6]*ba + b[7]*pba + b[8]*csi + b[9]*MinDBH;
        }

        else if( spp == 316 ) // RM
        {
            percent = b[10] + b[11]*ba + b[12]*pba + b[13]*csi + b[14]*MinDBH;

        }
        else if( spp == 97 || spp == 95 || spp == 94 ) // SPR: RS 97 BS 95 WS 94
        {
            percent = b[15] + b[16]*ba + b[17]*pba + b[18]*csi + b[19]*MinDBH;

        }
        else if( spp == 129 ) // WP
        {
            percent = b[20] + b[21]*ba + b[22]*pba + b[23]*csi + b[24]*MinDBH;

        }

        else if( spp == 9990 || spp == 746 || spp == 531 || spp == 318 ) // OH
        {
            percent = b[25] + b[26]*ba + b[27]*pba + b[28]*csi + b[29]*MinDBH;

        }

        else if( spp == 9991 ) // OS
        {
            percent = b[30] + b[31]*ba + b[32]*pba + b[33]*csi + b[34]*MinDBH;

        }

        percent = 1.0 / (1.0 + std::exp(-(percent)));

        species_group_percent[spp] += percent;

        total_percent += percent;
    }

    // normalize percentages (divide each group by total percentages)
    for( auto &[spp, sgp] : species_group_percent )
    {
        // normalize
        sgp /= total_percent;

        // allocated ingrowth
        sgp *= IPH;
    }

    // compute allocated ingrowth and create a new tree record for each plot-species combination
    for( auto &[spp, sba] : ba_spp )
    {
        // ingrowth for current species
        auto spp_ingrowth = species_group_percent[cross_walk.at(spp)] * sba / ba_grp_spp[cross_walk.at(spp)];

        // allocate ingrowth to plots
        for( auto &[plot, psba] : plot_species_ba )
        {
            // compute percentage of current species ba in this plot
            // this will be used to allocate species ingrowth
            auto plot_spp_percent = psba[spp] / sba;

            // create new tree
            if( plot_spp_percent > 0.0 )
            {
                TREE new_tree( plot, ++max_tree_id, spp, MinDBH, 0.0, spp_ingrowth*plot_spp_percent, 0.0, 0, 0 );
                trees.emplace_back( new_tree );
            }
        }
    }
}
 

// Ingrowth -- considers only annual growth cycles; returns ingrowth per hectare (IPH)
// Function of Li et al. (2011; CJFR 41, 2077-2089)
double STAND::ingrowth( INGROWTH_MODEL_TYPE type )

{    
    const std::array<double,7> gnls_a = { -0.2116, -0.0255,  -0.1396, -0.0054,   0.0433,   0.0409,   0.0 };
    const std::array<double,7> gnls_b = {  3.8982, -0.0257,  -0.3668,  0.0002,   0.0216,  -0.0514,   0.0 };
    const std::array<double,7> nlme_a = { -0.08217, 0.1113,  -1.2405, -0.2319,   0.03673, -0.7745,  -0.1301 };
    const std::array<double,7> nlme_b = {  2.8466, -0.03114, -0.2891,  0.003350, 0.2248,  -0.08223, -0.03548 };

    auto a = ( type == INGROWTH_MODEL_TYPE::GNLS ) ? gnls_a : nlme_a;
    auto b = ( type ==  INGROWTH_MODEL_TYPE::GNLS ) ? gnls_b : nlme_b;

    double link1 = a[0] + a[1]*ba + a[2]*(ba_hw/ba) + a[3]*(tph/1000.0) + a[4]*csi + a[5]*MinDBH + a[6]*qmd;

    double PI = (1.0/(1.0 + std::exp(-link1)));
    double eta = b[0] + b[1]*ba + b[2]*(ba_hw/ba) + b[3]*(tph/1000) + b[4]*csi + b[5]*(MinDBH) + b[6]*qmd;
    double IPH = std::exp(eta);

    if( cut_point == 0.0 )
        IPH *= PI;
    else
        IPH = PI >= cut_point ? IPH : 0.0;

    return(IPH);
}


// estimate annual diameter growth for the tree list 
void STAND::diameter_growth()
{
    try {
        for( auto &t : trees )
            t.dDBH(Region, csi, ba, percent_ba_removed, ba_pre_thin, qmd_ratio, thin_year, year, average_dbh_10_sw, topht, CDEF );

    } catch( ... ) {
        throw;
    }
}


// grow heights for the tree list one year
void STAND::height_growth()
{
    try {
        for( auto &t : trees )
            t.dHT( Region, csi, percent_ba_removed, ba_pre_thin, qmd_ratio, thin_year, year, average_dbh_10_sw, topht, CDEF );

    } catch( ... ) {
        throw;
    }
}

// estimate crown recession for the tree list one year
void STAND::crown_recession()
{
    try {
        for( auto &t : trees )
            t.dHCB( ccf, percent_ba_removed, ba_pre_thin, qmd_ratio, thin_year, year );

    } catch( ... ) {
        throw;
    }
}

// #### mortality modifiers ###
// SBW mortality modifier: stand level
double STAND::mort_sbw()
{
    double b1 = 0.0;
    double b2 = 0.0;
    double b3 = 0.0;
    double b4 = 0.0;
    double rat = 1.0;

    if( Region == "ME" )
    {
        b1 = -2.6380;
        b2 =  0.0114;
        b3 = -0.0076;
        b4 =  0.0074;
    }
    else if(Region == "NB")
    {
        b1 = -3.0893;
        b2 =  0.0071;
        b3 = -0.0037;
        b4 =  0.0000;
    }
    else
    {
        b1 = -2.6380;
        b2 =  0.0114;
        b3 = -0.0076;
        b4 =  0.0074;
    }

    try {
        double VOL = (topht/2.0)*ba;
        double aa = (1.0/(1.0 + std::exp(-b1)))*(1.0/(1.0 + std::exp(-(b3*VOL))));
        double bb = (1.0/(1.0 + std::exp(-b1)))*(1.0/(1.0 + std::exp(-(b2*CDEF*bf_ba + b3*VOL + b4*CDEF))));
        rat = (CDEF < 0.0) ? 1.0 : (aa > 0.0 ? bb/aa : 1.0);
    } catch( std::exception &e ) {
        std::cerr << "Exception in mort_sbw() for ba " << ba << ", bf_ba " << bf_ba << ", topht " << topht <<
                     ", CDEF " << CDEF << "\n" << e.what() << "\n";
        throw;
    }

    // result is mortality probability multiplier value >=1
    return rat;
}


// Thinning survival probability modifier: stand level
double STAND::mort_thin()
{
    double modifier = 1.0;

    if( thin_year >= 0 && thin_year <= year && percent_ba_removed > 0.0 && qmd_ratio > 0.0 && ba_pre_thin > 0.0 )
    {
        int time_since_thinning = year - thin_year;
        
        // parameter estimates 
        // double b30 =   -1.2402;
        // double b31 =  -24.5202;
        // double b32 =   -1.1302;
        // double b33 =    1.5884;
        double y30 =    8.3385;
        double y31 = -601.3096;
        double y32 =    0.5507;
        double y33 =    1.5798;
        
        try {
            modifier =  1.0 + std::exp(y30 + (y31/((100.0*(percent_ba_removed) + ba_pre_thin)+0.01))) * 
                        std::pow(y32, time_since_thinning) * std::pow(time_since_thinning,y33);
        } catch( std::exception &e ) {
            std::cerr << "Exception in mort_thin() for percent_ba_removed " << percent_ba_removed << 
                         ", ba_pre_thin " << ba_pre_thin << ", time_since_thinning " << time_since_thinning << "\n" << 
                         e.what() << "\n";
            throw;
        }
    }  

    return modifier;
}

// compute tree survival probability
void STAND::compute_survival_prob()
{
    try {
        for( auto &t : trees )
            t.survival_prob( Region, csi, ba, qmd, percent_ba_removed, ba_pre_thin, qmd_ratio, thin_year, year,
                             average_height_hw, average_height_sw, CDEF, use_sbw_mod, use_hw_mod, use_thin_mod );

    } catch( ... ) {
        throw;
    }
}

// apply growth and mortality to current trees for the stand
void STAND::apply_growth_mortality()
{
    for( auto &t : trees )
    {
        t.apply_growth_mortality();
    }
}

// compute mortality (change in trees per hectare)
void STAND::survival()
{
    try {
        // mortality modifiers
        double sbw_modifier = use_sbw_mod ? mort_sbw() : 1.0;
        double thin_modifier = use_thin_mod ? mort_thin() : 1.0;

        // estimate tree level mortality
        compute_survival_prob();

        size_t i = 0;
        try {
            // compute tree level mortality
            for( i = 0; i < trees.size(); i++ )
                trees[i].set_dtph( trees[i].tph * (1.0 - trees[i].get_survival())*sbw_modifier*thin_modifier );

        } catch( std::exception &e ) {
            std::cerr << "Exception in survival() for tree " << trees[i].tree_id << "\n" << e.what() << "\n";
            throw;
        }
            
    } catch( ... ) {
        throw;
    }
}

// search tree list for maximum tree id value
unsigned long STAND::find_max_tree_id()
{
    unsigned long t_max_tree_id = 0;
    for( auto &t : trees )
        if( t.tree_id > t_max_tree_id )
            t_max_tree_id = t.tree_id;

    return t_max_tree_id;
}

// Initialize a STAND object
void STAND::initialize()
{
    try {
        // expand the tree list if needed
        if( !expand_tree_list( 50.0  ) )
            throw std::runtime_error( "Error expanding tree list.\n" );

        max_tree_id = find_max_tree_id();

        compute_n_species();
        compute_ccf();    
        compute_ba_tph_bal(); 
        compute_ccfl();

        // impute missing heights
        auto r = Region == "ME" ? 0 : 1;
        for( auto &t : trees )
            t.ht_pred( ccf, r );

        compute_topht();
        predict_hcb();

        compute_tree_statistics();
        compute_sdi_rd();

        initialized = true;

    } catch( std::exception &e ) {
        std::cerr << "Exception in initialize()\n" << e.what() << "\n";
        throw;
    }
}

// grow a STAND object for n_years
// The stand's tree list is updated after execution
void STAND::grow( int n_years )
{
    try {
        // initialize stand and tree variables if needed
        if( !initialized )
            initialize();

        for( int i = 0; i < n_years; i++ )
        {
            if( use_ingrowth )
            {
                // compute total ingrowth 
                auto iph = ingrowth( INGROWTH_MODEL_TYPE::GNLS );

                // if ingrowth present, allocate it and create trees
                if( iph > 0.0 )
                {
                    // build species basal area map for ingrowth calculations
                    build_ba_spp_map();
                    ingrowth_composition( iph );
                    initialize();
                }
            }

            // estimate growth and mortality
            diameter_growth();
            height_growth();
            crown_recession();
            survival();

            // increment tree and crown dimensions
            apply_growth_mortality();

            // recompute stand and tree statistics
            compute_ccf();    
            compute_ba_tph_bal(); 
            compute_ccfl();
            compute_topht();
            compute_tree_statistics();
            compute_sdi_rd();

            year++;
        }

        // unexpand the tree list if needed
        if( !unexpand_tree_list() )
            throw std::runtime_error( "Error unexpanding tree list.\n" );
    } catch( ... ) {
        throw;
    }
}
