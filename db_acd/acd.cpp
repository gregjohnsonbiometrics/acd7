//
// Command line interface to ACD 7 with FIA DB Access
//
// (c) 2024 Greg Johnson
// Greg Johnson Biometrics LLC
// 05/21/2024
//

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <exception>
#include <vector>

#define OTL_ODBC      // Compile OTL 4.0/ODBC
#define OTL_CPP_20_ON
#define OTL_ODBC_LEGACY_RPC
#ifndef __linux__
	#define OTL_STRCAT_S(dest, dest_sz, src) strcat_s(dest, dest_sz, src)
	#define OTL_STRCPY_S(dest, dest_sz, src) strcpy_s(dest, dest_sz, src)
	#define OTL_STRNCPY_S(dest, dest_sz, src, count) strncpy_s(dest, dest_sz, src, count)
	#define OTL_SPRINTF_S sprintf_s
#else
	#define OTL_ODBC_UNIX // uncomment this line if UnixODBC is used
#endif

#include "otlv4.h"    // include the OTL 4.0 header file

#include "stand.hpp"

unsigned long get_ulong( std::istringstream &ss ) 
{
	std::string value;
	std::getline( ss, value, ',' );
	return std::stoul( value );
}

double get_double( std::istringstream &ss ) 
{
	std::string value;
	std::getline( ss, value, ',' );
	return std::stod( value );
}

int get_int( std::istringstream &ss ) 
{
	std::string value;
	std::getline( ss, value, ',' );
	return std::stoi( value );
}

bool get_bool( std::istringstream &ss ) 
{
	std::string value;
	std::getline( ss, value, ',' );
	return std::stoi( value ) > 0;
}

std::string get_string( std::istringstream &ss ) 
{
	std::string value;
	std::getline( ss, value, ',' );
	return value;
}

int main( int argc, char **argv )
{
	//////////////////////////////////////////////////////////////////////////////////////////////////////

	// check to see if proper number arguments present
	if( argc < 4 )
	{
		std::cerr << "Command line requires 3 parameters, " << argc-1 << " supplied\n";
		std::cerr << "Usage:\nacd.exe XX fia_db stand.csv\nwhere XX is number of years to project each tree list.";
		return -1;
	}

	// get number of years to project tree list
	int n_periods = 0;
	try { 
		n_periods = std::atoi( argv[1] );

	} catch( std::exception &e ) {
		std::cerr << "Error in specifying number of years to project tree list.\n";
		return -2;
	}

	// attempt to open FIA DB
	otl_connect db; // connect object

	try {
		otl_connect::otl_initialize(); // initialize ODBC environment
	} catch( std::exception &e ) {
		std::cerr << "could not initialize otl.\n";
		return -3;
	}

 	try {
		#ifdef __linux__ 
			std::string connect_string = std::string("DRIVER={SQLite3};Database=") + std::string(argv[2]) + 
		 							 std::string(";LongNames=0;Timeout=1000;NoTXN=0;SyncPragma=NORMAL;StepAPI=0;");
		#else
			std::string connect_string = std::string("DRIVER={SQLite3 ODBC Driver};Database=") + std::string(argv[2]) + 
									 std::string(";LongNames=0;Timeout=1000;NoTXN=0;SyncPragma=NORMAL;StepAPI=0;");	
		#endif								 
		db.rlogon( connect_string.c_str() ); 
	} catch(otl_exception& p) { // intercept OTL exceptions
		std::cerr << p.msg <<"\n"; // print out error message
		std::cerr << p.stm_text <<"\n"; // print out SQL that caused the error
		std::cerr << p.sqlstate <<"\n"; // print out SQLSTATE message
		std::cerr << p.var_info <<"\n"; // print out the variable that caused the error
		return -4;
 	}

	// attempt to open stand file
	std::string stand_filename{argv[3]};
	std::ifstream stand_file{stand_filename};
	if( !stand_file.is_open() )
	{
		std::cerr << "Did not find or could not open " << stand_filename << "\n";
		return -5;
	}

 	//////////////////////////////////////////////////////////////////////////////////////////////////////
  
	// read and parse stand information file
	// format:
	// stand_id, csi, cdef, use_sba, use_hw, use_thin, use_ingrowth, cut_point, MinDBH
	struct STAND_INFO {
		std::string region = "ME";
		std::string stand_id; // this is a FIA DB STAND_CN value
		int units = 1; // 0 = metric, 1 = imperial
		int year = 10; // arbitrary age
		double csi = 16.0;
		double elev = 0.0;
		double cdef = 0.0;
		bool use_sbw = false;
		bool use_hw = false;
		bool use_thin = false;
		bool use_ingrowth = false;
		double cut_point = 0.5; 
		double MinDBH = 3.0 / 2.54; // 3 cm
	};
 
	std::vector<STAND_INFO> stand_info;
	STAND_INFO s;

	// skip header record
	std::string line;
	std::getline(stand_file, line);
	std::getline(stand_file, line);

	try {
		while( !stand_file.eof() )
		{
			// parse stand info record
			std::istringstream ss(std::move(line));

			s.stand_id = get_string( ss );
			s.csi = get_double( ss );
			s.cdef = get_double( ss );
			s.use_sbw = get_bool( ss );
			s.use_hw = get_bool( ss );
			s.use_thin = get_bool( ss );
			s.use_ingrowth = get_bool( ss );
			s.cut_point = get_double( ss );
			s.MinDBH = get_double( ss );

			stand_info.emplace_back( s );

			std::getline(stand_file, line);
		}

	} catch( std::exception &e ) {
		std::cerr << "Could not read stand info record from " << stand_filename << "\n" << e.what() << "\n";
		return -41;
	}

	stand_file.close();

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// process each stand in the stand info file
	bool header_written = false;

	for( auto &s : stand_info )
	{
		const double ft_m =  (s.units == 0) ? 1.0 : 0.3048;
		const double in_cm = (s.units == 0) ? 1.0 : 2.54;
		const double ac_ha = (s.units == 0) ? 1.0 : 2.47105;	

		// read stand information from FIA DB
		std::string query = std::string("select BASAL_AREA_FACTOR, BRK_DBH, AGE, ELEVFT, SITE_INDEX FROM ") +
                            std::string("FVS_STANDINIT_PLOT where STAND_CN = \"") + s.stand_id + std::string("\"");

		double baf;
		double brk_dbh;
		double site;
	
		otl_stream i( 50, query.c_str(), db );
		if( i.eof() )
		{
			std::cerr << "Could not find STAND_CN = " << s.stand_id << "\n";
			continue;
		}

		i >> baf >> brk_dbh >> s.year >> s.elev >> site;
		site = (site > 1.0 && s.csi == 0.0) ? site*ft_m : s.csi;	
		baf = (baf == 0) ? 1.0 : 24.07219;

		// create STAND object
		STAND stand( s.region, s.year, site , s.elev*ft_m, s.cdef, 
					 s.use_sbw, s.use_hw, s.use_thin, s.use_ingrowth, s.cut_point, s.MinDBH );


		//////////////////////////////////////////////////////////////////////////////////////////////////////
		// read and parse tree list data
		char stand_id[30];
		double plot_id;
		double tree_id;
		char spp[5];
		double dbh;
		double ht;
		double expf;
		double cr;
		int form = 0;
		int risk = 0;

		query = std::string( "select STAND_CN, PLOT_ID, TREE_ID, TREE_COUNT, SPECIES, DIAMETER, HT, CRRATIO FROM " ) +
                std::string( "FVS_TREEINIT_PLOT where STAND_CN = \'" ) + s.stand_id + std::string("\'");

		try {
			otl_stream it( 50, query.c_str(), db );
			while( !it.eof() )
			{
				it >> stand_id >> plot_id >> tree_id >> expf >> spp >> dbh >> ht;
				if( it.is_null() )
					ht = 0.0;
				it >> cr;
				if( it.is_null() )
					cr = 0.0;

				expf *= (dbh < brk_dbh && brk_dbh != 999 ? 299.8611 : baf);

				auto ispp = std::atoi(spp);
				try {
					TREE t( static_cast<unsigned long>(plot_id), static_cast<unsigned long>(tree_id), ispp,  dbh*in_cm,  ht*ft_m,  expf*ac_ha, cr/100.0, form, risk );
					stand.trees.emplace_back( t );
				} catch( std::exception &e ) {
					std::cerr << e.what() << "\n";
					return -51;
				}
			}
		} catch(otl_exception& p) { // intercept OTL exceptions
			std::cerr << p.msg <<"\n"; // print out error message
			std::cerr << p.stm_text <<"\n"; // print out SQL that caused the error
			std::cerr << p.sqlstate <<"\n"; // print out SQLSTATE message
			std::cerr << p.var_info <<"\n"; // print out the variable that caused the error
			return 0;
		} catch( std::exception &e ) {
			std::cerr << "Could not read tree list data for STAND_CN " << s.stand_id << "\n" << e.what() << "\n";
			return -52;        
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////
		// grow the tree list
		try {
			stand.grow( n_periods );
		} catch ( std::exception &e ) {
			std::cerr << "Error growing tree list for " << s.stand_id << "\n" << e.what() << "\n";
			return -61;        
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////

		// write out grown tree list
		if( !header_written) 
		{
			std::cout << "stand_id, plot_id, tree_id, species, dbh, ht, tpa, cr, form, risk\n";
			header_written = true;
		}
		for( auto &t : stand.trees )
		{
			std::cout << s.stand_id << ", " << t.plot_id << ", " << t.tree_id << ", " << t.spp << ", " << 
						t.dbh/in_cm << ", " << t.ht/ft_m << ", " << t.tph/ac_ha << ", " << t.cr << ", " <<
						t.Form << ", " << t.Risk << "\n";
		}
	}

	return 0;
}
