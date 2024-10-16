//
// Command line interface to ACD 7

// (c) 2024 Greg Johnson
// Greg Johnson Biometrics LLC
// 05/07/2024
//

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <exception>
#include <vector>

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
	if( argc < 3 )
	{
		std::cerr << "Command line requires 2 parameters, " << argc-1 << " supplied\n";
		std::cerr << "Usage:\nacd.exe XX stand.csv\nwhere XX is number of years to project each tree list.";
		return -1;
	}

	// get number of years to project tree list
	int n_periods = 0;
	try{ 
		n_periods = std::atoi( argv[1] );

	} catch( std::exception &e ) {
		std::cerr << "Error in specifying number of years to project tree list.\n";
		return -2;
	}

	// attempt to open stand file
	std::string stand_filename{argv[2]};
	std::ifstream stand_file{stand_filename};
	if( !stand_file.is_open() )
	{
		std::cerr << "Did not find or could not open " << stand_filename << "\n";
		return -3;
	}

 	//////////////////////////////////////////////////////////////////////////////////////////////////////
  
	// read and parse stand information file
	// format:
	// stand_id, units, csi, elev, cdef, use_sba, use_hw, use_thin, use_ingrowth, cut_point, MinDBH
	struct STAND_INFO {
		std::string region = "ME";
		std::string stand_id;
		int units = 1; // 0 = metric, 1 = imperial
		int year = 0;
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

			s.region = get_string( ss );
			s.stand_id = get_string( ss );
			s.units = get_int( ss );
			s.year = get_int( ss );
			s.csi = get_double( ss );
			s.elev = get_double( ss );
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
		return -51;
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

		// create STAND object
		STAND stand( s.region, s.year, s.csi*ft_m, s.elev*ft_m, s.cdef, 
					 s.use_sbw, s.use_hw, s.use_thin, s.use_ingrowth, s.cut_point, s.MinDBH );

		// attempt to open tree list file
		std::string trees_filename = s.stand_id + ".csv";
		std::ifstream trees_file{trees_filename};
		if( !trees_file.is_open() )
		{
			std::cerr << "Did not find or could not open " << trees_filename << "\n";
			return -4;
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////
		// read and parse tree list file
		std::string stand_id;
		unsigned long plot_id;
		unsigned long tree_id;
		int spp;
		double dbh;
		double ht;
		double expf;
		double cr;
		int form_t;
		int risk_t;

		std::string header;
		int line_no = 1;
		try {
			// skip header and read first tree record
			std::getline(trees_file, header);
			std::getline( trees_file, line );

			while( !trees_file.eof() )
			{          
				std::istringstream ss(std::move(line));

				stand_id = get_string( ss );

				// check to make sure stand_id matches
				if( s.stand_id != stand_id )
				{
					std::cerr << "Stand ID in " << stand_filename << " does not match Stand ID in " << trees_filename << "\n";
					return -53;
				}

				plot_id = get_ulong( ss);
				tree_id = get_ulong( ss );
				spp = get_int( ss );
				dbh = get_double( ss );
				ht = get_double( ss );
				expf = get_double( ss );
				cr = get_double( ss );
				form_t = get_int( ss );
				risk_t = get_int( ss );

				TREE t( plot_id, tree_id, spp,  dbh*in_cm,  ht*ft_m,  expf*ac_ha, cr, form_t, risk_t );
				stand.trees.emplace_back( t );
				line_no++;

				// get next tree record
				std::getline( trees_file, line );
			}
		} catch( std::exception &e ) {
			std::cerr << "Could not read tree list data from " << trees_filename << " on line " << line_no << "\n" << e.what() << "\n";
			return -52;        
		}

		trees_file.close();

		//////////////////////////////////////////////////////////////////////////////////////////////////////
		// grow the tree list
		try {
			stand.grow( n_periods );
		} catch ( std::exception &e ) {
			std::cerr << "Error growing tree list from " << trees_filename << "\n" << e.what() << "\n";
			return -61;        
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////

		// write out grown tree list
		if( !header_written) 
		{
			std::cout << header << "\n";
			header_written = true;
		}
		for( auto &t : stand.trees )
		{
			std::cout << s.stand_id << ", " << t.plot_id << ", " << t.tree_id << ", " << t.spp << ", " << 
						t.dbh/in_cm << ", " << t.ht/ft_m << ", " << t.tph/ac_ha << ", " << t.cr << ", " <<
						t.Form << ", " << t.Risk << "\n";
		}
	}

	return 1;
}
