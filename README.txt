ACD v7
07/08/2024

==============================================

R Package
To install the R package install the acdR_0.76.zip using the local file install option in R or RStudio (a source package is in acdR_0.76.tar.gz). An example script using the AcdFiaBench.db is included in the file: acdR_example.R

==============================================

Command Line using CSV

A command line version of ACD v7 that processes csv files is included (acd.exe). To use it:

> acd.exe periods stand_info.csv > output.csv

where: 	* periods = number of years to project the tree list
       	* stand_info.csv contains the following fields (comma separated -- see the example test_stand_info.csv):
		region,stand_id,plot_id,units,year,csi,elev,cdef,use_sbw,use_hw,use_thin,use_ingrowth,cut_point,MinDBH

The program looks for a csv with the name of each stand_id in stand_info.csv. 

the tree_list.csv contains the following fields (comma separated -- see the example 0023200806050102903343.csv):
		Stand_ID, Plot_ID, Tree_ID, SPECIES, DBH, HT, tpa, cr, form, risk

Note that the header in each csv is required, but the text in the header is arbitrary. The field contents must be in the specified order.

acd.exe writes the grown tree list to stdout (the example above redirects it to "output.csv".

==============================================

Command Line Using FIA DB

A command line version of ACD v7 that processes stand conditions in FIA DB is included (db_acd.exe). To use it:

> db_acd.exe periods "path to FIA SQLite DB" stand_info.csv > output.csv

where: 	* periods = number of years to project the tree list
	* "path to FIA SQLite DB" for example C:/SQLite_FIADB_ME.db for the Maine FIA DB
       	* stand_info.csv contains the following fields (comma separated -- see the example db_stand_info.csv):
		STAND_CN,csi,cdef,use_sbw,use_hw,use_thin,use_ingrowth,cut_point,MinDBH

Note that the header in the db_stand_info csv is required, but the text in the header is arbitrary. The field contents must be in the specified order.

An Ubuntu 22.04 version is included (db_acd.linux). It requires SQLite3 to be installed and configured.

Future enhancements include an option to write the projected tree list to a table in an SQLite DB.