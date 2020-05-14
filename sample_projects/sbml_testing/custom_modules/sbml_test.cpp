/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./sbml_test.h"

int oxygen_substrate_idx; 
int glucose_substrate_idx; 
int energy_cell_idx; 
int ingest_oxy_cell_idx;
int ingest_glu_cell_idx;

#include "../intracellular/PhysiCell_intracellular.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

void setup_microenvironment( void )
{
	// domain parameters read from XML config file

	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
	// initialize BioFVM 
	initialize_microenvironment(); 	

	oxygen_substrate_idx = microenvironment.find_density_index( "oxygen" ); 
	glucose_substrate_idx = microenvironment.find_density_index( "glucose" ); 
	std::cout << "---------- setup_microenvironment() -----------\n";
	std::cout << "    oxygen_substrate_idx = " << oxygen_substrate_idx << std::endl;
	std::cout << "    glucose_substrate_idx = " << glucose_substrate_idx << std::endl;

	double oxy = 38.0;  // IC
	double oxy_del = 9.0;
	double glu = 32.0;  // IC
	double glu_del = 7.5;
	double x;

	int nregions = 5;
	// double xmin = -750.;   // NOTE! this depends on the config .xml !
	// double xmin = -150.;   // NOTE! this depends on the config .xml !
	double xmin = default_microenvironment_options.X_range[0]; 
	std::cout << "setup_microenvironment: xmin= " << xmin << std::endl;
	// double xdel = 1500./nregions;
	double xdel = (-2.0 * xmin)/nregions;
	std::cout << "setup_microenvironment: xdel= " << xdel << std::endl;

	// 5 regions across x: [-750:-450:-150:150:450:750]
	std::cout << "setup_microenvironment: num voxels= " << microenvironment.number_of_voxels() << std::endl;
	for( int n=0; n < microenvironment.number_of_voxels(); n++ )
	{
		// x coordinate of the nth voxel's center
		x = microenvironment.mesh.voxels[n].center[0];
		for( int iregion=1; iregion <= nregions; iregion++ )
		{
			float xmax = xmin + iregion*xdel;
			// std::cout << "setup_microenvironment: xmax= " << xmax<< std::endl;
			if (x < xmax)
			// if (x < (xmin + iregion*xdel))
			{
				microenvironment(n)[oxygen_substrate_idx] = oxy - (iregion-1)*oxy_del;
				microenvironment(n)[glucose_substrate_idx] = glu - (iregion-1)*glu_del;
				break;
			}
			// oxy -= iregion-5;
			// glu -= iregion-2;
		}
	}
	
	std::cout << "---------- finished setup_microenvironment() -----------\n";
	return; 
}

void create_cell_types( void )
{
	// set the random seed 
	// SeedRandom( parameters.ints("random_seed") );    // do in main.cpp
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	// cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.volume_update_function = NULL;
	// cell_defaults.functions.update_velocity = NULL;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	// cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.update_phenotype = 	energy_based_cell_phenotype; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	// do this BEFORE initialize_cell_def
	cell_defaults.phenotype.motility.migration_bias_direction = { 1.0, 0.0, 0.0 };  

	initialize_cell_definitions_from_pugixml(); 

	// ingest_glu_cell_idx = cell_defaults.custom_data.find_variable_index( "ingest_glu" ); 
	energy_cell_idx = cell_defaults.custom_data.find_variable_index( "energy" ); 
	// energy_cell_idx = cell_defaults->custom_data.find_variable_index( name );
	// energy_cell_idx = pCD->custom_data.find_variable_index( name );
	std::cout << "\n\n-------- create_cell_types():  energy_cell_idx = " << energy_cell_idx << std::endl;
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_tissue( void )
{
#ifdef LIBROADRUNNER
	// extern SBMLDocument_t *sbml_doc;
	rrc::RRVectorPtr vptr;
    rrc::RRCDataPtr result;  // start time, end time, and number of points
#endif

	Cell* pC;
	// float xval = -100.0;
	double xmin = default_microenvironment_options.X_range[0]; 
	double xval = xmin + 20.;

	std::cout << "\n---------- setup_tissue() -----------\n";

	std::cerr << "------------->>>>>  celltype1\n";
	pC = create_cell( get_cell_definition("celltype1") ); 
	double yval = -xmin / 2.0;
	pC->assign_position( xval, yval, 0.0 ); 

#ifdef LIBROADRUNNER
	std::cerr << "------------->>>>>  Creating rrHandle, loadSBML file\n\n";
	// std::cerr << "------------->>>>>  SBML file = " << cell_type_param.sbml_filename << std::endl;
	std::cerr << "------------->>>>>  SBML file = " << get_cell_definition("celltype1").sbml_filename << std::endl;
	rrc::RRHandle rrHandle = createRRInstance();
	// cell_defaults.phenotype.motility.persistence_time = parameters.doubles("persistence_time"); 
	// if (!rrc::loadSBML (rrHandle, "../Toy_Model_for_PhysiCell_1.xml")) {
	if (!rrc::loadSBML (rrHandle, get_cell_definition("celltype1").sbml_filename.c_str())) {
		std::cerr << "------------->>>>>  Error while loading SBML file  <-------------\n\n";
	// 	printf ("Error message: %s\n", getLastError());
	// 	getchar ();
	// 	exit (0);
	}
	pC->phenotype.molecular.model_rr = rrHandle;  // assign the intracellular model to each cell
	int r = rrc::getNumberOfReactions(rrHandle);
	int m = rrc::getNumberOfFloatingSpecies(rrHandle);
	int b = rrc::getNumberOfBoundarySpecies(rrHandle);
	int p = rrc::getNumberOfGlobalParameters(rrHandle);
	int c = rrc::getNumberOfCompartments(rrHandle);

	std::cerr << "Number of reactions = " << r << std::endl;
	std::cerr << "Number of floating species = " << m << std::endl;  // 4
	std::cerr << "Number of boundary species = " << b << std::endl;  // 0
	std::cerr << "Number of compartments = " << c << std::endl;  // 1

	std::cerr << "Floating species names:\n";
	std::cerr << "-----------------------\n";
	std::cerr << stringArrayToString(rrc::getFloatingSpeciesIds(rrHandle)) <<"\n"<< std::endl;

	vptr = rrc::getFloatingSpeciesConcentrations(rrHandle);
	std::cerr << vptr->Count << std::endl;
	for (int kdx=0; kdx<vptr->Count; kdx++)
		std::cerr << kdx << ") " << vptr->Data[kdx] << std::endl;
#endif

	//-------------------------
	std::cerr << "------------->>>>> celltype2\n";
	pC = create_cell( get_cell_definition("celltype2") ); 
	yval = -yval;
	pC->assign_position( xval, yval, 0.0 ); 

#ifdef LIBROADRUNNER
	std::cerr << "------------->>>>>  Creating rrHandle, loadSBML file\n\n";
	std::cerr << "------------->>>>>  SBML file = " << get_cell_definition("celltype2").sbml_filename << std::endl;
	rrHandle = createRRInstance();
	// if (!rrc::loadSBML (rrHandle, "../Toy_Model_for_PhysiCell_2.xml")) {
	if (!rrc::loadSBML (rrHandle, get_cell_definition("celltype2").sbml_filename.c_str())) {
		std::cerr << "------------->>>>>  Error while loading SBML file  <-------------\n\n";
		exit (0);
	}
	pC->phenotype.molecular.model_rr = rrHandle;  // assign the intracellular model to each cell
#endif

	return; 
}

void energy_based_cell_phenotype(Cell* pCell, Phenotype& phenotype , double dt)
{
	static int idx_glucose = 1;
	static int idx_oxygen = 3;

	std::cout << "------ energy_based_cell_phenotype ------, cell ID= " << pCell->ID << std::endl;

#ifdef LIBROADRUNNER
	rrc::RRVectorPtr vptr;
	rrc::RRCDataPtr result;  // start time, end time, and number of points
	// pC->phenotype.molecular.model_rr = rrHandle;  // assign the intracellular model to each cell
	vptr = rrc::getFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr);
	// std::cout << "--- before updating:" << std::endl;
	// for (int idx=0; idx<vptr->Count; idx++)
	// 	std::cout << idx << ", " << vptr->Data[idx] << std::endl;

	// // vptr->Data[idx_oxygen] += 0.1;
	// // rrc::setFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr, vptr);

	// vptr = rrc::getFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr);
	// // std::cout << vptr->Count << std::endl;
	// std::cout << "--- after updating oxygen:" << std::endl;
	// for (int idx=0; idx<vptr->Count; idx++)
	// 	std::cout << idx << ", " << vptr->Data[idx] << std::endl;

	int oxygen_i = microenvironment.find_density_index( "oxygen" ); 
	int glucose_i = microenvironment.find_density_index( "glucose" ); 
	// int energy_i = microenvironment.find_density_index( "energy" ); 
	int vi = microenvironment.nearest_voxel_index(pCell->position);
	double oxy_val = microenvironment(vi)[oxygen_i];
	double glucose_val = microenvironment(vi)[glucose_i];
	std::cout << "oxy_val at voxel of cell = " << oxy_val << std::endl;
	std::cout << "glucose_val at voxel of cell = " << glucose_val << std::endl;

	vptr->Data[idx_oxygen] = oxy_val;
	vptr->Data[idx_glucose] = glucose_val;
	rrc::setFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr, vptr);

	result = rrc::simulateEx (pCell->phenotype.molecular.model_rr, 0, 10, 10);  // start time, end time, and number of points
	int index = 0;
	// Print out column headers... typically time and species.
	for (int col = 0; col < result->CSize; col++)
	{
		// std::cout << result->ColumnHeaders[index++];
		std::cout << std::left << std::setw(15) << result->ColumnHeaders[index++];
		// if (col < result->CSize - 1)
		// {
		// 	// std::cout << "\t";
		// 	std::cout << "  ";
		// }
	}
	std::cout << "\n";

	index = 0;
	// Print out the data
	for (int row = 0; row < result->RSize; row++)
	{
		for (int col = 0; col < result->CSize; col++)
		{
			// std::cout << result->Data[index++];
			std::cout << std::left << std::setw(15) << result->Data[index++];
			// if (col < result->CSize -1)
			// {
			// 	// std::cout << "\t";
			// 	std::cout << "  ";
			// }
		}
		std::cout << "\n";
	}
	int idx = (result->RSize - 1) * result->CSize + 1;
	std::cout << "Saving last energy value (cell custom var) = " << result->Data[idx] << std::endl;
	pCell->custom_data[energy_cell_idx]  = result->Data[idx];
#endif

}

std::vector<std::string> energy_coloring_function( Cell* pCell )
{
	// color 0: cytoplasm fill 
	// color 1: outer outline 
	// color 2: nuclear fill 
	// color 3: nuclear outline 
	std::cout << "\n--- coloring fn: cell ID, energy = " << pCell->ID <<", "<< pCell->custom_data[energy_cell_idx] << std::endl; 
	std::cout << "-------- coloring fn():  energy_cell_idx = " << energy_cell_idx << std::endl;
	
	std::vector< std::string > output( 4, "white" ); 

	if (pCell->custom_data[energy_cell_idx] > 1.8)
		output[0] = "rgb(255,0,0)";
	else if (pCell->custom_data[energy_cell_idx] > 1.6)
		output[0] = "rgb(128,0,0)";
	else if (pCell->custom_data[energy_cell_idx] > 1.3)
		output[0] = "rgb(0,255,0)";
	else if (pCell->custom_data[energy_cell_idx] > 0.9)
		output[0] = "rgb(0,128,0)";
	else if (pCell->custom_data[energy_cell_idx] > 0.0)
		output[0] = "rgb(0,0,255)";
	else 
		output[0] = "rgb(0,0,128)";

	if (pCell->is_out_of_domain)
		output[0] = "rgb(128,128,128)";
	else if (!pCell->is_movable)
		output[0] = "rgb(0,0,0)";
	
	return output; 
}
// void celltype1_rule(Cell* pCell, Phenotype& phenotype , double dt)
// {
// 	std::cout << "------ celltype1_rule:";
// }
// void celltype2_rule(Cell* pCell, Phenotype& phenotype , double dt)
// {
// 	std::cout << "-- celltype2_rule:";
// }