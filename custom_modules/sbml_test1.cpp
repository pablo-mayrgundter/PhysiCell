
#include "./sbml_test1.h"

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
// #include <vector>
#include <string>

// extern Cell_Type_Parameters cell_type_param; 

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
	double xmin = -750.;
	int nregions = 5;
	double xdel = 1500./nregions;
	// 5 regions across x: [-750:-450:-150:150:450:750]
	std::cout << "setup_microenvironment: num voxels= " << microenvironment.number_of_voxels() << std::endl;
	for( int n=0; n < microenvironment.number_of_voxels(); n++ )
	{
		// x coordinate of the nth voxel's center
		x = microenvironment.mesh.voxels[n].center[0];
		for( int iregion=1; iregion <= nregions; iregion++ )
		{
			if (x < (xmin + iregion*xdel))
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
	// cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	cell_defaults.functions.volume_update_function = NULL;
	cell_defaults.functions.update_velocity = NULL;

	cell_defaults.functions.update_migration_bias = NULL; 
	// cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.update_phenotype = 	energy_based_cell_phenotype; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	// if( parameters.bools("predators_eat_prey") == true )
	// { get_cell_definition("predator").functions.custom_cell_rule = predator_hunting_function; }

	// if( parameters.bools("predators_cycle_if_big") == true )
	// { get_cell_definition("predator").functions.update_phenotype = predator_cycling_function; }

	// if( parameters.bools("prey_quorom_effect") == true )
	// { get_cell_definition("prey").functions.update_phenotype = prey_cycling_function; }
		
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
	// float xval = -600.0;
	// float yval = 1000.0;
	float xval = -100.0;

	std::cout << "\n---------- setup_tissue() -----------\n";

	// create just 3 cells, equally spaced in y; they'll migrate left-to-right
	// for (int idx=0; idx<3; idx++) {

		std::cerr << "------------->>>>>  celltype1\n";
		pC = create_cell( get_cell_definition("celltype1") ); 
		int yval = 100.;
		pC->assign_position( xval, yval, 0.0 ); 

		// pC->set_total_volume( pC->get_total_volume() * 3.0); 

		// Model_t *mm = SBMLDocument_getModel(sbml_doc);
		// std::cout << "mm =" << mm << std::endl;
		// pC->phenotype.molecular.molecular_model = mm;  // assign the intracellular model to each cell

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
		yval = -100.;
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

	// }
	return; 
}

// cell_defaults.functions.update_phenotype = energy_based_cell_phenotype; 
void energy_based_cell_phenotype(Cell* pCell, Phenotype& phenotype , double dt)
{
	static int idx_glucose = 1;
	static int idx_oxygen = 3;

	std::cout << "------ energy_based_cell_phenotype ------" << std::endl;

#ifdef LIBROADRUNNER
	rrc::RRVectorPtr vptr;
	rrc::RRCDataPtr result;  // start time, end time, and number of points
	// pC->phenotype.molecular.model_rr = rrHandle;  // assign the intracellular model to each cell
	vptr = rrc::getFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr);
	std::cout << "--- before updating:" << std::endl;
	for (int idx=0; idx<vptr->Count; idx++)
		std::cout << idx << ", " << vptr->Data[idx] << std::endl;

	// vptr->Data[idx_oxygen] += 0.1;
	// rrc::setFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr, vptr);

	vptr = rrc::getFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr);
	// std::cout << vptr->Count << std::endl;
	std::cout << "--- after updating oxygen:" << std::endl;
	for (int idx=0; idx<vptr->Count; idx++)
		std::cout << idx << ", " << vptr->Data[idx] << std::endl;

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
	
	std::vector< std::string > output( 4, "white" ); 

	std::cout << "--- coloring fn: cell ID, energy = " << pCell->ID <<", "<< pCell->custom_data[energy_cell_idx] << std::endl; 
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
void celltype1_rule(Cell* pCell, Phenotype& phenotype , double dt)
{
	std::cout << "------ celltype1_rule:";
}
void celltype2_rule(Cell* pCell, Phenotype& phenotype , double dt)
{
	std::cout << "-- celltype2_rule:";
}