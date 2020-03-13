#include "MultiCellWorld.h"

#include "../../Utilities/bitmap_image.hpp"


std::shared_ptr < ParameterLink<int>> MultiCellWorld::multiCellXPL =
Parameters::register_parameter("WORLD_MULTICELL_GEO-multiCellX", 5,
	"size in X of a multiCell");
std::shared_ptr < ParameterLink<int>> MultiCellWorld::multiCellYPL =
Parameters::register_parameter("WORLD_MULTICELL_GEO-multiCellY", 5,
	"size in Y of a multiCell");
std::shared_ptr < ParameterLink<double>> MultiCellWorld::cellResourceCollectionMinPL =
Parameters::register_parameter("WORLD_MULTICELL-cellResourceCollectionMin", .9,
	"resource amount a cell collects every update");
std::shared_ptr < ParameterLink<double>> MultiCellWorld::cellResourceCollectionMaxPL =
Parameters::register_parameter("WORLD_MULTICELL-cellResourceCollectionMax", 1.1,
	"resource amount a cell collects every update");
std::shared_ptr < ParameterLink<double>> MultiCellWorld::cellReproCostPL =
Parameters::register_parameter("WORLD_MULTICELL-cellReproCost", 10.0,
	"resource needed by a cell to reproduce");

std::shared_ptr < ParameterLink<int>> MultiCellWorld::worldXPL =
Parameters::register_parameter("WORLD_MULTICELL_GEO-worldX", 5,
	"size in X of world");
std::shared_ptr < ParameterLink<int>> MultiCellWorld::worldYPL =
Parameters::register_parameter("WORLD_MULTICELL_GEO-worldY", 5,
	"size in Y of world");

std::shared_ptr < ParameterLink<bool>> MultiCellWorld::randomGermPL =
Parameters::register_parameter("WORLD_MULTICELL-randomGerm", false,
	"when multi cells reproduce, should the germ cell or a random cell from the multi cell be used.");

std::shared_ptr < ParameterLink<bool>> MultiCellWorld::overwriteSelfOnFailPL =
Parameters::register_parameter("WORLD_MULTICELL-overwriteSelfOnFail", false,
	"when good cells reproduces and the target is not empty, place offspring over good cell.");

std::shared_ptr < ParameterLink<int>> MultiCellWorld::recordImageStepPL =
Parameters::register_parameter("WORLD_MULTICELL_IO-recordImageStep", -1,
	"record world images when update % recordImageStep == 0, if -1, do not record images.");

std::shared_ptr < ParameterLink<int>> MultiCellWorld::recordWorldStateStepPL =
Parameters::register_parameter("WORLD_MULTICELL_IO-recordWorldStateStep", -1,
	"record world state csv when update % recordWorldStateStep == 0, if -1, do not record status.");

std::shared_ptr < ParameterLink<double>> MultiCellWorld::initEvilPercentPL =
Parameters::register_parameter("WORLD_MULTICELL-initEvilPercent", -1.0,
	"use this to initalize first population genomes, if -1, seed genomes randomly, 0.0 = not evil, 1.0 = very evil");



////////////////////////////////////////////////////////////////////////
// spatial parameters
////////////////////////////////////////////////////////////////////////
std::shared_ptr < ParameterLink<bool>> MultiCellWorld::spatialWorldPL =
Parameters::register_parameter("WORLD_MULTICELL_SPATIAL-spatialWorld", false,
	"is the world spatial?");

std::shared_ptr < ParameterLink<bool>> MultiCellWorld::spatialMultiCellPL =
Parameters::register_parameter("WORLD_MULTICELL_SPATIAL-spatialMultiCell", false,
	"are the MultiCells spatial?");

std::shared_ptr < ParameterLink<int>> MultiCellWorld::spatialWorldDistPL =
Parameters::register_parameter("WORLD_MULTICELL_SPATIAL-spatialWorldDist", 1,
	"if world is spatial how far are offspring MultiCells distributed from parent on birth?");

std::shared_ptr < ParameterLink<int>> MultiCellWorld::spatialMultiCellDistPL =
Parameters::register_parameter("WORLD_MULTICELL_SPATIAL-spatialMultiCellDist", 1,
	"if MultiCells are spatial how far are offspring cells distributed from parent on birth?");

std::shared_ptr < ParameterLink<std::string>> MultiCellWorld::spatialWorldModelPL =
Parameters::register_parameter("WORLD_MULTICELL_SPATIAL-spatialWorldModel", (std::string)"square",
	"if world is spatial what model is used to place offspring MultiCells? (square or cardinal)");

std::shared_ptr < ParameterLink<std::string>> MultiCellWorld::spatialMultiCellModelPL =
Parameters::register_parameter("WORLD_MULTICELL_SPATIAL-spatialMultiCellModel", (std::string)"square",
	"if MultiCells are spatial what model is used to place offspring cells? (square or cardinal)");

std::shared_ptr < ParameterLink<std::string>> MultiCellWorld::spatialWorldEdgeRulePL =
Parameters::register_parameter("WORLD_MULTICELL_SPATIAL-spatialWorldEdgeRule", (std::string)"search",
	"if world is spatial what happens when offspring would be placed off the edge of the world? (fail, wrap or, search)");

std::shared_ptr < ParameterLink<std::string>> MultiCellWorld::spatialMultiCellEdgeRulePL =
Parameters::register_parameter("WORLD_MULTICELL_SPATIAL-spatialMultiCellEdgeRule", (std::string)"search",
	"if MultiCell is spatial happens when an offspring would be placed off the edge of the MultiCell? (fail, wrap or, search)");

std::shared_ptr < ParameterLink<double>> MultiCellWorld::evilCutoffMinPL =
Parameters::register_parameter("WORLD_MULTICELL-evilCutoffMin", .5,
	"");
std::shared_ptr < ParameterLink<double>> MultiCellWorld::evilCutoffMaxPL =
Parameters::register_parameter("WORLD_MULTICELL-evilCutoffMax", .5,
	"cell will be evil if cell alignment > evilCutoffMax.\ncell will be good if cell alignment < evilCutoffMin.\nif cell alignment is between evilCutoffMin and evilCutoffMax, cell will probabilisticly be evil");

std::shared_ptr < ParameterLink<double>> MultiCellWorld::alignmentMutationRatePL =
Parameters::register_parameter("WORLD_MULTICELL-alignmentMutationRate", .01,
	"per site mutation rate of alignment genome");

std::shared_ptr < ParameterLink<int>> MultiCellWorld::alignmentGenomeSizePL =
Parameters::register_parameter("WORLD_MULTICELL-alignmentGenomeSize", 10,
	"size of genome used to determin genome");

std::shared_ptr < ParameterLink<double>> MultiCellWorld::initWorldGenomeEvilPercentPL =
Parameters::register_parameter("WORLD_MULTICELL-initWorldGenomeEvilPercent", .5,
	"if useWorldGenome, an un-mutating genome (world genome) will be used as germ for all MultiCells.\n0.0 = not evil, 1.0 = very evil");
std::shared_ptr < ParameterLink<bool>> MultiCellWorld::useWorldGenomePL =
Parameters::register_parameter("WORLD_MULTICELL-useWorldGenome", false,
	"if useWorldGenome, an un-mutating genome (world genome) will be used as germ for all MultiCells.");

std::shared_ptr < ParameterLink<bool>> MultiCellWorld::suppressMultiCellOffspringPL =
Parameters::register_parameter("WORLD_MULTICELL-suppressMultiCellOffspring", false,
	"if true, MultiCells will overwrite self when they reproduce, but will not produce a sister MultiCell\n...used to measure MC fill times, without world level competition.");

void MultiCellWorld::saveWorldImage() {
	double cX = 1000;
	double cY = 1000;

	cartesian_canvas canvas(cX,cY);

	canvas.pen_color(0, 0, 0);
	canvas.fill_rectangle(-1 * (cX / 2), -1 * (cY / 2), cX/2, cY/2);

	double multiCellWidth = (cX / worldX);
	double multiCellHeight = (cY / worldY);

	double cellWidth = (multiCellWidth / multiCellX);
	double cellHeight = (multiCellHeight / multiCellY);

	double x1, x2, y1, y2;
	double multiCellX1, multiCellY1;
	for (int wY = 0; wY < worldY; wY++) {
		for (int wX = 0; wX < worldX; wX++) {
			if (!world(wX, wY).empty) {
				x1 = (wX * multiCellWidth) - (cX / 2);
				y1 = (wY * multiCellHeight) - (cY / 2);
				x2 = x1 + multiCellWidth;
				y2 = y1 + multiCellHeight;

				canvas.pen_color(255, 255, 255);
				canvas.fill_rectangle(x1,y1,x2,y2);
				//std::cout << "draw multi cell : " << x1 << ", " << y1 << ", " << x2 << ", " << y2 << std::endl;
				multiCellX1 = wX * multiCellWidth;
				multiCellY1 = wY * multiCellHeight;

				for (int y = 0; y < multiCellY; y++) {
					for (int x = 0; x < multiCellX; x++){
						if (!world(wX, wY).cells(x, y).empty) {
							x1 = multiCellX1 + x * cellWidth - (cX / 2);
							y1 = multiCellY1 + y * cellHeight - (cY / 2);
							x2 = x1 + cellWidth;
							y2 = y1 + cellHeight;
							if (world(wX, wY).cells(x, y).evil && world(wX, wY).cells(x, y).first) {
								canvas.pen_color(255, 0, 255);
							}
							else if (world(wX, wY).cells(x, y).murder) {
								canvas.pen_color(125, 0, 0);
							}
							else if (world(wX, wY).cells(x, y).evil) {
								canvas.pen_color(255, 0, 0);
							}
							else if (world(wX, wY).cells(x, y).first) {
								canvas.pen_color(0, 0, 255);
							}
							else {
								canvas.pen_color(0, 255, 0);
							}
							canvas.fill_rectangle(x1, y1, x2, y2);
							//std::cout << "  draw cell : " << x1 << ", " << y1 << ", " << x2 << ", " << y2 << std::endl;

						}
					}
				}
			}
		}
	}

	// add a grid
	canvas.pen_color(0, 0, 0);
	double pos = -1 * (cX / 2);
	while (pos < (cX / 2)) {
		canvas.fill_rectangle(pos-1, -1 * (cY / 2), pos, (cY / 2));
		pos += multiCellWidth;
	}
	pos = -1 * (cY / 2);
	while (pos < (cY / 2)) {
		canvas.fill_rectangle(-1 * (cX / 2), pos-1, (cX / 2), pos);
		pos += multiCellHeight;
	}
	canvas.image().save_image("output_"+std::to_string(100000+Global::update)+".bmp");
}

MultiCellWorld::MultiCellWorld(std::shared_ptr<ParametersTable> PT_)
	: AbstractWorld(PT_) {


	std::cout << "in MultiCellWorld constructor :: " << std::endl;
	multiCellX = multiCellXPL->get(PT);
	multiCellY = multiCellYPL->get(PT);
	cellResourceCollectionMin = cellResourceCollectionMinPL->get(PT);
	cellResourceCollectionMax = cellResourceCollectionMaxPL->get(PT);

	worldX = worldXPL->get(PT);
	worldY = worldYPL->get(PT);

	cellReproCost = cellReproCostPL->get(PT);

	randomGerm = randomGermPL->get(PT);

	overwriteSelfOnFail = overwriteSelfOnFailPL->get(PT);

	recordImageStep = recordImageStepPL->get(PT);
	recordWorldStateStep = recordWorldStateStepPL->get(PT);

	initEvilPercent = initEvilPercentPL->get(PT);




	spatialWorld=spatialWorldPL->get(PT);
	spatialMultiCell= spatialMultiCellPL->get(PT);
	spatialWorldDist = spatialWorldDistPL->get(PT);
	spatialMultiCellDist = spatialMultiCellDistPL->get(PT);

	if (spatialWorldModelPL->get(PT) == "cardinal") {
		spatialWorldModel = "cardinal";
	}
	else if (spatialWorldModelPL->get(PT) == "square") {
		spatialWorldModel = "square";
	}
	else {
		std::cout << "  in MultiCellWorld constructor :: parameter spatialWorldModel has invalid value \"" << spatialWorldModelPL->get(PT) << "\". Please correct and retry.\n  exiting..." << std::endl;
		exit(1);
	}

	if (spatialMultiCellModelPL->get(PT) == "cardinal") {
		spatialMultiCellModel = "cardinal";
	}
	else if (spatialMultiCellModelPL->get(PT) == "square") {
		spatialMultiCellModel = "square";
	}
	else {
		std::cout << "  in MultiCellWorld constructor :: parameter spatialMultiCellModel has invalid value \"" << spatialMultiCellModelPL->get(PT) << "\". Please correct and retry.\n  exiting..." << std::endl;
		exit(1);
	}

	if (spatialWorldEdgeRulePL->get(PT) == "fail") {
		spatialWorldEdgeRule = "fail";
	}
	else if (spatialWorldEdgeRulePL->get(PT) == "wrap") {
		spatialWorldEdgeRule = "wrap";
	}
	else if (spatialWorldEdgeRulePL->get(PT) == "search") {
		spatialWorldEdgeRule = "search";
	}
	else {
		std::cout << "  in MultiCellWorld constructor :: parameter spatialWorldEdgeRule has invalid value \"" << spatialWorldEdgeRulePL->get(PT) << "\". Please correct and retry.\n  exiting..." << std::endl;
		exit(1);
	}

	if (spatialMultiCellEdgeRulePL->get(PT) == "fail") {
		spatialMultiCellEdgeRule = "fail";
	}
	else if (spatialMultiCellEdgeRulePL->get(PT) == "wrap") {
		spatialMultiCellEdgeRule = "wrap";
	}
	else if (spatialMultiCellEdgeRulePL->get(PT) == "search") {
		spatialMultiCellEdgeRule = "search";
	}
	else {
		std::cout << "  in MultiCellWorld constructor :: parameter spatialMultiCellEdgeRule has invalid value \"" << spatialMultiCellEdgeRulePL->get(PT) << "\". Please correct and retry.\n  exiting..." << std::endl;
		exit(1);
	}

	evilCutoff.first = evilCutoffMinPL->get(PT);
	evilCutoff.second = evilCutoffMaxPL->get(PT);
	if (evilCutoff.first > evilCutoff.second) {
		std::cout << "  in MultiCellWorld constructor :: parameter evilCutoffMin > evilCutoffMax. Please correct and retry.\n  exiting..." << std::endl;
		exit(1);
	}
	if (evilCutoff.first == evilCutoff.second) {
		evilCutoff.second = -1;
	}

	alignmentMutationRate = alignmentMutationRatePL->get(PT);
	alignmentGenomeSize = alignmentGenomeSizePL->get(PT);

	// make lists of world and multicell locations
	// these will be randomized every generation to
	// allow the order cells are visted in to be randomized
	for (int yy = 0; yy < worldY; yy++) {
		for (int xx = 0; xx < worldY; xx++) {
			worldOrder.push_back({ xx,yy });
		}
	}
	for (int yy = 0; yy < multiCellY; yy++) {
		for (int xx = 0; xx < multiCellX; xx++) {
			multiCellOrder.push_back({ xx,yy });
		}
	}


	initWorldGenomeEvilPercent = initWorldGenomeEvilPercentPL->get(PT);
	useWorldGenome = useWorldGenomePL->get(PT);

	worldGenome.resize(alignmentGenomeSize);
	double s = worldGenome.size() * initWorldGenomeEvilPercent;
	for (double i = 0; i < worldGenome.size(); i++) {
		if (i < s) {
			worldGenome[i] = 1;
		}
		else {
			worldGenome[i] = 0;
		}
	}
	std::cout << "Using WorldGneome: " << useWorldGenome << "   World Genome has been initialized to: ";
	for (auto v : worldGenome) {
		std::cout << v << " ";
	}
	std::cout << std::endl;

	suppressMultiCellOffspring = suppressMultiCellOffspringPL->get(PT);

	// columns to be added to ave file
	popFileColumns.clear();
	//popFileColumns.push_back("germ");
	//popFileColumns.push_back("evil");
	//popFileColumns.push_back("murder");
	//popFileColumns.push_back("alignment");
	//}

}


void MultiCellWorld::evaluate(std::map<std::string, std::shared_ptr<Group>>& groups, int analyze, int visualize, int debug) {

	std::cout << "  Starting evaluation in MultiCellWorld." << std::endl;

	world.reset(worldX, worldY);

	if ((worldX * worldY) < groups["root::"]->population.size()) {
		std::cout << "  worldX * worldY < inital population size. Either decrease inital population size or increase world size.\n  exiting..." << std::endl;
		exit(1);
	}
	auto tempPop = groups["root::"]->population;

	for (auto org : tempPop) {
		// place this multi cell somewhere
		int x, y;
		do {
			x = Random::getIndex(worldX);
			y = Random::getIndex(worldY);
		} while (!world(x, y).empty);

		// make a new MultiCell form org
		std::vector<bool> genome;
		
		initGenome(genome);

		if (useWorldGenome) {
			genome = worldGenome;
		}

		birthMultiCell(groups, org, genome, x, y); // -1,-1 = no parent

		// place org in kill list (we don't need it any more
		killList.insert(org);
	}

	std::cout << "done setting up initial MultiCells..." << std::endl;
	std::cout << "pop size: " << groups["root::"]->population.size() << std::endl;

	// run the evaluation
	do {
		// every generation, randomize the order we will visit both world locations and multiCell locations
		std::shuffle(worldOrder.begin(), worldOrder.end(), Random::Generator());
		std::shuffle(multiCellOrder.begin(), multiCellOrder.end(), Random::Generator());

		for (auto worldLocation : worldOrder) {
			int wX = worldLocation.first;
			int wY = worldLocation.second;

			if (!world(wX, wY).empty) {
				for (auto MCLocation : multiCellOrder) {
					int x = MCLocation.first;
					int y = MCLocation.second;
					if (!world(wX, wY).cells(x, y).empty && world(wX, wY).cells(x, y).org->timeOfBirth < Global::update) { // if this is a cell not empty and not brand new

						// gain resource
						world(wX, wY).cells(x, y).resource += Random::getDouble(cellResourceCollectionMin, cellResourceCollectionMax);

						if (world(wX, wY).cells(x, y).resource > cellReproCost) { // cell will try to make an offspring
							// attempt to make offspring
							birthCell(groups, wX, wY, x, y);

							//int numOffspring = Random::getBinomial(4, .33);
							//for (int i = 0; i < numOffspring; i++) {
							//	birthCell(groups, wX, wY, x, y);
							//}

						}
					}
				}
			}

			// are all cells full in this MC?
			// assume this MC is full
			world(wX, wY).full = true;
			// now look at each cell, if any are empty, the MC is not full
			for (int y = 0; y < world(wX, wY).cells.y(); y++) {
				for (int x = 0; x < world(wX, wY).cells.x(); x++) {
					if (world(wX, wY).cells(x, y).empty) {
						world(wX, wY).full = false; // if any cell is empty, this MC is not full
					}
				}
			}

			if (world(wX, wY).full) { // this is here so we can change the rule in the future...
				world(wX, wY).canProduce = true;
			}
		}

		// we are about to make new MCs, so randomize the order that we visit world locations
		std::shuffle(worldOrder.begin(), worldOrder.end(), Random::Generator());

		for (auto worldLocation : worldOrder) {
			int wX = worldLocation.first;
			int wY = worldLocation.second;
			// if cell canProduce, add mutiCell
			if (world(wX, wY).canProduce) {
				if (!world(wX, wY).empty && world(wX, wY).germ->timeOfBirth < Global::update) { // do not look at this MC if it is new
					
																								// this MC will reproduce this will destroy this MC (will be overwritten by self) and another MC
					// save some info about this MC
					addMCFate(world(wX, wY), 1); // 1 = death by reproducing over self

					// put all cells in this MC in kill list
					auto parentGenome = world(wX, wY).genome;
					killList.insert(world(wX, wY).germ);

					//if random germ, place a random cell in germ, so it will be used to make new MCs
					if (randomGerm) {
						int randX, randY;
						do {
							randX = Random::getIndex(multiCellX);
							randY = Random::getIndex(multiCellY);
						} while (world(wX, wY).cells(randX, randY).empty);

						world(wX, wY).germ = world(wX, wY).cells(randX, randY).org;
						// use this cells genome
						parentGenome = world(wX, wY).cells(randX, randY).genome;
					}

					// kill all cells in this MC
					for (int y = 0; y < multiCellX; y++) {
						for (int x = 0; x < multiCellY; x++) {
							if (!world(wX, wY).cells(x, y).empty) { // if this is empty...
								killList.insert(world(wX, wY).cells(x, y).org);
							}
						}
					}

					// make a new MC over the current MC
					std::shared_ptr<Organism> germ_org = world(wX, wY).germ->makeMutatedOffspringFrom(world(wX, wY).germ);
					
					if (useWorldGenome) {
						birthMultiCell(groups, germ_org, worldGenome, wX, wY);
					}
					else {
						birthMultiCell(groups, germ_org, parentGenome, wX, wY);
					}

					if (!suppressMultiCellOffspring) {
						// make a new MC at a random world location (not here)
						int newX, newY;
						pickGridLoc(true, wX, wY, newX, newY); // true = forWorld
						if (!world(newX, newY).empty) {

							// we have selected an MC that is not empty, it's about this be killed!
							//save some info about this MC before it gets destroyed
							addMCFate(world(newX, newY), 0); // 0 = death by overwrite from anther MC

							// kill everything in the target location if not empty
							if (!world(newX, newY).empty) {
								killList.insert(world(newX, newY).germ);
								for (int y = 0; y < multiCellX; y++) {
									for (int x = 0; x < multiCellY; x++) {
										if (!world(newX, newY).cells(x, y).empty) { // if this is not empty...
											killList.insert(world(newX, newY).cells(x, y).org);
										}
									}
								}
							}
						}
						germ_org = world(wX, wY).germ->makeMutatedOffspringFrom(world(wX, wY).germ);
						if (useWorldGenome) {
							birthMultiCell(groups, germ_org, worldGenome, newX, newX);
						}
						else {
							birthMultiCell(groups, germ_org, parentGenome, newX, newX);
						}
					}
				}
			}
		}

		// scan all cells to get MC and cell counts
		int MC_count = 0;
		int cellCount = 0;
		for (int wY = 0; wY < worldY; wY++) {
			for (int wX = 0; wX < worldX; wX++) {
				if (!world(wX, wY).empty) {
					MC_count++;
					for (int y = 0; y < world(wX, wY).cells.y(); y++) {
						for (int x = 0; x < world(wX, wY).cells.x(); x++) {
							if (!world(wX, wY).cells(x, y).empty) { // if this is empty...
								cellCount++;
							}
						}
					}
				}
			}
		}

		// save image of world
		if (recordImageStep > 0 && Global::update % recordImageStep == 0) {
			saveWorldImage();
		}

		// save some world details to record csv file
		if (recordWorldStateStep > 0 && Global::update % recordWorldStateStep == 0) {
			saveMCFateFile(); // save fates files

			double state_count_MC = 0; // this will also be count_all_germ
			double state_count_all_evil_germs = 0;
			double state_count_all_cells = 0;
			double state_count_all_evil_cells = 0;

			double state_count_MC_cells = 0;
			double state_count_MC_evil_cell = 0;

			std::vector<double> all_cell_alginments;
			std::vector<double> all_germ_alginments;
			std::vector<double> evil_rogue_rates; // % evil cell in MC with good germ
			std::vector<double> good_rogue_rates; // % good cell in MC with evil germ


			std::string outputLine = "";
			for (int wY = 0; wY < worldY; wY++) {
				for (int wX = 0; wX < worldX; wX++) {
					//if (!world(wX, wY).empty && world(wX, wY).germ->timeOfBirth != Global::update) {
					if (!world(wX, wY).empty) {
						state_count_MC++;
						bool evilGerm = world(wX, wY).evil;
						state_count_all_evil_germs += evilGerm;
						all_germ_alginments.push_back(world(wX, wY).alignment);
						state_count_MC_cells = 0;
						state_count_MC_evil_cell = 0;
						for (int y = 0; y < world(wX, wY).cells.y(); y++) {
							for (int x = 0; x < world(wX, wY).cells.x(); x++) {
								if (!world(wX, wY).cells(x, y).empty) { // if this is not empty...
									state_count_all_cells++;
									state_count_MC_cells++;
									if (world(wX, wY).cells(x, y).evil) {
										state_count_all_evil_cells++;
										state_count_MC_evil_cell++;
									}
									all_cell_alginments.push_back(world(wX, wY).cells(x, y).alignment);
								}
							}
						}
						if (!evilGerm) {
							evil_rogue_rates.push_back(state_count_MC_evil_cell / state_count_MC_cells);
						}
						else {
							good_rogue_rates.push_back(1.0 - (state_count_MC_evil_cell / state_count_MC_cells) );
						}
					}
				}
			}
			std::string headerString = "update,MC,cells,evilGermRate,evilCellRate,aveGermAlignments,aveCellAlignments,aveEvilRogueRate,aveGoodRogueRate";
			std::string dataString = std::to_string(Global::update) + ",";;
			dataString += std::to_string(state_count_MC) + ",";
			dataString += std::to_string(state_count_all_cells) + ",";
			dataString += std::to_string(state_count_all_evil_germs / state_count_MC) + ",";
			dataString += std::to_string(state_count_all_evil_cells / state_count_all_cells) + ",";

			dataString += std::to_string(std::accumulate(
				all_germ_alginments.begin(),
				all_germ_alginments.end(), 0.0)
				/ static_cast<double>(all_germ_alginments.size())) + ",";
			dataString += std::to_string(std::accumulate(
				all_cell_alginments.begin(),
				all_cell_alginments.end(), 0.0)
				/ static_cast<double>(all_cell_alginments.size())) + ",";

			dataString += std::to_string(std::accumulate(
				evil_rogue_rates.begin(),
				evil_rogue_rates.end(), 0.0)
				/ static_cast<double>(evil_rogue_rates.size())) + ",";
			dataString += std::to_string(std::accumulate(
				good_rogue_rates.begin(),
				good_rogue_rates.end(), 0.0)
				/ static_cast<double>(good_rogue_rates.size()));
			FileManager::writeToFile("MultiCellReport.csv", dataString, headerString);
		}

		// record population data
		groups["root::"]->archive();

		// decide which cells survive
		std::vector<std::shared_ptr<Organism>> newPopulation;
		for (size_t i = 0; i < groups["root::"]->population.size(); i++) {
			if (killList.find(groups["root::"]->population[i]) == killList.end()) { // if not in kill list
				newPopulation.push_back(groups["root::"]->population[i]);           // move into new population
			}
		}

		groups["root::"]->population = newPopulation;

		// now kill everyone in kill list so that the memory can be reused
		for (auto org : killList) {
			org->kill();
		}
		killList.clear();

		// print a short accounting to terminal
		std::cout << "finished update: " << Global::update << "  pop size: " << groups["root::"]->population.size() << "  MC count: " << MC_count << "  cell count: " << cellCount << std::endl;
		Global::update++;
	} while (!groups["root::"]->archivist->finished_);
	std::cout << "finished run!" << std::endl;
}

std::unordered_map<std::string, std::unordered_set<std::string>>
MultiCellWorld::requiredGroups() {
	//return { {"root::",{"G:root::"}} };
	return { {"root::",{}} };
	// this world needs orgs, with a genome.
}



