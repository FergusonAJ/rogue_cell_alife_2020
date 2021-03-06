//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

#pragma once

#include "../AbstractWorld.h"
#include "../../Genome/CircularGenome/CircularGenome.h"

#include <cstdlib>
#include <thread>
#include <iostream>
#include <iomanip>

#include <vector>

class MultiCellWorld : public AbstractWorld {

public:

	static std::shared_ptr < ParameterLink<int>> multiCellXPL;
	static std::shared_ptr < ParameterLink<int>> multiCellYPL;
	int multiCellX;
	int multiCellY;
	static std::shared_ptr < ParameterLink<double>> cellResourceCollectionMinPL;
	static std::shared_ptr < ParameterLink<double>> cellResourceCollectionMaxPL;
	double cellResourceCollectionMin;
	double cellResourceCollectionMax;
	static std::shared_ptr < ParameterLink<double>> cellReproCostPL;
	double cellReproCost;

	static std::shared_ptr < ParameterLink<int>> worldXPL;
	static std::shared_ptr < ParameterLink<int>> worldYPL;
	int worldX;
	int worldY;

	static std::shared_ptr < ParameterLink<bool>> randomGermPL;
	bool randomGerm;
	static std::shared_ptr < ParameterLink<bool>> overwriteSelfOnFailPL;
	bool overwriteSelfOnFail;

	static std::shared_ptr < ParameterLink<double>> initEvilPercentPL;
	double initEvilPercent;

	static std::shared_ptr < ParameterLink<int>> recordImageStepPL;
	int recordImageStep;
	static std::shared_ptr < ParameterLink<int>> recordWorldStateStepPL;
	int recordWorldStateStep;


	static std::shared_ptr < ParameterLink<bool>> spatialWorldPL;
	bool spatialWorld;
	static std::shared_ptr < ParameterLink<bool>> spatialMultiCellPL;
	bool spatialMultiCell;
	static std::shared_ptr < ParameterLink<int>> spatialWorldDistPL;
	int spatialWorldDist;
	static std::shared_ptr < ParameterLink<int>> spatialMultiCellDistPL;
	int spatialMultiCellDist;
	static std::shared_ptr < ParameterLink<std::string>> spatialWorldModelPL;
	std::string spatialWorldModel;
	static std::shared_ptr < ParameterLink<std::string>> spatialMultiCellModelPL;
	std::string spatialMultiCellModel;
	static std::shared_ptr < ParameterLink<std::string>> spatialWorldEdgeRulePL;
	std::string spatialWorldEdgeRule;
	static std::shared_ptr < ParameterLink<std::string>> spatialMultiCellEdgeRulePL;
	std::string spatialMultiCellEdgeRule;

	static std::shared_ptr < ParameterLink<double>> evilCutoffMinPL;
	static std::shared_ptr < ParameterLink<double>> evilCutoffMaxPL;
	std::pair<double,double> evilCutoff;
	
	static std::shared_ptr < ParameterLink<double>> alignmentMutationRatePL;
	double alignmentMutationRate;
	static std::shared_ptr < ParameterLink<int>> alignmentGenomeSizePL;
	int alignmentGenomeSize;

	static std::shared_ptr < ParameterLink<double>> initWorldGenomeEvilPercentPL;
	double initWorldGenomeEvilPercent;
	static std::shared_ptr < ParameterLink<bool>> useWorldGenomePL;
	bool useWorldGenome = 0;
	std::vector<bool> worldGenome;

	static std::shared_ptr < ParameterLink<bool>> suppressMultiCellOffspringPL;
	bool suppressMultiCellOffspring;

	std::vector<std::pair<int,int>> worldOrder, multiCellOrder;

	// Vector2d wraps a vector<T> and provides (x,y) style access
	// no error checking is provided for out of range errors
	// internally this class uses R(ow) and C(olumn) (i.e. how the data is stored in the data vector)
	// the user sees (x,y) where x = column, y = row
	template <typename T> class Vector2d {
		std::vector<T> data;
		int R, C;

		// get index into data vector for a given x,y
		inline int getIndex(int r, int c) { return (r * C) + c; }

	public:
		Vector2d() {
			R = 0;
			C = 0;
		}
		// construct a vector of size x * y
		Vector2d(int x, int y) : R(y), C(x) { data.resize(R * C); }

		Vector2d(int x, int y, T value) : R(y), C(x) { data.resize(R * C, value); }

		void reset(int x, int y) {
			R = y;
			C = x;
			data.clear();
			data.resize(R * C);
		}

		void reset(int x, int y, T value) {
			R = y;
			C = x;
			data.clear();
			data.resize(R * C, value);
		}

		// overwrite this classes data (vector<T>) with data coppied from newData
		void assign(std::vector<T> newData) {
			if ((int)newData.size() != R * C) {
				std::cout << "  ERROR :: in Vector2d::assign() vector provided does not "
					"fit. provided vector is size "
					<< newData.size() << " but Rows(" << R << ") * Columns(" << C
					<< ") == " << R * C << ". Exitting." << std::endl;
				exit(1);
			}
			data = newData;
		}

		// provides access to value x,y can be l-value or r-value (i.e. used for
		// lookup of assignment)
		T& operator()(int x, int y) { return data[getIndex(y, x)]; }

		T& operator()(double x, double y) {
			return data[getIndex((int)(y), (int)(x))];
		}

		T& operator()(std::pair<int, int> loc) {
			return data[getIndex(loc.second, loc.first)];
		}

		T& operator()(std::pair<double, double> loc) {
			return data[getIndex((int)(loc.second), (int)(loc.first))];
		}

		// show the contents of this Vector2d with index values, and x,y values
		void show() {
			for (int r = 0; r < R; r++) {
				for (int c = 0; c < C; c++) {
					std::cout << getIndex(r, c) << " : " << c << "," << r << " : "
						<< data[getIndex(r, c)] << "\n";
				}
			}
		}

		// show the contents of this Vector2d in a grid
		void showGrid(int precision = -1) {
			if (precision < 0) {
				for (int r = 0; r < R; r++) {
					for (int c = 0; c < C; c++) {
						std::cout << data[getIndex(r, c)] << " ";
					}
					std::cout << "\n";
				}
			}
			else {
				for (int r = 0; r < R; r++) {
					for (int c = 0; c < C; c++) {
						if (data[getIndex(r, c)] == 0) {
							std::cout << std::setfill(' ') << std::setw((precision * 2) + 2)
								<< " ";
						}
						else {
							std::cout << std::setfill(' ') << std::setw((precision * 2) + 1)
								<< std::fixed << std::setprecision(precision)
								<< data[getIndex(r, c)] << " ";
						}
					}
					std::cout << "\n";
				}
			}
		}
		int x() { return C; }

		int y() { return R; }
	};

	class Cell { // a location in a MC
	public:
		std::shared_ptr<Organism> org;
		double resource = 0;
		bool empty = true;    // is there a cell in this location
		bool murder = false;  // set true if this cell replaced a cell on birth
		bool first = false;   // set true if this cell is the first cell in an MC

		double alignment = 0; // derived from genome
		bool evil = false;    // set true if this cell is evil

		std::vector<bool> genome;
	};

	class MultiCell {
	public:
		std::shared_ptr<Organism> germ;
		Vector2d<Cell> cells;
		bool empty = true;
		bool full = false;
		bool canProduce = false;

		double alignment = 0; // derived from germ genome
		bool evil = 0;

		int goodCellBirthCounts = 0;
		int evilCellBirthCounts = 0; // includes murder
		int murderCellBirthCounts = 0;
		double cellAlignmentSum = 0.0;

		std::vector<bool> genome;
	};

	int MCReportID = 0;
	std::string MCReportHeader = "reportID,timeOfBirth,ageAtDeath,goodCellBirths,evilCellBirths,murderCellBirths,cellAlignmentAve,germAlignment,reproDeath,overwriteDeath";
	std::string MCReport = "";

	// call when an MC is being destroyed to record some data about that MC. reproDeath = 1 if this death is from a MC reproducing over itself and 0 if overwritten by some other MC
	void addMCFate(MultiCell& MC,int reproDeath) {
		if (MCReport.size() != 0) { // there is already data in MCReport, add a new line to the last line
			MCReport += "\n";
		}
		MCReport += std::to_string(MCReportID++) + ",";
		MCReport += std::to_string(MC.germ->timeOfBirth) + ",";
		MCReport += std::to_string(Global::update-MC.germ->timeOfBirth) + ",";
		MCReport += std::to_string(MC.goodCellBirthCounts) + ",";
		MCReport += std::to_string(MC.evilCellBirthCounts) + ",";
		MCReport += std::to_string(MC.murderCellBirthCounts) + ",";
		MCReport += std::to_string(MC.cellAlignmentSum/static_cast<double>(MC.goodCellBirthCounts + MC.evilCellBirthCounts)) + ",";
		MCReport += std::to_string(MC.alignment) + ",";
		MCReport += std::to_string(reproDeath) + ",";
		MCReport += std::to_string(1 - reproDeath);
	}

	// write all MC reports to file
	void saveMCFateFile() {
		if (MCReport.size() > 0) { // if the string is not empty
			FileManager::writeToFile("MC_Fates.csv", MCReport, MCReportHeader);
			MCReport = ""; // clear out the string
		}
	}

	double getAlignment(Cell& cell) {
		return(static_cast<double>(std::accumulate(cell.genome.begin(), cell.genome.end(), 0.0)) / static_cast<double>(cell.genome.size()));
	}

	double getAlignment(MultiCell& multiCell) {
		/*
		std::cout << "MC genome: ";
		for (auto v : multiCell.genome) {
			std::cout << v;
		}
		std::cout << "  " << static_cast<double>(std::accumulate(multiCell.genome.begin(), multiCell.genome.end(), 0.0))
			<< " / " << static_cast<double>(multiCell.genome.size()) << " = "
			<< static_cast<double>(std::accumulate(multiCell.genome.begin(), multiCell.genome.end(), 0.0)) / static_cast<double>(multiCell.genome.size())
			<< std::endl;
			*/

		return(static_cast<double>(std::accumulate(multiCell.genome.begin(), multiCell.genome.end(), 0.0)) / static_cast<double>(multiCell.genome.size()));
	}

	double getEvil(double alignment) {
		if (evilCutoff.second == -1) {
			return(alignment > evilCutoff.first);
		}
		else {
			if (alignment > evilCutoff.second) {
				return(1);
			}
			else if (alignment < evilCutoff.first) {
				return(0);
			}
			return(Random::P((alignment - evilCutoff.first) / (evilCutoff.second - evilCutoff.first)));
		}
	}


	Vector2d<MultiCell> world;
	std::set<std::shared_ptr<Organism>> killList;

	// pick a location form world that is not parentX, parentY...
	// place location in targetX, targetY
	// if forWorld, select a world location, if !forWorld, select an MC location
	void pickGridLoc(bool forWorld, int parentX, int parentY, int& targetX, int& targetY) {

		int sizeX = (forWorld) ? worldX : multiCellX;
		int sizeY = (forWorld) ? worldY : multiCellY;

		auto spatialModel = (forWorld) ? spatialWorldModel : spatialMultiCellModel;
		auto spatialEdgeRule = (forWorld) ? spatialWorldEdgeRule : spatialMultiCellEdgeRule;

		if (spatialWorld == 0) { // world is well mixed, pick random location that is not parent location
			do { // keep picking until we get a good one!
				targetX = Random::getIndex(sizeX);
				targetY = Random::getIndex(sizeY);
			} while (targetX == parentX && targetY == parentY);
		}
		else { // world is spatial
			if (spatialModel == "square") {
				bool pickIsValid = false;
				while (!pickIsValid) {
					do { // keep picking until we get a good one!
						targetX = parentX + Random::getInt(-1 * spatialWorldDist, spatialWorldDist);
						targetY = parentY + Random::getInt(-1 * spatialWorldDist, spatialWorldDist);
					} while (targetX == parentX && targetY == parentY);
					if (targetX < 0 || targetX >= sizeX || targetY < 0 || targetY >= sizeY) {
						if (spatialEdgeRule == "fail") {
							targetX = -1;
							targetY = -1;
							pickIsValid = true;
						}
						else if (spatialEdgeRule == "wrap") {
							targetX = loopMod(targetX, sizeX);
							targetY = loopMod(targetY, sizeY);
							pickIsValid = true;
						}
						// else, option is "search", do nothing, pickIsValid is false, so we will pick again.
					}
					else {
						pickIsValid = true;
					}
				}
			}
			else if (spatialModel == "cardinal") {
				bool pickIsValid = false;
				while (!pickIsValid) {

					int direction = Random::getIndex(4);

					if (direction == 0) {
						targetX = parentX + Random::getInt(1, spatialWorldDist);
						targetY = parentY;
					}
					else if (direction == 1) {
						targetX = parentX;
						targetY = parentY + Random::getInt(1, spatialWorldDist);
					}
					else if (direction == 2) {
						targetX = parentX - Random::getInt(1, spatialWorldDist);
						targetY = parentY;
					}
					else if (direction == 3) {
						targetX = parentX;
						targetY = parentY - Random::getInt(1, spatialWorldDist);
					}
					if (targetX < 0 || targetX >= sizeX || targetY < 0 || targetY >= sizeY) {
						if (spatialEdgeRule == "fail") {
							targetX = -1;
							targetY = -1;
							pickIsValid = true;
						}
						else if (spatialEdgeRule == "wrap") {
							targetX = loopMod(targetX, sizeX);
							targetY = loopMod(targetY, sizeY);
							pickIsValid = true;
						}
						// else, option is "search", do nothing, pickIsValid is false, so we will pick again.
					}
					else {
						pickIsValid = true;
					}
				}
			}
		}
		if (targetX == parentX && targetY == parentY) {
			pickGridLoc(forWorld, parentX, parentY, targetX, targetY);
		}
	}

	void initGenome(std::vector<bool>& genome) {
		genome.resize(alignmentGenomeSize);
		if (initEvilPercent < 0) {
			for (double i = 0; i < genome.size(); i++) {
				genome[i] = Random::getIndex(2);
			}
		}
		else {
			double s = genome.size() * initEvilPercent;
			for (double i = 0; i < genome.size(); i++) {
				if (i < s) {
					genome[i] = 1;
				}
				else {
					genome[i] = 0;
				}
			}
		}
	}

	void mutateGenome(std::vector<bool>& genome) {
		int numMutations = Random::getBinomial(alignmentGenomeSize, alignmentMutationRate);
		for (int s = 0; s < numMutations; s++) {
			genome[Random::getIndex(alignmentGenomeSize)] = Random::getIndex(2);
		}
	}
	
	// make a new cell in an MC given a world location, the parent location within the MC. groups are needed to add germ and first cell to population
	// note that this function is not used when the first cell is added to an MC since it's parent is the germ and has no MC location.
	void birthCell(std::map<std::string, std::shared_ptr<Group>>& groups, int parentWorldX, int parentWorldY, int parentCellX, int parentCellY) {		

		// parent spends resources, even if they did not produce an offspring
		world(parentWorldX, parentWorldY).cells(parentCellX, parentCellY).resource = 0;

		// pick a location for offspring
		int offspringCellX, offspringCellY;
		pickGridLoc(false, parentCellX, parentCellY, offspringCellX, offspringCellY);

		// record if target cell location is currently empty
		bool targetEmpty = world(parentWorldX, parentWorldY).cells(offspringCellX, offspringCellY).empty;

		// if target location is empty or parent is evil or overwriteSelfOnFail an offspring will be made
		if (targetEmpty || world(parentWorldX, parentWorldY).cells(parentCellX, parentCellY).evil || overwriteSelfOnFail) {
			
			// get a link to the parent org
			auto parentOrg = world(parentWorldX, parentWorldY).cells(parentCellX, parentCellY).org;

			// if target is not empty and the parent is evil
			if (!targetEmpty && world(parentWorldX, parentWorldY).cells(parentCellX, parentCellY).evil) {
				// kill what's in the target cell now...
				killList.insert(world(parentWorldX, parentWorldY).cells(offspringCellX, offspringCellY).org);
			}
			
			else if (!targetEmpty) { // ... and a good parent
				// parent will oversrite self set target location to parent location
				offspringCellX = parentCellX;
				offspringCellY = parentCellY;
				killList.insert(parentOrg); // kill what's here now... i.e. the parent
			}

			// make a new org from parent in the new location
			world(parentWorldX, parentWorldY).cells(offspringCellX, offspringCellY).org = parentOrg->makeMutatedOffspringFrom(parentOrg);
			world(parentWorldX, parentWorldY).cells(offspringCellX, offspringCellY).genome = world(parentWorldX, parentWorldY).cells(parentCellX, parentCellY).genome;
			
			mutateGenome(world(parentWorldX, parentWorldY).cells(offspringCellX, offspringCellY).genome);

			// add new org to population
			groups["root::"]->population.push_back(world(parentWorldX, parentWorldY).cells(offspringCellX, offspringCellY).org);

			// reset resource on offspring
			world(parentWorldX, parentWorldY).cells(offspringCellX, offspringCellY).resource = 0;

			// set flags on cell
			world(parentWorldX, parentWorldY).cells(offspringCellX, offspringCellY).alignment = getAlignment(world(parentWorldX, parentWorldY).cells(offspringCellX, offspringCellY));
			world(parentWorldX, parentWorldY).cells(offspringCellX, offspringCellY).evil = getEvil(world(parentWorldX, parentWorldY).cells(offspringCellX, offspringCellY).alignment);
			world(parentWorldX, parentWorldY).cells(offspringCellX, offspringCellY).murder = world(parentWorldX, parentWorldY).cells(parentCellX,parentCellY).evil && !targetEmpty;
			world(parentWorldX, parentWorldY).cells(offspringCellX, offspringCellY).first = false;
			world(parentWorldX, parentWorldY).cells(offspringCellX, offspringCellY).empty = false;

			// add tracking data to MC
			world(parentWorldX, parentWorldY).goodCellBirthCounts += !world(parentWorldX, parentWorldY).cells(offspringCellX, offspringCellY).evil;
			world(parentWorldX, parentWorldY).evilCellBirthCounts += world(parentWorldX, parentWorldY).cells(offspringCellX, offspringCellY).evil;
			world(parentWorldX, parentWorldY).murderCellBirthCounts += world(parentWorldX, parentWorldY).cells(offspringCellX, offspringCellY).murder;
			world(parentWorldX, parentWorldY).cellAlignmentSum += world(parentWorldX, parentWorldY).cells(offspringCellX, offspringCellY).alignment;
		}
	}

	// make a new mutliCell given a parent org, a genome, a world location, groups are needed to add germ and first cell to population
	// note, the first cell is created directly here as opposed to using birthCell because the parent (i.e. the germ) does not have a
	// location in the MC
	void birthMultiCell(std::map<std::string, std::shared_ptr<Group>>& groups, std::shared_ptr<Organism> parent,
		std::vector<bool> genome_, int newWorldX, int newWorldY) {
		
		world(newWorldX, newWorldY).cells.reset(multiCellX, multiCellY);
		
		world(newWorldX, newWorldY).germ = parent;
		groups["root::"]->population.push_back(world(newWorldX, newWorldY).germ);
		world(newWorldX, newWorldY).genome = genome_;

		if (!useWorldGenome) {
			mutateGenome(world(newWorldX, newWorldY).genome);
		}
		world(newWorldX, newWorldY).alignment = getAlignment(world(newWorldX, newWorldY));
		world(newWorldX, newWorldY).evil = getEvil(world(newWorldX, newWorldY).alignment);
		
		world(newWorldX, newWorldY).empty = false;
		world(newWorldX, newWorldY).full = false;
		world(newWorldX, newWorldY).canProduce = false;

		// pick location for first cell
		int offspringCellX = Random::getIndex(multiCellX);
		int offspringCellY = Random::getIndex(multiCellY);

		// make first cell from germ
		world(newWorldX, newWorldY).cells(offspringCellX, offspringCellY).org = world(newWorldX, newWorldY).germ->makeMutatedOffspringFrom(world(newWorldX, newWorldY).germ);
		world(newWorldX, newWorldY).cells(offspringCellX, offspringCellY).genome = world(newWorldX, newWorldY).genome;
		
		mutateGenome(world(newWorldX, newWorldY).cells(offspringCellX, offspringCellY).genome);

		// add new org to population
		groups["root::"]->population.push_back(world(newWorldX, newWorldY).cells(offspringCellX, offspringCellY).org);

		world(newWorldX, newWorldY).cells(offspringCellX, offspringCellY).alignment = getAlignment(world(newWorldX, newWorldY).cells(offspringCellX, offspringCellY));
		world(newWorldX, newWorldY).cells(offspringCellX, offspringCellY).evil = getEvil(world(newWorldX, newWorldY).cells(offspringCellX, offspringCellY).alignment);
		world(newWorldX, newWorldY).cells(offspringCellX, offspringCellY).murder = false;
		world(newWorldX, newWorldY).cells(offspringCellX, offspringCellY).first = true;
		world(newWorldX, newWorldY).cells(offspringCellX, offspringCellY).empty = false;

		// add tracking data to MC, this is the first cell, so use '=' not '+='
		world(newWorldX, newWorldY).goodCellBirthCounts = !world(newWorldX, newWorldY).cells(offspringCellX, offspringCellY).evil;
		world(newWorldX, newWorldY).evilCellBirthCounts = world(newWorldX, newWorldY).cells(offspringCellX, offspringCellY).evil;
		world(newWorldX, newWorldY).murderCellBirthCounts = world(newWorldX, newWorldY).cells(offspringCellX, offspringCellY).murder;
		world(newWorldX, newWorldY).cellAlignmentSum = world(newWorldX, newWorldY).cells(offspringCellX, offspringCellY).alignment;
	}

	void displayMultiCell(MultiCell& MC) {
		std::cout << "types in MultiCell:" << std::endl;
		for (int y = 0; y < MC.cells.y(); y++) {
			for (int x = 0; x < MC.cells.x(); x++) {
				if (MC.cells(x, y).empty) {
					std::cout << "    ";
				}
				else if (MC.cells(x, y).first) {
					std::cout << "f   ";
				}
				else if (MC.cells(x, y).murder) {
					std::cout << "m   ";
				}
				else if (MC.cells(x, y).evil) {
					std::cout << "e   ";
				}
				else{
					std::cout << "g   ";
				}
			}
			std::cout << "  |" << std::endl;
		}
		std::cout << std::endl;
	}

	void displayMultiCellResource(MultiCell& MC) {
		double MC_total = 0;
		std::cout << "resource values in MultiCell:" << std::endl;
		for (int y = 0; y < MC.cells.y(); y++) {
			for (int x = 0; x < MC.cells.x(); x++) {
				if (MC.cells(x, y).empty) {
					std::cout << "-   ";
				}
				else {
					std::cout << MC.cells(x, y).resource << " ";
					MC_total += MC.cells(x, y).resource;
				}
			}
			std::cout << " |" << std::endl;
		}
		std::cout << "total for MC: " << MC_total << std::endl;
	}

	void saveWorldImage();

	double vectAve(std::vector<double> vect) {
		return std::accumulate(vect.begin(), vect.end(), 0.0) / static_cast<double>(vect.size());
	}

	MultiCellWorld(std::shared_ptr<ParametersTable> PT_ = nullptr);
	virtual ~MultiCellWorld() = default;

	void evaluate(std::map<std::string, std::shared_ptr<Group>>& groups,
		int analyze, int visualize, int debug);

	virtual std::unordered_map<std::string, std::unordered_set<std::string>>
		requiredGroups() override;
};

