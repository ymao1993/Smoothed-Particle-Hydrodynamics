#ifndef SPATIALGRID_H
#define SPATIALGRID_H

#include <vector>
#include "SPHCommon.hpp"

/* 
 * SpatialGrid discretizes the space into equally-spaced grids
 * each grid stores the indeices of the particle in order to
 * speed up the neighborhood particle finding. 
 *
 * Che-Yuan Liang Apr. 2016 
 */


class SpatialGrid
{
protected:
	const double base_x;
	const double base_y;
	const double base_z;

	const double len_x;
	const double len_y;
	const double len_z;

	const double grid_h;
	const double radius;

	const int dim_x;
	const int dim_y;
	const int dim_z;
	const int dim_xy;

	const unsigned int arraySize;
	SPHSim::Grid* grid;

	inline int particleToIndex(const SPHSim::SPHParticle &p)
	{
		const double px = p.position.x;
		const double py = p.position.y;
		const double pz = p.position.z;

		return floor((px-base_x)/grid_h) + floor((py-base_y)/grid_h) * dim_x + floor((pz-base_z)/grid_h) * dim_x * dim_y;
	}

public:
	SpatialGrid(double base_x, double base_y, double base_z,
				double length_x, double length_y, double length_z, double h, double r);
	~SpatialGrid();
	
	void addParticle(const SPHSim::SPHParticle &p, int i);
	void clearGrid();
	void findNeighbors(SPHSim::SPHParticle &p_i, const std::vector<SPHSim::SPHParticle> &particles);
	inline SPHSim::vec3 indexToPosition(int i)
	{
		return SPHSim::vec3(base_x + (i % dim_x) * grid_h, base_y + ((i % dim_xy) / dim_x) * grid_h, base_z + (i / dim_xy) * grid_h);
	}

};

#endif