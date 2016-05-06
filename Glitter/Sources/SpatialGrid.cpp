#include <cstdlib>
#include <iostream>
#include <cmath>

#include "SpatialGrid.hpp"
#include "SPHMath.hpp"


static inline bool compareDistance(const SPHSim::SPHParticle &p1, const SPHSim::SPHParticle &p2, double distance)
{
	SPHSim::vec3 v = p1.position - p2.position;
	return v.norm() < distance;
}


// setup the boundary by extending radius on the both end automatically
SpatialGrid::SpatialGrid(double bx, double by, double bz, double lx, double ly, double lz, double h, double r)
: base_x(bx - 2*r), base_y(by - 2*r), base_z(bz - 2*r), 
  len_x(lx + 4*r), len_y(ly + 4*r), len_z(lz + 4*r), grid_h(h), radius(r),
  dim_x(ceil(len_x / grid_h)), dim_y(ceil(len_y / grid_h)), dim_z(ceil(len_z / grid_h)),
  dim_xy(dim_x * dim_y),
  arraySize(dim_x * dim_y * dim_z),
  grid(new SPHSim::Grid[arraySize])
{
}

void SpatialGrid::clearGrid()
{
	for (int i = 0; i < arraySize; i++)
	{
		grid[i].particleIndices.clear();
	}
}

void SpatialGrid::addParticle(const SPHSim::SPHParticle &p_i, int i)
{
	const int gridIndex = particleToIndex(p_i);

	if (gridIndex >= arraySize || gridIndex < 0) return;

	grid[gridIndex].particleIndices.push_back(i);
}

void SpatialGrid::findNeighbors(SPHSim::SPHParticle &p_i,
						 	   const std::vector<SPHSim::SPHParticle> &particles) 
{
	p_i.neighbors.clear();

	const int thisGridIndex = particleToIndex(p_i);
	// check out all neighboring grid cells
	for (int dz = -1; dz <= 1; dz++) 
	{
		for (int dy = -1; dy <= 1; dy++) 
		{
			for (int dx = -1; dx <= 1; dx++) 	
			{
				const int neighborGridIndex = thisGridIndex + dx + dy * dim_x + dz * dim_xy;
				if (neighborGridIndex >= arraySize || neighborGridIndex < 0) continue;
				if (grid[neighborGridIndex].particleIndices.size() == 0) continue;

				for (int j : grid[neighborGridIndex].particleIndices) 
				{
					if (compareDistance(p_i, particles[j], KERNEL_H)) 
					{
						p_i.neighbors.push_back(j);
					}
				}

			}
		}
	}
}

SpatialGrid::~SpatialGrid()
{
	printf("deleting grid\n");
	delete [] grid;
}
