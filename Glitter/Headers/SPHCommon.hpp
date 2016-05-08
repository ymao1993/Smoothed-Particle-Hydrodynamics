#ifndef SPHCOMMOM_H
#define SPHCOMMOM_H

#include <vector>
#include "SPHMath.hpp"

#define KERNEL_H 		(1.1f)
#define K 				(1e3)
#define RESET_DENSITY   (1.2f) 
#define G 				(9.8f)
#define MIU				(0.1f)
#define BOUNDARYDAMPING (0.9f)
#define DAMPING 		(0.5f)
#define TENSION_COEF	(1.f)
	
namespace SPHSim
{
	/*
	 * SPH Particle struct
	 **/

	struct SPHParticle
	{
		float mass;
		float pressure;
		float density;
		vec3 velocity; 
		vec3 position;
		vec3 f; //total forces
		std::vector<int> neighbors;
		SPHParticle(float mass, vec3 position):mass(mass),position(position){}
	};

	struct UnitCube
	{
		bool vertexFlag[8];
		double vertexValue[8];

		SPHSim::vec3 cubePos;

		inline SPHSim::vec3 addOffset(double dx, double dy, double dz) 
		{
			return cubePos + SPHSim::vec3(dx, dy, dz);
		}
	};

	struct Grid
	{
		bool flag;
		double value;
		SPHSim::vec3 position;
		std::vector<int> particleIndices;
	};

}
#endif