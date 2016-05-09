#ifndef SPHSIMULATOR_H
#define SPHSIMULATOR_H

#include "SPHCommon.hpp"

/**
 * Fluid simulator based on Smoothed Particle Hydrodynamics(SPH)
 *
 * References:
 *
 * Müller M, Charypar D, Gross M. Particle-based fluid simulation for i
 * nteractive applications[C]//Proceedings of the 2003 ACM SIGGRAPH/Eurographics 
 * symposium on Computer animation. Eurographics Association, 2003: 154-159.
 * (http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.121.844&rep=rep1&type=pdf) 
 *
 * Lecture notes from CMU 15-467 (by Stelian Coros):
 * (http://www.cs.cmu.edu/~scoros/cs15467-s16/lectures/10-fluids1.pdf)
 * (http://www.cs.cmu.edu/~scoros/cs15467-s16/lectures/11-fluids2.pdf)
 *
 * Author: Yu Mao
 * Contact: ymao1@andrew.cmu.edu
 */

namespace SPHSim
{

	struct BoxDef
	{
		double spanX;
		double spanY;
		double spanZ;
	};

	/**
	 * SPH configuration paramters
	 */
	struct SPHConfig
	{
		double kernel_h; //compact support of the smooth function
		double k; //(density - density0) * k = pressure
		double density0; 
		double g; //gravity constant
		double miu;
		double boundarydamping; //0~1
		double damping;
		double tensionCoe; //tension coefficient

		BoxDef box;
	};

	class SPHSimulator
	{
	public:
		/**
		 * constructor
		 */
		SPHSimulator(SPHConfig conf);

		/**
		 * destructor
		 */
		~SPHSimulator();

		/**
		 * setup the simulator at the first time
		 */
		void setup();

		/**
		 * reset the all the particle states
		 */
		void reset();

		/**
		 * update
		 */
		void update(float delta);

		/**
		 * get data for rendering
		 */
		void getData(float** position, int& vertexNum);

		/**
		 * setGradient direction
		 */
		void setGravityDirection(float x, float y, float z)
		{
			gravityDir.x = x;
			gravityDir.y = y;
			gravityDir.z = z;
			gravityDir.normalize();
		}

	private:
		SPHConfig conf;

		std::vector<SPHParticle> particles;
		std::vector<int> indices;
		vec3 gravityDir;

		/*  <Device memory>  */
		// position
		double* device_position;
		// velocity
		double* device_velocity;
		// force
		double* device_force;
		// mass
		float* device_mass;
		// density
		float* device_density;
		// presure
		float* device_pressure;
		// neighbor
		int* device_neighbors;
		// bucket
		Bucket* device_bucket;
		// index
		int* device_index;
		// bucket index
		int* device_bucket_index;
		// bucket_mask
		int* device_bucket_mask;
		// exclusive scan
		int* device_scan;
		// compressed bucket
		int* device_compressed_bucket;

		

		/*  <Host memory>  */
		// position
		double* host_position;
		// velocity
		double* host_velocity;	
		// force
		double* host_force;
		// mass
		float* host_mass;
		// neighbors
		int* host_neighbors;

		// numOfParticle
		int N;
		int dim_x;
		int dim_y;
		int dim_z; 
		int arraySize;



	};
}


#endif
