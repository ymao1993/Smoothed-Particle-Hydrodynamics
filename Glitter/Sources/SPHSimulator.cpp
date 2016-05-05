#include "SPHSimulator.hpp"
#include "SPHCommon.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>


namespace SPHSim
{
	SPHSimulator::SPHSimulator(SPHConfig config)
	:
	conf(config),
	particleGrid(-config.box.spanX/2, -config.box.spanY/2, -config.box.spanZ/2, 
				  config.box.spanX, config.box.spanY, config.box.spanZ, config.kernel_h, 
				  config.kernel_h)
	{
		setGravityDirection(0, -1, 0);
	}

	SPHSimulator::~SPHSimulator()
	{

	}

	void SPHSimulator::setup()
	{
		double mass = 1.f;
		double start_x = -conf.box.spanX/3;
		double start_y = -conf.box.spanY/3;
		double start_z = -conf.box.spanZ/3;
		double end_x   =  conf.box.spanX/3;
		double end_y   =  conf.box.spanY/3;
		double end_z   =  conf.box.spanZ/3;

		double f = 0.9;

		for (float x = start_x; x < end_x; x += conf.kernel_h * f)
			for (float y = start_y; y < end_y; y += conf.kernel_h * f)
				for (float z = start_z; z < end_z; z += conf.kernel_h * f)
					particles.push_back(SPHParticle(mass, vec3(x,y,z)));

		return;
	}

	void SPHSimulator::reset()
	{
		return;
	}

	inline float SPHSimulator::wPoly6(vec3 rvec, float h)
	{
		float r = rvec.norm();
		float result = 0;
		if(r<=h)
		{
			result = (315.f/(64.f * PI * pow(h,9))) * pow(pow(h,2) - pow(r,2),3);
		}
		return result;
	}

	inline vec3 SPHSimulator::wPoly6Gradient(vec3 rvec, float h)
	{
		float r = rvec.norm();
		vec3 result;
		if(r<=h)
		{
			result = -945.0f / (32.0f * PI * std::pow(h,9)) * std::pow(std::pow(h,2) - std::pow(r,2),2) * rvec;
		}
		return result;
	}

	inline float SPHSimulator::wPoly6Laplacian(vec3 rvec, float h)
	{
		float r = rvec.norm();
		float result = 0;
		if(r<=h)
		{
			result = 945.0f / (32.0f * PI * std::pow(h,9)) * (std::pow(h,2) - std::pow(r,2)) * (7.0f * std::pow(r,2) - 3.0f * std::pow(h,2));
		}
		return result;
	}

	inline vec3 SPHSimulator::wSpikyGradient(vec3 rvec, float h)
	{
		float r = rvec.norm();
		vec3 result;
		if(r<=h && r > EPSILON)
		{
			result = ((- 45.f * std::pow(h-r,2)) / (r * PI * std::pow(h,6))) * rvec;
		}
		return result;
	}

	inline float SPHSimulator::wViscosityLaplacian(vec3 rvec, float h)
	{
		float r = rvec.norm();
		float result = 0;
		if(r<=h && r>EPSILON)
		{
			result = 45.f / (PI * std::pow(h,6)) * (h - r);
		}
		return result;
	}

	void SPHSimulator::update(float delta)
	{
		//TODO: for all particles, find neighbourhoods
		// Spatial Grid Neibor Finding
		particleGrid.clearGrid();

		for (int i = 0; i < particles.size(); i++)
		{
			particleGrid.addParticle(particles[i], i);
		}
        for (int i = 0; i < particles.size(); i++)
		{
			particleGrid.findNeighbors(particles[i], particles);
		}

        for(int i=0; i<particles.size(); i++)
		{
			float density = 0;;
			for(int j: particles[i].neighbors)
			{
				vec3 rvec = particles[i].position - particles[j].position;
				float w = wPoly6(rvec, conf.kernel_h);
				density += particles[j].mass * w;
			}
			float pressure = conf.k * (density - conf.density0); //max(xxx,0)??
			particles[i].density = density;
			particles[i].pressure = pressure;
		}

		//TODO: for all particles, compute forces
        for(int i=0; i<particles.size(); i++)
		{
			vec3 fpressure, fviscosity, ftension, fgravity;
			float color_laplacian = 0;
			vec3 color_gradient;
			for(auto j: particles[i].neighbors)
			{
				vec3 rvec = particles[i].position - particles[j].position;

				//pressure force
				vec3 wspikyGradient = wSpikyGradient(rvec, conf.kernel_h);
				fpressure  -= (particles[j].mass / particles[j].density) * ((particles[i].pressure + particles[j].pressure) / 2.f) * wspikyGradient;

				//viscosity force
				float wviscosityLaplacian = wViscosityLaplacian(rvec, conf.kernel_h);
				fviscosity += conf.miu * (particles[j].mass / particles[j].density) * (particles[j].velocity - particles[i].velocity) * wviscosityLaplacian;

				//compute gradient of color field
				vec3 wpoly6gradient = wPoly6Gradient(rvec, conf.kernel_h);
				color_gradient += (particles[j].mass / particles[j].density) * wpoly6gradient;

				//compute laplacian of color filed
				float wpoly6Laplacian = wPoly6Laplacian(rvec, conf.kernel_h);
				color_laplacian += (particles[j].mass / particles[j].density) * wpoly6Laplacian;
			}
			fgravity = gravityDir * conf.g * particles[i].density;

			if(color_gradient.norm() > EPSILON)
			{
				ftension =  -conf.tensionCoe * color_laplacian * color_gradient.unit();
			}

			vec3 f = fpressure + fviscosity + fgravity + ftension;
			particles[i].f = f;
		}

		//TODO: for all particles, update velocity and position
		//using symplectic euler scheme
        for(int i=0; i<particles.size(); i++)
		{
			vec3 accel = (particles[i].f / particles[i].density);
			accel -= conf.damping * particles[i].velocity; //the faster the velocity, the quicker the damping
			particles[i].velocity += delta * accel;
			particles[i].position += delta * particles[i].velocity;

			//handle boundary condition
			vec3 min(-conf.box.spanX/2,-conf.box.spanY/2,-conf.box.spanZ/2);
			vec3 max( conf.box.spanX/2, conf.box.spanY/2, conf.box.spanZ/2);
			box_clamp_and_reflect(particles[i].position, particles[i].velocity ,min, max, conf.boundarydamping);
		}
	}

	void SPHSimulator::getData(float** position, int& vertexNum)
	{
		//free the space if needed
		if(*position != NULL)
			free(*position);

		//allocate new memory
		*position = (float*)malloc(sizeof(float) * particles.size() * 3);
		vertexNum = particles.size();

		for(int i=0; i<vertexNum; i++)
		{
			(*position)[3*i+0] = particles[i].position.x;
			(*position)[3*i+1] = particles[i].position.y;
			(*position)[3*i+2] = particles[i].position.z;
		}
		return;
	}

}
