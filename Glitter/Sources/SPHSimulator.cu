#include "SPHSimulator.cuh"
#include "SPHCommon.hpp"
#include <thrust/sort.h>
#include <thrust/device_ptr.h>

#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
//#include <omp.h>

#define SPH_TEST
#ifdef SPH_TEST
#define WIDTH		30
#define HEIGHT		30
#define DEPTH		30
#define SCALE		1.0f
#define STARTX		1.f
#define STARTY		1.f
#define STARTZ 		15.0f
#endif

#define px(id) (d_position[3*id])
#define py(id) (d_position[3*id+1])
#define pz(id) (d_position[3*id+2])
#define mass(id) (d_mass[id])
#define density(id) (d_density[id])
#define pressure(id) (d_pressure[id])

struct cuda_vec3
{
	double x;
	double y;
	double z;
};

struct cuda_fvec3
{
	float x;
	float y;
	float z;
};

// Device code
__device__
 float wPoly6(cuda_vec3 rvec, float h)
{
	float r = sqrt(rvec.x*rvec.x + rvec.y*rvec.y + rvec.z*rvec.z);
	float result = 0;
	if(r<=h)
	{
		result = (315.f/(64.f * SPHSim::PI * pow(h,9.f))) * pow(pow(h,2.f) - pow(r,2.f),3.f);
	}
	return result;
}

__device__ 
cuda_vec3 wPoly6Gradient(cuda_vec3 rvec, float h)
{
	float r = sqrt(rvec.x*rvec.x + rvec.y*rvec.y + rvec.z*rvec.z);
	cuda_vec3 result = {0, 0, 0};
	float c = -945.0f / (32.0f * SPHSim::PI * pow(h,9.f)) * pow(pow(h,2.f) - pow(r,2.f),2.f);
	if(r<=h)
	{
		result.x = c * rvec.x;
		result.y = c * rvec.y;
		result.z = c * rvec.z;
	}
	return result;
}

__device__
float wPoly6Laplacian(cuda_vec3 rvec, float h)
{
	float r = sqrt(rvec.x*rvec.x + rvec.y*rvec.y + rvec.z*rvec.z);
	float result = 0;
	if(r<=h)
	{
		result = 945.0f / (32.0f * SPHSim::PI * pow(h,9.f)) * (pow(h,2.f) - pow(r,2.f)) * (7.0f * pow(r,2.f) - 3.0f * pow(h,2.f));
	}
	return result;
}

__device__
cuda_vec3 wSpikyGradient(cuda_vec3 rvec, float h)
{
	float r = sqrt(rvec.x*rvec.x + rvec.y*rvec.y + rvec.z*rvec.z);
	cuda_vec3 result = {0, 0, 0};
	float c = ((- 45.f * pow(h-r,2.f)) / (r * SPHSim::PI * pow(h,6.f)));
	if((r<=h) && (r > SPHSim::EPSILON))
	{
		result.x = c * rvec.x;
		result.y = c * rvec.y;
		result.z = c * rvec.z;
	}
	return result;
}

__device__
float wViscosityLaplacian(cuda_vec3 rvec, float h)
{
	float r = sqrt(rvec.x*rvec.x + rvec.y*rvec.y + rvec.z*rvec.z);
	float result = 0;
	if(r<=h && r>SPHSim::EPSILON)
	{
		result = 45.f / (SPHSim::PI * pow(h,6.f)) * (h - r);
	}
	return result;
}

__device__
int particleToIndex(cuda_vec3 v, double kernel_h, int dim_x, int dim_y)
{
	return (int)floor(v.x/kernel_h) + (int)floor(v.y/kernel_h) * dim_x + (int)floor(v.z/kernel_h) * dim_x * dim_y;
}

__global__
void kernelMakeBlock(int* d_compressed_block, int* d_mask, int* d_scan, int arraySize)
{
	int idx = threadIdx.x + ((gridDim.x * blockIdx.y) + blockIdx.x) * blockDim.x;
	if (idx >= arraySize) return;
	if (d_mask[idx] == 1)
	{
		d_compressed_block[d_scan[idx]] = idx;
	}
}

__global__
void kernelSplitBlock(SPHSim::Bucket* d_bucket, int* mask, int arraySize)
{
	int idx = threadIdx.x + ((gridDim.x * blockIdx.y) + blockIdx.x) * blockDim.x;
	// printf("threadx: %d thready: %d threadz: %d\n", threadIdx.x, threadIdx.y, threadIdx.z);
	// printf("blockidx: %d, blockidxy: %d, blockIdxz: %d\n", blockIdx.x, blockIdx.y, blockIdx.z);
	if (idx >= arraySize) return;

	if (d_bucket[idx].size == 0)
	{
		mask[idx] = 0;
	}
	else
	{
		mask[idx] = 1;
	}	
}
__global__
void kernelComputeBucket(SPHSim::Bucket* d_bucket, int* d_bucket_index, int N, int arraySize)
{
	 int idx = threadIdx.x + ((gridDim.x * blockIdx.y) + blockIdx.x) * blockDim.x;
    // guard
    if ((idx >= 0) && (idx < N) && (d_bucket_index[idx] < arraySize) && (d_bucket_index[idx] >= 0))
    {
        atomicAdd(&d_bucket[d_bucket_index[idx]].size, 1);
        atomicCAS(&d_bucket[d_bucket_index[idx]].startIndex, 0, idx);
        atomicMin(&d_bucket[d_bucket_index[idx]].startIndex, idx);
    }
}

__global__
void kernelComputeIndex(double* d_position, int* d_index, int* d_bucket_index, int N,
									 float kernel_h, int dim_x, int dim_y)
{
	int id = threadIdx.x + ((gridDim.x * blockIdx.y) + blockIdx.x)*blockDim.x;
	if (id >= N) return;

	d_index[id] = id;
	cuda_vec3 v = {px(id), py(id), pz(id)};
	int index = particleToIndex(v, kernel_h, dim_x, dim_y);
	d_bucket_index[id] = index;
	// printf("particle: %d, position: %f, %f, %f index: %d kernel: %f, dim_x: %d, dim_y: %d N: %d\n", id, v.x, v.y, v.z, index, kernel_h, dim_x, dim_y, N);	
}

__global__
void kernelSortParticle(double* d_position, double* d_sorted_position, double* d_velocity, double* d_sorted_velocity, int* d_index, int N)
{
	int id = threadIdx.x + ((gridDim.x * blockIdx.y) + blockIdx.x)*blockDim.x;

	if (id >= N) return;
	d_sorted_position[3*id] = d_position[3*d_index[id]];
	d_sorted_position[3*id+1] = d_position[3*d_index[id]+1];
	d_sorted_position[3*id+2] = d_position[3*d_index[id]+2];

	d_sorted_velocity[3*id] = d_velocity[3*d_index[id]];
	d_sorted_velocity[3*id+1] = d_velocity[3*d_index[id]+1];
	d_sorted_velocity[3*id+2] = d_velocity[3*d_index[id]+2];

}

__global__
void kernelUpdateAndHandle(double* d_position, double* d_velocity, double* d_force, float* d_density,
									 int N, SPHSim::SPHConfig conf, double delta, double spanX, double spanY, double spanZ)
{
	int id = threadIdx.x + ((gridDim.x * blockIdx.y) + blockIdx.x)*blockDim.x;
	if (id >= N) return;
	
	// update position & velocity
	cuda_vec3 accel = {d_force[3*id] / d_density[id], d_force[3*id+1] / d_density[id], d_force[3*id+2] / d_density[id]};
	accel.x -= conf.damping * d_velocity[3*id];
	accel.y -= conf.damping * d_velocity[3*id+1];
	accel.z -= conf.damping * d_velocity[3*id+2];
	d_velocity[3*id] += delta * accel.x;
	d_velocity[3*id+1] += delta * accel.y;
	d_velocity[3*id+2] += delta * accel.z;
	d_position[3*id] += delta * d_velocity[3*id];
	d_position[3*id+1] += delta * d_velocity[3*id+1];
	d_position[3*id+2] += delta * d_velocity[3*id+2];

	// collision handle

	//handle boundary condition
	cuda_vec3 min = {0, 0, 0};
	cuda_vec3 max = {spanX, spanY, spanZ};
	if(d_position[3*id] < min.x) { d_position[3*id] = min.x; d_velocity[3*id] *= -(1.f-conf.damping); }
	if(d_position[3*id] > max.x) { d_position[3*id] = max.x; d_velocity[3*id] *= -(1.f-conf.damping); }
	if(d_position[3*id+1] < min.y) { d_position[3*id+1] = min.y; d_velocity[3*id+1] *= -(1.f-conf.damping); }
	if(d_position[3*id+1] > max.y) { d_position[3*id+1] = max.y; d_velocity[3*id+1] *= -(1.f-conf.damping); }
	if(d_position[3*id+2] < min.z) { d_position[3*id+2] = min.z; d_velocity[3*id+2] *= -(1.f-conf.damping); }
	if(d_position[3*id+2] > max.z) { d_position[3*id+2] = max.z; d_velocity[3*id+2] *= -(1.f-conf.damping); }

	
}

// for all particles, compute density and pressure
__global__ 
void kernelComputeDensity(double* d_position, SPHSim::Bucket* device_bucket, int* d_compressed_bucket, float* d_density, float* d_mass,
								float* d_pressure, int N, int numNonZeroBucket, SPHSim::SPHConfig conf, int dim_x, int dim_y, int arraySize)
{

	if (blockIdx.x >= numNonZeroBucket) return;
	if (threadIdx.x >= 27) return;

	int bucketId = d_compressed_bucket[blockIdx.x];
	int neighborBucketId;
	int count = threadIdx.x;

	for (int k = -1; k <= 1; k++)
	{
		for (int j = -1; j <= 1; j++)
		{
			for (int i = -1; i <= 1; i++)
			{
				if (--count < 0)
				{
					// printf("threadIdx: %d got %d, %d, %d\n", threadIdx.x, i, j, k);
					neighborBucketId = bucketId+i+j*dim_x+k*dim_x*dim_y;
					goto run;	
				}
			}
		}
	}

	run:
	for (int m = 0; m < device_bucket[bucketId].size; m++)
	{
		int particleId = device_bucket[bucketId].startIndex + m;
		// printf("particle: %d nieghborId: %d\n", particleId, neighborBucketId);
		__shared__ 
		float density;
		density = 0;
		__shared__
		double mx, my, mz, nx, ny, nz;
		mx = px(particleId);
		my = py(particleId);
		mz = pz(particleId);

		if (neighborBucketId >= arraySize || neighborBucketId < 0) continue; 
		for (int l = 0; l < device_bucket[neighborBucketId].size; l++)
		{
			int neighborId = device_bucket[neighborBucketId].startIndex + l;
			nx = px(neighborId);
			ny = py(neighborId);
			nz = pz(neighborId);
			// printf("neighbor particle %d: %f, %f, %f\n", l, px(neighborId), py(neighborId), pz(neighborId));
			cuda_vec3 v = {mx - nx, my - ny, mz - nz};
			
			float w = wPoly6(v, conf.kernel_h);
			atomicAdd(&density, mass(neighborId) * wPoly6(v, conf.kernel_h));
		}
		__syncthreads();
		float pressure = conf.k * (density - conf.density0); //max(xxx,0)??
		d_density[particleId] = density;
		d_pressure[particleId] = pressure;
		// printf("particle %d presure: %f density: %f\n", particleId, pressure, density);
	}
}
// compute force
__global__
void kernelComputeForce(double* d_position, double* d_velocity, SPHSim::Bucket* device_bucket, int* d_compressed_bucket, double* d_force, 
		float* d_density, float* d_mass, float* d_pressure, int N, int numNonZeroBucket, SPHSim::SPHConfig conf, int dim_x, int dim_y, int arraySize, 
		float gx, float gy, float gz)
{
	if (blockIdx.x >= numNonZeroBucket) return;
	if (threadIdx.x >= 27) return;

	int bucketId = d_compressed_bucket[blockIdx.x];
	int neighborBucketId;
	int count = threadIdx.x;

	for (int k = -1; k <= 1; k++)
	{
		for (int j = -1; j <= 1; j++)
		{
			for (int i = -1; i <= 1; i++)
			{
				if (--count < 0)
				{
					// printf("threadIdx: %d got %d, %d, %d\n", threadIdx.x, i, j, k);
					neighborBucketId = bucketId+i+j*dim_x+k*dim_x*dim_y;
					goto run;	
				}
			}
		}
	}

	run:
	for (int m = 0; m < device_bucket[bucketId].size; m++)
	{
		int particleId = device_bucket[bucketId].startIndex + m;

		__shared__ cuda_fvec3 fpressure, fviscosity, ftension, fgravity; 
		fpressure = {0, 0, 0};
		fviscosity = {0, 0, 0};
		ftension = {0, 0, 0};
		fgravity = {0, 0, 0};

		__shared__ float color_laplacian;
		color_laplacian = 0;
		__shared__ cuda_fvec3 color_gradient;
		color_gradient = {0, 0, 0};
		if (neighborBucketId >= arraySize || neighborBucketId < 0) continue; 
		for (int l = 0; l < device_bucket[neighborBucketId].size; l++)
		{
			int neighborId = device_bucket[neighborBucketId].startIndex + l;
			cuda_vec3 rvec = {px(particleId) - px(neighborId), py(particleId) - py(neighborId), pz(particleId) - pz(neighborId)};

			// pressure force
			cuda_vec3 wspikyGradient = wSpikyGradient(rvec, conf.kernel_h);
			float c = (mass(neighborId) / density(neighborId)) * ((pressure(particleId) + pressure(neighborId)) / 2.f);
			atomicAdd(&fpressure.x, -c * wspikyGradient.x);
			atomicAdd(&fpressure.y, -c * wspikyGradient.y);
			atomicAdd(&fpressure.z, -c * wspikyGradient.z);

			// viscosity force
			double wviscosityLaplacian = wViscosityLaplacian(rvec, conf.kernel_h);
			c = conf.miu * (mass(neighborId) / density(neighborId)) * wviscosityLaplacian;
			atomicAdd(&fviscosity.x, c * (d_velocity[3*neighborId] - d_velocity[3*particleId]));
			atomicAdd(&fviscosity.y, c * (d_velocity[3*neighborId+1] - d_velocity[3*particleId+1]));
			atomicAdd(&fviscosity.z, c * (d_velocity[3*neighborId+2] - d_velocity[3*particleId+2]));
			// printf("%f\n", c );

			// compute gradient of color field
			cuda_vec3 wpoly6gradient = wPoly6Gradient(rvec, conf.kernel_h);
			c = (mass(neighborId) / density(neighborId));
			atomicAdd(&color_gradient.x, c * wpoly6gradient.x);
			atomicAdd(&color_gradient.y, c * wpoly6gradient.y);
			atomicAdd(&color_gradient.z, c * wpoly6gradient.z);

			// compute laplacian of color filed
			double wpoly6Laplacian = wPoly6Laplacian(rvec, conf.kernel_h);
			atomicAdd(&color_laplacian, (mass(neighborId) / density(neighborId)) * wpoly6Laplacian);
		}
		__syncthreads();

		// gravity force
		float tmp = (float)conf.g * density(particleId);
		fgravity = {gx*tmp, gy*tmp, gz*tmp};

		float norm = sqrt(color_gradient.x*color_gradient.x + color_gradient.y*color_gradient.y + color_gradient.z*color_gradient.z);
		if(norm > SPHSim::EPSILON)
		{
			float inv_norm = 1.f / norm;
			float c = -conf.tensionCoe * color_laplacian * inv_norm;
			ftension.x = c * color_gradient.x; 
			ftension.y = c * color_gradient.y;
			ftension.z = c * color_gradient.z;
		}
		d_force[3*particleId] = fpressure.x + fviscosity.x + fgravity.x + ftension.x;
		d_force[3*particleId+1] = fpressure.y + fviscosity.y + fgravity.y + ftension.y;
		d_force[3*particleId+2] = fpressure.z + fviscosity.y + fgravity.z + ftension.z;

	}
}
 
// Host code
namespace SPHSim
{
	SPHSimulator::SPHSimulator(SPHConfig config)
	:
	conf(config)
	{
  		dim_x = ceil(conf.box.spanX / conf.kernel_h);
  		dim_y = ceil(conf.box.spanY / conf.kernel_h);
  	    dim_z = ceil(conf.box.spanZ / conf.kernel_h);
  		arraySize = dim_x * dim_y * dim_z;

        // printf("SpatialGrid initialized:\n");
		// printf("arraySize: %d\n", arraySize);
		// printf("grid_h: %f\n", conf.kernel_h);

		gravityDir = vec3(0.f,-1.f,0.f);
	}

	SPHSimulator::~SPHSimulator()
	{
		delete [] host_position;
		delete [] host_velocity;
		delete [] host_force;
		delete [] host_mass;

		cudaFree( device_position );
		cudaFree( device_velocity );
		cudaFree( device_force );
		cudaFree( device_mass );
		cudaFree( device_pressure );
		cudaFree( device_density );
		cudaFree( device_index );
		cudaFree( device_bucket_index );
		cudaFree( device_bucket_mask );
		cudaFree( device_scan );
	}

	void SPHSimulator::setup()
	{
#ifdef SPH_TEST

		double mass = 1.f;
		double start_x =  conf.box.spanX/5;
		double start_y =  conf.box.spanY/5;
		double start_z =  conf.box.spanZ/5;
		double end_x   =  (conf.box.spanX/5) * 4;
		double end_y   =  (conf.box.spanY/5) * 4;
		double end_z   =  (conf.box.spanZ/5) * 4;

		double f = 0.9;

		for (float x = start_x; x < end_x; x += conf.kernel_h * f)
			for (float y = start_y; y < end_y; y += conf.kernel_h * f)
				for (float z = start_z; z < end_z; z += conf.kernel_h * f)
					particles.push_back(SPHParticle(mass, vec3(x,y,z)));

		std::cout << "particle number:" << particles.size() << std::endl;


#endif
		N = particles.size();
		/*  <Device memory>  */
		// position
		cudaMalloc((void **)&device_position, sizeof(double) * N * 3); 
		// velocity
		cudaMalloc((void **)&device_velocity, sizeof(double) * N * 3);
		// force
		cudaMalloc((void **)&device_force, sizeof(double) * N * 3);
		// mass
		cudaMalloc((void **)&device_mass, sizeof(float) * N);
		// density
		cudaMalloc((void **)&device_density, sizeof(float) * N);
		// presure
		cudaMalloc((void **)&device_pressure, sizeof(float) * N);
		// index
		cudaMalloc((void **)&device_index, sizeof(int) * N);
		// bucket index
		cudaMalloc((void **)&device_bucket_index, sizeof(int) * N);
		// bucket
		cudaMalloc((void **)&device_bucket, sizeof(Bucket) * arraySize);
		// bucket mask
		cudaMalloc((void** )&device_bucket_mask, sizeof(int) * arraySize);
		// exclusive scan
		cudaMalloc((void** )&device_scan, sizeof(int) * arraySize);

		/*  <Host memory>  */
		// position
		host_position = new double[N * 3];
		// velocity	
		host_velocity = new double[N * 3];
		// force
		host_force = new double[N * 3];
		// mass
		host_mass = new float[N];

		return;
	}

	void SPHSimulator::reset()
	{
		return;
	}

	void SPHSimulator::update(float delta)
	{


//TODO: for all particles, find neighborhoods
		/*
		>>>> copy host particle attributes to device
		*/
		//double t1_cpHtoD = omp_get_wtime();

		for (int i = 0; i < N; i++)
		{
			host_position[3 * i    ] = particles[i].position.x;
			host_position[3 * i + 1] = particles[i].position.y;
			host_position[3 * i + 2] = particles[i].position.z;
			host_velocity[3 * i    ] = particles[i].velocity.x;
			host_velocity[3 * i + 1] = particles[i].velocity.y;
			host_velocity[3 * i + 2] = particles[i].velocity.z;
			host_mass[i] = particles[i].mass;
		}
		cudaMemcpy( device_position, host_position, sizeof(double) * N * 3, cudaMemcpyHostToDevice );
		cudaMemcpy( device_velocity, host_velocity, sizeof(double) * N * 3, cudaMemcpyHostToDevice );
		cudaMemcpy( device_mass, host_mass, sizeof(float) * N, cudaMemcpyHostToDevice);

		// double t2_cpHtoD = omp_get_wtime();
		// std::cout << "Delta cp host to device : " 
		// 		<< (t2_cpHtoD - t1_cpHtoD) * 1000
		// 		<< " ms" << std::endl;

		dim3 threadsPerBlock(32, 1, 1);
		dim3 numBlocks(N/32+1, 1, 1);		
		// double t1_kernel = omp_get_wtime();

		kernelComputeIndex<<<numBlocks, threadsPerBlock>>>(device_position, device_index, device_bucket_index, N,
			 												conf.kernel_h, dim_x, dim_y);
		cudaDeviceSynchronize();
		
		thrust::sort_by_key(thrust::device_pointer_cast<int>(device_bucket_index),
							thrust::device_pointer_cast<int>(device_bucket_index) + N,
							thrust::device_pointer_cast<int>(device_index));
		cudaDeviceSynchronize();

		double* device_sorted_position;
		double* device_sorted_velocity;
		cudaMalloc((void** )&device_sorted_position, sizeof(double) * N * 3);
		cudaMalloc((void** )&device_sorted_velocity, sizeof(double) * N * 3);
		kernelSortParticle<<<numBlocks, threadsPerBlock>>>(device_position, device_sorted_position, device_velocity, device_sorted_velocity, device_index, N);
		cudaDeviceSynchronize();
		cudaMemcpy(device_position, device_sorted_position, sizeof(double) * N * 3, cudaMemcpyDeviceToDevice);
		cudaMemcpy(device_velocity, device_sorted_velocity, sizeof(double) * N * 3, cudaMemcpyDeviceToDevice);
		cudaFree(device_sorted_position);
		cudaFree(device_sorted_velocity);

		cudaMemset(device_bucket, 0, sizeof(Bucket) * arraySize);
		kernelComputeBucket<<<numBlocks, threadsPerBlock>>>(device_bucket, device_bucket_index, N, arraySize);
		cudaDeviceSynchronize();

		// make bucket
		dim3 threadsPerBlock2(256, 1, 1);
		dim3 numBlocks2((arraySize/256)+1, 1, 1);		
		
		kernelSplitBlock<<<numBlocks2, threadsPerBlock2>>>(device_bucket, device_bucket_mask, arraySize);
		cudaDeviceSynchronize();

		
		thrust::exclusive_scan(thrust::device_pointer_cast<int>(device_bucket_mask), 
							   thrust::device_pointer_cast<int>(device_bucket_mask) + arraySize,
							   thrust::device_pointer_cast<int>(device_scan));
		cudaDeviceSynchronize();

		int numNonZeroBucket;
		cudaMemcpy((int*)&numNonZeroBucket, device_scan+arraySize-1, sizeof(int), cudaMemcpyDeviceToHost);
		cudaMalloc((void** )&device_compressed_bucket, sizeof(int) * numNonZeroBucket);
		kernelMakeBlock<<<numBlocks2, threadsPerBlock2>>>(device_compressed_bucket, device_bucket_mask, device_scan, arraySize);
		cudaDeviceSynchronize();


		int* host_compressed_bucket = new int[numNonZeroBucket];
		cudaMemcpy(host_compressed_bucket, device_compressed_bucket, sizeof(int) * numNonZeroBucket, cudaMemcpyDeviceToHost);

		// double t1_density = omp_get_wtime();

		dim3 threadsPerBlock3(32, 1, 1);
		dim3 numBlocks3(numNonZeroBucket, 1, 1);		
		kernelComputeDensity<<<numBlocks3, threadsPerBlock3>>>(device_position, device_bucket, device_compressed_bucket,
												device_density, device_mass, device_pressure, N, numNonZeroBucket, conf, dim_x, dim_y, arraySize);
		cudaDeviceSynchronize();

		// double t2_density = omp_get_wtime();
		// std::cout << "Delta density : " 
		// 		<< (t2_density - t1_density) * 1000
		// 		<< " ms" << std::endl;

		//double t1_force = omp_get_wtime();
		kernelComputeForce<<<numBlocks3, threadsPerBlock3>>>(device_position, device_velocity, device_bucket, device_compressed_bucket,
									device_force, device_density, device_mass, device_pressure, N, numNonZeroBucket, conf, dim_x, dim_y, arraySize,
									gravityDir.x, gravityDir.y, gravityDir.z);
		cudaDeviceSynchronize();
		// double t2_force = omp_get_wtime();
		// std::cout << "Delta force : " 
		// 		<< (t2_force - t1_force) * 1000
		// 		<< " ms" << std::endl;

		kernelUpdateAndHandle<<<numBlocks, threadsPerBlock>>>(device_position, device_velocity, device_force,
																					device_density, N, conf, delta, conf.box.spanX, conf.box.spanY, conf.box.spanZ);
		cudaDeviceSynchronize();

		// double t2_kernel = omp_get_wtime();
		// std::cout << "Delta kernel : "
		// 		  << (t2_kernel-t1_kernel) * 1000
		// 		  << " ms" << std::endl;

		cudaFree(device_compressed_bucket);
		/*
		>>>> copy devices particle attributes to host
		*/
		// double t1_cpDtoH = omp_get_wtime();
		cudaMemcpy( host_position, device_position, sizeof(double) * N * 3, cudaMemcpyDeviceToHost );
		cudaMemcpy( host_velocity, device_velocity, sizeof(double) * N * 3, cudaMemcpyDeviceToHost );

		for (int i = 0; i < particles.size(); i++)
		{
			particles[i].position.x = host_position[3 * i];
			particles[i].position.y = host_position[3 * i + 1];
			particles[i].position.z = host_position[3 * i + 2];
			particles[i].velocity.x = host_velocity[3 * i];
			particles[i].velocity.y = host_velocity[3 * i + 1];
			particles[i].velocity.z = host_velocity[3 * i + 2];

		 }
		// double t2_cpDtoH = omp_get_wtime();
		// std::cout << "Delta copy device to host : "
		// 		<< (t2_cpDtoH - t1_cpDtoH) * 1000
		// 		<< " ms" << std::endl;
	}

	void SPHSimulator::getData(float** position, int& vertexNum)
	{
		//allocate new memory
   		if(*position == NULL)
		   *position = (float*)malloc(sizeof(float) * particles.size() * 3);
		vertexNum = 0;

		for(int i=0; i<particles.size(); i++)
		{
			if(true)
			{
				(*position)[3*vertexNum+0] = particles[i].position.x;
				(*position)[3*vertexNum+1] = particles[i].position.y;
				(*position)[3*vertexNum+2] = particles[i].position.z;
				vertexNum ++;
			}
		}
		return;
	}

}
