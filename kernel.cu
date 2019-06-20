#include <cuda_runtime.h>
#include <cuda.h>
#include <math_functions.h>

#include "kernel.h"

__device__ float tempParticle1[NUM_OF_DIMENSIONS];
__device__ float tempParticle2[NUM_OF_DIMENSIONS];


__device__ float cuda_calculate_F(float x,float y, bool b_[Nb_x][Nb_y], struct point2d r[Nb_x][Nb_y], float phi[Nb_x][Nb_y])
{
  float d = 0;
  for (int i = 0; i < Nb_x; i++)
  {
    for (int j = 0; j < Nb_y; j++)
    {
      if (b_[i][j] == 0) d = d + phi[i][j]*sqrt((x-r[i][j].x)*(x-r[i][j].x) + (y-r[i][j].y)*(y-r[i][j].y));
    }
  } 
  return d;
}

__device__ float cuda_calc_d(float x1,float y1,float x2,float y2)
{
    return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

__device__ bool cuda_check_if_grid_is_inside_R(float x,float y, struct point2d r)
{
    float rx1 = r.x + res/2;
    float rx2 = r.x + res/2;
    float rx3 = r.x - res/2;
    float rx4 = r.x - res/2;

    float ry1 = r.y - res/2;
    float ry2 = r.y + res/2;
    float ry3 = r.y + res/2;
    float ry4 = r.y - res/2;
    if ( (cuda_calc_d(x,y,rx1,ry1) < R) && (cuda_calc_d(x,y,rx2,ry2) < R) && (cuda_calc_d(x,y,rx3,ry3) < R) && (cuda_calc_d(x,y,rx4,ry4) < R) )
    {
        return 1;
    }
    else 
    {
        return 0;
    }
}

__device__ void cuda_forward_euler(float *x, float *y, float *psi, float v_w, float psi_w, float u_phi[], float u_v[], float Dt)
{   
    float dx, dy, v_g[N+1], chi[N+1], dchi;

    v_g[0] = sqrt( pow(u_v[0]*cos(psi[0])+v_w*cos(psi_w), 2) + pow(u_v[0]*sin(psi[0])+v_w*sin(psi_w), 2) );

    chi[0] = atan2( ( u_v[0]*sin(psi[0]) + v_w*sin(psi_w) ) , ( u_v[0]*cos(psi[0]) + v_w*cos(psi_w) ));
    if (chi[0] < 0) chi[0] = chi[0] + 2*pi;

  //FORWARD EULER METHOD TO INTEGRATE THE STATES OF AGENT (Eq. 02) (Eq. 04)
  for (int k = 0; k < N; k++)
  {
    //x
    dx = v_g[k]*cos(chi[k]);
    x[k+1] = x[k] + dx*Dt;
    //y
    dy = v_g[k]*sin(chi[k]);
    y[k+1] = y[k] + dy*Dt;
    //chi
    dchi = g*tan(u_phi[k])*cos(chi[k] - psi[k])/v_g[k];
    chi[k+1] = chi[k] + dchi*Dt;
    if (chi[k+1]>2*pi) chi[k+1] = chi[k+1]-2*pi;
    //psi
    psi[k+1] = chi[k+1] - asin(v_w/u_v[k+1]*sin(psi_w - chi[k+1]));
    //v_g
    v_g[k+1] = sqrt( pow(u_v[k+1]*cos(psi[k+1])+v_w*cos(psi_w), 2) + pow(u_v[k+1]*sin(psi[k+1])+v_w*sin(psi_w), 2) );
  }
}

__device__ float cuda_cost_function(int index, float controls[], struct agent agents[], struct point2d r[Nb_x][Nb_y], float phi[Nb_x][Nb_y])
{
    for (int k = 0; k < N; k++) agents[index].u_phi[k] = controls[k];
    for (int k = 0; k < N; k++) agents[index].u_v[k] = controls[k+N];
    cuda_forward_euler(agents[index].x, agents[index].y, agents[index].psi, agents[index].v_w, agents[index].psi_w, agents[index].u_phi, agents[index].u_v, Dt_MPC);
    
    //EVALUATE COST FUNCTION
    float F = 0;
    float total_cost = 0.0;
    float sum_phi_y = 0;
    float d[I];

    for (int i = 0; i < I; i++)
    {
        if (i != index) 
        {
          for (int j = 0; j < Nb_x; j++)
          {
            for (int k = 0; k < Nb_y; k++)
            {
              agents[index].B[j][k] = agents[index].B[j][k] || agents[i].b[j][k];
            }
          }
        }
    } 

    //anti colision for future steps. current step doesn't matter
    for (int k = 1; k < N; k++)
    {
      for (int i = 0; i < I; i++) 
      {
        d[i] = sqrt(pow(agents[index].x[k] - agents[i].x[k], 2) + pow(agents[index].y[k] - agents[i].y[k], 2));
        if ((d[i] < r_c) && (i != index)) return 99999999999999999999999999999999999999999999999999999999999999999999.9;
      }
    }

    for (int k = 0; k < (N+1); k++)
    { 
        //Lagrangian term
        float phi_y = 0;
        //check if 50x50 cells are visited
        int x_idx = (int)agents[index].x[k]/res;
        int y_idx = (int)agents[index].y[k]/res;
        for (int i = (x_idx - R/res); i <= (x_idx + R/res); i++)
        {
            for (int j = (y_idx - R/res); j <= (y_idx + R/res); j++)
            {
                //printf("x=%f,y=%f \n",x1[k],y1[k]);
                if (cuda_check_if_grid_is_inside_R(agents[index].x[k],agents[index].y[k],r[i][j]))
                {
                    agents[index].B[i][j] = 1;
                }
            }
        }
        for (int i = 0; i < Nb_x; i++)
        {
            for (int j = 0; j < Nb_y; j++)
            {
                phi_y = phi_y + phi[i][j]*agents[index].B[i][j];
            }
        }
        sum_phi_y = sum_phi_y + phi_y;
    }

    F = cuda_calculate_F(agents[index].x[N],agents[index].y[N],agents[index].B,r,phi); //r_x[Nb], r_y[Nb]: position of center of cells; b[Nb] is the binary of cell visit
    total_cost = F - sum_phi_y; //Eq 13
    return total_cost;
}

__global__ void kernelUpdateParticle(float *positions, float *velocities, 
                                     float *pBests, float *gBest, float r1, 
                                     float r2)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if(i >= NUM_OF_PARTICLES * NUM_OF_DIMENSIONS)
        return;
    
    float rp = r1;
    float rg = r2;

    velocities[i] = OMEGA * velocities[i] + c1 * rp * (pBests[i] - positions[i])
            + c2 * rg * (gBest[i % NUM_OF_DIMENSIONS] - positions[i]);

    // Update posisi particle
    positions[i] += velocities[i];

    if (((i % (N*2)) < N) && (positions[i] < (phi_min))) positions[i] = phi_min;
    if (((i % (N*2)) < N) && (positions[i] > (phi_max))) positions[i] = phi_max;

    if (((i % (N*2)) >= N) && (positions[i] < v_min)) positions[i] = v_min;
    if (((i % (N*2)) >= N) && (positions[i] > v_max)) positions[i] = v_max;
}
//ok
__global__ void kernelUpdatePBest(int index, float *positions, float *pBests, float* gBest, struct agent agents[], struct point2d r[Nb_x][Nb_y], float phi[Nb_x][Nb_y])
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if(i >= NUM_OF_PARTICLES * NUM_OF_DIMENSIONS || i % NUM_OF_DIMENSIONS != 0)
        return;

    for (int j = 0; j < NUM_OF_DIMENSIONS; j++)
    {
        tempParticle1[j] = positions[i + j];
        tempParticle2[j] = pBests[i + j];
    }

    if (cuda_cost_function(index,tempParticle1,agents,r,phi) < cuda_cost_function(index,tempParticle2,agents,r,phi))
    {
        for (int k = 0; k < NUM_OF_DIMENSIONS; k++)
            pBests[i + k] = positions[i + k];       
    }
}

extern "C" void cuda_pso(int index, float positions[], float velocities[], float pBests[], 
                         float *gBest, struct agent agents[], struct point2d r[Nb_x][Nb_y], float phi[Nb_x][Nb_y])
{
    int size = NUM_OF_PARTICLES * NUM_OF_DIMENSIONS;
    float *devPos;
    float *devVel;
    float *devPBest;
    float *devGBest;
    float temp[NUM_OF_DIMENSIONS];
    // Memory allocation
    cudaMalloc((void**)&devPos, sizeof(float) * size);
    cudaMalloc((void**)&devVel, sizeof(float) * size);
    cudaMalloc((void**)&devPBest, sizeof(float) * size);
    cudaMalloc((void**)&devGBest, sizeof(float) * NUM_OF_DIMENSIONS);
    // Thread & Block number
    int threadsNum = 256;
    int blocksNum = NUM_OF_PARTICLES / threadsNum;
    // Copy particle datas from host to device
    cudaMemcpy(devPos, positions, sizeof(float) * size, cudaMemcpyHostToDevice);
    cudaMemcpy(devVel, velocities, sizeof(float) * size, 
               cudaMemcpyHostToDevice);
    cudaMemcpy(devPBest, pBests, sizeof(float) * size, cudaMemcpyHostToDevice);
    cudaMemcpy(devGBest, gBest, sizeof(float) * NUM_OF_DIMENSIONS, 
               cudaMemcpyHostToDevice);

    // PSO main function
    for (int iter = 0; iter < MAX_ITER; iter++)
    {     
        
        // Update position and velocity
        kernelUpdateParticle<<<blocksNum, threadsNum>>>(devPos, devVel, 
                                                        devPBest, devGBest, 
                                                        getRandomClamped(), 
                                                        getRandomClamped());
        // Update pBest
        kernelUpdatePBest<<<blocksNum, threadsNum>>>(index, devPos, devPBest, 
                                                     devGBest,agents,r,phi);
        // Update gBest
        cudaMemcpy(pBests, devPBest, 
                   sizeof(float) * NUM_OF_PARTICLES * NUM_OF_DIMENSIONS, 
                   cudaMemcpyDeviceToHost);
        for(int i = 0; i < size; i += NUM_OF_DIMENSIONS)
        {
            for(int k = 0; k < NUM_OF_DIMENSIONS; k++)
                temp[k] = pBests[i + k];
        
            if (cost_function(index,temp,agents,r,phi) < cost_function(index, gBest,agents,r,phi))
            {
                for (int k = 0; k < NUM_OF_DIMENSIONS; k++)
                    gBest[k] = temp[k];
            }
        }
        cudaMemcpy(devGBest, gBest, sizeof(float) * NUM_OF_DIMENSIONS, 
                   cudaMemcpyHostToDevice);
    }
    // Retrieve particle datas from device to host
    cudaMemcpy(positions, devPos, sizeof(float) * size, cudaMemcpyDeviceToHost);
    cudaMemcpy(velocities, devVel, sizeof(float) * size, 
               cudaMemcpyDeviceToHost);
    cudaMemcpy(pBests, devPBest, sizeof(float) * size, cudaMemcpyDeviceToHost);
    cudaMemcpy(gBest, devGBest, sizeof(float) * NUM_OF_DIMENSIONS, 
               cudaMemcpyDeviceToHost); 
    // cleanup
    cudaFree(devPos);
    cudaFree(devVel);
    cudaFree(devPBest);
    cudaFree(devGBest);
}
