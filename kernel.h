const int 		I = 3; 						//Number of agents

//Limits for the search and rescue area [m]
const float x_min = 0; 
const float x_max = 5000;
const float y_min = 0;
const float y_max = 5000;

const float res = 100; //resolution of the grids [m]
const float R = 200.0; //sensor radius [m]

const int Nb_x = (x_max-x_min)/res; //number of grids in x axis: 50
const int Nb_y = (y_max-y_min)/res; //number of grids in y axis: 50

const float 	horizon = 25.0; 				//time horizon [s]
const int 		N = 25; 						//N (discretization steps)
const int 		N_max = 30;						//N max (discretization steps) - IMC message
const float 	Dt_MPC = horizon/(float)N; 		//delta time [s]

const int 		NUM_OF_PARTICLES = 1024;	//number of PSO particles
const int 		NUM_OF_DIMENSIONS = N*2;	//number of PSO dimensions
const int 		MAX_ITER = 40;				//number of iterations
const float 	c1 = 1;						//PSO local coefficient
const float 	c2 = 2;						//PSO global coefficient
const float 	OMEGA = 1;					//inertia weight

const float 	pi = 3.14159265359;			//the constant PI
const float 	g = 9.81; 					//gravity acceleration [m/s^2]

const float 	a = 10.0;  					//for Eq. 11
const float 	b = 1.0;    				//for Eq. 11
const float 	c = 10.0;    				//for Eq. 11

const float 	r_c = 100.0;    			//anti colision

const float 	v_min = 10.0; 				//min airspeed [m/s]
const float 	v_max = 22.0; 				//max airspeed [m/s]
const float 	phi_min = -45*pi/180; 		//min roll angle [rad]
const float 	phi_max = 45*pi/180;  		//max roll angle [rad]

const float 	dv_a = 30.0; 				//airspeed rate (acceleration) [m/s2] (abs)
const float 	dphi = 150*pi/180; 			//min roll angle rate [rad/s]

//X8 aerodynamic model
const float rho             = 1.225;    //Air density 
const float S               = 0.74;     //Wing surface - Kristoffer says 0.75 in his thesis, but I think 0.65 (effective) is more realistic
const float W               = 32.96;      //Aircraft weight (Newtons)
const float eta_p_bldc      = 0.5;      //Propulsion efficiency - shotuld be adjusted to v_br 
const float k_               = 1.0 / (pi * 3.9 * 0.9);

//Reference latitude and longitude of the MPC-PSO NED Frame
const double rlat = 61.561599*pi/180;            //reference latitude
const double rlon = 4.839016*pi/180;             //reference longitude


//Declaring here as static because of segmentation fault (core dumped)
// Particle
static float positions[NUM_OF_PARTICLES * NUM_OF_DIMENSIONS];
static float velocities[NUM_OF_PARTICLES * NUM_OF_DIMENSIONS];
static float pBests[NUM_OF_PARTICLES * NUM_OF_DIMENSIONS];
// gBest
static float gBest[NUM_OF_DIMENSIONS];
static float temp[NUM_OF_DIMENSIONS];

static float phi_out[Nb_x][Nb_y];


struct point2d 
{
    float x;
    float y;
};

static struct point2d r_out[Nb_x][Nb_y];             //x,y position of grid 


struct agent {
	int index;
	double lat;
	double lon;
	float height;
	float z;
	float u_phi[N_max] = { 0 };
	float u_v[N_max] = { 0 };
	float x[N_max+1];
	float y[N_max+1];
	float psi[N_max+1];
	float v_a0;
	float phi0;
	float v_w = 9.899494937;
	float psi_w = 0.785398163;
	int type = 0;	//type (0 UAV, 1 ASV, 2 GS)
	float constraints[4]; //dphi_min, dphi_max, dv_min, dv_max
	bool B[Nb_x][Nb_y] = {{0}}; //if cell was already visited
    bool b[Nb_x][Nb_y] = {{0}}; //planned cells to visit
};

float 	cost_function(int index, float controls[], struct agent agents[], struct point2d r[Nb_x][Nb_y], float phi[Nb_x][Nb_y]);
void 	forward_euler(float *x, float *y, float *psi, float v_w, float psi_w, float u_phi[], float u_v[], float Dt);
float 	getRandom(float low, float high);
float 	getRandomClamped();
void 	PSO(int index, struct agent *agents, struct point2d r[Nb_x][Nb_y], float phi[Nb_x][Nb_y]);

float 	calculate_F(float x,float y, bool b_[Nb_x][Nb_y], struct point2d r[Nb_x][Nb_y], float phi[Nb_x][Nb_y]);
bool 	check_if_grid_is_inside_R(float x, float y, struct point2d r);
float 	calc_d(float x1,float y1,float x2,float y2);

extern "C" void cuda_pso(int index, float positions[], float velocities[], float pBests[], 
                       float *gBest, struct agent agents[], struct point2d r[Nb_x][Nb_y], float phi[Nb_x][Nb_y]);
