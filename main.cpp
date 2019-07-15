	#include <string>
	#include <sstream>
	#include <cstdlib>
	#include <random>
	#include <cmath>
	#include "iostream"
	#include "Grid/grid_dist_id.hpp"
	#include "data_type/aggregate.hpp"
	#include "Decomposition/CartDecomposition.hpp"
	#include "timer.hpp"

	// Definition of molecules indexes (for concentration fiels)

const size_t I = 0;
const size_t R = 1;
const size_t RA2 = 2;
const size_t A = 3;
const size_t RA = 4;

const size_t Parameters = 5;

	// Global markers
const size_t DIFFUSION_NOISE = 1;
const size_t EXTRINSIC_NOISE = 1;

	// Function evolve_time
void evolve_time(  grid_dist_id<2,  double, aggregate< double, double, double, double, double, double[23]>> & g_dist_read,   grid_dist_id<2,  double, aggregate< double, double, double, double, double, double[23]>> & g_dist_write, const  double T, const  double sT, openfpm::vector< double> & x5, std::mt19937_64 & en, std::normal_distribution< double> & NDist, double stats_vec[],size_t Ncells, int j,double dAe, double Induction)
{
	Vcluster<> & v_cl = create_vcluster();

	  //Init means
	double tot_A_partial = 0.0;
	double x1mean = 0.0;
	double x1var = 0.0;
	double x2mean = 0.0;
	double x2var = 0.0;
	double x3mean = 0.0;
	double x3var = 0.0;
	double x4mean = 0.0;
	double x4var = 0.0;

	double x1 = 0.0;
	double x2 = 0.0;
	double x3 = 0.0;
	double x4 = 0.0;
	size_t primera = 0;

	  // Create iterator for the distributed grid
	  // Iterator is like a for loop, but more intelligent and efficient.
	  // It loops over all the elements of the grid even if they are in different processors.
	auto dome = g_dist_read.getDomainIterator();
	  // Iterate the iterator
	while (dome.isNext())
	{
	    //  std::cout << "In iterator dome" << std::endl;

	    // Get the actual position from the iterator in the subdomain
		auto key = dome.get();
	    // Get the global position, from the local one.
		auto key_g = g_dist_read.getGKey(key);

	    // Copy parameters from the distributed grid variable Parameters into the named variables

	    // I parameters
		double dI = g_dist_read.template get<Parameters>(key)[0];
		double pI = g_dist_read.template get<Parameters>(key)[1];
		double kI = g_dist_read.template get<Parameters>(key)[2];
		double pN_I = g_dist_read.template get<Parameters>(key)[3];
		double dmI = g_dist_read.template get<Parameters>(key)[4];
		double kdLux = g_dist_read.template get<Parameters>(key)[5];
		double alphaI = g_dist_read.template get<Parameters>(key)[6];
	    // R Parameters
		double dR = g_dist_read.template get<Parameters>(key)[7];
		double pR = g_dist_read.template get<Parameters>(key)[8];
		double cR =  g_dist_read.template get<Parameters>(key)[9];
		double dmR = g_dist_read.template get<Parameters>(key)[10];
		double k_1 =  g_dist_read.template get<Parameters>(key)[11];
		double kd1 =  g_dist_read.template get<Parameters>(key)[12];
	    // Dimer parameters
		double k_2 = g_dist_read.template get<Parameters>(key)[13];
		double kd2 = g_dist_read.template get<Parameters>(key)[14];
		double dRA2 = g_dist_read.template get<Parameters>(key)[15];
	    // A parameters
		double kA =  g_dist_read.template get<Parameters>(key)[16];
		double dA = g_dist_read.template get<Parameters>(key)[17];
		double D =  g_dist_read.template get<Parameters>(key)[18];
		double Vcell =  g_dist_read.template get<Parameters>(key)[19];
		double Vext =  g_dist_read.template get<Parameters>(key)[20];

	    // RA parameter
		double dRA = g_dist_read.template get<Parameters>(key)[21];


	    // Propensities
	    // I 
		double X11 = dI * g_dist_read.template get<I>(key);
		double c1 = (pI * kI * pN_I)/dmI;
		double c2 = 1/(kdLux + g_dist_read.template get<RA2>(key));
		double X12 = c1*c2*(kdLux + alphaI * g_dist_read.template get<RA2>(key));

	    // R
		double  X21 = dR * g_dist_read.template get<R>(key);
		double  X22 = pR*cR/dmR;
		double  X23 = k_1 * g_dist_read.template get<RA>(key);
		double  X24 = ( k_1 / kd1 ) * g_dist_read.template get<R>(key) * g_dist_read.template get<A>(key);

	    // RA2
		double X31 = (k_2 + dRA2)* g_dist_read.template get<RA2>(key);
		double X32 = k_2 * g_dist_read.template get<RA>(key) * g_dist_read.template get<RA>(key) / kd2;

	    // A
		double X41 = X23;
		double X42 = kA * g_dist_read.template get<I>(key);
		double X43 = dA * g_dist_read.template get<A>(key);
		double X44 = D * (Vcell/Vext * x5.last() - g_dist_read.template get<A>(key));
		double X45 = X24;

	    // Noises
		double x1_noise1 = NDist(en);
		double x1_noise2 = NDist(en);

		double x2_noise1 = NDist(en);
		double x2_noise2 = NDist(en);
		double x2_noise3 = NDist(en);
		double x2_noise4 = NDist(en);

		double x3_noise1 = NDist(en);
		double x3_noise2 = NDist(en);

		double x4_noise1 = x2_noise3;
		double x4_noise2 = NDist(en);
		double x4_noise3 = NDist(en);
		double x4_noise_difu = NDist(en);
		double x4_noise5 = x2_noise4;
		if (DIFFUSION_NOISE == 0){
			x4_noise1 = NDist(en);
			x4_noise5 = NDist(en);
			x4_noise_difu = 0;
		}

	    // Deterministic and stochastic terms with stoicheometry included
		double x1_det = T*( -X11 + X12 );
		double x1_sto = sT*( -sqrt(std::abs(X11))*x1_noise1 + sqrt(std::abs(X12))*x1_noise2 );

		double x2_det = T*( -X21 + X22 + X23 - X24  );
		double x2_sto = sT*( -sqrt(std::abs(X21))*x2_noise1 + sqrt(std::abs(X22))*x2_noise2 + sqrt(std::abs(X23))*x2_noise3 - sqrt(std::abs(X24))*x2_noise4 );

		double x3_det = T*( -X31 + X32 );
		double x3_sto = sT*( -sqrt(std::abs(X31))*x3_noise1 + sqrt(std::abs(X32))*x3_noise2 );

		double x4_det = T*( X41 + X42 - X43 + X44 - X45 );
		double x4_sto = sT*( sqrt(std::abs(X41))*x4_noise1 + sqrt(std::abs(X42))*x4_noise2 - sqrt(std::abs(X43))*x4_noise3 + sqrt(std::abs(X44))*x4_noise_difu - sqrt(std::abs(X45))*x4_noise5 );

	    // Compute new value
		g_dist_write.template get<I>(key) = g_dist_read.template get<I>(key) + x1_det + x1_sto;
		g_dist_write.template get<R>(key) = g_dist_read.template get<R>(key) + x2_det + x2_sto;
		g_dist_write.template get<RA2>(key) = g_dist_read.template get<RA2>(key) + x3_det + x3_sto;
		g_dist_write.template get<A>(key) = g_dist_read.template get<A>(key) + x4_det + x4_sto;

		if ( g_dist_write.template get<I>(key) < 0 ) {
			g_dist_write.template get<I>(key) = 0;
		}
		if ( g_dist_write.template get<R>(key) < 0 ) {
			g_dist_write.template get<R>(key) = 0;
		}
		if ( g_dist_write.template get<RA2>(key) < 0 ) {
			g_dist_write.template get<RA2>(key) = 0;
		}
		if ( g_dist_write.template get<A>(key) < 0 ) {
			g_dist_write.template get<A>(key) = 0;
		}

	    // Algebraic constraints
		double c6 = 2*k_2*kd1*g_dist_write.template get<RA2>(key) + k_1 * g_dist_write.template get<R>(key) * g_dist_write.template get<A>(key);
		double c7 = 8*k_2*c6;
		double c8 = k_1 + dRA;
		double c9 = (kd1*kd2)*c8*c8;
		double c10 = (kd2*c8)/(4*k_2);
		double c11 = c7/c9 + 1;
	    g_dist_write.template get<RA>(key)  = c10*(sqrt(std::abs(c11)) - 1);  //% This is R.A

	    // Write the partial value of Ae from the present cell.
	    tot_A_partial += T* (-1)* X44 + sT * (-1) * sqrt(std::abs(X44))*x4_noise_difu;

	    // Add to the mean the term corresponding to the present cell.
	    x1mean += g_dist_write.template get<I>(key)/ Ncells;
	    x2mean += g_dist_write.template get<R>(key)/ Ncells;
	    x3mean += g_dist_write.template get<RA2>(key)/ Ncells;
	    x4mean += g_dist_write.template get<A>(key)/ Ncells;

	    // Increment domain iterator.
	    ++dome;
	}

	  // Excecute the sum of the means and A_partial over all the processors
	v_cl.sum(  tot_A_partial );
	v_cl.sum(  x1mean );
	v_cl.sum(  x2mean );
	v_cl.sum(  x3mean );
	v_cl.sum(  x4mean );
	v_cl.execute();

	  // Create iterator for the distributed grid
	auto dom_var = g_dist_read.getDomainIterator();
	while (dom_var.isNext())
	{
	    // Get the actual position from the iterator in the subdomain
		auto key = dom_var.get();

	    //Calculate the variance incrementaly by adding the term corresponding to the i-th cell (with the iterator)
		x1var += (g_dist_write.template get<I>(key)-x1mean)*(g_dist_write.template get<I>(key)-x1mean)/(Ncells-1);
		x2var += (g_dist_write.template get<R>(key)-x2mean)*(g_dist_write.template get<R>(key)-x2mean)/(Ncells-1);
		x3var += (g_dist_write.template get<RA2>(key)-x3mean)*(g_dist_write.template get<RA2>(key)-x3mean)/(Ncells-1);
		x4var += (g_dist_write.template get<A>(key)-x4mean)*(g_dist_write.template get<A>(key)-x4mean)/(Ncells-1);
	    // Increment domain iterator.
		++dom_var;
	}
	  // Excecute the sum of the variances over all the processors
	v_cl.sum(  x1var );
	v_cl.sum(  x2var );
	v_cl.sum(  x3var );
	v_cl.sum(  x4var );
	v_cl.execute();

	  // Put the mean and variances of x's  into the stats_vec array
	stats_vec[0] = x1mean;
	stats_vec[1] = x1var;
	stats_vec[2] = x2mean;
	stats_vec[3] = x2var;
	stats_vec[4] = x3mean;
	stats_vec[5] = x3var;
	stats_vec[6] = x4mean;
	stats_vec[7] = x4var;

	  // Noise for x5
	double x5_noise = NDist(en);
	  // Stochastic part of x5, x5.last is the previous value of x5
	double x5_sto = sT*(-sqrt( std::abs(dAe * x5.last()) )*x5_noise);
	  //Calculate the new x5 value

	double x5to_add = x5.last() + T*(- dAe * x5.last() ) + x5_sto + tot_A_partial;

	  //Add the las calculated value to x5 vector
	x5.add( x5to_add );

	  // End of evolve_time function
}

void input_data(grid_dist_id<2,  double, aggregate< double, double, double, double, double, double[23]>> & g1,   grid_dist_id<2,  double, aggregate< double, double, double, double, double, double[23]>> & g2, std::mt19937_64 & en, std::normal_distribution< double> & NDist, size_t Ncells, double *dAe_p,double enoise_sigma)
{
	  // Function to get the parametrs from the file param.dat, and write them into the distributed grid, in each cell.

	  // Read from file
	std::ifstream ifs_param ("param.dat", std::ifstream::in);
	char delim;
	double parameters[23] = {};

	ifs_param >> parameters[0];
	for (int j = 1 ; j < 23 ; j++)
	{
		ifs_param >> delim;
		ifs_param >> parameters[j];
	}

	ifs_param.close();

	  // Create iterator for the distributed grid
	auto dom_init = g1.getDomainIterator();
	  // Iterate
	while (dom_init.isNext())
	{
	    // Get the actual position from the iterator in the subdomain
		auto key = dom_init.get();
	    // Initialize the grid
		g1.template get<I>(key) = 0.0;
		g1.template get<R>(key) = 0.0;
		g1.template get<RA2>(key) = 0.0;
		g1.template get<A>(key) = 0.0;
		g1.template get<RA>(key) = 0.0;

		for (int l = 0; l < 19; l++)
		{
			double aux_param = (EXTRINSIC_NOISE*enoise_sigma*NDist(en) + 1)*parameters[l];

			g1.template get<Parameters>(key)[l] = aux_param;
			g2.template get<Parameters>(key)[l] = aux_param;

		}
		g1.template get<Parameters>(key)[19] = parameters[19];
		g2.template get<Parameters>(key)[19] = parameters[19];

		g1.template get<Parameters>(key)[20] = parameters[20];
		g2.template get<Parameters>(key)[20] = parameters[20];

		double aux_param = (EXTRINSIC_NOISE*enoise_sigma*NDist(en) + 1)*parameters[21];
		g1.template get<Parameters>(key)[21] = aux_param;
		g2.template get<Parameters>(key)[21] = aux_param;

	    // Increment the iterator
		++dom_init;
	}
	  // Save the external dAe into the corresponding variable.
	*dAe_p =  parameters[22];
	  // End of function input_data
}


int main(int argc, char* argv[])
{
	  // Initialize the library 
	openfpm_init(&argc,&argv);
	Vcluster<> & v_cl = create_vcluster();

	  // The input parameters though the param file, is the 1th argument
	  // Reading file names from argv
	char param_file_name[40] ={};

	std::strcpy(argv[1],param_file_name);

	  // The number of cells, is the 2th argument
	std::istringstream iss(argv[2]);
	int val;

	if (!(iss >> val)) std::cerr << "Invalid number " << argv[2] << '\n';

	  // The variance of the noise, is the 3th argument
	std::istringstream iss3(argv[3]);
	double val3;
	if (!(iss3 >> val3)) std::cerr << "Invalid number " << argv[3] << '\n';

	  // The initial condition for A extrnal is the fourth argument
	std::istringstream iss4(argv[4]);
	double val4;
	if (!(iss4 >> val4)) std::cerr << "Invalid number " << argv[4] << '\n';

	  // Stochastic simulation FLAG is the fifth argument
	std::istringstream iss5(argv[5]);
	size_t STOCHASTIC;
	if (!(iss5 >> STOCHASTIC)) std::cerr << "Invalid number " << argv[5] << '\n';

	 // Writing histogram as an output (of the last point) FLAG is the sixth argument
	std::istringstream iss6(argv[6]);
	size_t WRITE_OUTPUT_H;
	if (!(iss6 >> WRITE_OUTPUT_H)) std::cerr << "Invalid number " << argv[6] << '\n';

	 // Writing temporal means and std dev as an output FLAG is the sixth argument
	std::istringstream iss7(argv[7]);
	size_t WRITE_OUTPUT_T;
	if (!(iss7 >> WRITE_OUTPUT_T)) std::cerr << "Invalid number " << argv[7] << '\n';

	  // Writing temporal of all cells as an output FLAG is the sixth argument
	std::istringstream iss8(argv[8]);
	size_t WRITE_OUTPUT_T_TOTAL;
	if (!(iss8 >> WRITE_OUTPUT_T_TOTAL)) std::cerr << "Invalid number " << argv[8] << '\n';

	  // Random number generator
	size_t seed1 = v_cl.getProcessUnitID();
	srand (time(NULL));
	size_t seed0 = rand();
	srand (seed0);
	size_t seed2 = getuid() + rand();

	seed2+=seed1<<16;
	seed1+=seed2<<11;
	seed2+=((signed int)seed1)>>7;
	seed1^=((signed int)seed2)>>3;
	seed2*=0xA5366B4D;
	seed2^=seed2>>10;
	seed2^=((signed int)seed2)>>19;
	seed1+=seed2^0x6d2d4e11;
	seed1 =0x79dedea3*(seed1^(((signed int)seed1)>>14));
	seed2 =(seed1 + seed2) ^ (((signed int)seed1)>>8);
	size_t MTseed=0xABCB96F7 + (seed2>>1);

	  std::mt19937_64 engine(MTseed);  //Use the global key as a seed for the PRNG
	  std::normal_distribution <double> NDist(0,1); // double Normal distribution with mean 0 and std 1

	  //Creation of the distributed grid domain
	  Box<2, double> domain({0.0,0.0},{1.0,1.0});
	  size_t a = 10;
	  size_t b = 24;
	  size_t Ncells = val;

      // Incremental determination of the size of the 2D domain
	  if (Ncells > 240)
	  	a = 50;
	  if (Ncells > 1200)
	  	b = 48;
	  if (Ncells > 2400)
	  	a = 100;
	  if (Ncells > 4800)
	  	b = 120;

	  size_t sz[2] {a,b};

	  // Ghost
	  Ghost<2, double> g(0.01);

	  // Vector for x5 (as normal variable, not a distriduted one)
	  openfpm::vector<double> x5;
	  // Vector for the statistics, mean and variance (as normal variable, not a distriduted one)
	  openfpm::vector<openfpm::vector<double>> stats;
	  // Vector for the histogram of I at final values
	  openfpm::vector<double> I_hist;
	    // Vector for the for the 4 states and all the cells (at time t) The first one is I
	  openfpm::vector<double> States_time_t;
	  openfpm::vector<double> States_time_t_R;
	  openfpm::vector<double> States_time_t_RA2;
	  openfpm::vector<double> States_time_t_A;

	  // Create a distributed grid in 2D
	  grid_dist_id<2,  double, aggregate< double, double, double, double, double, double[23]>> g1(sz,domain,g);
	  grid_dist_id<2,  double, aggregate< double, double, double, double, double, double[23]>>  g2(g1.getDecomposition(),sz,g);

	  // Initialize the grids and put the values of the changing parameters
	  double dAee = 0;
	  input_data(g1, g2, engine, NDist,  Ncells, &dAee,val3);

	  // Initialization of x5
	  x5.add(val4);

	  //
	  //  Now we start iterating in time to integrate the reaction equations ///
	  //

	  //Simulation parameters
	  double T = 25e-3;     // Simulation sampling time 
	  double sT = STOCHASTIC*sqrt(T);   // Euler-Maruyana
	  double Tsim = 800;     // sec total simulation time

	  // Time initialization for the simulation
	  // Convert timeline in number of simulation steps
	  double N = floor(Tsim/T);       // number of simulation steps round minus infinity
	  Tsim = N*T;
	  double t = 0;

	  // Define the size of stats 5 means, 4 variances and the time = 10
	  const int Size_stats = 10;

	  // Define the size of stats_vec 4 means, 4 variances
	  double stats_vec[8];


	  // Allocating memory for vectors: stats and hist.
	  stats.resize(N);
	  for (int k = 0; k < N; k++)
	  {
	  	stats.get(k).resize(Size_stats);
	  }

	  stats.get(0).get(0) = t;
	  for (int l = 1; l < 9; l++)
	  {
	  	stats.get(0).get(l) = 0;
	  }
	  stats.get(0).get(9) = x5.last();

	  // Initilize histo counter
	  size_t j_hist = 1;
	  // Initialize the cell counter every time (for States temporal writing)
	  size_t j_histE = 1;

	  // Now we start a loop in time to solve the differential equations
	  for (int j = 0; j < (N-1); j++)
	  {

	    // J increments by 1 in each time step. Always read the first argument and and write in second. But change the positions in memory swaping between 1 and 2.

	  	if (j%2 == 0)
	  	{
	  		evolve_time(g1, g2, T, sT, x5, engine, NDist, stats_vec, Ncells,j,dAee,val4);
	  	}
	  	else
	  	{
	  		evolve_time(g2, g1, T, sT, x5, engine, NDist, stats_vec, Ncells,j,dAee,val4);
	  	}

	    // Write the output if needed (WRITE_OUTPUT_T_TOTAL=1) of the whole grid (all the cells) every 10 time steps to save space.
	  	if (WRITE_OUTPUT_T_TOTAL && j%10 == 0)
	  	{
	      // Allocating memory for vector States_time_t
	  		States_time_t.resize(0);
	  		States_time_t_R.resize(0);
	  		States_time_t_RA2.resize(0);
	  		States_time_t_A.resize(0);

	      // Here put the code to write all the variables into a text file, and append every time step
		        // Save the last time point histogram
	      // Create iterator for g1
	  		auto dom_hist = g1.getDomainIterator();
	  		size_t contad = 0;

	      //Iterate
	  		while (dom_hist.isNext())
	  		{
	        // Get the actual position from the iterator in the subdomain
	  			auto key = dom_hist.get();
	        // Add the value to I_hist
	  			States_time_t.add(g1.template get<I>(key));
	  			States_time_t_R.add(g1.template get<R>(key));
	  			States_time_t_RA2.add(g1.template get<RA2>(key));
	  			States_time_t_A.add(g1.template get<A>(key));

	        //Increment the counter
	  			contad++;
	        //Increment the iterator
	  			++dom_hist;
	  		}

	  		openfpm::vector< double> States_col;
	  		openfpm::vector< double> States_colR;
	  		openfpm::vector< double> States_colRA;
	  		openfpm::vector< double> States_colA;

	      // Saving all the values into one processor
	  		v_cl.SGather(States_time_t,States_col,0);
	  		v_cl.SGather(States_time_t_R,States_colR,0);
	  		v_cl.SGather(States_time_t_RA2,States_colRA,0);
	  		v_cl.SGather(States_time_t_A,States_colA,0);

		  // Write the files only from Processot number 0
	  		if (v_cl.getProcessUnitID()==0 )
	  		{
	  			std::ofstream ofs_state_I ( "States_I.dat" , std::ofstream::out|std::ofstream::app);
	  			std::ofstream ofs_state_R ( "States_R.dat" , std::ofstream::out|std::ofstream::app);
	  			std::ofstream ofs_state_RA2 ( "States_RA2.dat" , std::ofstream::out|std::ofstream::app);
	  			std::ofstream ofs_state_A ( "States_A.dat" , std::ofstream::out|std::ofstream::app);

	        // The time j
	  			ofs_state_I << j*T;
	  			ofs_state_R << j*T;
	  			ofs_state_RA2 << j*T;
	  			ofs_state_A << j*T;

	        // Value of the states for each cell
	  			for (int m = 0 ; m < (States_col.size()-1) ; m++)
	  			{
	  				size_t numb = m+1;
	  				ofs_state_I << ", " << States_col.get(m);
	  				ofs_state_R << ", " << States_colR.get(m);
	  				ofs_state_RA2 << ", " << States_colRA.get(m);
	  				ofs_state_A << ", " << States_colA.get(m);
	  			}

	  			ofs_state_I << std::endl;
	  			ofs_state_R << std::endl;
	  			ofs_state_RA2 << std::endl;
	  			ofs_state_A << std::endl;

	  			ofs_state_I.close();
	  			ofs_state_R.close();
	  			ofs_state_RA2.close();
	  			ofs_state_A.close();
	  		}
	  	}

	   // To obtain the long-term histogram of the variable of interest (in this case I)
	   if (j==N-2 && WRITE_OUTPUT_H ) //This is the last point.
	   {
	      // Save the last time point histogram
	      // Create iterator for g1
	   	auto dom_hist = g1.getDomainIterator();
	   	size_t contad = 0;
	   	I_hist.resize(0);
	      //Iterate
	   	while (dom_hist.isNext())
	   	{
	        // Get the actual position from the iterator in the subdomain
	   		auto key = dom_hist.get();
	        // Add the value to I_hist
	   		I_hist.add(g1.template get<I>(key));

	        //Increment the counter
	   		contad++;
	        //Increment the iterator
	   		++dom_hist;
	   	}
	   	openfpm::vector<double> Lux_col;

	      // Saving i histogram
	   	v_cl.SGather(I_hist,Lux_col,0);


	   	if (v_cl.getProcessUnitID()==0 )
	   	{
	   		std::ofstream ofs_hist ( "I_hist.00" + std::to_string(j_hist), std::ofstream::out);

	   		ofs_hist << "cell, I" << std::endl;

	        //  ofs_var << counter << ", " << dR[n] << ", " << cR[m] << ", ";
	   		for (int jjjj = 0 ; jjjj < Lux_col.size() ; jjjj++)
	   		{
	   			size_t numb = jjjj+1;
	   			ofs_hist << numb << ", " << Lux_col.get(jjjj) << std::endl;
	   		}
	   		ofs_hist.close();
	   		j_hist++;

	   	}
	   }

	    //Get means and variances for each time from evolve_time
	   stats.get(j+1).get(0) = t;
	   for (int l = 1; l < 9; l++)
	   {
	   	stats.get(j+1).get(l) = stats_vec[l-1];
	   }

	   stats.get(j+1).get(9) = x5.last();

	    //Increment the time
	   t += T;
	}

	// Outputs of population statistics
	if (v_cl.getProcessUnitID()==0 )
	{
	    // Saving the population statistics though time in file media_temporal.txt
		if (WRITE_OUTPUT_T == 1)
		{
			std::ofstream ofs ("media_temporal.txt", std::ofstream::out);

		 // std::cout << "Saving to file... " <<  std::endl;
			ofs << "t, mx1, vx1, mx2, vx2, mx3, vx3, mx4, vx4, mx5, vx5" << std::endl;
	      //      0   1    2    3    4    5    6    7    8    9
			for (int k = 0; k < N; k++)
			{
				if (k%100 == 0)
				{
					for (int j = 0 ; j < 9 ; j++)
					{
						ofs << stats.get(k).get(j) << ", ";
					}
					ofs << stats.get(k).get(9) << std::endl;
				}
			}

			ofs.close();
	    } //End of temporal writing if



	    // Calculation of the long-term population statistics
	    // Calculating total variances
	    double means[5] {};
	    double mean_vars[5] {};
	    double var_means[5] {};
	    double varis[5] {};
	    double sdts[5] {};
	    // Taking only the last two-thirds of the signal (i.e. steady state)
	    size_t Ni = ceil(N*2/3);
	    for (int k = Ni; k < N; k++)
	    {
	      // Means of means and means of variances for x1 to x4
	    	for (int j = 0 ; j < 4 ; j++)
	    	{
	    		means[j] += stats.get(k).get(2*j+1)/(N-Ni);
	    		mean_vars[j] += stats.get(k).get(2*j+2)/(N-Ni);
	    	}
	      // Mean of means x5
	    	means[4] += stats.get(k).get(9)/(N-Ni);
	    	mean_vars[4] = 0.0;
	    }


	    for (int k = Ni; k < N; k++)
	    {
	      // Variance of means  for x1 to x5
	    	for (int j = 0 ; j < 5 ; j++)
	    	{
	    		var_means[j] += (stats.get(k).get(2*j+1)-means[j])*(stats.get(k).get(2*j+1)-means[j])/(N-Ni-1);
	    	}
	    }
	    // Total variance
	    for (int j = 0 ; j < 5 ; j++)
	    {
	    	varis[j] =  var_means[j] + mean_vars[j];
	    	sdts[j] = sqrt(varis[j]);
	    }

	    // Saving variances

	    std::ofstream ofs_var ("output.dat", std::ofstream::out);

	    ofs_var << "mx1, std1, mx2, std2, mx3, std3, mx4, std4, mx5, std5" << std::endl;
	    //      0   1    2    3    4    5    6    7    8    9
	    //  ofs_var << counter << ", " << dR[n] << ", " << cR[m] << ", ";
	    for (int j = 0 ; j < 4 ; j++)
	    {
	    	ofs_var << means[j] << ", " << sdts[j] << ", ";
	    }
	    ofs_var << means[4] << ", " << sdts[4] << std::endl;
	    ofs_var.close();
	} 

	  // Deinitialize the library
	openfpm_finalize();
}
