#ifndef __EVOLUTION_H_
#define __EVOLUTION_H_

#include "../Algorithm.h"
#include "dmoeaclass.h"

class CMOEAD:public Algorithm
{

public:
	enum DecomFun { _TCHE1, _TCHE2, _NBI1, _NBI2, _NBI3 };
	CMOEAD(ParamMap &v);
	virtual ~CMOEAD();


	void init_uniformweight();               // initialize the weights for subproblems
	void init_neighbourhood();               // calculate the neighbourhood of each subproblem
	void init_population();                  // initialize the population


	void update_reference(CMOEADInd &ind);                 // update ideal point which is used in Tchebycheff or NBI method
	void update_problem(CMOEADInd &ind, int id, int type); // update current solutions in the neighbourhood

	void diffevolution();                                      // DE-based recombination
	void matingselection(vector<int> &list, int cid, int size, int type);  // select mating parents

	// execute MOEAD
	ReturnFlag run_();

    void calc_distance();

	void save_front(char savefilename[1024]);       // save the pareto front into files
	void save_ps(char savefilename[1024]);

	double fitnessfunction(vector <double> &y_obj, MyVector &namda);
	double dist_vector(vector <double> &vec1, vector <double> &vec2);

    vector <CSUB>       population;
	vector <vector<double>>  ps;
	vector <CSUB>       offspring;
	vector <int>        array;
	CMOEADInd           *ind_arr;

	double              distance;                   // generational distance from PF to solutions found
	int                 popsize;

	void operator=(const CMOEAD &moea);
	vector <double> idealpoint;
	int pops;
	int	   niche = 20,
		limit = 2,
		unit = 33;

	DecomFun m_decomFunction= _TCHE1;
	

	double realb = 0.9;     // probability of selecting mating parents from neighborhood
};

void minfastsort(vector<double> &x, vector<int> &idx, int n, int m);
#endif