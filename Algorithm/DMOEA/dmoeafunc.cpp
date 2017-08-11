#include "dmoeafunc.h"

#include "../../Measure/mMultiObj.h"
#include "recombination.h"

CMOEAD::CMOEAD(ParamMap &v) :Algorithm(-1, "MOEAD"),pops(v[param_popSize])
{

	ind_arr = new CMOEADInd[Global::msp_global->mp_problem->getNumObj()];

	// initialize ideal point
	for (int n = 0; n<Global::msp_global->mp_problem->getNumObj(); n++)
	{
		idealpoint.push_back(1.0e+30);
		ind_arr[n].rnd_init();
		ind_arr[n].obj_eval();
	}

	//load the pareto front
	Global::msp_global->mp_problem->getObjGlobalOpt(ps);
	if (Global::msp_global->mp_problem->getOptType() == MAX_OPT) {
		for (auto &i : ps) {
			for (auto &j : i) j *= -1;
		}
	}
}

CMOEAD::~CMOEAD()
{
	idealpoint.clear();
	delete[] ind_arr;
}


void CMOEAD::init_population()
{

	for (int i = 0; i<population.size(); i++)
	{
		population[i].indiv.rnd_init();
		population[i].indiv.obj_eval();
		update_reference(population[i].indiv);
		//nfes++;
	}
}

void CMOEAD::operator=(const CMOEAD &alg)
{

	population = alg.population;
	ps = alg.ps;
	ind_arr = alg.ind_arr;
	offspring = alg.offspring;
	distance = alg.distance;
	popsize = alg.popsize;
}


// createt the weight vectors with uniform distribution
void CMOEAD::init_uniformweight()
{
	if (Global::msp_global->mp_problem->getNumObj() == 2)
	{
		//vector<CMOEADInd> ws;
		//loadpfront("F6Weight500.dat",ws);
		//pops = 500;

		for (int n = 0; n<pops; n++)
		{
			CSUB sub;
			double a = 1.0*n / (pops - 1);
			sub.namda.push_back(a);
			sub.namda.push_back(1 - a);

			//load weight vectors from file
			//sub.namda.push_back(ws[n].y_obj[0]);
			//sub.namda.push_back(ws[n].y_obj[1]);

			population.push_back(sub);
		}
		popsize = pops;
	}
	else
	{
		for (int i = 0; i <= unit; i++)
		{
			for (int j = 0; j <= unit; j++)
			{
				if (i + j <= unit)
				{
					CSUB sub;
					sub.array.push_back(i);
					sub.array.push_back(j);
					sub.array.push_back(unit - i - j);
					for (int k = 0; k<sub.array.size(); k++)
						sub.namda.push_back(1.0*sub.array[k] / unit);
					population.push_back(sub);
				}
			}
		}

		popsize = population.size();
		pops = popsize;
	}
}

void CMOEAD::init_neighbourhood()
{

	vector<double> dis(population.size());
	vector<int> index(population.size());

	for (int i = 0; i<population.size(); i++)
	{
		// calculate the distances based on weight vectors
		for (int j = 0; j<population.size(); j++)
		{
			//x[j]    = dist_vector(population[i].namda,population[j].namda);
			dis[j] = population[i].namda.getDis(population[j].namda);
			index[j] = j;
		}

		// find 'niche' nearest neighboring subproblems
		minfastsort(dis, index,population.size(),niche);
		
		for (int k = 0; k<niche; k++)
		{
			population[i].table.push_back(index[k]);
		}

	}
}

void CMOEAD::update_problem(CMOEADInd &indiv, int id, int type)
{
	// indiv: child solution
	// id:   the id of current subproblem
	// type: update solutions in - neighborhood (1) or whole population (otherwise)
	int size, time = 0;
	if (type == 1)	size = population[id].table.size();
	else        size = population.size();
	//int *perm = new int[size];
	//random_permutation(perm, size);
	vector<int> perm(size);
	Global::msp_global->initializeRandomArray<vector<int> >(perm, size);
	for (int i = 0; i<size; i++)
	{
		int k;
		if (type == 1) k = population[id].table[perm[i]];
		else        k = perm[i];

		// calculate the values of objective function regarding the current subproblem
		double f1, f2;
		f1 = fitnessfunction(population[k].indiv.y_obj, population[k].namda);
		f2 = fitnessfunction(indiv.y_obj, population[k].namda);
		if (f2<f1)
		{
			population[k].indiv = indiv;
			time++;
		}
		// the maximal number of solutions updated is not allowed to exceed 'limit'
		if (time >= limit)
		{
			return;
		}
	}
	//delete[] perm;
}

void CMOEAD::update_reference(CMOEADInd &ind)
{
	//ind: child solution
	for (int n = 0; n<Global::msp_global->mp_problem->getNumObj(); n++)
	{
		if (ind.y_obj[n]<idealpoint[n])
		{
			idealpoint[n] = ind.y_obj[n];
			ind_arr[n] = ind;
		}
	}
}

void CMOEAD::matingselection(vector<int> &list, int cid, int size, int type) {
	// list : the set of the indexes of selected mating parents
	// cid  : the id of current subproblem
	// size : the number of selected mating parents
	// type : 1 - neighborhood; otherwise - whole population
	int ss = population[cid].table.size(), r, p;
	while (list.size()<size)
	{
		if (type == 1) {
			r = int(ss*Global::msp_global->mp_uniformAlg->Next());
			p = population[cid].table[r];
		}
		else
			p = int(population.size()*Global::msp_global->mp_uniformAlg->Next());

		bool flag = true;
		for (int i = 0; i<list.size(); i++)
		{
			if (list[i] == p) // p is in the list
			{
				flag = false;
				break;
			}
		}

		if (flag) list.push_back(p);
	}
}

void CMOEAD::diffevolution()
{
	pops = population.size();

	vector<int> perm(pops);
	Global::msp_global->initializeRandomArray<vector<int> >(perm, pops);

	for (int i = 0; i<pops; i++)
	{
		int n = perm[i];
		// or int n = i;
		int type;
		double rnd = Global::msp_global->mp_uniformAlg->Next();

		// mating selection based on probability
		if (rnd<realb)    type = 1;   // neighborhood
		else             type = 2;   // whole population

									 // select the indexes of mating parents
		vector<int> p;
		matingselection(p, n, 2, type);  // neighborhood selection

										 // produce a child solution
		CMOEADInd child;
		diff_evo_xover2(population[n].indiv, population[p[0]].indiv, population[p[1]].indiv, child);

		// apply polynomial mutation
		realmutation(child, 1.0 / Global::msp_global->mp_problem->getNumDim());

		// evaluate the child solution
		child.obj_eval();

		// update the reference points and other solutions in the neighborhood or the whole population
		update_reference(child);
		update_problem(child, n, type);

		p.clear(); 	//nfes++; 
	}

}


ReturnFlag CMOEAD::run_()
{


	vector<double> gd;
	char filename[1024];

	// initialization 
	int gen = 1;
	//nfes      = 0;
	init_uniformweight();
	init_neighbourhood();
	init_population();
	 calc_distance();

	//gd.push_back(0);  gd.push_back(distance);  // id and igd value

	vector<vector<double>*> pobj;
	for (int i = 0; i<popsize; i++)
		pobj.push_back(&population[i].indiv.y_obj);

#ifdef OFEC_CONSOLE
	if (mMultiObj::getMultiObj() && Global::msp_global->mp_problem->isGlobalOptKnown())
		mMultiObj::getMultiObj()->recordDistance(Global::msp_global.get(), Global::msp_global->m_runId, pobj);
#endif

	// evolution
	while (!ifTerminating())
	{

		diffevolution();

		//int dd = int(max_gen/25.0);
		int dd = Global::g_arg[param_sampleFre];
		// calculate igd-values
		if(gen%dd==0)
		{
		calc_distance();
		//cout<<"gen = "<<gen<<"  gd = "<<distance<<"  "<<endl;
		//gd.push_back(int(1.0*gen/dd)); gd.push_back(distance);
		#ifdef OFEC_CONSOLE
					if (mMultiObj::getMultiObj() && Global::msp_global->mp_problem->isGlobalOptKnown())
						mMultiObj::getMultiObj()->recordDistance(Global::msp_global.get(), Global::msp_global->m_runId, pobj);
		#endif
				}		
				gen++;
	}

#ifdef OFEC_CONSOLE
	if (mMultiObj::getMultiObj()) {
		mMultiObj::getMultiObj()->reInitialize(Global::msp_global.get(), popsize);
		for (int i = 0; i < popsize; i++) {
			if (Global::msp_global->mp_problem->getOptType() == MAX_OPT)
				for (int n = 0; n < Global::msp_global->mp_problem->getNumObj(); n++) population[i].indiv.y_obj[n] *= -1;
			mMultiObj::getMultiObj()->record(Global::msp_global.get(), i, population[i].indiv.y_obj, population[i].indiv.x_var);
		}
			
	}
#endif

	population.clear();
	ps.clear();
	std::cout<<" Outcome of the "<<Global::msp_global->m_runId<<"th run:  distance= "<<distance<<" nfes = "<<Global::msp_global->mp_problem->getEvaluations()<<endl;
	return Return_Terminate;
}


void CMOEAD::save_front(char saveFilename[1024])
{
	std::fstream fout;
	fout.open(saveFilename, std::ios::out);
	for (int n = 0; n<popsize; n++)
	{
		for (int k = 0; k<Global::msp_global->mp_problem->getNumObj(); k++)
			fout << population[n].indiv.y_obj[k] << "  ";
		fout << "\n";
	}
	fout.close();
}

void CMOEAD::save_ps(char saveFilename[1024])
{
	std::fstream fout;
	fout.open(saveFilename, std::ios::out);
	for (int n = 0; n<popsize; n++)
	{
		for (int k = 0; k<Global::msp_global->mp_problem->getNumDim(); k++)
			fout << population[n].indiv.x_var[k] << "  ";
		fout << "\n";
	}
	fout.close();
}


void CMOEAD::calc_distance()
{
	distance = 0;
	for (int i = 0; i<ps.size(); i++)
	{
		double min_d = 1.0e+10;
		for (int j = 0; j<population.size(); j++)
		{
			double d = dist_vector(ps[i], population[j].indiv.y_obj);
			if (d<min_d)  min_d = d;
		}
		distance += min_d;
	}
	distance /= ps.size();
}

double CMOEAD::fitnessfunction(vector <double> &y_obj, MyVector &namda)
{
	// Chebycheff Scalarizing Function
	double fitness = 0;
	int nobj = Global::msp_global->mp_problem->getNumObj();
	if (m_decomFunction == _TCHE1)
	{
		double max_fun = -1.0e+30;
		for (int n = 0; n<nobj; n++)
		{
			//double diff = fabs(y_obj[n] - idealpoint[n] + scale[n]);
			//double diff = fabs(y_obj[n] - idealpoint[n] + 0.05);
			double diff = fabs(y_obj[n] - idealpoint[n]);
			//double diff = fabs(y_obj[n] - 0);
			double feval;
			if (namda[n] == 0)
				feval = 0.0001*diff;
			else
				feval = diff*namda[n];
			if (feval>max_fun) max_fun = feval;

		}
		fitness = max_fun;
	}

	if (m_decomFunction == _TCHE2)
	{
		// reference point in the CHIM
		vector<double> scale(nobj, 1);
		double max_fun = -1.0e+30;
		for (int n = 0; n<nobj; n++)
		{
			double diff = (y_obj[n] - idealpoint[n]) / scale[n]; // error: scale not initialized in the original code
			double feval;
			if (namda[n] == 0)
				feval = 0.0001*diff;
			else
				feval = diff*namda[n];
			if (feval>max_fun) max_fun = feval;

		}
		fitness = max_fun;
	}



	// CHIM + Tchebycheff
	// CHIM is not available in 3 objectives
	if (m_decomFunction == _NBI1) {

		// quasi normal direction
		MyVector norm;
		for (int i = 0; i<nobj; i++)
		{
			norm.push_back(0.0);
			for (int j = 0; j<nobj; j++) {
				norm[i] += -ind_arr[j].y_obj[i];
			}
		}

		// normalization
		norm.normalize();


		// reference point in the CHIM
		vector <double> base;
		for (int i = 0; i<nobj; i++)
		{
			double tp2 = 0;
			for (int j = 0; j<nobj; j++)
				tp2 += ind_arr[j].y_obj[i] * namda[j];
			base.push_back(tp2);
		}

		// Tchebycheff function
		double max_fun = -1.0e+30;
		for (int n = 0; n<nobj; n++)
		{
			double diff = y_obj[n] - base[n];
			double feval = -diff*norm[n];
			if (feval>max_fun) max_fun = feval;

		}
		fitness = max_fun;
	}

	//* Boundary intersection approach
	//* reference point is chosen as the ideal point
	//* the direction is independent of CHIM
	if (m_decomFunction == _NBI2)
	{

		namda.normalize();

		// penalty method 
		// temporary vectors NBI method

		MyVector realA(nobj);
		MyVector realB(nobj);

		// difference beween current point and reference point
		for (int n = 0; n<nobj; n++)
			realA[n] = (y_obj[n] - idealpoint[n]);

		// distance along the search direction norm

		double d1 = fabs(realA*namda);
		//double d1 = prod_vector(realA,norm);

		// distance to the search direction norm
		for (int n = 0; n<nobj; n++)
			realB[n] = (y_obj[n] - (idealpoint[n] + d1*namda[n]));

		double d2 = realB.length();

		fitness = (d1 + 5 * d2);

		//t2 = clock();
		//total_sec+=(t2 - t1);
	}

	// NBI method
	if (m_decomFunction == _NBI3) {

		// quasi normal direction
		MyVector norm;
		for (int i = 0; i<nobj; i++)
		{
			norm.push_back(0.0);
			for (int j = 0; j<nobj; j++) {
				norm[i] += -ind_arr[j].y_obj[i];
			}
		}

		// normalization
		norm.normalize();


		// reference point in the CHIM
		vector <double> base;
		for (int i = 0; i<nobj; i++)
		{
			double tp2 = 0;
			for (int j = 0; j<nobj; j++)
				tp2 += ind_arr[j].y_obj[i] * namda[j];
			base.push_back(tp2);
		}

		// penalty method 
		// temporary vectors NBI method
		MyVector realA;
		MyVector realB;

		// difference beween current point and reference point
		for (int n = 0; n<nobj; n++)
			realA.push_back(y_obj[n] - base[n]);

		// distance along the search direction norm
		double d1 = realA*norm;

		// distance to the search direction norm
		for (int n = 0; n<nobj; n++)
			realB.push_back(y_obj[n] - (base[n] + d1*norm[n]));
		double d2 = realB.length();

		fitness = -d1 + 2 * d2;
	}


	return fitness;
}

double CMOEAD::dist_vector(vector <double> &vec1, vector <double> &vec2)
{
	int dim = vec1.size();
	double sum = 0;
	for (int n = 0; n<dim; n++)
		sum += (vec1[n] - vec2[n])*(vec1[n] - vec2[n]);
	return sqrt(sum);
}

void minfastsort(vector<double> &x, vector<int> &idx, int n, int m)
{
	for (int i = 0; i<m; i++)
	{
		for (int j = i + 1; j<n; j++)
			if (x[i]>x[j])
			{
				double temp = x[i];
				x[i] = x[j];
				x[j] = temp;
				int id = idx[i];
				idx[i] = idx[j];
				idx[j] = id;
			}
	}
}