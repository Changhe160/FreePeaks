#ifndef __INDIVIDUAL_H_
#define __INDIVIDUAL_H_


#include "../../Utility/myVector.h"

class CMOEADInd{
public:
	CMOEADInd();
	virtual ~CMOEADInd();

	vector <double> x_var;
	vector <double> y_obj;

	void   rnd_init();
	void   obj_eval();

    bool   operator<(const CMOEADInd &ind2);
	bool   operator<<(const CMOEADInd &ind2);
    bool   operator==(const CMOEADInd &ind2);
    void   operator=(const CMOEADInd &ind2);

	void show_objective();
	void show_variable();

	int    rank;

};


class CSUB
{
public:
	CSUB();
	virtual ~CSUB();

	void show();

	CMOEADInd       indiv;     // best solution
	vector <int>    array;     // lattice point in a simplex
							   //vector <double> namda;     // weight vector
	MyVector namda;
	vector <int>    table;     // neighbourhood table

	double          density, fitness;

	void  operator=(const CSUB &sub2);
};


#endif