
#include "dmoeaclass.h"
#include "../../Global/global.h"
#include "../../Problem/ContinuousProblem.h"


CMOEADInd::CMOEADInd()
{
	for (int i = 0; i<Global::msp_global->mp_problem->getNumDim(); i++)
		x_var.push_back(0.0);
	for (int n = 0; n<Global::msp_global->mp_problem->getNumObj(); n++)
		y_obj.push_back(0.0);
	rank = 0;
}

CMOEADInd::~CMOEADInd()
{

}

void CMOEADInd::rnd_init()
{
	double lowBound, uppBound;

	for (int n = 0; n < Global::msp_global->mp_problem->getNumDim(); n++) {
		lowBound = CAST_PROBLEM_CONT->getSearchRange().getDomain(n).m_lower;
		uppBound = CAST_PROBLEM_CONT->getSearchRange().getDomain(n).m_upper;
		x_var[n] = lowBound + Global::msp_global->mp_uniformAlg->Next()*(uppBound - lowBound);
	}


}

void CMOEADInd::obj_eval()
{
	Solution<CodeVReal> s;
	s.data().m_x = x_var;
	s.evaluate();
	if (Global::msp_global->mp_problem->getOptType() == MIN_OPT)     y_obj = s.data().m_obj;
	else for (int n = 0; n<Global::msp_global->mp_problem->getNumObj(); n++) y_obj[n] = -s.data().m_obj[n];

}


void CMOEADInd::show_objective()
{
	for (int n = 0; n<Global::msp_global->mp_problem->getNumObj(); n++)
		printf("%f ", y_obj[n]);
	printf("\n");
}

void CMOEADInd::show_variable()
{
	for (int n = 0; n<Global::msp_global->mp_problem->getNumDim(); n++)
		printf("%f ", x_var[n]);
	printf("\n");
}

void CMOEADInd::operator=(const CMOEADInd &ind2)
{
	x_var = ind2.x_var;
	y_obj = ind2.y_obj;
	rank = ind2.rank;
}

bool CMOEADInd::operator<(const CMOEADInd &ind2)
{
	bool dominated = true;
	for (int n = 0; n<Global::msp_global->mp_problem->getNumObj(); n++)
	{
		if (ind2.y_obj[n]<y_obj[n]) return false;
	}
	if (ind2.y_obj == y_obj) return false;
	return dominated;
}


bool CMOEADInd::operator<<(const CMOEADInd &ind2)
{
	bool dominated = true;
	for (int n = 0; n<Global::msp_global->mp_problem->getNumObj(); n++)
	{
		if (ind2.y_obj[n]<y_obj[n] - 0.0001) return false;
	}
	if (ind2.y_obj == y_obj) return false;
	return dominated;
}

bool CMOEADInd::operator==(const CMOEADInd &ind2)
{
	if (ind2.y_obj == y_obj) return true;
	else return false;
}

CSUB::CSUB()
{
}

CSUB::~CSUB() {
}

void CSUB::show()
{
	for (int n = 0; n<namda.size(); n++) {
		printf("%f ", namda[n]);
	}
	printf("\n");
}

void CSUB::operator=(const CSUB &sub2) {
	indiv = sub2.indiv;
	array = sub2.array;
	table = sub2.table;
	namda = sub2.namda;
}
