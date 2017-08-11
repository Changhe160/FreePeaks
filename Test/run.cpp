#include "run.h"
#include "test.h"

Run::Run(int runId,ParamMap &v):m_runId(runId){
	if(!Global::msp_global.get()){
		#ifdef OFEC_CONSOLE
		Global::msp_global.reset(new Global(m_runId,1./7,(m_runId+1.)/((int)(MAX_NUM_RUN)+1.)));
		#endif
		
		Global::msp_global->mp_problem.reset(Global::ms_classFactory.constructProblem((v[param_proName]))(v));
		
		#ifdef OFEC_CONSOLE
		if(Global::msp_global->mp_problem->isProTag(SOP)){
			if(Global::msp_global->mp_problem->isProTag(CONT)){
				if(CAST_PROBLEM_CONT->getGOpt().flagGloObj()) Global::g_arg[param_gOptFlag]=true;
				else  Global::g_arg[param_gOptFlag]=false;
			}

			if(mSingleObj::getSingleObj()==nullptr){
				if(Global::msp_global->mp_problem->isProTag(DOP))
					mSingleObjDyn::initialize(Global::g_arg); 
				else
					mSingleObj::initialize(Global::g_arg); 
				mSingleObj::getSingleObj()->setFileName(Global::g_arg);
			}
			
			vector<double> gOpt;
			if(mSingleObj::getSingleObj() &&Global::msp_global->mp_problem->getObjGlobalOpt(gOpt)){
				mSingleObj::getSingleObj()->addGOpt(m_runId,gOpt[0]);
			}
			
			if(m_runId==0){
					mSingleObj::getSingleObj()->setAccuracy(Global::msp_global->mp_problem->getAccuracy());
					mSingleObj::getSingleObj()->setCompareType(Global::msp_global->mp_problem->getOptType());
					mSingleObj::getSingleObj()->setProParameter(Global::msp_global->mp_problem->m_proPar);
					
			}
			
		}
		if(Global::msp_global->mp_problem->isProTag(MMP)){
			if(mMultiModal::getPopInfor()==nullptr){
					mMultiModal::initialize(Global::g_arg);
					mMultiModal::getPopInfor()->setFileName(Global::g_arg);	
			}
		}
		if(Global::msp_global->mp_problem->isProTag(MOP)){
			if(mMultiObj::getMultiObj()==nullptr){	
				if(Global::msp_global->mp_problem->isProTag(DOP))
					mMultiObj::initialize(Global::msp_global.get(),Global::g_arg);
				else mMultiObj::initialize(Global::msp_global.get(),Global::g_arg); 
				mMultiObj::getMultiObj()->setFileName(Global::g_arg);
			}
		}
		
		#endif

		#ifndef OFEC_PROBLEM_DEBUG
			Global::msp_global->mp_algorithm.reset(Global::ms_classFactory.constructAlgorithm((v[param_algName]))(v));
		#endif
		#ifdef OFEC_CONSOLE
			if (mSingleObj::getSingleObj() && Global::msp_global->mp_algorithm)	mSingleObj::getSingleObj()->setAlgParameter(Global::msp_global->mp_algorithm->m_algPar);
		#endif
	}

}
Run::~Run(){
	Global::msp_global.reset();	
}
ReturnFlag Run::go(){
	return Global::msp_global->mp_algorithm->run();
}
