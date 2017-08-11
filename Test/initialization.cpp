#include "initialization.h"
#include "../Global/global.h"
#include "types.h"


template<typename List, typename Map, typename Base, int n >
struct Rigister {
	static void regist(const vector<Info> & info, Map &map,  STRING2ID &msi) {
		for (int i = 0; i<info.size(); i++)
		{
			if (info[i].ClassIndex == n)
			{
				Global::registerItem(msi, info[i].Name);
				map.insert(make_pair(info[i].Name, &createObject<Base, typename Loki::TL::TypeAt<List, n>::Result>));
			}
		}
		Rigister<List, Map,Base,n - 1>::regist(info,map,msi);
	}
};

template<typename List, typename Map, typename Base>
struct Rigister<List,Map,Base,-1> {
	static void regist(const vector<Info> & info, Map &map,  STRING2ID &msi) { }
};

void registerProNdAlg(){
	
	//register problems and algorithms here
	vector<Info> proInfo, algInfo,termInfo;
	LoadData(proInfo,algInfo, termInfo);
	
	
	Rigister<ProList, classFactory::ClassMapProblem,Problem, Loki::TL::Length<ProList>::value - 1>::regist(proInfo, Global::ms_classFactory.m_theMapProblem,Global::msm_pro);
	Rigister<AlgList, classFactory::ClassMapAlgorithm, Algorithm, Loki::TL::Length<AlgList>::value - 1>::regist(algInfo, Global::ms_classFactory.m_theMapAlgorithm,Global::msm_alg);
	Rigister<TermList, classFactory::ClassMapTermination, Termination, Loki::TL::Length<TermList>::value - 1>::regist(termInfo, Global::ms_classFactory.m_theMapTermination,Global::msm_term);


	//bind algorithms and problems here
	for(int i=0;i<algInfo.size();i++)
	{
		for(int j=0;j<proInfo.size();j++)
		{
			for (set<string>::iterator iter = proInfo[j].Type.begin();;)
			{
				if (!algInfo[i].Type.count(*iter))
					break;
				if (++iter == proInfo[j].Type.end())
				{
					Global::registerAlg4Pro(make_pair(algInfo[i].Name, proInfo[j].Name));
					break;
				}
			}			
		}
	}
	// check if the algorithm and the problem match
	if (Global::msm_alg4pro.find(make_pair(Global::g_arg[param_algName], Global::g_arg[param_proName])) == Global::msm_alg4pro.end()) {
		throw myException("The algorithm is not for solving the problem");
	}

}

void setAddtParameters(){
	ParamMap &v=Global::g_arg;
	v[param_algId]=Global::msm_alg[(v[param_algName])]; 
	v[param_proId]=Global::msm_pro[(v[param_proName])];

	if(v.find(param_solutionValidationMode)==v.end())	v[param_solutionValidationMode]=VALIDATION_SETTOBOUND;
	if(v.find(param_populationInitialMethod)==v.end())	v[param_populationInitialMethod]=0;
	
	/*if (gGetProblemName(Global::msm_pro[Global::g_arg[param_proName]]).find("_FreePeak_D_") != string::npos){
			#ifdef OFEC_CONSOLE
			v[param_maxEvals] = v[param_numChange] * v[param_changeFre];
			#endif
	}*/
	if (Global::g_arg.find(param_dataFile1) == Global::g_arg.end() ||
		Global::g_arg.find(param_dataFile2) == Global::g_arg.end() ||
		Global::g_arg.find(param_dataFile3) == Global::g_arg.end()) {
		FreePeak::setup(v[param_proName]);
	}
	
	v[param_evalCountFlag]=true;	 

}

void LoadData(vector<Info> & proInfo,vector<Info> & algInfo, vector<Info> &termInfo)
{
	ifstream infile;
	ostringstream os;
	os<<Global::g_arg[param_workingDir]<<"Test/data/problem.txt";
	infile.open(os.str());
	if(!infile)		throw myException("load problem file fail in the @LoadData()");
	string line;
	while(getline(infile,line))
	{	
		string::size_type pos=string::npos;
		if((pos=line.find('\r'))&&pos!=string::npos){  line.erase(pos,1);}
		if(line.compare("#begin")==0)	break;
	}
	
	//load problem data
	string word;
	while(getline(infile,line))
	{
		Info temp;
		istringstream stream(line);
		stream>>temp.Name;
		stream>>temp.ClassIndex;
		while(stream>>word)
			temp.Type.insert(word);
		proInfo.push_back(move(temp));
	}
	infile.close();
	infile.clear();
	os.str("");

	os<<Global::g_arg[param_workingDir]<<"Test/data/algorithm.txt";
	infile.open(os.str());
	if(!infile) throw myException("load algorithm file fail in the @LoadData()");
	while(getline(infile,line))
	{	
		string::size_type pos=string::npos;
		if((pos=line.find('\r'))&&pos!=string::npos){  line.erase(pos,1);}
		if(line.compare("#begin")==0)	break;
	}
	
	//load algorithm data
	while(getline(infile,line))
	{
		Info temp;
		istringstream stream(line);
		stream>>temp.Name; 
		stream>>temp.ClassIndex;
		while(stream>>word)
			temp.Type.insert(word);
		algInfo.push_back(move(temp));
	}
	infile.close();
	infile.clear();
	os.str("");

	os << Global::g_arg[param_workingDir] << "Test/data/termination.txt";
	infile.open(os.str());
	if (!infile) throw myException("load algorithm file fail in the @LoadData()");
	while (getline(infile, line))
	{
		string::size_type pos = string::npos;
		if ((pos = line.find('\r')) && pos != string::npos) { line.erase(pos, 1); }
		if (line.compare("#begin") == 0)	break;
	}

	//load termination data
	while (getline(infile, line))
	{
		Info temp;
		istringstream stream(line);
		stream >> temp.Name;
		stream >> temp.ClassIndex;		
		termInfo.push_back(move(temp));
	}
	infile.close();
	infile.clear();
}
void setGlobalParameters(int argn,char *argv[]){
	//format of argument list from main: param_proName=FUN_Sphere param_numDim=5  param_numPeak=10 
	Global::registerParamter(); //register all system parameters first
	string letter;
	for(auto i='a';i<='z';i++) {
		letter+=i; 
		letter+=i-32;
	}
	string remove="\r\t\n\b\v";
	letter+="\\/:";
	for(int i=1;i<argn;i++){
		string par=argv[i];
		while(size_t pos=par.find_first_of(remove)){
			if(string::npos==pos) break;
			par.erase(par.begin()+pos);
		}
		size_t pos=par.find('=');
		if(pos==string::npos)  throw myException("Invalid argument:gSetGlobalParameters()");
		string value=par.substr(pos+1,par.size()-1),name=par.substr(0,pos);
		name.insert(0,"param_");
		if(Global::msm_param.find(name)==Global::msm_param.end()) throw myException("Invalid argument:gSetGlobalParameters()");
		if(value.compare("true")==0){
			Global::g_arg[Global::msm_param[name]]=true;
		}else if(value.compare("false")==0){
			Global::g_arg[Global::msm_param[name]]=false;
		}else if(value.find_first_of(letter)!=string::npos){
			if(value.size()==1) 
				Global::g_arg[Global::msm_param[name]]=value[0];
			else	
				Global::g_arg[Global::msm_param[name]]= value;
		}else if(string::npos!=par.find('.',pos+1)){
			double val=atof(value.c_str());
			Global::g_arg[Global::msm_param[name]]=val;
		}else{
			int val=atoi(value.c_str());
			Global::g_arg[Global::msm_param[name]]=val;
		}
	}
#ifdef OFEC_CONSOLE
	if(Global::g_arg.find(param_sampleFre)==Global::g_arg.end()) Global::g_arg[param_sampleFre]=2; 
	if(Global::g_arg.find(param_workingDir)==Global::g_arg.end())  Global::g_arg[param_workingDir]=string("./");
#endif
	Global::ms_curAlgId=Global::msm_alg[(Global::g_arg[param_algName])];
	Global::ms_curProId=Global::msm_pro[(Global::g_arg[param_proName])];

}

