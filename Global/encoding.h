#ifndef ENCONDING_H
#define ENCONDING_H

#include "../Utility/include.h"
#include "../Utility/TypeVar/typeVar.h"
#include "../Utility/TypeList/Typelist.h"

enum Encoding{Code_VReal=0,Code_VInt,Code_VBinary,Code_VVar};

struct VirtualEncoding{
	vector<double> m_obj;
	VirtualEncoding(size_t s=0):m_obj(s){}
	VirtualEncoding(const vector<double> &o) :m_obj(o){}
	virtual ~VirtualEncoding()=0;
	double getObjDistance(const vector<double>&rhs, DistanceMode mode = DIS_MANHATTAN)const{
		if (this->m_obj.size() != rhs.size()){
			throw myException("the number of objetives must be the same @ double VirtualEncoding::getObjDistance()");
		}
		double dis = 0;
		switch (mode){
			case DIS_EUCLIDEAN:
				for (size_t i = 0; i<this->m_obj.size(); ++i){
					dis += (this->m_obj[i] - rhs[i])*(this->m_obj[i] - rhs[i]);
				}
				dis=sqrt(dis);
				break;
			case DIS_MANHATTAN:
				for (size_t i = 0; i<this->m_obj.size(); ++i){
					dis += fabs(this->m_obj[i] - rhs[i]);
				}
				break;			
			case DIS_HAMMING:			
				for (size_t i = 0; i<this->m_obj.size(); ++i){
					if(this->m_obj[i] !=rhs[i]) dis+=1;
				}
				break;	
		}		
		return dis;
	}
};

inline VirtualEncoding::~VirtualEncoding(){ }

struct CodeVReal:public VirtualEncoding
{
	vector<double> m_x;
	CodeVReal() = default;
	CodeVReal(size_t s1,size_t s2):VirtualEncoding(s2),m_x(s1){}
	CodeVReal(int s1,int s2):VirtualEncoding(s2),m_x(s1){}
	CodeVReal(const vector<double> &x, const vector<double> &o) :VirtualEncoding(o), m_x(x){}
	CodeVReal(const vector<double> &x, const int s) :VirtualEncoding(s), m_x(x){}
	template<typename Itor>
	CodeVReal(Itor begin,Itor end):m_x(begin,end){}
	double operator[]( int i)const{
		return m_x[i];
	}
	double& operator[]( int i){
		return m_x[i];
	}
};

#endif