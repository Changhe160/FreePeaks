#ifndef TYPES_ALG_PRO_H
#define TYPES_ALG_PRO_H
#include "../Utility/TypeList/Typelist.h"
#include "test.h"


typedef LOKI_TYPELIST_2(CMAES, CMOEAD) AlgList;

typedef LOKI_TYPELIST_5(FOnePeak, FFreePeak_OnePeak, FFreePeak_D_OnePeak, FFreePeak_M_OnePeak, FreePeak_D_M_OnePeak) ProList;



#endif

