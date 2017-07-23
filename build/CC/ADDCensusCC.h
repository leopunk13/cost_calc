#pragma once
#include "../CommFunc.h"
#include "../CCMethod.h"

#define CENCUS_WND 5
#define CENCUS_BIT 80
#define W1 0.5
//#define W1 1.5
//#define W2 8
#define t1 0.35
#define t2 0.65

class ADDCensusCC :
	public CCMethod
{
public:
	ADDCensusCC(void)
	{
		printf("\n\t\tADDCensus for Cost Computation");
	}
	~ADDCensusCC(void){}
public:
	void buildCV(const Mat& lImg, const Mat& rImg, const int maxDis, Mat* costVol);
#ifdef COMPUTE_RIGHT
	void buildRightCV(const Mat& lImg, const Mat& rImg, const int maxDis, Mat* rCostVol);
#endif
};
