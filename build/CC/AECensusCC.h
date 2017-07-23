#pragma once
#include "../CommFunc.h"
#include "../CCMethod.h"

#define AECENCUS_WND 9
#define AECENCUS_BIT 8
#define WEIGHT_M 0.2
#define WEIGHT_GM 0.1
#define WEIGHT_RM 0.9
#define W1 0.25
#define t1 0.35
#define t2 0.65
#define PI 3.14159265
class AECencusCC :
	public CCMethod
{
public:
	AECencusCC(void)
	{
		printf("\n\t\tAECencus for Cost Computation");
	}
	~AECencusCC(void){}
public:
	void buildCV(const Mat& lImg, const Mat& rImg, const int maxDis, Mat* costVol);
#ifdef COMPUTE_RIGHT
	void buildRightCV(const Mat& lImg, const Mat& rImg, const int maxDis, Mat* rCostVol);
#endif
};
