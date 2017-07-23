#pragma once
#include "../CommFunc.h"
#include "../PPMethod.h"

#define MED_SZ 19
#define SIG_CLR 0.1
#define SIG_DIS 9
#define MaxDisparity 256
#define Tc 8
#define Ts 10
//
// Weight-Median Post-processing
//
class WMPP :
	public PPMethod
{
public:
	WMPP(void) 
	{
		printf( "\n\t\tWeight-Median Post-processing" );
	}
	~WMPP(void) {}
public:
	void postProcess( const Mat& lImg, const Mat& rImg, const int maxDis, const int disSc,
		Mat& lDis, Mat& rDis, Mat& lSeg, Mat& lChk);
};
