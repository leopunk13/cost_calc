///////////////////////////////////////////////////////
// File: CommonFunc
// Desc: Common function + Hearder files
//
// Author: Zhang Kang
// Date: 2013/09/06
///////////////////////////////////////////////////////
#pragma  once
#define DOUBLE_MAX 1e10
#define COMPUTE_RIGHT

#include<opencv2/opencv.hpp>
#include<string>
#include<iostream>
#include<bitset>
using namespace std;
using namespace cv;

//
// Opencv Lib 2.4.6
//
#ifdef _DEBUG
#pragma comment( lib, "opencv_calib3d300d.lib" )
#pragma comment( lib, "opencv_core300d.lib" )
#pragma comment( lib, "opencv_features2d300d.lib" )
#pragma comment( lib, "opencv_flann300d.lib" )
#pragma comment( lib, "opencv_highgui300d.lib" )
#pragma comment( lib, "opencv_imgproc300d.lib" )
#pragma comment( lib, "opencv_ml300d.lib" )
#pragma comment( lib, "opencv_objdetect300d.lib" )
#pragma comment( lib, "opencv_photo300d.lib" )
#pragma comment( lib, "opencv_stitching300d.lib" )
#pragma comment( lib, "opencv_superres300d.lib" )
#pragma comment( lib, "opencv_ts300d.lib" )
#pragma comment( lib, "opencv_video300d.lib" )
#pragma comment( lib, "opencv_videostab300d.lib" )
#else
#pragma comment( lib, "opencv_calib3d300.lib" )
#pragma comment( lib, "opencv_core300.lib" )
#pragma comment( lib, "opencv_features2d300.lib" )
#pragma comment( lib, "opencv_flann300.lib" )
#pragma comment( lib, "opencv_highgui300.lib" )
#pragma comment( lib, "opencv_imgproc300.lib" )
#pragma comment( lib, "opencv_ml300.lib" )
#pragma comment( lib, "opencv_objdetect300.lib" )
#pragma comment( lib, "opencv_photo300.lib" )
#pragma comment( lib, "opencv_stitching300.lib" )
#pragma comment( lib, "opencv_superres300.lib" )
#pragma comment( lib, "opencv_ts300.lib" )
#pragma comment( lib, "opencv_video300.lib" )
#pragma comment( lib, "opencv_videostab300.lib" )
#endif

// output matrix
template<class T>
void PrintMat( const Mat& mat )
{
	int rows = mat.rows;
	int cols = mat.cols;
	printf( "\n%d x %d Matrix\n", rows, cols );
	for( int r = 0; r < rows; r ++ ) {
		for( int c = 0; c < cols; c ++  ) {
			cout << mat.at<T>( r, c ) << "\t";
		}
		printf( "\n" );
	}
	printf( "\n" );
}
