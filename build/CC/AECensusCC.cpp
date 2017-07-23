#include "AECensusCC.h"

inline double cg_myCostGrd(double* lC, double* rC)
{
	double clrDiff = 0;
	// three color
	for (int c = 0; c < 3; c++) {
		double temp = fabs(lC[c] - rC[c]);
		clrDiff += temp;
	}
	
	return clrDiff;
}

inline Mat calcBDIP(Mat img, int blockSize)
{
	int nr = img.rows;
	int nc = img.cols;
	int radius = blockSize / 2;
	Mat imgBlueChannels, imgGreenChannels, imgRedChannels;
	Mat BDIP_RGB;
	vector<Mat> splitchannels;
	vector<Mat> mergechannels;
	img.convertTo(img, CV_32F);
	cvtColor(img, img, CV_RGB2BGR);
	img.convertTo(img, CV_32F, 255);

	split(img, splitchannels);
	imgBlueChannels = splitchannels.at(0);
	imgGreenChannels = splitchannels.at(1);
	imgRedChannels = splitchannels.at(2);
	Mat BDIP_R = imgRedChannels.clone();
	Mat BDIP_G = imgGreenChannels.clone();
	Mat BDIP_B = imgBlueChannels.clone();
	
	for (int y = 0; y < nr; y++)
	{
		for (int x = 0; x < nc; x++)
		{
			float maxR = 0.0f, maxG = 0.0f, maxB = 0.0f;
			float sumR = 0.0f, sumG = 0.0f, sumB = 0.0f;

			for (int m = y - radius; m <= y + radius; m++)
			{
				for (int n = x - radius; n <= x + radius; n++)
				{
					if (m < 0 || n < 0 || n >= nc || m >= nr)continue;
					sumR += imgRedChannels.ptr<float>(m)[n];
					sumG += imgGreenChannels.ptr<float>(m)[n];
					sumB += imgBlueChannels.ptr<float>(m)[n];
					if (maxR < imgRedChannels.ptr<float>(m)[n])
						maxR = imgRedChannels.ptr<float>(m)[n];
					if (maxG < imgGreenChannels.ptr<float>(m)[n])
						maxG = imgGreenChannels.ptr<float>(m)[n];
					if (maxB < imgBlueChannels.ptr<float>(m)[n])
						maxB = imgBlueChannels.ptr<float>(m)[n];
				}
			}
			if (maxR == 0)
				maxR = 1;
			if (maxG == 0)
				maxG = 1;
			if (maxB == 0)
				maxB = 1;
			BDIP_R.ptr<float>(y)[x] = blockSize * blockSize - (sumR / maxR);
			BDIP_G.ptr<float>(y)[x] = blockSize * blockSize - (sumG / maxG);
			BDIP_B.ptr<float>(y)[x] = blockSize * blockSize - (sumB / maxB);
		}
	}
	mergechannels.push_back(BDIP_B);
	mergechannels.push_back(BDIP_G);
	mergechannels.push_back(BDIP_R);
	merge(mergechannels, BDIP_RGB);
	return BDIP_RGB;
}

inline Mat calcDE(Mat img, int windowSize)
{
	int nr = img.rows;
	int nc = img.cols;
	int radius = windowSize / 2;
	Mat gimg;
	Mat tmp;
	img.convertTo(tmp, CV_32F);
	cvtColor(tmp, gimg, CV_RGB2GRAY);
	gimg.convertTo(gimg, CV_8U, 255);
	Mat DE = gimg.clone();
	double temp1 = 0;
	int temp2 = 0;
	const int hist_sz = 256;//0��255��һ��256���Ҷ�ֵ    
	double hist[hist_sz];
	memset(hist, 0, sizeof(hist));
	for (int i = 0; i < 256; i++)
	{
		hist[i] = 0.0f;
	}

	for (int y = 0; y < nr; y++)
	{
		for (int x = 0; x < nc; x++)
		{
			double sum_DE = 0.0f;
			for (int i = 0; i < 256; i++)
			{
				hist[i] = 0.0f;
			}

			for (int m = y - radius; m <= y + radius; m++)
			{
				for (int n = x - radius; n <= x + radius; n++)
				{
					if (m < 0 || n < 0 || n >= nc || m >= nr)continue;
					temp1 = gimg.ptr<uchar>(m)[n];
					temp2 = int(temp1);//������ת��  
					hist[temp2]++; //����ʵ����hist�д洢���Ҷ�ֵ���ֵĴ���
				}
			}
			for (int i = 0; i < 256; i++)
			{
				if (hist[i] != 0)
					sum_DE += -hist[i] * (log(hist[i]) / log(2));
			}
			//printf("%d,%d,		 %f\n",x,y, -sum_DE);
			DE.ptr<uchar>(y)[x] = static_cast<uchar>(sum_DE);
		}
	}
	return DE;
}

//inline Mat MaskFilter(Mat img, int winsize, float gamma_r, float gamma_s)
//{
//	int radius = winsize / 2;
//	int nr = img.rows;
//	int nc = img.cols;
//	Mat result = img.clone();
//	//for (int it = 0; it < 3;it++)
//	//{
//	for (int y = 0; y < nr; y++)
//	{
//
//		for (int x = 0; x < nc; x++)
//		{
//			int upY = min(nr, y + radius);
//			int upX = min(nc, x + radius);
//			int downY = max(1, y - radius);
//			int downX = max(1, x - radius);
//			int medianVal = -1;
//			int c = 0;
//			int wc = 0;
//			double rangediff = 0.0f;
//			double spatialdiff = 0.0f;
//			double  weight = 0.0f;
//			double weightsum = 0.0f;
//			double weightsum0 = 0.0f;
//			double weightsum1 = 0.0f;
//			double weightsum2 = 0.0f;
//			//double *weightedhist = new double[255];
//			for (int i = downY; i < upY; i++)
//			{
//				int maxDispVal = 0;
//				for (int j = downX; j < upX; j++)
//				{
//					/*int dispVals = dispimgPtr[j];
//					if (maxDispVal < dispVals)
//					maxDispVal = dispVals;*/
//					if (i< 0 || j < 0 || i>nr || j>nc)continue;
//					//uchar * datadisp = img.ptr<uchar>(i);
//					//compute weight
//					//rangediff = sqrt(static_cast<double>(pow(leftimg.ptr<uchar>(y)[x]-leftimg.ptr<uchar>(i)[j], 2)));
//					rangediff = sqrt(static_cast<double>(pow((img.at<Vec3d>(y, x)[0] - img.at<Vec3d>(i, j)[0]), 2) + pow((img.at<Vec3d>(y, x)[1] - img.at<Vec3d>(i, j)[1]), 2) + pow((img.at<Vec3d>(y, x)[2] - img.at<Vec3d>(i, j)[2]), 2)));
//					spatialdiff = static_cast<double>(sqrt(pow(x - i, 2) + pow(y - j, 2)));
//					weight += exp(-(pow(rangediff, 2) / (2 * pow(gamma_r, 2))) - (pow(spatialdiff, 2) / (2 * pow(gamma_s, 2))));
//					c++;
//				}
//			}
//			for (int i = downY; i < upY; i++)
//			{
//				for (int j = downX; j < upX; j++)
//				{
//					if (i< 0 || j < 0 || i>nr || j>nc)continue;
//					double * datadisp = img.ptr<double>(i);
//					//rangediff = sqrt(static_cast<double>(pow(leftimg.ptr<uchar>(y)[x]-leftimg.ptr<uchar>(i)[j], 2)));
//					rangediff = sqrt(static_cast<double>(pow((img.at<Vec3d>(y, x)[0] - img.at<Vec3d>(i, j)[0]), 2) + pow((img.at<Vec3d>(y, x)[1] - img.at<Vec3d>(i, j)[1]), 2) + pow((img.at<Vec3d>(y, x)[2] - img.at<Vec3d>(i, j)[2]), 2)));
//					spatialdiff = static_cast<double>(sqrt(pow(x - i, 2) + pow(y - j, 2)));
//					if (img.channels() == 1)
//					{
//						weightsum += (1 / weight)*exp(-(pow(rangediff, 2) / (2 * pow(gamma_r, 2))) - (pow(spatialdiff, 2) / (2 * pow(gamma_s, 2))))*datadisp[j];
//					}
//					if (img.channels() == 3)
//					{
//						weightsum0 += (1 / weight)*exp(-(pow(rangediff, 2) / (2 * pow(gamma_r, 2))) - (pow(spatialdiff, 2) / (2 * pow(gamma_s, 2))))*img.at<Vec3d>(i, j)[0];
//						weightsum1 += (1 / weight)*exp(-(pow(rangediff, 2) / (2 * pow(gamma_r, 2))) - (pow(spatialdiff, 2) / (2 * pow(gamma_s, 2))))*img.at<Vec3d>(i, j)[1];
//						weightsum2 += (1 / weight)*exp(-(pow(rangediff, 2) / (2 * pow(gamma_r, 2))) - (pow(spatialdiff, 2) / (2 * pow(gamma_s, 2))))*img.at<Vec3d>(i, j)[2];
//
//					}
//
//				}
//			}
//			if (result.channels() == 1)
//				result.ptr<double>(y)[x] = weightsum;
//			if (result.channels() == 3)
//			{
//				result.at<Vec3d>(y, x)[0] = weightsum0;
//				result.at<Vec3d>(y, x)[1] = weightsum1;
//				result.at<Vec3d>(y, x)[2] = weightsum2;
//			}
//		}
//
//	}
//
//	//}
//	return result;
//}

void AECencusCC::buildCV(const Mat& lImg, const Mat& rImg, const int maxDis, Mat* costVol)
{
	CV_Assert(lImg.type() == CV_64FC3 && rImg.type() == CV_64FC3);

	int hei = lImg.rows;
	int wid = lImg.cols;
	Mat lGray, rGray;

	/*Mat leftLab, rightLab;
	Mat left_img = imread("im0.png", -1);
	Mat right_img = imread("im1.png", -1);
	cvtColor(left_img, leftLab, CV_BGR2Lab);
	cvtColor(right_img, rightLab, CV_BGR2Lab);*/

	Mat ltmp = lImg.clone();
	Mat leftLab, rightLab;
	ltmp.convertTo(ltmp, CV_32F, 255);
	cvtColor(ltmp, leftLab, CV_RGB2Lab);
	Mat rtmp = rImg.clone();
	rtmp.convertTo(rtmp, CV_32F, 255);
	cvtColor(rtmp, rightLab, CV_RGB2Lab);
	//Mat GM0 = MaskFilter(lImg, 9, 60, 50);
	//Mat GM1 = MaskFilter(rImg, 9, 60, 50);

	Mat lImgBlueChannels, rImgBlueChannels;// GM0BlueChannels, GM1BlueChannels;
	Mat lImgGreenChannels, rImgGreenChannels;// GM0GreenChannels, GM1GreenChannels;
	Mat lImgRedChannels, rImgRedChannels;// GM0RedChannels, GM1RedChannels;

	Mat lImgRGBGradx = lImg.clone();
	Mat rImgRGBGradx = rImg.clone();
	Mat lImgRGBGrady = lImg.clone();
	Mat rImgRGBGrady = rImg.clone();
	//Mat GM0RGBGradx = GM0.clone();
	//Mat GM0RGBGrady = GM0.clone();
	//Mat GM1RGBGradx = GM1.clone();
	//Mat GM1RGBGrady = GM1.clone();
	vector<Mat> splitLImgChannels, splitRImgChannels;// splitGM0channels, splitGM1channels;
	split(lImg, splitLImgChannels);
	split(rImg, splitRImgChannels);
	//split(GM0, splitGM0channels);
	//split(GM1, splitGM1channels);
	lImgBlueChannels = splitLImgChannels.at(0);
	lImgGreenChannels = splitLImgChannels.at(1);
	lImgRedChannels = splitLImgChannels.at(2);
	rImgBlueChannels = splitRImgChannels.at(0);
	rImgGreenChannels = splitRImgChannels.at(1);
	rImgRedChannels = splitRImgChannels.at(2);
	//GM0BlueChannels = splitGM0channels.at(0);
	//GM0GreenChannels = splitGM0channels.at(1);
	//GM0RedChannels = splitGM0channels.at(2);
	//GM1BlueChannels = splitGM1channels.at(0);
	//GM1GreenChannels = splitGM1channels.at(1);
	//GM1RedChannels = splitGM1channels.at(2);
	/*double horizontal_k[3][3] = { { -1, 0, 1 }, { -2, 0, 2 }, { -1, 0, 1 } };
	double vertical_k[3][3] = { { -1, -2, -1 }, { 0, 0, 0 }, { 1, 2, 1 } };
	Mat maskx = Mat(3, 3, CV_32FC1, horizontal_k);
	Mat masky = Mat(3, 3, CV_32FC1, vertical_k);*/
	Mat gradLImgBluex, gradRImgBluex;// gradGM0Bluex, gradGM1Bluex;
	gradLImgBluex = lImgBlueChannels.clone();
	gradRImgBluex = rImgBlueChannels.clone();
	//gradGM0Bluex = GM0BlueChannels.clone();
	//gradGM1Bluex = GM1BlueChannels.clone();
	Mat gradLImgGreenx, gradRImgGreenx;// gradGM0Greenx, gradGM1Greenx;
	gradLImgGreenx = lImgGreenChannels.clone();
	gradRImgGreenx = rImgGreenChannels.clone();
	//gradGM0Greenx = GM0GreenChannels.clone();
	//gradGM1Greenx = GM1GreenChannels.clone();
	Mat gradLImgRedx, gradRImgRedx;// gradGM0Redx, gradGM1Redx;
	gradLImgRedx = lImgRedChannels.clone();
	gradRImgRedx = rImgRedChannels.clone();
	//gradGM0Redx = GM0RedChannels.clone();
	//gradGM1Redx = GM1RedChannels.clone();
	Mat gradLImgBluey, gradRImgBluey;// gradGM0Bluey, gradGM1Bluey;
	gradLImgBluey = lImgBlueChannels.clone();
	gradRImgBluey = rImgBlueChannels.clone();
	//gradGM0Bluey = GM0BlueChannels.clone();
	//gradGM1Bluey = GM1BlueChannels.clone();
	Mat gradLImgGreeny, gradRImgGreeny;// gradGM0Greeny, gradGM1Greeny;
	gradLImgGreeny = lImgGreenChannels.clone();
	gradRImgGreeny = rImgGreenChannels.clone();
	//gradGM0Greeny = GM0GreenChannels.clone();
	//gradGM1Greeny = GM1GreenChannels.clone();
	Mat gradLImgRedy, gradRImgRedy;// gradGM0Redy, gradGM1Redy;
	gradLImgRedy = lImgRedChannels.clone();
	gradRImgRedy = rImgRedChannels.clone();
	//gradGM0Redy = GM0RedChannels.clone();
	//gradGM1Redy = GM1RedChannels.clone();
	

	for (int y = 0; y < hei; y++)
	{
		for (int x = 0; x < wid; x++)
		{
			if (x - 1 < 0 || x + 1 >= wid || y - 1 < 0 || y + 1 >= hei)
			{
				gradLImgBluex.ptr<double>(y)[x] = 0;
				gradLImgGreenx.ptr<double>(y)[x] = 0;
				gradLImgRedx.ptr<double>(y)[x] = 0;
				gradRImgBluey.ptr<double>(y)[x] = 0;
				gradRImgGreeny.ptr<double>(y)[x] = 0;
				gradRImgRedy.ptr<double>(y)[x] = 0;
				gradLImgBluex.ptr<double>(y)[x] = 0;
				gradLImgGreenx.ptr<double>(y)[x] = 0;
				gradLImgRedx.ptr<double>(y)[x] = 0;
				gradRImgBluey.ptr<double>(y)[x] = 0;
				gradRImgGreeny.ptr<double>(y)[x] = 0;
				gradRImgRedy.ptr<double>(y)[x] = 0;
				//gradGM0Bluex.ptr<double>(y)[x] = 0;
				//gradGM0Greenx.ptr<double>(y)[x] = 0;
				//gradGM0Redx.ptr<double>(y)[x] = 0;
				//gradGM0Bluey.ptr<double>(y)[x] = 0;
				//gradGM0Greeny.ptr<double>(y)[x] = 0;
				//gradGM0Redy.ptr<double>(y)[x] = 0;
				//gradGM1Bluex.ptr<double>(y)[x] = 0;
				//gradGM1Greenx.ptr<double>(y)[x] = 0;
				//gradGM1Redx.ptr<double>(y)[x] = 0;
				//gradGM1Bluey.ptr<double>(y)[x] = 0;
				//gradGM1Greeny.ptr<double>(y)[x] = 0;
				//gradGM1Redy.ptr<double>(y)[x] = 0;
			}
			else{
				gradLImgBluex.ptr<double>(y)[x] = lImgBlueChannels.ptr<double>(y)[x-1] - lImgBlueChannels.ptr<double>(y)[x + 1];
				gradLImgGreenx.ptr<double>(y)[x] = lImgGreenChannels.ptr<double>(y)[x - 1] - lImgGreenChannels.ptr<double>(y)[x + 1];
				gradLImgRedx.ptr<double>(y)[x] = lImgRedChannels.ptr<double>(y)[x - 1] - lImgRedChannels.ptr<double>(y)[x + 1];
				gradLImgBluey.ptr<double>(y)[x] = lImgBlueChannels.ptr<double>(y - 1)[x] - lImgBlueChannels.ptr<double>(y + 1)[x];
				gradLImgGreeny.ptr<double>(y)[x] = lImgGreenChannels.ptr<double>(y - 1)[x] - lImgGreenChannels.ptr<double>(y + 1)[x];
				gradLImgRedy.ptr<double>(y)[x] = lImgRedChannels.ptr<double>(y - 1)[x] - lImgRedChannels.ptr<double>(y + 1)[x];
				gradRImgBluex.ptr<double>(y)[x] = rImgBlueChannels.ptr<double>(y)[x-1] - rImgBlueChannels.ptr<double>(y)[x + 1];
				gradRImgGreenx.ptr<double>(y)[x] = rImgGreenChannels.ptr<double>(y)[x-1] - rImgGreenChannels.ptr<double>(y)[x + 1];
				gradRImgRedx.ptr<double>(y)[x] = rImgRedChannels.ptr<double>(y)[x-1] - rImgRedChannels.ptr<double>(y)[x + 1];
				gradRImgBluey.ptr<double>(y)[x] = rImgBlueChannels.ptr<double>(y - 1)[x] - rImgBlueChannels.ptr<double>(y + 1)[x];
				gradRImgGreeny.ptr<double>(y)[x] = rImgGreenChannels.ptr<double>(y - 1)[x] - rImgGreenChannels.ptr<double>(y + 1)[x];
				gradRImgRedy.ptr<double>(y)[x] = rImgRedChannels.ptr<double>(y - 1)[x] - rImgRedChannels.ptr<double>(y + 1)[x];
			}

		}
	}

	Mat BDIP0 = calcBDIP(lImg, 7);
	Mat BDIP1 = calcBDIP(rImg, 7);
	Mat DE0 = calcDE(lImg, 7);
	Mat DE1 = calcDE(rImg, 7);
	
	vector<Mat> splitBDIP0channels, splitBDIP1channels;
	Mat BDIP0BlueChannels, BDIP0GreenChannels, BDIP0RedChannels;
	Mat BDIP1BlueChannels, BDIP1GreenChannels, BDIP1RedChannels;
	split(BDIP0, splitBDIP0channels);
	split(BDIP1, splitBDIP1channels);
	BDIP0BlueChannels = splitBDIP0channels.at(0);
	BDIP0GreenChannels = splitBDIP0channels.at(1);
	BDIP0RedChannels = splitBDIP0channels.at(2);
	BDIP1BlueChannels = splitBDIP1channels.at(0);
	BDIP1GreenChannels = splitBDIP1channels.at(1);
	BDIP1RedChannels = splitBDIP1channels.at(2);

	//for (int d = 0; d < maxDis; d++) {
	//	for (int y = 0; y < hei; y++)
	//	{
	//		
	//		double* lData = (double*)lImg.ptr<double>(y);
	//		double* rData = (double*)rImg.ptr<double>(y);
	//		double* cost = (double*)costVol[d].ptr<double>(y);
	//		for (int x = 0; x < wid; x++)
	//		{
	//			if (x - d >= 0)
	//			{
	//				double cAD = fabs(static_cast<double>(leftLab.at<Vec3f>(y, x)[0] - rightLab.at<Vec3f>(y, x - d)[0])) + fabs(static_cast<double>(leftLab.at<Vec3f>(y, x)[1] - rightLab.at<Vec3f>(y, x - d)[1])) + fabs(static_cast<double>(leftLab.at<Vec3f>(y, x)[2] - rightLab.at<Vec3f>(y, x - d)[2]));
	//				//double cADgx = 1 / 3 * (fabs(gradLImgBluex.ptr<double>(y)[x] - gradRImgBluex.ptr<double>(y)[x - d]) + fabs(gradLImgGreenx.ptr<double>(y)[x] - gradRImgGreenx.ptr<double>(y)[x - d])
	//				//	+ fabs(gradLImgRedx.ptr<double>(y)[x] - gradRImgRedx.ptr<double>(y)[x - d]));
	//				/*+ fabs(static_cast<double>(w2*gradGM0Bluex.ptr<uchar>(m)[n + d] - w2*gradGM1Bluex.ptr<uchar>(m)[n]))
	//				+ fabs(static_cast<double>(w2*gradGM0Greenx.ptr<uchar>(m)[n + d] - w2*gradGM1Greenx.ptr<uchar>(m)[n])) + fabs(static_cast<double>(w2*gradGM0Redx.ptr<uchar>(m)[n + d] - w2*gradGM1Redx.ptr<uchar>(m)[n]));*/
	//				//double cADgy = 1 / 3 * (fabs(gradLImgBluey.ptr<double>(y)[x] - gradRImgBluey.ptr<double>(y)[x - d]) + fabs(gradLImgGreeny.ptr<double>(y)[x] - gradRImgGreeny.ptr<double>(y)[x - d])
	//				//	+ fabs(gradLImgRedy.ptr<double>(y)[x] - gradRImgRedy.ptr<double>(y)[x - d]));
	//				/*+ fabs(static_cast<double>(w2*gradGM0Bluey.ptr<uchar>(m)[n + d] - w2*gradGM1Bluey.ptr<uchar>(m)[n]))
	//				+ fabs(static_cast<double>(w2*gradGM0Greeny.ptr<uchar>(m)[n + d] - w2*gradGM1Greeny.ptr<uchar>(m)[n])) + fabs(static_cast<double>(w2*gradGM0Redy.ptr<uchar>(m)[n + d] - w2*gradGM1Redy.ptr<uchar>(m)[n]));*/
	//				//cost[x] = cAD * WEIGHT_AD;// +cADgx * WEIGHT_GRDX + cADgy * WEIGHT_GRDY;
	//			}
	//		}
	//	}
	//}
	/*
	Mat tmp;
	lImg.convertTo(tmp, CV_32F);
	cvtColor(tmp, lGray, CV_RGB2GRAY);
	lGray.convertTo(lGray, CV_8U, 255);
	rImg.convertTo(tmp, CV_32F);
	cvtColor(tmp, rGray, CV_RGB2GRAY);
	rGray.convertTo(rGray, CV_8U, 255);

	int w0[3][3] = { { -1, 0, 1 }, { -2, 0, 2 }, { -1, 0, 1 } };//�Ƕ�Ϊ0�ȵ�Ȩֵ
	int w45[3][3] = { { 2, 1, 0 }, { 1, 0, -1 }, { 0, -1, -2 } };//�Ƕ�Ϊ45��
	int w90[3][3] = { { -1, -2, -1 }, { 0, 0, 0 }, { 1, 2, 1 } };//�Ƕ�Ϊ90��
	int w135[3][3] = { { 0, 1, 2 }, { -1, 0, 1 }, { -2, -1, 0 } };//�Ƕ�Ϊ135��
	int *ML = new int[8];
	int *MR = new int[8];
	int *WDL = new int[8];
	int *WDR = new int[8];
	int *DL = new int[8];
	int *DR = new int[8];

	bitset<AECENCUS_BIT>* lCode = new bitset<AECENCUS_BIT>[wid * hei];
	bitset<AECENCUS_BIT>* rCode = new bitset<AECENCUS_BIT>[wid * hei];
	bitset<AECENCUS_BIT>* wLCode0 = new bitset<AECENCUS_BIT>[wid * hei];
	bitset<AECENCUS_BIT>* wRCode0 = new bitset<AECENCUS_BIT>[wid * hei];
	bitset<AECENCUS_BIT>* wLCode45 = new bitset<AECENCUS_BIT>[wid * hei];
	bitset<AECENCUS_BIT>* wRCode45 = new bitset<AECENCUS_BIT>[wid * hei];
	bitset<AECENCUS_BIT>* wLCode90 = new bitset<AECENCUS_BIT>[wid * hei];
	bitset<AECENCUS_BIT>* wRCode90 = new bitset<AECENCUS_BIT>[wid * hei];
	bitset<AECENCUS_BIT>* wLCode135 = new bitset<AECENCUS_BIT>[wid * hei];
	bitset<AECENCUS_BIT>* wRCode135 = new bitset<AECENCUS_BIT>[wid * hei];
	bitset<AECENCUS_BIT>* pLCode = lCode;
	bitset<AECENCUS_BIT>* pRCode = rCode;
	bitset<AECENCUS_BIT>* pwLCode0 = wLCode0;
	bitset<AECENCUS_BIT>* pwRCode0 = wRCode0;
	bitset<AECENCUS_BIT>* pwLCode45 = wLCode45;
	bitset<AECENCUS_BIT>* pwRCode45 = wRCode45;
	bitset<AECENCUS_BIT>* pwLCode90 = wLCode90;
	bitset<AECENCUS_BIT>* pwRCode90 = wRCode90;
	bitset<AECENCUS_BIT>* pwLCode135 = wLCode135;
	bitset<AECENCUS_BIT>* pwRCode135 = wRCode135;
	for (int y = 0; y < hei; y++)
	{
		uchar* pLData = (uchar*)(lGray.ptr<uchar>(y));
		uchar* pRData = (uchar*)(rGray.ptr<uchar>(y));
		for (int x = 0; x < wid; x++)
		{
			
			ML[0] = 1 / 8 * (lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x]
				+ lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y)[x]);
			ML[1] = 1 / 8 * (lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y)[x]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0]);
			ML[2] = 1 / 8 * (lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ lGray.ptr<uchar>(y)[x]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]);
			ML[3] = 1 / 8 * (lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x]
				+ lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ lGray.ptr<uchar>(y)[x]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x]);
			ML[4] = 1 / 8 * (lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ lGray.ptr<uchar>(y)[x]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]);
			ML[5] = 1 / 8 * (lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y)[x]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 2 >= 0 ? x - 2 : (x - 2 >= 0 ? wid - 1 : wid - 2)]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x]);
			ML[6] = 1 / 8 * (lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y)[x]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0]);
			ML[7] = 1 / 8 * (lGray.ptr<uchar>(y)[x]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]);

			MR[0] = 1 / 8 * (rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x]
				+ rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y)[x]);
			MR[1] = 1 / 8 * (rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y)[x]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0]);
			MR[2] = 1 / 8 * (rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ rGray.ptr<uchar>(y)[x]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]);
			MR[3] = 1 / 8 * (rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x]
				+ rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ rGray.ptr<uchar>(y)[x]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x]);
			MR[4] = 1 / 8 * (rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ rGray.ptr<uchar>(y)[x]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]);
			MR[5] = 1 / 8 * (rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y)[x]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 2 >= 0 ? x - 2 : (x - 2 >= 0 ? wid - 1 : wid - 2)]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x]);
			MR[6] = 1 / 8 * (rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y)[x]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0]);
			MR[7] = 1 / 8 * (rGray.ptr<uchar>(y)[x]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]);
			for (int i = 0; i < 8; i++)
			{
				(*pLCode)[i] = (pLData[x] > ML[i]);
				(*pRCode)[i] = (pRData[x] > MR[i]);
			}
			pLCode++;
			pRCode++;
		}
	}
		
	for (int y = 0; y < hei; y++)
	{
		for (int x = 0; x < wid; x++)
		{
			WDL[0] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w0[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w0[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w0[1][2]
				+ lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[2][0]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : hei - 1] * w0[2][1]
				+ lGray.ptr<uchar>(y)[x] * w0[2][2]);
			WDL[1] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w0[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w0[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w0[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w0[1][2]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[2][0]
				+ lGray.ptr<uchar>(y)[x] * w0[2][1]
				+ lGray.ptr<uchar>(y)[x + 1 <wid ? x + 1 : 0] * w0[2][2]);
			WDL[2] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w0[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w0[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w0[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[1][2]
				+ lGray.ptr<uchar>(y)[x] * w0[2][0]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w0[2][1]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[2][2]);
			WDL[3] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[0][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[0][1]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w0[0][2]
				+ lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[1][0]
				+ lGray.ptr<uchar>(y)[x] * w0[1][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[2][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 2] * w0[2][1]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w0[2][2]);
			WDL[4] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w0[0][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w0[0][1]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[0][2]
				+ lGray.ptr<uchar>(y)[x] * w0[1][0]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[1][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w0[2][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w0[2][1]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[2][2]);
			WDL[5] = abs(lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[0][0]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[0][1]
				+ lGray.ptr<uchar>(y)[x] * w0[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w0[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w0[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w0[2][2]);
			WDL[6] = abs(lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[0][0]
				+ lGray.ptr<uchar>(y)[x] * w0[0][1]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w0[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w0[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w0[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w0[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w0[2][2]);
			WDL[7] = abs(lGray.ptr<uchar>(y)[x] * w0[0][0]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w0[0][1]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w0[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w0[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w0[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[2][2]);

			DL[0] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[1] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] - lGray.ptr<uchar>(y)[x]);
			DL[2] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);
			DL[3] = abs(lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[4] = abs(lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);
			DL[5] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[6] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] - lGray.ptr<uchar>(y)[x]);
			DL[7] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);

			WDR[0] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w0[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w0[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w0[1][2]
				+ rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[2][0]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : hei - 1] * w0[2][1]
				+ rGray.ptr<uchar>(y)[x] * w0[2][2]);
			WDR[1] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w0[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w0[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w0[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w0[1][2]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[2][0]
				+ rGray.ptr<uchar>(y)[x] * w0[2][1]
				+ rGray.ptr<uchar>(y)[x + 1 <wid ? x + 1 : 0] * w0[2][2]);
			WDR[2] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w0[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w0[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w0[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[1][2]
				+ rGray.ptr<uchar>(y)[x] * w0[2][0]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w0[2][1]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[2][2]);
			WDR[3] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[0][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[0][1]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w0[0][2]
				+ rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[1][0]
				+ rGray.ptr<uchar>(y)[x] * w0[1][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[2][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 2] * w0[2][1]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w0[2][2]);
			WDR[4] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w0[0][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w0[0][1]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[0][2]
				+ rGray.ptr<uchar>(y)[x] * w0[1][0]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[1][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w0[2][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w0[2][1]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[2][2]);
			WDR[5] = abs(rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[0][0]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[0][1]
				+ rGray.ptr<uchar>(y)[x] * w0[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w0[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w0[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w0[2][2]);
			WDR[6] = abs(rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[0][0]
				+ rGray.ptr<uchar>(y)[x] * w0[0][1]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w0[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w0[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w0[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w0[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w0[2][2]);
			WDR[7] = abs(rGray.ptr<uchar>(y)[x] * w0[0][0]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w0[0][1]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w0[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w0[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w0[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[2][2]);

			DR[0] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[1] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] - rGray.ptr<uchar>(y)[x]);
			DR[2] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);
			DR[3] = abs(rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[4] = abs(rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);
			DR[5] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[6] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] - rGray.ptr<uchar>(y)[x]);
			DR[7] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);

			for (int i = 0; i < 8; i++)
			{
				(*pwLCode0)[i] = (DL[i] > WDL[i]);
				(*pwRCode0)[i] = (DR[i] > WDR[i]);
			}
			pwLCode0++;
			pwRCode0++;
		}
	}
	
	for (int y = 0; y < hei; y++)
	{
		for (int x = 0; x < wid; x++)
		{
			WDL[0] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w45[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w45[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w45[1][2]
				+ lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[2][0]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : hei - 1] * w45[2][1]
				+ lGray.ptr<uchar>(y)[x] * w45[2][2]);
			WDL[1] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w45[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w45[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w45[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w45[1][2]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[2][0]
				+ lGray.ptr<uchar>(y)[x] * w45[2][1]
				+ lGray.ptr<uchar>(y)[x + 1 <wid ? x + 1 : 0] * w45[2][2]);
			WDL[2] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w45[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w45[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w45[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[1][2]
				+ lGray.ptr<uchar>(y)[x] * w45[2][0]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w45[2][1]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[2][2]);
			WDL[3] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[0][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[0][1]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w45[0][2]
				+ lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[1][0]
				+ lGray.ptr<uchar>(y)[x] * w45[1][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[2][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 2] * w45[2][1]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w45[2][2]);
			WDL[4] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w45[0][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w45[0][1]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[0][2]
				+ lGray.ptr<uchar>(y)[x] * w45[1][0]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[1][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w45[2][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w45[2][1]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[2][2]);
			WDL[5] = abs(lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[0][0]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[0][1]
				+ lGray.ptr<uchar>(y)[x] * w45[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w45[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w45[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w45[2][2]);
			WDL[6] = abs(lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[0][0]
				+ lGray.ptr<uchar>(y)[x] * w45[0][1]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w45[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w45[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w45[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w45[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w45[2][2]);
			WDL[7] = abs(lGray.ptr<uchar>(y)[x] * w45[0][0]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w45[0][1]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w45[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w45[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w45[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[2][2]);

			DL[0] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[1] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] - lGray.ptr<uchar>(y)[x]);
			DL[2] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);
			DL[3] = abs(lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[4] = abs(lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);
			DL[5] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[6] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] - lGray.ptr<uchar>(y)[x]);
			DL[7] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);

			WDR[0] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w45[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w45[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w45[1][2]
				+ rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[2][0]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : hei - 1] * w45[2][1]
				+ rGray.ptr<uchar>(y)[x] * w45[2][2]);
			WDR[1] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w45[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w45[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w45[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w45[1][2]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[2][0]
				+ rGray.ptr<uchar>(y)[x] * w45[2][1]
				+ rGray.ptr<uchar>(y)[x + 1 <wid ? x + 1 : 0] * w45[2][2]);
			WDR[2] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w45[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w45[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w45[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[1][2]
				+ rGray.ptr<uchar>(y)[x] * w45[2][0]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w45[2][1]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[2][2]);
			WDR[3] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[0][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[0][1]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w45[0][2]
				+ rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[1][0]
				+ rGray.ptr<uchar>(y)[x] * w45[1][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[2][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 2] * w45[2][1]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w45[2][2]);
			WDR[4] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w45[0][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w45[0][1]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[0][2]
				+ rGray.ptr<uchar>(y)[x] * w45[1][0]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[1][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w45[2][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w45[2][1]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[2][2]);
			WDR[5] = abs(rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[0][0]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[0][1]
				+ rGray.ptr<uchar>(y)[x] * w45[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w45[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w45[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w45[2][2]);
			WDR[6] = abs(rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[0][0]
				+ rGray.ptr<uchar>(y)[x] * w45[0][1]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w45[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w45[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w45[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w45[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w45[2][2]);
			WDR[7] = abs(rGray.ptr<uchar>(y)[x] * w45[0][0]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w45[0][1]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w45[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w45[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w45[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[2][2]);

			DR[0] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[1] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] - rGray.ptr<uchar>(y)[x]);
			DR[2] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);
			DR[3] = abs(rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[4] = abs(rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);
			DR[5] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[6] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] - rGray.ptr<uchar>(y)[x]);
			DR[7] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);

			for (int i = 0; i < 8; i++)
			{
				(*pwLCode45)[i] = (DL[i] > WDL[i]);
				(*pwRCode45)[i] = (DR[i] > WDR[i]);
			}
			pwLCode45++;
			pwRCode45++;
		}
	}

	for (int y = 0; y < hei; y++)
	{
		for (int x = 0; x < wid; x++)
		{
			WDL[0] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w90[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w90[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w90[1][2]
				+ lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[2][0]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : hei - 1] * w90[2][1]
				+ lGray.ptr<uchar>(y)[x] * w90[2][2]);
			WDL[1] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w90[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w90[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w90[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w90[1][2]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[2][0]
				+ lGray.ptr<uchar>(y)[x] * w90[2][1]
				+ lGray.ptr<uchar>(y)[x + 1 <wid ? x + 1 : 0] * w90[2][2]);
			WDL[2] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w90[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w90[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w90[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[1][2]
				+ lGray.ptr<uchar>(y)[x] * w90[2][0]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w90[2][1]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[2][2]);
			WDL[3] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[0][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[0][1]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w90[0][2]
				+ lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[1][0]
				+ lGray.ptr<uchar>(y)[x] * w90[1][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[2][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 2] * w90[2][1]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w90[2][2]);
			WDL[4] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w90[0][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w90[0][1]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[0][2]
				+ lGray.ptr<uchar>(y)[x] * w90[1][0]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[1][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w90[2][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w90[2][1]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[2][2]);
			WDL[5] = abs(lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[0][0]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[0][1]
				+ lGray.ptr<uchar>(y)[x] * w90[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w90[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w90[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w90[2][2]);
			WDL[6] = abs(lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[0][0]
				+ lGray.ptr<uchar>(y)[x] * w90[0][1]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w90[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w90[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w90[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w90[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w90[2][2]);
			WDL[7] = abs(lGray.ptr<uchar>(y)[x] * w90[0][0]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w90[0][1]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w90[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w90[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w90[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[2][2]);

			DL[0] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[1] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] - lGray.ptr<uchar>(y)[x]);
			DL[2] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);
			DL[3] = abs(lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[4] = abs(lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);
			DL[5] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[6] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] - lGray.ptr<uchar>(y)[x]);
			DL[7] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);

			WDR[0] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w90[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w90[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w90[1][2]
				+ rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[2][0]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : hei - 1] * w90[2][1]
				+ rGray.ptr<uchar>(y)[x] * w90[2][2]);
			WDR[1] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w90[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w90[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w90[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w90[1][2]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[2][0]
				+ rGray.ptr<uchar>(y)[x] * w90[2][1]
				+ rGray.ptr<uchar>(y)[x + 1 <wid ? x + 1 : 0] * w90[2][2]);
			WDR[2] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w90[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w90[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w90[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[1][2]
				+ rGray.ptr<uchar>(y)[x] * w90[2][0]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w90[2][1]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[2][2]);
			WDR[3] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[0][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[0][1]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w90[0][2]
				+ rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[1][0]
				+ rGray.ptr<uchar>(y)[x] * w90[1][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[2][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 2] * w90[2][1]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w90[2][2]);
			WDR[4] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w90[0][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w90[0][1]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[0][2]
				+ rGray.ptr<uchar>(y)[x] * w90[1][0]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[1][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w90[2][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w90[2][1]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[2][2]);
			WDR[5] = abs(rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[0][0]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[0][1]
				+ rGray.ptr<uchar>(y)[x] * w90[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w90[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w90[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w90[2][2]);
			WDR[6] = abs(rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[0][0]
				+ rGray.ptr<uchar>(y)[x] * w90[0][1]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w90[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w90[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w90[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w90[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w90[2][2]);
			WDR[7] = abs(rGray.ptr<uchar>(y)[x] * w90[0][0]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w90[0][1]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w90[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w90[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w90[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[2][2]);

			DR[0] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[1] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] - rGray.ptr<uchar>(y)[x]);
			DR[2] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);
			DR[3] = abs(rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[4] = abs(rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);
			DR[5] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[6] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] - rGray.ptr<uchar>(y)[x]);
			DR[7] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);

			for (int i = 0; i < 8; i++)
			{
				(*pwLCode90)[i] = (DL[i] > WDL[i]);
				(*pwRCode90)[i] = (DR[i] > WDR[i]);
			}
			pwLCode90++;
			pwRCode90++;
		}
	}

	for (int y = 0; y < hei; y++)
	{
		for (int x = 0; x < wid; x++)
		{
			WDL[0] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w135[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w135[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w135[1][2]
				+ lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[2][0]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : hei - 1] * w135[2][1]
				+ lGray.ptr<uchar>(y)[x] * w135[2][2]);
			WDL[1] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w135[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w135[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w135[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w135[1][2]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[2][0]
				+ lGray.ptr<uchar>(y)[x] * w135[2][1]
				+ lGray.ptr<uchar>(y)[x + 1 <wid ? x + 1 : 0] * w135[2][2]);
			WDL[2] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w135[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w135[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w135[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[1][2]
				+ lGray.ptr<uchar>(y)[x] * w135[2][0]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w135[2][1]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[2][2]);
			WDL[3] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[0][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[0][1]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w135[0][2]
				+ lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[1][0]
				+ lGray.ptr<uchar>(y)[x] * w135[1][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[2][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 2] * w135[2][1]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w135[2][2]);
			WDL[4] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w135[0][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w135[0][1]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[0][2]
				+ lGray.ptr<uchar>(y)[x] * w135[1][0]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[1][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w135[2][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w135[2][1]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[2][2]);
			WDL[5] = abs(lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[0][0]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[0][1]
				+ lGray.ptr<uchar>(y)[x] * w135[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w135[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w135[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w135[2][2]);
			WDL[6] = abs(lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[0][0]
				+ lGray.ptr<uchar>(y)[x] * w135[0][1]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w135[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w135[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w135[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w135[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w135[2][2]);
			WDL[7] = abs(lGray.ptr<uchar>(y)[x] * w135[0][0]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w135[0][1]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w135[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w135[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w135[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[2][2]);

			DL[0] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[1] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] - lGray.ptr<uchar>(y)[x]);
			DL[2] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);
			DL[3] = abs(lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[4] = abs(lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);
			DL[5] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[6] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] - lGray.ptr<uchar>(y)[x]);
			DL[7] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);

			WDR[0] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w135[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w135[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w135[1][2]
				+ rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[2][0]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : hei - 1] * w135[2][1]
				+ rGray.ptr<uchar>(y)[x] * w135[2][2]);
			WDR[1] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w135[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w135[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w135[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w135[1][2]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[2][0]
				+ rGray.ptr<uchar>(y)[x] * w135[2][1]
				+ rGray.ptr<uchar>(y)[x + 1 <wid ? x + 1 : 0] * w135[2][2]);
			WDR[2] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w135[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w135[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w135[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[1][2]
				+ rGray.ptr<uchar>(y)[x] * w135[2][0]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w135[2][1]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[2][2]);
			WDR[3] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[0][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[0][1]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w135[0][2]
				+ rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[1][0]
				+ rGray.ptr<uchar>(y)[x] * w135[1][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[2][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 2] * w135[2][1]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w135[2][2]);
			WDR[4] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w135[0][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w135[0][1]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[0][2]
				+ rGray.ptr<uchar>(y)[x] * w135[1][0]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[1][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w135[2][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w135[2][1]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[2][2]);
			WDR[5] = abs(rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[0][0]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[0][1]
				+ rGray.ptr<uchar>(y)[x] * w135[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w135[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w135[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w135[2][2]);
			WDR[6] = abs(rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[0][0]
				+ rGray.ptr<uchar>(y)[x] * w135[0][1]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w135[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w135[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w135[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w135[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w135[2][2]);
			WDR[7] = abs(rGray.ptr<uchar>(y)[x] * w135[0][0]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w135[0][1]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w135[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w135[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w135[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[2][2]);

			DR[0] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[1] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] - rGray.ptr<uchar>(y)[x]);
			DR[2] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);
			DR[3] = abs(rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[4] = abs(rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);
			DR[5] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[6] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] - rGray.ptr<uchar>(y)[x]);
			DR[7] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);

			for (int i = 0; i < 8; i++)
			{
				(*pwLCode135)[i] = (DL[i] > WDL[i]);
				(*pwRCode135)[i] = (DR[i] > WDR[i]);
			}
			pwLCode135++;
			pwRCode135++;
		}
	}

	// build cost volume
	bitset<AECENCUS_BIT> lB, wlB0, wlB45, wlB90, wlB135;
	bitset<AECENCUS_BIT> rB, wrB0, wrB45, wrB90, wrB135;
	pLCode = lCode;
	pwLCode0 = wLCode0;
	pwLCode45 = wLCode45;
	pwLCode90 = wLCode90;
	pwLCode135 = wLCode135;
	for (int y = 0; y < hei; y++) {
		int index = y * wid;
		for (int x = 0; x < wid; x++) {
			lB = *pLCode;
			wlB0 = *pwLCode0;
			wlB45 = *pwLCode45;
			wlB90 = *pwLCode90;
			wlB135 = *pwLCode135;
			for (int d = 0; d < maxDis; d++) {
				double* cost = (double*)costVol[d].ptr<double>(y);
				double tmpCost1 = AECENCUS_BIT, tmpCost2 = AECENCUS_BIT;
				if (x - d >= 0) {
					rB = rCode[index + x - d];
					wrB0 = wRCode0[index + x - d];
					wrB45 = wRCode45[index + x - d];
					wrB90 = wRCode90[index + x - d];
					wrB135 = wRCode135[index + x - d];
					tmpCost1 = (lB ^ rB).count();
					tmpCost2 = (wlB0 ^ wrB0).count() + (wlB45 ^ wrB45).count() + (wlB90 ^ wrB90).count() + (wlB135 ^ wrB135).count();
				}
				//cost[x] = WEIGHT_M *min(tmpCost1,7.2) + (1 - WEIGHT_M) * min(tmpCost2,28.8);
			}
			pLCode++;
			pwLCode0++;
			pwLCode45++;
			pwLCode90++;
			pwLCode135++;
		}
	}
	delete[] lCode;
	delete[] rCode;
	delete[] wLCode0;
	delete[] wRCode0;
	delete[] wLCode45;
	delete[] wRCode45;
	delete[] wLCode90;
	delete[] wRCode90;
	delete[] wLCode135;
	delete[] wRCode135;
*/
	for (int y = 0; y < hei; y++)
	{
		for (int x = 0; x < wid; x++)
		{
			for (int d = 0; d < maxDis; d++)
			{
				double* cost = (double*)costVol[d].ptr<double>(y);
				if (x - d >= 0) {
					double w_pq = 0.0f;
					double w_p_q_ = 0.0f;
					double wq_DE1 = 0.0f;
					double wq_DE2 = 0.0f;
					double w_q_DE1 = 0.0f;
					double w_q_DE2 = 0.0f;
					double maxq_DE = 0.0f;
					if (y - 1 >= 0 && y + 1 < hei)
					{
						if (x - 1 >= 0)
						{
							if (maxq_DE < DE0.ptr<uchar>(y - 1)[x - 1])
								maxq_DE = DE0.ptr<uchar>(y - 1)[x - 1];
							if (maxq_DE < DE0.ptr<uchar>(y)[x - 1])
								maxq_DE = DE0.ptr<uchar>(y)[x - 1];
							if (maxq_DE < DE0.ptr<uchar>(y + 1)[x - 1])
								maxq_DE = DE0.ptr<uchar>(y + 1)[x - 1];
						}
						if (maxq_DE < DE0.ptr<uchar>(y - 1)[x])
							maxq_DE = DE0.ptr<uchar>(y - 1)[x];
						if (maxq_DE < DE0.ptr<uchar>(y)[x])
							maxq_DE = DE0.ptr<uchar>(y)[x];
						if (maxq_DE < DE0.ptr<uchar>(y + 1)[x])
							maxq_DE = DE0.ptr<uchar>(y + 1)[x];
						if (x + 1 < wid)
						{
							if (maxq_DE < DE0.ptr<uchar>(y - 1)[x + 1])
								maxq_DE = DE0.ptr<uchar>(y - 1)[x + 1];
							if (maxq_DE < DE0.ptr<uchar>(y)[x + 1])
								maxq_DE = DE0.ptr<uchar>(y)[x + 1];
							if (maxq_DE < DE0.ptr<uchar>(y + 1)[x + 1])
								maxq_DE = DE0.ptr<uchar>(y + 1)[x + 1];
						}
					}
					else if (y - 1 < 0)
					{
						if (x - 1 >= 0)
						{
							if (maxq_DE < DE0.ptr<uchar>(y)[x - 1])
								maxq_DE = DE0.ptr<uchar>(y)[x - 1];
							if (maxq_DE < DE0.ptr<uchar>(y + 1)[x - 1])
								maxq_DE = DE0.ptr<uchar>(y + 1)[x - 1];
						}
						if (maxq_DE < DE0.ptr<uchar>(y)[x])
							maxq_DE = DE0.ptr<uchar>(y)[x];
						if (maxq_DE < DE0.ptr<uchar>(y + 1)[x])
							maxq_DE = DE0.ptr<uchar>(y + 1)[x];
						if (x + 1 < wid)
						{
							if (maxq_DE < DE0.ptr<uchar>(y)[x + 1])
								maxq_DE = DE0.ptr<uchar>(y)[x + 1];
							if (maxq_DE < DE0.ptr<uchar>(y + 1)[x + 1])
								maxq_DE = DE0.ptr<uchar>(y + 1)[x + 1];
						}
					}
					else
					{
						if (x - 1 >= 0)
						{
							if (maxq_DE < DE0.ptr<uchar>(y - 1)[x - 1])
								maxq_DE = DE0.ptr<uchar>(y - 1)[x - 1];
							if (maxq_DE < DE0.ptr<uchar>(y)[x - 1])
								maxq_DE = DE0.ptr<uchar>(y)[x - 1];
						}
						if (maxq_DE < DE0.ptr<uchar>(y - 1)[x])
							maxq_DE = DE0.ptr<uchar>(y - 1)[x];
						if (maxq_DE < DE0.ptr<uchar>(y)[x])
							maxq_DE = DE0.ptr<uchar>(y)[x];
						if (x + 1 < wid)
						{
							if (maxq_DE < DE0.ptr<uchar>(y - 1)[x + 1])
								maxq_DE = DE0.ptr<uchar>(y - 1)[x + 1];
							if (maxq_DE < DE0.ptr<uchar>(y)[x + 1])
								maxq_DE = DE0.ptr<uchar>(y)[x + 1];
						}
					}
					if (maxq_DE == 0)
					{
						maxq_DE = 1;
					}
					wq_DE1 = 1 - static_cast<double>(DE0.ptr<uchar>(y)[x]) / maxq_DE;
					wq_DE2 = static_cast<double>(DE0.ptr<uchar>(y)[x]) / maxq_DE;
					double max_q_DE = 0.0f;
					if (y - 1 >= 0 && y + 1 < hei)
					{
						if (x - d - 1 >= 0)
						{
							if (max_q_DE < DE1.ptr<uchar>(y - 1)[x - d - 1])
								max_q_DE = DE1.ptr<uchar>(y - 1)[x - d - 1];
							if (max_q_DE < DE1.ptr<uchar>(y)[x - d - 1])
								max_q_DE = DE1.ptr<uchar>(y)[x - d - 1];
							if (max_q_DE < DE1.ptr<uchar>(y + 1)[x - d - 1])
								max_q_DE = DE1.ptr<uchar>(y + 1)[x - d - 1];
						}
						if (max_q_DE < DE1.ptr<uchar>(y - 1)[x - d])
							max_q_DE = DE1.ptr<uchar>(y - 1)[x - d];
						if (max_q_DE < DE1.ptr<uchar>(y)[x - d])
							max_q_DE = DE1.ptr<uchar>(y)[x - d];
						if (max_q_DE < DE1.ptr<uchar>(y + 1)[x - d])
							max_q_DE = DE1.ptr<uchar>(y + 1)[x - d];
						if (x - d + 1 < wid)
						{
							if (max_q_DE < DE1.ptr<uchar>(y - 1)[x - d + 1])
								max_q_DE = DE1.ptr<uchar>(y - 1)[x - d + 1];
							if (max_q_DE < DE1.ptr<uchar>(y)[x - d + 1])
								max_q_DE = DE1.ptr<uchar>(y)[x - d + 1];
							if (max_q_DE < DE1.ptr<uchar>(y + 1)[x - d + 1])
								max_q_DE = DE1.ptr<uchar>(y + 1)[x - d + 1];
						}
					}
					else if (y - 1 < 0)
					{
						if (x - d - 1 >= 0)
						{
							if (max_q_DE < DE1.ptr<uchar>(y)[x - d - 1])
								max_q_DE = DE1.ptr<uchar>(y)[x - d - 1];
							if (max_q_DE < DE1.ptr<uchar>(y + 1)[x - d - 1])
								max_q_DE = DE1.ptr<uchar>(y + 1)[x - d - 1];
						}
						if (max_q_DE < DE1.ptr<uchar>(y)[x - d])
							max_q_DE = DE1.ptr<uchar>(y)[x - d];
						if (max_q_DE < DE1.ptr<uchar>(y + 1)[x - d])
							max_q_DE = DE1.ptr<uchar>(y + 1)[x - d];
						if (x - d + 1 < wid)
						{
							if (max_q_DE < DE1.ptr<uchar>(y)[x - d + 1])
								max_q_DE = DE1.ptr<uchar>(y)[x - d + 1];
							if (max_q_DE < DE1.ptr<uchar>(y + 1)[x - d + 1])
								max_q_DE = DE1.ptr<uchar>(y + 1)[x - d + 1];
						}
					}
					else
					{
						if (x - d - 1 >= 0)
						{
							if (max_q_DE < DE1.ptr<uchar>(y - 1)[x - d - 1])
								max_q_DE = DE1.ptr<uchar>(y - 1)[x - d - 1];
							if (max_q_DE < DE1.ptr<uchar>(y)[x - d - 1])
								max_q_DE = DE1.ptr<uchar>(y)[x - d - 1];
						}
						if (max_q_DE < DE1.ptr<uchar>(y - 1)[x - d])
							max_q_DE = DE1.ptr<uchar>(y - 1)[x - d];
						if (max_q_DE < DE1.ptr<uchar>(y)[x - d])
							max_q_DE = DE1.ptr<uchar>(y)[x - d];
						if (x - d + 1 < wid)
						{
							if (max_q_DE < DE1.ptr<uchar>(y - 1)[x - d + 1])
								max_q_DE = DE1.ptr<uchar>(y - 1)[x - d + 1];
							if (max_q_DE < DE1.ptr<uchar>(y)[x - d + 1])
								max_q_DE = DE1.ptr<uchar>(y)[x - d + 1];
						}
					}
					if (max_q_DE == 0)
					{
						max_q_DE = 1;
					}
					w_q_DE1 = 1 - static_cast<double>(DE1.ptr<uchar>(y)[x - d]) / max_q_DE;
					w_q_DE2 = static_cast<double>(DE1.ptr<uchar>(y)[x - d]) / max_q_DE;
					const int dim_sz = 17;//��17��ά��    
					double dim_q[dim_sz], dim_q_[dim_sz];
					memset(dim_q, 0, sizeof(dim_q));
					memset(dim_q_, 0, sizeof(dim_q_));
					for (int i = 0; i < 17; i++)
					{
						dim_q[i] = 0.0f;
						dim_q_[i] = 0.0f;
					}
					dim_q[0] = wq_DE1 * static_cast<double>(leftLab.at<Vec3f>(y, x)[0]);
					dim_q_[0] =w_q_DE1 *  static_cast<double>( rightLab.at<Vec3f>(y, x - d)[0]);
					dim_q[1] =w_q_DE1 * static_cast<double>(leftLab.at<Vec3f>(y, x)[1]);
					dim_q_[1] =w_q_DE1 * static_cast<double>(rightLab.at<Vec3f>(y, x - d)[1]);
					dim_q[2] =wq_DE1 * static_cast<double>(leftLab.at<Vec3f>(y, x)[2]);
					dim_q_[2] =w_q_DE1 * static_cast<double>(rightLab.at<Vec3f>(y, x - d)[2]);

					if (wq_DE2 >= t1)
					{
						dim_q[3] = wq_DE2 * (gradLImgRedx.ptr<double>(y)[x]);// *WEIGHT_RM +gradGM0Redx.ptr<double>(y)[x] * WEIGHT_GM);
						dim_q[4] = wq_DE2 * (gradLImgGreenx.ptr<double>(y)[x]);// *WEIGHT_RM + gradGM0Greenx.ptr<double>(y)[x] * WEIGHT_GM);
						dim_q[5] = wq_DE2 * (gradLImgBluex.ptr<double>(y)[x]);// *WEIGHT_RM + gradGM0Bluex.ptr<double>(y)[x] * WEIGHT_GM);
						dim_q[6] = wq_DE2 * (gradLImgRedy.ptr<double>(y)[x]);// *WEIGHT_RM + gradGM0Redy.ptr<double>(y)[x] * WEIGHT_GM);
						dim_q[7] = wq_DE2 * (gradLImgGreeny.ptr<double>(y)[x]);// *WEIGHT_RM + gradGM0Greeny.ptr<double>(y)[x] * WEIGHT_GM);
						dim_q[8] = wq_DE2 * (gradLImgBluey.ptr<double>(y)[x]);// *WEIGHT_RM + gradGM0Bluey.ptr<double>(y)[x] * WEIGHT_GM);
						dim_q[9] = wq_DE2 * atan2(dim_q[6], dim_q[3]);
						dim_q[10] = wq_DE2 * atan2(dim_q[7], dim_q[4]);
						dim_q[11] = wq_DE2 * atan2(dim_q[8], dim_q[5]);
						dim_q[12] = dim_q[3] + dim_q[4] + dim_q[5] + dim_q[6] + dim_q[7] + dim_q[8];
						dim_q[13] = dim_q[9] + dim_q[10] + dim_q[11];
					}
					if (w_q_DE2 >= t1)
					{
						dim_q_[3] = w_q_DE2 * (gradRImgRedx.ptr<double>(y)[x - d]);// *WEIGHT_RM + gradGM1Redx.ptr<double>(y)[x - d] * WEIGHT_GM);
						dim_q_[4] = w_q_DE2 * (gradRImgGreenx.ptr<double>(y)[x - d]);// *WEIGHT_RM + gradGM1Greenx.ptr<double>(y)[x - d] * WEIGHT_GM);
						dim_q_[5] = w_q_DE2 * (gradRImgBluex.ptr<double>(y)[x - d]);// *WEIGHT_RM + gradGM1Bluex.ptr<double>(y)[x - d] * WEIGHT_GM);
						dim_q_[6] = w_q_DE2 * (gradRImgRedy.ptr<double>(y)[x - d]);// *WEIGHT_RM + gradGM1Redy.ptr<double>(y)[x - d] * WEIGHT_GM);
						dim_q_[7] = w_q_DE2 * (gradRImgGreeny.ptr<double>(y)[x - d]);// *WEIGHT_RM + gradGM1Greeny.ptr<double>(y)[x - d] * WEIGHT_GM);
						dim_q_[8] = w_q_DE2 * (gradRImgBluey.ptr<double>(y)[x - d]);// *WEIGHT_RM + gradGM1Bluey.ptr<double>(y)[x - d] * WEIGHT_GM);
						dim_q_[9] = w_q_DE2 * atan2(dim_q_[6], dim_q_[3]);
						dim_q_[10] = w_q_DE2 * atan2(dim_q_[7], dim_q_[4]);
						dim_q_[11] = w_q_DE2 * atan2(dim_q_[8], dim_q_[5]);
						dim_q_[12] = dim_q_[3] + dim_q_[4] + dim_q_[5] + dim_q_[6] + dim_q_[7] + dim_q_[8];
						dim_q_[13] = dim_q_[9] + dim_q_[10] + dim_q_[11];
					}
					if (wq_DE2 < t2)
					{
						dim_q[9] = wq_DE2 * atan2((gradLImgRedy.ptr<double>(y)[x]), (gradLImgRedx.ptr<double>(y)[x]));
						dim_q[10] = wq_DE2 * atan2((gradLImgGreeny.ptr<double>(y)[x]), (gradLImgGreenx.ptr<double>(y)[x]));
						dim_q[11] = wq_DE2 * atan2((gradLImgBluey.ptr<double>(y)[x]), (gradLImgBluex.ptr<double>(y)[x]));
						dim_q[12] = wq_DE2 * ((gradLImgRedx.ptr<double>(y)[x]) + (gradLImgRedy.ptr<double>(y)[x]) + (gradLImgGreenx.ptr<double>(y)[x])
							+ (gradLImgGreeny.ptr<double>(y)[x]) + (gradLImgBluex.ptr<double>(y)[x]) + (gradLImgBluey.ptr<double>(y)[x]));
						dim_q[13] = dim_q[9] + dim_q[10] + dim_q[11];
						dim_q[14] = static_cast<double>(wq_DE2 * BDIP0RedChannels.ptr<float>(y)[x]);
						dim_q[15] = static_cast<double>(wq_DE2 * BDIP0GreenChannels.ptr<float>(y)[x]);
						dim_q[16] = static_cast<double>(wq_DE2 * BDIP0BlueChannels.ptr<float>(y)[x]);
					}
					if (w_q_DE2 < t2)
					{
						dim_q_[9] = w_q_DE2 * atan2((gradRImgRedy.ptr<double>(y)[x - d]), (gradRImgRedx.ptr<double>(y)[x - d]));
						dim_q_[10] = w_q_DE2 * atan2((gradRImgGreeny.ptr<double>(y)[x - d]), (gradRImgGreenx.ptr<double>(y)[x - d] ));
						dim_q_[11] = w_q_DE2 * atan2((gradRImgBluey.ptr<double>(y)[x - d]), (gradRImgBluex.ptr<double>(y)[x - d]));
						dim_q_[12] = w_q_DE2 * ((gradRImgRedx.ptr<double>(y)[x - d]) + (gradRImgRedy.ptr<double>(y)[x - d]) + (gradRImgGreenx.ptr<double>(y)[x - d])
							+ (gradRImgGreeny.ptr<double>(y)[x - d]) + (gradRImgBluex.ptr<double>(y)[x - d]) + (gradRImgBluey.ptr<double>(y)[x - d]));
						dim_q_[13] = dim_q_[9] + dim_q_[10] + dim_q_[11];
						dim_q_[14] = static_cast<double>(w_q_DE2 * BDIP1RedChannels.ptr<float>(y)[x - d]);
						dim_q_[15] = static_cast<double>(w_q_DE2 * BDIP1GreenChannels.ptr<float>(y)[x - d]);
						dim_q_[16] = static_cast<double>(w_q_DE2 * BDIP1BlueChannels.ptr<float>(y)[x - d]);
					}
					double Dp_sum = 0.0f;
					for (int i = 0; i < 17; i++)
					{
						Dp_sum += fabs(dim_q[i] - dim_q_[i]);
					}
					//cost[x] = Dp_sum * W1 + cost[x] * (1 - W1);
					cost[x] = Dp_sum;
				}
			}
		}
	}
}

#ifdef COMPUTE_RIGHT
void AECencusCC::buildRightCV(const Mat& lImg, const Mat& rImg, const int maxDis, Mat* rCostVol)
{
	CV_Assert(lImg.type() == CV_64FC3 && rImg.type() == CV_64FC3);

	int hei = lImg.rows;
	int wid = lImg.cols;
	Mat lGray, rGray;
	
	Mat ltmp = lImg.clone();
	Mat leftLab, rightLab;
	ltmp.convertTo(ltmp, CV_32F, 255);
	cvtColor(ltmp, leftLab, CV_RGB2Lab);
	Mat rtmp = rImg.clone();
	rtmp.convertTo(rtmp, CV_32F, 255);
	cvtColor(rtmp, rightLab, CV_RGB2Lab);

	//Mat leftLab, rightLab;
	//Mat left_img = imread("im0.png", -1);
	//Mat right_img = imread("im1.png", -1);
	//cvtColor(left_img, leftLab, CV_BGR2Lab);
	//cvtColor(right_img, rightLab, CV_BGR2Lab);
	//Mat GM0 = MaskFilter(lImg, 9, 60, 50);
	//Mat GM1 = MaskFilter(rImg, 9, 60, 50);
	Mat lImgBlueChannels, rImgBlueChannels;// GM0BlueChannels, GM1BlueChannels;
	Mat lImgGreenChannels, rImgGreenChannels;// GM0GreenChannels, GM1GreenChannels;
	Mat lImgRedChannels, rImgRedChannels;// GM0RedChannels, GM1RedChannels;

	Mat lImgRGBGradx = lImg.clone();
	Mat rImgRGBGradx = rImg.clone();
	Mat lImgRGBGrady = lImg.clone();
	Mat rImgRGBGrady = rImg.clone();
	//Mat GM0RGBGradx = GM0.clone();
	//Mat GM0RGBGrady = GM0.clone();
	//Mat GM1RGBGradx = GM1.clone();
	//Mat GM1RGBGrady = GM1.clone();
	vector<Mat> splitLImgChannels, splitRImgChannels;// splitGM0channels, splitGM1channels;
	split(lImg, splitLImgChannels);
	split(rImg, splitRImgChannels);
	//split(GM0, splitGM0channels);
	//split(GM1, splitGM1channels);
	lImgBlueChannels = splitLImgChannels.at(0);
	lImgGreenChannels = splitLImgChannels.at(1);
	lImgRedChannels = splitLImgChannels.at(2);
	rImgBlueChannels = splitRImgChannels.at(0);
	rImgGreenChannels = splitRImgChannels.at(1);
	rImgRedChannels = splitRImgChannels.at(2);
	//GM0BlueChannels = splitGM0channels.at(0);
	//GM0GreenChannels = splitGM0channels.at(1);
	//GM0RedChannels = splitGM0channels.at(2);
	//GM1BlueChannels = splitGM1channels.at(0);
	//GM1GreenChannels = splitGM1channels.at(1);
	//GM1RedChannels = splitGM1channels.at(2);
	//double horizontal_k[3][3] = { { -1, 0, 1 }, { -2, 0, 2 }, { -1, 0, 1 } };
	//double vertical_k[3][3] = { { -1, -2, -1 }, { 0, 0, 0 }, { 1, 2, 1 } };
	//Mat maskx = Mat(3, 3, CV_32FC1, horizontal_k);
	//Mat masky = Mat(3, 3, CV_32FC1, vertical_k);
	Mat gradLImgBluex, gradRImgBluex;// gradGM0Bluex, gradGM1Bluex;
	gradLImgBluex = lImgBlueChannels.clone();
	gradRImgBluex = rImgBlueChannels.clone();
	//gradGM0Bluex = GM0BlueChannels.clone();
	//gradGM1Bluex = GM1BlueChannels.clone();
	Mat gradLImgGreenx, gradRImgGreenx;// gradGM0Greenx, gradGM1Greenx;
	gradLImgGreenx = lImgGreenChannels.clone();
	gradRImgGreenx = rImgGreenChannels.clone();
	//gradGM0Greenx = GM0GreenChannels.clone();
	//gradGM1Greenx = GM1GreenChannels.clone();
	Mat gradLImgRedx, gradRImgRedx;// gradGM0Redx, gradGM1Redx;
	gradLImgRedx = lImgRedChannels.clone();
	gradRImgRedx = rImgRedChannels.clone();
	//gradGM0Redx = GM0RedChannels.clone();
	//gradGM1Redx = GM1RedChannels.clone();
	Mat gradLImgBluey, gradRImgBluey;// gradGM0Bluey, gradGM1Bluey;
	gradLImgBluey = lImgBlueChannels.clone();
	gradRImgBluey = rImgBlueChannels.clone();
	//gradGM0Bluey = GM0BlueChannels.clone();
	//gradGM1Bluey = GM1BlueChannels.clone();
	Mat gradLImgGreeny, gradRImgGreeny;// gradGM0Greeny, gradGM1Greeny;
	gradLImgGreeny = lImgGreenChannels.clone();
	gradRImgGreeny = rImgGreenChannels.clone();
	//gradGM0Greeny = GM0GreenChannels.clone();
	//gradGM1Greeny = GM1GreenChannels.clone();
	Mat gradLImgRedy, gradRImgRedy;// gradGM0Redy, gradGM1Redy;
	gradLImgRedy = lImgRedChannels.clone();
	gradRImgRedy = rImgRedChannels.clone();
	//gradGM0Redy = GM0RedChannels.clone();
	//gradGM1Redy = GM1RedChannels.clone();

	for (int y = 0; y < hei; y++)
	{
		for (int x = 0; x < wid; x++)
		{
			if (x - 1 < 0 || x + 1 >= wid || y - 1 < 0 || y + 1 >= hei)
			{
				gradLImgBluex.ptr<double>(y)[x] = 0;
				gradLImgGreenx.ptr<double>(y)[x] = 0;
				gradLImgRedx.ptr<double>(y)[x] = 0;
				gradRImgBluey.ptr<double>(y)[x] = 0;
				gradRImgGreeny.ptr<double>(y)[x] = 0;
				gradRImgRedy.ptr<double>(y)[x] = 0;
				gradLImgBluex.ptr<double>(y)[x] = 0;
				gradLImgGreenx.ptr<double>(y)[x] = 0;
				gradLImgRedx.ptr<double>(y)[x] = 0;
				gradRImgBluey.ptr<double>(y)[x] = 0;
				gradRImgGreeny.ptr<double>(y)[x] = 0;
				gradRImgRedy.ptr<double>(y)[x] = 0;
				//gradGM0Bluex.ptr<double>(y)[x] = 0;
				//gradGM0Greenx.ptr<double>(y)[x] = 0;
				//gradGM0Redx.ptr<double>(y)[x] = 0;
				//gradGM0Bluey.ptr<double>(y)[x] = 0;
				//gradGM0Greeny.ptr<double>(y)[x] = 0;
				//gradGM0Redy.ptr<double>(y)[x] = 0;
				//gradGM1Bluex.ptr<double>(y)[x] = 0;
				//gradGM1Greenx.ptr<double>(y)[x] = 0;
				//gradGM1Redx.ptr<double>(y)[x] = 0;
				//gradGM1Bluey.ptr<double>(y)[x] = 0;
				//gradGM1Greeny.ptr<double>(y)[x] = 0;
				//gradGM1Redy.ptr<double>(y)[x] = 0;
			}
			else{
				gradLImgBluex.ptr<double>(y)[x] = lImgBlueChannels.ptr<double>(y)[x-1] - lImgBlueChannels.ptr<double>(y)[x + 1];
				gradLImgGreenx.ptr<double>(y)[x] = lImgGreenChannels.ptr<double>(y)[x - 1] - lImgGreenChannels.ptr<double>(y)[x + 1];
				gradLImgRedx.ptr<double>(y)[x] = lImgRedChannels.ptr<double>(y)[x - 1] - lImgRedChannels.ptr<double>(y)[x + 1];
				gradLImgBluey.ptr<double>(y)[x] = lImgBlueChannels.ptr<double>(y - 1)[x] - lImgBlueChannels.ptr<double>(y + 1)[x];
				gradLImgGreeny.ptr<double>(y)[x] = lImgGreenChannels.ptr<double>(y - 1)[x] - lImgGreenChannels.ptr<double>(y + 1)[x];
				gradLImgRedy.ptr<double>(y)[x] = lImgRedChannels.ptr<double>(y - 1)[x] - lImgRedChannels.ptr<double>(y + 1)[x];
				gradRImgBluex.ptr<double>(y)[x] = rImgBlueChannels.ptr<double>(y)[x-1] - rImgBlueChannels.ptr<double>(y)[x + 1];
				gradRImgGreenx.ptr<double>(y)[x] = rImgGreenChannels.ptr<double>(y)[x-1] - rImgGreenChannels.ptr<double>(y)[x + 1];
				gradRImgRedx.ptr<double>(y)[x] = rImgRedChannels.ptr<double>(y)[x-1] - rImgRedChannels.ptr<double>(y)[x + 1];
				gradRImgBluey.ptr<double>(y)[x] = rImgBlueChannels.ptr<double>(y - 1)[x] - rImgBlueChannels.ptr<double>(y + 1)[x];
				gradRImgGreeny.ptr<double>(y)[x] = rImgGreenChannels.ptr<double>(y - 1)[x] - rImgGreenChannels.ptr<double>(y + 1)[x];
				gradRImgRedy.ptr<double>(y)[x] = rImgRedChannels.ptr<double>(y - 1)[x] - rImgRedChannels.ptr<double>(y + 1)[x];
			
			}

		}
	}

	Mat BDIP0 = calcBDIP(lImg, 7);
	Mat BDIP1 = calcBDIP(rImg, 7);
	Mat DE0 = calcDE(lImg, 7);
	Mat DE1 = calcDE(rImg, 7);

	vector<Mat> splitBDIP0channels, splitBDIP1channels;
	Mat BDIP0BlueChannels, BDIP0GreenChannels, BDIP0RedChannels;
	Mat BDIP1BlueChannels, BDIP1GreenChannels, BDIP1RedChannels;
	split(BDIP0, splitBDIP0channels);
	split(BDIP1, splitBDIP1channels);
	BDIP0BlueChannels = splitBDIP0channels.at(0);
	BDIP0GreenChannels = splitBDIP0channels.at(1);
	BDIP0RedChannels = splitBDIP0channels.at(2);
	BDIP1BlueChannels = splitBDIP1channels.at(0);
	BDIP1GreenChannels = splitBDIP1channels.at(1);
	BDIP1RedChannels = splitBDIP1channels.at(2);

	//for (int d = 0; d < maxDis; d++) {
	//	for (int y = 0; y < hei; y++)
	//	{

	//		double* lData = (double*)lImg.ptr<double>(y);
	//		double* rData = (double*)rImg.ptr<double>(y);
	//		double* cost = (double*)costVol[d].ptr<double>(y);
	//		for (int x = 0; x < wid; x++)
	//		{
	//			if (x + d < wid)
	//			{
	//				double cAD = fabs(static_cast<double>(leftLab.at<Vec3f>(y, x + d)[0] - rightLab.at<Vec3f>(y, x)[0])) + fabs(static_cast<double>(leftLab.at<Vec3f>(y, x + d)[1] - rightLab.at<Vec3f>(y, x)[1])) + fabs(static_cast<double>(leftLab.at<Vec3f>(y, x + d)[2] - rightLab.at<Vec3f>(y, x)[2]));
	//				//double cADgx = 1 / 3 * (fabs(gradLImgBluex.ptr<double>(y)[x] - gradRImgBluex.ptr<double>(y)[x - d]) + fabs(gradLImgGreenx.ptr<double>(y)[x] - gradRImgGreenx.ptr<double>(y)[x - d])
	//				//	+ fabs(gradLImgRedx.ptr<double>(y)[x] - gradRImgRedx.ptr<double>(y)[x - d]));
	//				/*+ fabs(static_cast<double>(w2*gradGM0Bluex.ptr<uchar>(m)[n + d] - w2*gradGM1Bluex.ptr<uchar>(m)[n]))
	//				+ fabs(static_cast<double>(w2*gradGM0Greenx.ptr<uchar>(m)[n + d] - w2*gradGM1Greenx.ptr<uchar>(m)[n])) + fabs(static_cast<double>(w2*gradGM0Redx.ptr<uchar>(m)[n + d] - w2*gradGM1Redx.ptr<uchar>(m)[n]));*/
	//				//double cADgy = 1 / 3 * (fabs(gradLImgBluey.ptr<double>(y)[x] - gradRImgBluey.ptr<double>(y)[x - d]) + fabs(gradLImgGreeny.ptr<double>(y)[x] - gradRImgGreeny.ptr<double>(y)[x - d])
	//				//	+ fabs(gradLImgRedy.ptr<double>(y)[x] - gradRImgRedy.ptr<double>(y)[x - d]));
	//				/*+ fabs(static_cast<double>(w2*gradGM0Bluey.ptr<uchar>(m)[n + d] - w2*gradGM1Bluey.ptr<uchar>(m)[n]))
	//				+ fabs(static_cast<double>(w2*gradGM0Greeny.ptr<uchar>(m)[n + d] - w2*gradGM1Greeny.ptr<uchar>(m)[n])) + fabs(static_cast<double>(w2*gradGM0Redy.ptr<uchar>(m)[n + d] - w2*gradGM1Redy.ptr<uchar>(m)[n]));*/
	//				//cost[x] = cAD * WEIGHT_AD;// +cADgx * WEIGHT_GRDX + cADgy * WEIGHT_GRDY;
	//			}
	//		}
	//	}
	//}
/*
	Mat tmp;
	lImg.convertTo(tmp, CV_32F);
	cvtColor(tmp, lGray, CV_RGB2GRAY);
	lGray.convertTo(lGray, CV_8U, 255);
	rImg.convertTo(tmp, CV_32F);
	cvtColor(tmp, rGray, CV_RGB2GRAY);
	rGray.convertTo(rGray, CV_8U, 255);

	int w0[3][3] = { { -1, 0, 1 }, { -2, 0, 2 }, { -1, 0, 1 } };//�Ƕ�Ϊ0�ȵ�Ȩֵ
	int w45[3][3] = { { 2, 1, 0 }, { 1, 0, -1 }, { 0, -1, -2 } };//�Ƕ�Ϊ45��
	int w90[3][3] = { { -1, -2, -1 }, { 0, 0, 0 }, { 1, 2, 1 } };//�Ƕ�Ϊ90��
	int w135[3][3] = { { 0, 1, 2 }, { -1, 0, 1 }, { -2, -1, 0 } };//�Ƕ�Ϊ135��
	int *ML = new int[8];
	int *MR = new int[8];
	int *WDL = new int[8];
	int *WDR = new int[8];
	int *DL = new int[8];
	int *DR = new int[8];

	bitset<AECENCUS_BIT>* lCode = new bitset<AECENCUS_BIT>[wid * hei];
	bitset<AECENCUS_BIT>* rCode = new bitset<AECENCUS_BIT>[wid * hei];
	bitset<AECENCUS_BIT>* wLCode0 = new bitset<AECENCUS_BIT>[wid * hei];
	bitset<AECENCUS_BIT>* wRCode0 = new bitset<AECENCUS_BIT>[wid * hei];
	bitset<AECENCUS_BIT>* wLCode45 = new bitset<AECENCUS_BIT>[wid * hei];
	bitset<AECENCUS_BIT>* wRCode45 = new bitset<AECENCUS_BIT>[wid * hei];
	bitset<AECENCUS_BIT>* wLCode90 = new bitset<AECENCUS_BIT>[wid * hei];
	bitset<AECENCUS_BIT>* wRCode90 = new bitset<AECENCUS_BIT>[wid * hei];
	bitset<AECENCUS_BIT>* wLCode135 = new bitset<AECENCUS_BIT>[wid * hei];
	bitset<AECENCUS_BIT>* wRCode135 = new bitset<AECENCUS_BIT>[wid * hei];
	bitset<AECENCUS_BIT>* pLCode = lCode;
	bitset<AECENCUS_BIT>* pRCode = rCode;
	bitset<AECENCUS_BIT>* pwLCode0 = wLCode0;
	bitset<AECENCUS_BIT>* pwRCode0 = wRCode0;
	bitset<AECENCUS_BIT>* pwLCode45 = wLCode45;
	bitset<AECENCUS_BIT>* pwRCode45 = wRCode45;
	bitset<AECENCUS_BIT>* pwLCode90 = wLCode90;
	bitset<AECENCUS_BIT>* pwRCode90 = wRCode90;
	bitset<AECENCUS_BIT>* pwLCode135 = wLCode135;
	bitset<AECENCUS_BIT>* pwRCode135 = wRCode135;
	for (int y = 0; y < hei; y++)
	{
		uchar* pLData = (uchar*)(lGray.ptr<uchar>(y));
		uchar* pRData = (uchar*)(rGray.ptr<uchar>(y));
		for (int x = 0; x < wid; x++)
		{

			ML[0] = 1 / 8 * (lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x]
				+ lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y)[x]);
			ML[1] = 1 / 8 * (lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y)[x]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0]);
			ML[2] = 1 / 8 * (lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ lGray.ptr<uchar>(y)[x]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]);
			ML[3] = 1 / 8 * (lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x]
				+ lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ lGray.ptr<uchar>(y)[x]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x]);
			ML[4] = 1 / 8 * (lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ lGray.ptr<uchar>(y)[x]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]);
			ML[5] = 1 / 8 * (lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y)[x]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 2 >= 0 ? x - 2 : (x - 2 >= 0 ? wid - 1 : wid - 2)]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x]);
			ML[6] = 1 / 8 * (lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y)[x]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0]);
			ML[7] = 1 / 8 * (lGray.ptr<uchar>(y)[x]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]);

			MR[0] = 1 / 8 * (rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x]
				+ rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y)[x]);
			MR[1] = 1 / 8 * (rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y)[x]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0]);
			MR[2] = 1 / 8 * (rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ rGray.ptr<uchar>(y)[x]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]);
			MR[3] = 1 / 8 * (rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x]
				+ rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ rGray.ptr<uchar>(y)[x]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x]);
			MR[4] = 1 / 8 * (rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ rGray.ptr<uchar>(y)[x]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]);
			MR[5] = 1 / 8 * (rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y)[x]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 2 >= 0 ? x - 2 : (x - 2 >= 0 ? wid - 1 : wid - 2)]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x]);
			MR[6] = 1 / 8 * (rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y)[x]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0]);
			MR[7] = 1 / 8 * (rGray.ptr<uchar>(y)[x]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)]);
			for (int i = 0; i < 8; i++)
			{
				(*pLCode)[i] = (pLData[x] > ML[i]);
				(*pRCode)[i] = (pRData[x] > MR[i]);
			}
			pLCode++;
			pRCode++;
		}
	}

	for (int y = 0; y < hei; y++)
	{
		for (int x = 0; x < wid; x++)
		{
			WDL[0] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w0[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w0[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w0[1][2]
				+ lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[2][0]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : hei - 1] * w0[2][1]
				+ lGray.ptr<uchar>(y)[x] * w0[2][2]);
			WDL[1] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w0[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w0[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w0[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w0[1][2]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[2][0]
				+ lGray.ptr<uchar>(y)[x] * w0[2][1]
				+ lGray.ptr<uchar>(y)[x + 1 <wid ? x + 1 : 0] * w0[2][2]);
			WDL[2] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w0[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w0[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w0[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[1][2]
				+ lGray.ptr<uchar>(y)[x] * w0[2][0]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w0[2][1]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[2][2]);
			WDL[3] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[0][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[0][1]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w0[0][2]
				+ lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[1][0]
				+ lGray.ptr<uchar>(y)[x] * w0[1][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[2][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 2] * w0[2][1]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w0[2][2]);
			WDL[4] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w0[0][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w0[0][1]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[0][2]
				+ lGray.ptr<uchar>(y)[x] * w0[1][0]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[1][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w0[2][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w0[2][1]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[2][2]);
			WDL[5] = abs(lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[0][0]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[0][1]
				+ lGray.ptr<uchar>(y)[x] * w0[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w0[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w0[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w0[2][2]);
			WDL[6] = abs(lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[0][0]
				+ lGray.ptr<uchar>(y)[x] * w0[0][1]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w0[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w0[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w0[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w0[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w0[2][2]);
			WDL[7] = abs(lGray.ptr<uchar>(y)[x] * w0[0][0]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w0[0][1]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w0[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w0[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w0[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[2][2]);

			DL[0] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[1] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] - lGray.ptr<uchar>(y)[x]);
			DL[2] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);
			DL[3] = abs(lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[4] = abs(lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);
			DL[5] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[6] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] - lGray.ptr<uchar>(y)[x]);
			DL[7] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);

			WDR[0] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w0[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w0[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w0[1][2]
				+ rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[2][0]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : hei - 1] * w0[2][1]
				+ rGray.ptr<uchar>(y)[x] * w0[2][2]);
			WDR[1] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w0[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w0[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w0[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w0[1][2]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[2][0]
				+ rGray.ptr<uchar>(y)[x] * w0[2][1]
				+ rGray.ptr<uchar>(y)[x + 1 <wid ? x + 1 : 0] * w0[2][2]);
			WDR[2] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w0[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w0[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w0[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[1][2]
				+ rGray.ptr<uchar>(y)[x] * w0[2][0]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w0[2][1]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[2][2]);
			WDR[3] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[0][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[0][1]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w0[0][2]
				+ rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[1][0]
				+ rGray.ptr<uchar>(y)[x] * w0[1][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[2][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 2] * w0[2][1]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w0[2][2]);
			WDR[4] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w0[0][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w0[0][1]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[0][2]
				+ rGray.ptr<uchar>(y)[x] * w0[1][0]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[1][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w0[2][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w0[2][1]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[2][2]);
			WDR[5] = abs(rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[0][0]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[0][1]
				+ rGray.ptr<uchar>(y)[x] * w0[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w0[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w0[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w0[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w0[2][2]);
			WDR[6] = abs(rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[0][0]
				+ rGray.ptr<uchar>(y)[x] * w0[0][1]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w0[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] * w0[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w0[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w0[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w0[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w0[2][2]);
			WDR[7] = abs(rGray.ptr<uchar>(y)[x] * w0[0][0]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w0[0][1]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w0[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w0[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w0[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w0[2][2]);

			DR[0] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[1] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] - rGray.ptr<uchar>(y)[x]);
			DR[2] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);
			DR[3] = abs(rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[4] = abs(rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);
			DR[5] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[6] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] - rGray.ptr<uchar>(y)[x]);
			DR[7] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);

			for (int i = 0; i < 8; i++)
			{
				(*pwLCode0)[i] = (DL[i] > WDL[i]);
				(*pwRCode0)[i] = (DR[i] > WDR[i]);
			}
			pwLCode0++;
			pwRCode0++;
		}
	}

	for (int y = 0; y < hei; y++)
	{
		for (int x = 0; x < wid; x++)
		{
			WDL[0] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w45[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w45[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w45[1][2]
				+ lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[2][0]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : hei - 1] * w45[2][1]
				+ lGray.ptr<uchar>(y)[x] * w45[2][2]);
			WDL[1] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w45[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w45[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w45[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w45[1][2]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[2][0]
				+ lGray.ptr<uchar>(y)[x] * w45[2][1]
				+ lGray.ptr<uchar>(y)[x + 1 <wid ? x + 1 : 0] * w45[2][2]);
			WDL[2] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w45[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w45[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w45[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[1][2]
				+ lGray.ptr<uchar>(y)[x] * w45[2][0]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w45[2][1]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[2][2]);
			WDL[3] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[0][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[0][1]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w45[0][2]
				+ lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[1][0]
				+ lGray.ptr<uchar>(y)[x] * w45[1][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[2][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 2] * w45[2][1]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w45[2][2]);
			WDL[4] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w45[0][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w45[0][1]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[0][2]
				+ lGray.ptr<uchar>(y)[x] * w45[1][0]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[1][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w45[2][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w45[2][1]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[2][2]);
			WDL[5] = abs(lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[0][0]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[0][1]
				+ lGray.ptr<uchar>(y)[x] * w45[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w45[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w45[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w45[2][2]);
			WDL[6] = abs(lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[0][0]
				+ lGray.ptr<uchar>(y)[x] * w45[0][1]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w45[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w45[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w45[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w45[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w45[2][2]);
			WDL[7] = abs(lGray.ptr<uchar>(y)[x] * w45[0][0]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w45[0][1]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w45[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w45[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w45[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[2][2]);

			DL[0] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[1] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] - lGray.ptr<uchar>(y)[x]);
			DL[2] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);
			DL[3] = abs(lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[4] = abs(lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);
			DL[5] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[6] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] - lGray.ptr<uchar>(y)[x]);
			DL[7] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);

			WDR[0] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w45[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w45[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w45[1][2]
				+ rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[2][0]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : hei - 1] * w45[2][1]
				+ rGray.ptr<uchar>(y)[x] * w45[2][2]);
			WDR[1] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w45[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w45[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w45[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w45[1][2]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[2][0]
				+ rGray.ptr<uchar>(y)[x] * w45[2][1]
				+ rGray.ptr<uchar>(y)[x + 1 <wid ? x + 1 : 0] * w45[2][2]);
			WDR[2] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w45[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w45[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w45[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[1][2]
				+ rGray.ptr<uchar>(y)[x] * w45[2][0]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w45[2][1]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[2][2]);
			WDR[3] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[0][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[0][1]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w45[0][2]
				+ rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[1][0]
				+ rGray.ptr<uchar>(y)[x] * w45[1][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[2][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 2] * w45[2][1]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w45[2][2]);
			WDR[4] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w45[0][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w45[0][1]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[0][2]
				+ rGray.ptr<uchar>(y)[x] * w45[1][0]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[1][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w45[2][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w45[2][1]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[2][2]);
			WDR[5] = abs(rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[0][0]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[0][1]
				+ rGray.ptr<uchar>(y)[x] * w45[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w45[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w45[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w45[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w45[2][2]);
			WDR[6] = abs(rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[0][0]
				+ rGray.ptr<uchar>(y)[x] * w45[0][1]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w45[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] * w45[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w45[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w45[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w45[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w45[2][2]);
			WDR[7] = abs(rGray.ptr<uchar>(y)[x] * w45[0][0]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w45[0][1]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w45[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w45[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w45[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w45[2][2]);

			DR[0] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[1] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] - rGray.ptr<uchar>(y)[x]);
			DR[2] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);
			DR[3] = abs(rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[4] = abs(rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);
			DR[5] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[6] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] - rGray.ptr<uchar>(y)[x]);
			DR[7] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);

			for (int i = 0; i < 8; i++)
			{
				(*pwLCode45)[i] = (DL[i] > WDL[i]);
				(*pwRCode45)[i] = (DR[i] > WDR[i]);
			}
			pwLCode45++;
			pwRCode45++;
		}
	}

	for (int y = 0; y < hei; y++)
	{
		for (int x = 0; x < wid; x++)
		{
			WDL[0] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w90[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w90[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w90[1][2]
				+ lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[2][0]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : hei - 1] * w90[2][1]
				+ lGray.ptr<uchar>(y)[x] * w90[2][2]);
			WDL[1] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w90[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w90[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w90[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w90[1][2]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[2][0]
				+ lGray.ptr<uchar>(y)[x] * w90[2][1]
				+ lGray.ptr<uchar>(y)[x + 1 <wid ? x + 1 : 0] * w90[2][2]);
			WDL[2] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w90[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w90[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w90[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[1][2]
				+ lGray.ptr<uchar>(y)[x] * w90[2][0]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w90[2][1]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[2][2]);
			WDL[3] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[0][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[0][1]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w90[0][2]
				+ lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[1][0]
				+ lGray.ptr<uchar>(y)[x] * w90[1][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[2][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 2] * w90[2][1]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w90[2][2]);
			WDL[4] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w90[0][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w90[0][1]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[0][2]
				+ lGray.ptr<uchar>(y)[x] * w90[1][0]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[1][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w90[2][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w90[2][1]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[2][2]);
			WDL[5] = abs(lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[0][0]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[0][1]
				+ lGray.ptr<uchar>(y)[x] * w90[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w90[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w90[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w90[2][2]);
			WDL[6] = abs(lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[0][0]
				+ lGray.ptr<uchar>(y)[x] * w90[0][1]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w90[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w90[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w90[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w90[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w90[2][2]);
			WDL[7] = abs(lGray.ptr<uchar>(y)[x] * w90[0][0]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w90[0][1]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w90[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w90[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w90[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[2][2]);

			DL[0] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[1] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] - lGray.ptr<uchar>(y)[x]);
			DL[2] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);
			DL[3] = abs(lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[4] = abs(lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);
			DL[5] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[6] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] - lGray.ptr<uchar>(y)[x]);
			DL[7] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);

			WDR[0] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w90[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w90[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w90[1][2]
				+ rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[2][0]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : hei - 1] * w90[2][1]
				+ rGray.ptr<uchar>(y)[x] * w90[2][2]);
			WDR[1] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w90[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w90[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w90[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w90[1][2]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[2][0]
				+ rGray.ptr<uchar>(y)[x] * w90[2][1]
				+ rGray.ptr<uchar>(y)[x + 1 <wid ? x + 1 : 0] * w90[2][2]);
			WDR[2] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w90[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w90[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w90[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[1][2]
				+ rGray.ptr<uchar>(y)[x] * w90[2][0]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w90[2][1]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[2][2]);
			WDR[3] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[0][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[0][1]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w90[0][2]
				+ rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[1][0]
				+ rGray.ptr<uchar>(y)[x] * w90[1][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[2][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 2] * w90[2][1]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w90[2][2]);
			WDR[4] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w90[0][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w90[0][1]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[0][2]
				+ rGray.ptr<uchar>(y)[x] * w90[1][0]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[1][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w90[2][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w90[2][1]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[2][2]);
			WDR[5] = abs(rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[0][0]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[0][1]
				+ rGray.ptr<uchar>(y)[x] * w90[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w90[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w90[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w90[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w90[2][2]);
			WDR[6] = abs(rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[0][0]
				+ rGray.ptr<uchar>(y)[x] * w90[0][1]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w90[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] * w90[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w90[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w90[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w90[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w90[2][2]);
			WDR[7] = abs(rGray.ptr<uchar>(y)[x] * w90[0][0]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w90[0][1]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w90[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w90[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w90[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w90[2][2]);

			DR[0] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[1] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] - rGray.ptr<uchar>(y)[x]);
			DR[2] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);
			DR[3] = abs(rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[4] = abs(rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);
			DR[5] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[6] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] - rGray.ptr<uchar>(y)[x]);
			DR[7] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);

			for (int i = 0; i < 8; i++)
			{
				(*pwLCode90)[i] = (DL[i] > WDL[i]);
				(*pwRCode90)[i] = (DR[i] > WDR[i]);
			}
			pwLCode90++;
			pwRCode90++;
		}
	}

	for (int y = 0; y < hei; y++)
	{
		for (int x = 0; x < wid; x++)
		{
			WDL[0] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w135[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w135[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w135[1][2]
				+ lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[2][0]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : hei - 1] * w135[2][1]
				+ lGray.ptr<uchar>(y)[x] * w135[2][2]);
			WDL[1] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w135[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w135[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w135[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w135[1][2]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[2][0]
				+ lGray.ptr<uchar>(y)[x] * w135[2][1]
				+ lGray.ptr<uchar>(y)[x + 1 <wid ? x + 1 : 0] * w135[2][2]);
			WDL[2] = abs(lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w135[0][0]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w135[0][1]
				+ lGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[0][2]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w135[1][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[1][2]
				+ lGray.ptr<uchar>(y)[x] * w135[2][0]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w135[2][1]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[2][2]);
			WDL[3] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[0][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[0][1]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w135[0][2]
				+ lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[1][0]
				+ lGray.ptr<uchar>(y)[x] * w135[1][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[2][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 2] * w135[2][1]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w135[2][2]);
			WDL[4] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w135[0][0]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w135[0][1]
				+ lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[0][2]
				+ lGray.ptr<uchar>(y)[x] * w135[1][0]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[1][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w135[2][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w135[2][1]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[2][2]);
			WDL[5] = abs(lGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[0][0]
				+ lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[0][1]
				+ lGray.ptr<uchar>(y)[x] * w135[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w135[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w135[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w135[2][2]);
			WDL[6] = abs(lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[0][0]
				+ lGray.ptr<uchar>(y)[x] * w135[0][1]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w135[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w135[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w135[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w135[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w135[2][2]);
			WDL[7] = abs(lGray.ptr<uchar>(y)[x] * w135[0][0]
				+ lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w135[0][1]
				+ lGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[0][2]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w135[1][0]
				+ lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[1][2]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w135[2][0]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w135[2][1]
				+ lGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[2][2]);

			DL[0] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[1] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] - lGray.ptr<uchar>(y)[x]);
			DL[2] = abs(lGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);
			DL[3] = abs(lGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[4] = abs(lGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);
			DL[5] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] - lGray.ptr<uchar>(y)[x]);
			DL[6] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] - lGray.ptr<uchar>(y)[x]);
			DL[7] = abs(lGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] - lGray.ptr<uchar>(y)[x]);

			WDR[0] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w135[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w135[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w135[1][2]
				+ rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[2][0]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : hei - 1] * w135[2][1]
				+ rGray.ptr<uchar>(y)[x] * w135[2][2]);
			WDR[1] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x - 1 >= 0 ? x - 1 : wid - 1] * w135[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w135[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w135[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w135[1][2]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[2][0]
				+ rGray.ptr<uchar>(y)[x] * w135[2][1]
				+ rGray.ptr<uchar>(y)[x + 1 <wid ? x + 1 : 0] * w135[2][2]);
			WDR[2] = abs(rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x] * w135[0][0]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 1 < wid ? x + 1 : 0] * w135[0][1]
				+ rGray.ptr<uchar>(y - 2 >= 0 ? y - 2 : (y - 1 >= 0 ? hei - 1 : hei - 2))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[0][2]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w135[1][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[1][2]
				+ rGray.ptr<uchar>(y)[x] * w135[2][0]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w135[2][1]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[2][2]);
			WDR[3] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[0][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[0][1]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w135[0][2]
				+ rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[1][0]
				+ rGray.ptr<uchar>(y)[x] * w135[1][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[2][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 2] * w135[2][1]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w135[2][2]);
			WDR[4] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] * w135[0][0]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] * w135[0][1]
				+ rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[0][2]
				+ rGray.ptr<uchar>(y)[x] * w135[1][0]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[1][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w135[2][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w135[2][1]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[2][2]);
			WDR[5] = abs(rGray.ptr<uchar>(y)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[0][0]
				+ rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[0][1]
				+ rGray.ptr<uchar>(y)[x] * w135[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w135[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 2 >= 0 ? x - 2 : (x - 1 >= 0 ? wid - 1 : wid - 2)] * w135[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w135[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w135[2][2]);
			WDR[6] = abs(rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[0][0]
				+ rGray.ptr<uchar>(y)[x] * w135[0][1]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w135[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] * w135[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] * w135[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x - 1 >= 0 ? x - 1 : wid - 1] * w135[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w135[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w135[2][2]);
			WDR[7] = abs(rGray.ptr<uchar>(y)[x] * w135[0][0]
				+ rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] * w135[0][1]
				+ rGray.ptr<uchar>(y)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[0][2]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] * w135[1][0]
				+ rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[1][2]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x] * w135[2][0]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 1 < wid ? x + 1 : 0] * w135[2][1]
				+ rGray.ptr<uchar>(y + 2 < hei ? y + 2 : (y + 1 < hei ? 0 : 1))[x + 2 < wid ? x + 2 : (x + 1 < wid ? 0 : 1)] * w135[2][2]);

			DR[0] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[1] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x] - rGray.ptr<uchar>(y)[x]);
			DR[2] = abs(rGray.ptr<uchar>(y - 1 >= 0 ? y - 1 : hei - 1)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);
			DR[3] = abs(rGray.ptr<uchar>(y)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[4] = abs(rGray.ptr<uchar>(y)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);
			DR[5] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x - 1 >= 0 ? x - 1 : wid - 1] - rGray.ptr<uchar>(y)[x]);
			DR[6] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x] - rGray.ptr<uchar>(y)[x]);
			DR[7] = abs(rGray.ptr<uchar>(y + 1 < hei ? y + 1 : 0)[x + 1 < wid ? x + 1 : 0] - rGray.ptr<uchar>(y)[x]);

			for (int i = 0; i < 8; i++)
			{
				(*pwLCode135)[i] = (DL[i] > WDL[i]);
				(*pwRCode135)[i] = (DR[i] > WDR[i]);
			}
			pwLCode135++;
			pwRCode135++;
		}
	}

	// build cost volume
	bitset<AECENCUS_BIT> lB, wlB0, wlB45, wlB90, wlB135;
	bitset<AECENCUS_BIT> rB, wrB0, wrB45, wrB90, wrB135;
	pRCode = rCode;
	pwRCode0 = wRCode0;
	pwRCode45 = wRCode45;
	pwRCode90 = wRCode90;
	pwRCode135 = wRCode135;
	for (int y = 0; y < hei; y++) {
		int index = y * wid;
		for (int x = 0; x < wid; x++) {
			rB = *pRCode;
			wrB0 = *pwRCode0;
			wrB45 = *pwRCode45;
			wrB90 = *pwRCode90;
			wrB135 = *pwRCode135;
			for (int d = 0; d < maxDis; d++) {
				double* cost = (double*)rCostVol[d].ptr<double>(y);
				double tmpCost1 = AECENCUS_BIT, tmpCost2 = AECENCUS_BIT;
				if (x + d < wid) {
					lB = lCode[index + x + d];
					wlB0 = wLCode0[index + x + d];
					wlB45 = wLCode45[index + x + d];
					wlB90 = wLCode90[index + x + d];
					wlB135 = wLCode135[index + x + d];
					tmpCost1 = (lB ^ rB).count();
					tmpCost2 = (wlB0 ^ wrB0).count() + (wlB45 ^ wrB45).count() + (wlB90 ^ wrB90).count() + (wlB135 ^ wrB135).count();
				}
				//cost[x] = WEIGHT_M * min(tmpCost1,7.2) + (1 - WEIGHT_M) * min(tmpCost2,28.8);
			}
			pRCode++;
			pwRCode0++;
			pwRCode45++;
			pwRCode90++;
			pwRCode135++;
		}
	}
	delete[] lCode;
	delete[] rCode;
	delete[] wLCode0;
	delete[] wRCode0;
	delete[] wLCode45;
	delete[] wRCode45;
	delete[] wLCode90;
	delete[] wRCode90;
	delete[] wLCode135;
	delete[] wRCode135;
*/
	for (int y = 0; y < hei; y++)
	{
		
		for (int x = 0; x < wid; x++)
		{
			for (int d = 0; d < maxDis; d++)
			{
				double* cost = (double*)rCostVol[d].ptr<double>(y);
				if (x + d < wid) {
					double w_pq = 0.0f;
					double w_p_q_ = 0.0f;
					double wq_DE1 = 0.0f;
					double wq_DE2 = 0.0f;
					double w_q_DE1 = 0.0f;
					double w_q_DE2 = 0.0f;
					double maxq_DE = 0.0f;
					if (y - 1 >= 0 && y + 1 < hei)
					{
						if (x + d - 1 >= 0)
						{
							if (maxq_DE < DE0.ptr<uchar>(y - 1)[x + d - 1])
								maxq_DE = DE0.ptr<uchar>(y - 1)[x + d - 1];
							if (maxq_DE < DE0.ptr<uchar>(y)[x + d - 1])
								maxq_DE = DE0.ptr<uchar>(y)[x + d - 1];
							if (maxq_DE < DE0.ptr<uchar>(y + 1)[x + d - 1])
								maxq_DE = (DE0.ptr<uchar>(y + 1)[x + d - 1]);
						}
						if (maxq_DE < (DE0.ptr<uchar>(y - 1)[x + d]))
							maxq_DE = (DE0.ptr<uchar>(y - 1)[x + d]);
						if (maxq_DE < (DE0.ptr<uchar>(y)[x + d]))
							maxq_DE = (DE0.ptr<uchar>(y)[x + d]);
						if (maxq_DE < (DE0.ptr<uchar>(y + 1)[x + d]))
							maxq_DE = (DE0.ptr<uchar>(y + 1)[x + d]);
						if (x + d + 1 < wid)
						{
							if (maxq_DE < (DE0.ptr<uchar>(y - 1)[x + d + 1]))
								maxq_DE = (DE0.ptr<uchar>(y - 1)[x + d + 1]);
							if (maxq_DE < (DE0.ptr<uchar>(y)[x + d + 1]))
								maxq_DE = (DE0.ptr<uchar>(y)[x + d + 1]);
							if (maxq_DE < (DE0.ptr<uchar>(y + 1)[x + d + 1]))
								maxq_DE = (DE0.ptr<uchar>(y + 1)[x + d + 1]);
						}
					}
					else if (y - 1 < 0)
					{
						if (x + d - 1 >= 0)
						{
							if (maxq_DE < (DE0.ptr<uchar>(y)[x + d - 1]))
								maxq_DE = (DE0.ptr<uchar>(y)[x + d - 1]);
							if (maxq_DE < (DE0.ptr<uchar>(y + 1)[x + d - 1]))
								maxq_DE = (DE0.ptr<uchar>(y + 1)[x + d - 1]);
						}
						if (maxq_DE < (DE0.ptr<uchar>(y)[x + d]))
							maxq_DE = (DE0.ptr<uchar>(y)[x + d]);
						if (maxq_DE < (DE0.ptr<uchar>(y + 1)[x + d]))
							maxq_DE = (DE0.ptr<uchar>(y + 1)[x + d]);
						if (x + d + 1 < wid)
						{
							if (maxq_DE < (DE0.ptr<uchar>(y)[x + d + 1]))
								maxq_DE = (DE0.ptr<uchar>(y)[x + d + 1]);
							if (maxq_DE < (DE0.ptr<uchar>(y + 1)[x + d + 1]))
								maxq_DE = (DE0.ptr<uchar>(y + 1)[x + d + 1]);
						}
					}
					else
					{
						if (x + d - 1 >= 0)
						{
							if (maxq_DE < (DE0.ptr<uchar>(y - 1)[x + d - 1]))
								maxq_DE = (DE0.ptr<uchar>(y - 1)[x + d - 1]);
							if (maxq_DE < (DE0.ptr<uchar>(y)[x + d - 1]))
								maxq_DE = (DE0.ptr<uchar>(y)[x + d - 1]);
						}
						if (maxq_DE < (DE0.ptr<uchar>(y - 1)[x + d]))
							maxq_DE = (DE0.ptr<uchar>(y - 1)[x + d]);
						if (maxq_DE < (DE0.ptr<uchar>(y)[x + d]))
							maxq_DE = (DE0.ptr<uchar>(y)[x + d]);
						if (x + d + 1 < wid)
						{
							if (maxq_DE < (DE0.ptr<uchar>(y - 1)[x + d + 1]))
								maxq_DE = (DE0.ptr<uchar>(y - 1)[x + d + 1]);
							if (maxq_DE < (DE0.ptr<uchar>(y)[x + d + 1]))
								maxq_DE = (DE0.ptr<uchar>(y)[x + d + 1]);
						}
					}
					if (maxq_DE == 0)
					{
						maxq_DE = 1;
					}
					wq_DE1 = 1 - static_cast<double>(DE0.ptr<uchar>(y)[x + d]) / maxq_DE;
					wq_DE2 = static_cast<double>(DE0.ptr<uchar>(y)[x + d]) / maxq_DE;
					double max_q_DE = 0.0f;
					if (y - 1 >= 0 && y + 1 < hei)
					{
						if (x - 1 >= 0)
						{
							if (max_q_DE < (DE1.ptr<uchar>(y - 1)[x - 1]))
								max_q_DE = (DE1.ptr<uchar>(y - 1)[x - 1]);
							if (max_q_DE < (DE1.ptr<uchar>(y)[x - 1]))
								max_q_DE = (DE1.ptr<uchar>(y)[x - 1]);
							if (max_q_DE < (DE1.ptr<uchar>(y + 1)[x - 1]))
								max_q_DE = (DE1.ptr<uchar>(y + 1)[x - 1]);
						}
						if (max_q_DE < (DE1.ptr<uchar>(y - 1)[x]))
							max_q_DE = (DE1.ptr<uchar>(y - 1)[x]);
						if (max_q_DE < (DE1.ptr<uchar>(y)[x]))
							max_q_DE = (DE1.ptr<uchar>(y)[x]);
						if (max_q_DE < (DE1.ptr<uchar>(y + 1)[x]))
							max_q_DE = (DE1.ptr<uchar>(y + 1)[x]);
						if (x + 1 < wid)
						{
							if (max_q_DE < (DE1.ptr<uchar>(y - 1)[x + 1]))
								max_q_DE = (DE1.ptr<uchar>(y - 1)[x + 1]);
							if (max_q_DE < (DE1.ptr<uchar>(y)[x + 1]))
								max_q_DE = (DE1.ptr<uchar>(y)[x + 1]);
							if (max_q_DE < (DE1.ptr<uchar>(y + 1)[x + 1]))
								max_q_DE = (DE1.ptr<uchar>(y + 1)[x + 1]);
						}
					}
					else if (y - 1 < 0)
					{
						if (x - 1 >= 0)
						{
							if (max_q_DE < (DE1.ptr<uchar>(y)[x - 1]))
								max_q_DE = (DE1.ptr<uchar>(y)[x - 1]);
							if (max_q_DE < (DE1.ptr<uchar>(y + 1)[x - 1]))
								max_q_DE = (DE1.ptr<uchar>(y + 1)[x - 1]);
						}
						if (max_q_DE < (DE1.ptr<uchar>(y)[x]))
							max_q_DE = (DE1.ptr<uchar>(y)[x]);
						if (max_q_DE < (DE1.ptr<uchar>(y + 1)[x]))
							max_q_DE = (DE1.ptr<uchar>(y + 1)[x]);
						if (x + 1 < wid)
						{
							if (max_q_DE < (DE1.ptr<uchar>(y)[x + 1]))
								max_q_DE = (DE1.ptr<uchar>(y)[x + 1]);
							if (max_q_DE < (DE1.ptr<uchar>(y + 1)[x + 1]))
								max_q_DE = (DE1.ptr<uchar>(y + 1)[x + 1]);
						}
					}
					else
					{
						if (x - 1 >= 0)
						{
							if (max_q_DE < (DE1.ptr<uchar>(y - 1)[x - 1]))
								max_q_DE = (DE1.ptr<uchar>(y - 1)[x - 1]);
							if (max_q_DE < (DE1.ptr<uchar>(y)[x - 1]))
								max_q_DE = (DE1.ptr<uchar>(y)[x - 1]);
						}
						if (max_q_DE < (DE1.ptr<uchar>(y - 1)[x]))
							max_q_DE = (DE1.ptr<uchar>(y - 1)[x]);
						if (max_q_DE < (DE1.ptr<uchar>(y)[x]))
							max_q_DE = (DE1.ptr<uchar>(y)[x]);
						if (x + 1 < wid)
						{
							if (max_q_DE < (DE1.ptr<uchar>(y - 1)[x + 1]))
								max_q_DE = (DE1.ptr<uchar>(y - 1)[x + 1]);
							if (max_q_DE < (DE1.ptr<uchar>(y)[x + 1]))
								max_q_DE = (DE1.ptr<uchar>(y)[x + 1]);
						}
					}
					if (max_q_DE == 0)
					{
						max_q_DE = 1;
					}
					w_q_DE1 = 1 - static_cast<double>(DE1.ptr<uchar>(y)[x]) / max_q_DE;
					w_q_DE2 = static_cast<double>(DE1.ptr<uchar>(y)[x]) / max_q_DE;
					const int dim_sz = 17;//��17��ά��    
					double dim_q[dim_sz], dim_q_[dim_sz];
					memset(dim_q, 0, sizeof(dim_q));
					memset(dim_q_, 0, sizeof(dim_q_));
					for (int i = 0; i < 17; i++)
					{
						dim_q[i] = 0.0f;
						dim_q_[i] = 0.0f;
					}
					dim_q[0] =wq_DE1 * static_cast<double>(leftLab.at<Vec3f>(y, x + d)[0]);
					dim_q_[0] = w_q_DE1 * static_cast<double>(rightLab.at<Vec3f>(y, x)[0]);
					dim_q[1] = wq_DE1 * static_cast<double>(leftLab.at<Vec3f>(y, x + d)[1]);
					dim_q_[1] =w_q_DE1 * static_cast<double>(rightLab.at<Vec3f>(y, x)[1]);
					dim_q[2] =wq_DE1 * static_cast<double>(leftLab.at<Vec3f>(y, x + d)[2]);
					dim_q_[2] =w_q_DE1 * static_cast<double>(rightLab.at<Vec3f>(y, x)[2]);

					if (wq_DE2 >= t1)
					{
						dim_q[3] = wq_DE2 * (gradLImgRedx.ptr<double>(y)[x + d]);// *WEIGHT_RM + gradGM0Redx.ptr<double>(y)[x + d] * WEIGHT_GM);
						dim_q[4] = wq_DE2 * (gradLImgGreenx.ptr<double>(y)[x + d]);// *WEIGHT_RM + gradGM0Greenx.ptr<double>(y)[x + d] * WEIGHT_GM);
						dim_q[5] = wq_DE2 * (gradLImgBluex.ptr<double>(y)[x + d]);// *WEIGHT_RM + gradGM0Bluex.ptr<double>(y)[x + d] * WEIGHT_GM);
						dim_q[6] = wq_DE2 * (gradLImgRedy.ptr<double>(y)[x + d]);// *WEIGHT_RM + gradGM0Redy.ptr<double>(y)[x + d] * WEIGHT_GM);
						dim_q[7] = wq_DE2 * (gradLImgGreeny.ptr<double>(y)[x + d]);// *WEIGHT_RM + gradGM0Greeny.ptr<double>(y)[x + d] * WEIGHT_GM);
						dim_q[8] = wq_DE2 * (gradLImgBluey.ptr<double>(y)[x + d]);// *WEIGHT_RM + gradGM0Bluey.ptr<double>(y)[x + d] * WEIGHT_GM);
						dim_q[9] = wq_DE2 * atan2(dim_q[6], dim_q[3]);
						dim_q[10] = wq_DE2 * atan2(dim_q[7], dim_q[4]);
						dim_q[11] = wq_DE2 * atan2(dim_q[8], dim_q[5]);
						dim_q[12] = dim_q[3] + dim_q[4] + dim_q[5] + dim_q[6] + dim_q[7] + dim_q[8];
						dim_q[13] = dim_q[9] + dim_q[10] + dim_q[11];
					}
					if (w_q_DE2 >= t1)
					{
						dim_q_[3] = w_q_DE2 * (gradRImgRedx.ptr<double>(y)[x]);// *WEIGHT_RM + gradGM1Redx.ptr<double>(y)[x] * WEIGHT_GM);
						dim_q_[4] = w_q_DE2 * (gradRImgGreenx.ptr<double>(y)[x]);// *WEIGHT_RM + gradGM1Greenx.ptr<double>(y)[x] * WEIGHT_GM);
						dim_q_[5] = w_q_DE2 * (gradRImgBluex.ptr<double>(y)[x]);// *WEIGHT_RM + gradGM1Bluex.ptr<double>(y)[x] * WEIGHT_GM);
						dim_q_[6] = w_q_DE2 * (gradRImgRedy.ptr<double>(y)[x]);// *WEIGHT_RM + gradGM1Redy.ptr<double>(y)[x] * WEIGHT_GM);
						dim_q_[7] = w_q_DE2 * (gradRImgGreeny.ptr<double>(y)[x]);// *WEIGHT_RM + gradGM1Greeny.ptr<double>(y)[x] * WEIGHT_GM);
						dim_q_[8] = w_q_DE2 * (gradRImgBluey.ptr<double>(y)[x]);// *WEIGHT_RM + gradGM1Bluey.ptr<double>(y)[x] * WEIGHT_GM);
						dim_q_[9] = w_q_DE2 * atan2(dim_q_[6], dim_q_[3]);
						dim_q_[10] = w_q_DE2 * atan2(dim_q_[7], dim_q_[4]);
						dim_q_[11] = w_q_DE2 * atan2(dim_q_[8], dim_q_[5]);
						dim_q_[12] = dim_q_[3] + dim_q_[4] + dim_q_[5] + dim_q_[6] + dim_q_[7] + dim_q_[8];
						dim_q_[13] = dim_q_[9] + dim_q_[10] + dim_q_[11];
					}
					if (wq_DE2 < t2)
					{
						dim_q[9] = wq_DE2 * atan2((gradLImgRedy.ptr<double>(y)[x + d]), (gradLImgRedx.ptr<double>(y)[x + d]));
						dim_q[10] = wq_DE2 * atan2((gradLImgGreeny.ptr<double>(y)[x + d]), (gradLImgGreenx.ptr<double>(y)[x + d]));
						dim_q[11] = wq_DE2 * atan2((gradLImgBluey.ptr<double>(y)[x + d]), (gradLImgBluex.ptr<double>(y)[x + d]));
						dim_q[12] = wq_DE2 * ((gradLImgRedx.ptr<double>(y)[x + d]) + (gradLImgRedy.ptr<double>(y)[x + d]) + (gradLImgGreenx.ptr<double>(y)[x + d])
							+ (gradLImgGreeny.ptr<double>(y)[x + d]) + (gradLImgBluex.ptr<double>(y)[x + d]) + (gradLImgBluey.ptr<double>(y)[x + d]));
						dim_q[13] = dim_q[9] + dim_q[10] + dim_q[11];
						dim_q[14] = static_cast<double>(wq_DE2 * BDIP0RedChannels.ptr<float>(y)[x + d]);
						dim_q[15] = static_cast<double>(wq_DE2 * BDIP0GreenChannels.ptr<float>(y)[x + d]);
						dim_q[16] = static_cast<double>(wq_DE2 * BDIP0BlueChannels.ptr<float>(y)[x + d]);
					}
					if (w_q_DE2 < t2)
					{
						dim_q_[9] = w_q_DE2 * atan2((gradRImgRedy.ptr<double>(y)[x]), (gradRImgRedx.ptr<double>(y)[x]));
						dim_q_[10] = w_q_DE2 * atan2((gradRImgGreeny.ptr<double>(y)[x]), (gradRImgGreenx.ptr<double>(y)[x]));
						dim_q_[11] = w_q_DE2 * atan2((gradRImgBluey.ptr<double>(y)[x]), (gradRImgBluex.ptr<double>(y)[x]));
						dim_q_[12] = w_q_DE2 * ((gradRImgRedx.ptr<double>(y)[x]) + (gradRImgRedy.ptr<double>(y)[x]) + (gradRImgGreenx.ptr<double>(y)[x])
							+ (gradRImgGreeny.ptr<double>(y)[x]) + (gradRImgBluex.ptr<double>(y)[x]) + (gradRImgBluey.ptr<double>(y)[x]));
						dim_q_[13] = dim_q_[9] + dim_q_[10] + dim_q_[11];
						dim_q_[14] = static_cast<double>(w_q_DE2 * BDIP1RedChannels.ptr<float>(y)[x]);
						dim_q_[15] = static_cast<double>(w_q_DE2 * BDIP1GreenChannels.ptr<float>(y)[x]);
						dim_q_[16] = static_cast<double>(w_q_DE2 * BDIP1BlueChannels.ptr<float>(y)[x]);
					}
					double Dp_sum = 0.0f;
					for (int i = 0; i < 17; i++)
					{
						Dp_sum += fabs(dim_q[i] - dim_q_[i]);
					}
					//cost[x] = cost[x] * W1 + Dp_sum * (1 - W1);
					cost[x] = Dp_sum;
				}
			}
		}
	}
}
#endif
