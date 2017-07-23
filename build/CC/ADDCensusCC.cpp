#include "ADDCensusCC.h"
#include <fstream>
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

void ADDCensusCC::buildCV( const Mat& lImg, const Mat& rImg, const int maxDis, Mat* costVol )
{
	// for TAD + Grd input image must be CV_64FC3
	CV_Assert( lImg.type() == CV_64FC3 && rImg.type() == CV_64FC3 );

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


	Mat lImgBlueChannels, rImgBlueChannels;// GM0BlueChannels, GM1BlueChannels;
	Mat lImgGreenChannels, rImgGreenChannels;// GM0GreenChannels, GM1GreenChannels;
	Mat lImgRedChannels, rImgRedChannels;// GM0RedChannels, GM1RedChannels;

	Mat lImgRGBGradx = lImg.clone();
	Mat rImgRGBGradx = rImg.clone();
	Mat lImgRGBGrady = lImg.clone();
	Mat rImgRGBGrady = rImg.clone();

	vector<Mat> splitLImgChannels, splitRImgChannels;// splitGM0channels, splitGM1channels;
	split(lImg, splitLImgChannels);
    split(rImg, splitRImgChannels);
	lImgBlueChannels = splitLImgChannels.at(0);
	lImgGreenChannels = splitLImgChannels.at(1);
	lImgRedChannels = splitLImgChannels.at(2);
	rImgBlueChannels = splitRImgChannels.at(0);
	rImgGreenChannels = splitRImgChannels.at(1);
	rImgRedChannels = splitRImgChannels.at(2);

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

	Mat tmp;
	lImg.convertTo( tmp, CV_32F );
	cvtColor( tmp, lGray, CV_RGB2GRAY );
	lGray.convertTo( lGray, CV_8U, 255 );
	rImg.convertTo( tmp, CV_32F );
	cvtColor( tmp, rGray, CV_RGB2GRAY );
	rGray.convertTo( rGray, CV_8U, 255 );
	// prepare binary code 
	int H_WD = CENCUS_WND / 2;
	bitset<CENCUS_BIT>* lCode = new bitset<CENCUS_BIT>[ wid * hei ];
	bitset<CENCUS_BIT>* rCode = new bitset<CENCUS_BIT>[ wid * hei ];
	bitset<CENCUS_BIT>* pLCode = lCode;
	bitset<CENCUS_BIT>* pRCode = rCode;
	for( int y = 0; y < hei; y ++ ) {
		uchar* pLData = ( uchar* ) ( lGray.ptr<uchar>( y ) );
		uchar* pRData = ( uchar* ) ( rGray.ptr<uchar>( y ) );
		for( int x = 0; x < wid; x ++ ) {
			int bitCnt = 0;
			for( int wy = - H_WD; wy <= H_WD; wy ++ ) {
				int qy = ( y + wy + hei ) % hei;
				uchar* qLData = ( uchar* ) ( lGray.ptr<uchar>( qy ) );
				uchar* qRData = ( uchar* ) ( rGray.ptr<uchar>( qy ) );
				for( int wx = - H_WD; wx <= H_WD; wx ++ ) {
					if( wy != 0 || wx != 0 ) {
						int qx = ( x + wx + wid ) % wid;
						( *pLCode )[ bitCnt ] = ( pLData[ x ] > qLData[ qx ] );
						( *pRCode )[ bitCnt ] = ( pRData[ x ] > qRData[ qx ] );
						bitCnt ++;
					}

				}
			}
			pLCode ++;
			pRCode ++;
		}
	}
	// build cost volume
	bitset<CENCUS_BIT> lB;
	bitset<CENCUS_BIT> rB;
	pLCode = lCode;
	for( int y = 0; y < hei; y ++ ) {
		int index = y * wid;
		for( int x = 0; x < wid; x ++ ) {
			lB = *pLCode;
			for( int d = 0; d < maxDis; d ++ ) {
				double* cost   = ( double* ) costVol[ d ].ptr<double>( y );
				cost[ x ] = CENCUS_BIT;
				if( x - d >= 0 ) {
					rB = rCode[ index + x - d ];
					cost[ x ] = ( lB ^ rB ).count();
				}

			}
			pLCode ++;
		}
	}
	delete [] lCode;
	delete [] rCode;

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
					cost[x] = cost[x] * W1 + Dp_sum * (1 - W1);
					//cost[x] = (1 - exp(-cost[x]/W1))  + (1 - exp(-Dp_sum/W2));
				}
			}
		}
	}
}
#ifdef COMPUTE_RIGHT
void ADDCensusCC::buildRightCV( const Mat& lImg, const Mat& rImg, const int maxDis, Mat* rCostVol )
{
	// for TAD + Grd input image must be CV_64FC3
	CV_Assert( lImg.type() == CV_64FC3 && rImg.type() == CV_64FC3 );
        //ofstream fout1, fout2;
	//fout1.open("cost_1.txt");
	//fout2.open("cost_2.txt");
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


	Mat tmp;
	lImg.convertTo( tmp, CV_32F );
	cvtColor( tmp, lGray, CV_RGB2GRAY );
	lGray.convertTo( lGray, CV_8U, 255 );
	rImg.convertTo( tmp, CV_32F );
	cvtColor( tmp, rGray, CV_RGB2GRAY );
	rGray.convertTo( rGray, CV_8U, 255 );
	// prepare binary code 
	int H_WD = CENCUS_WND / 2;
	bitset<CENCUS_BIT>* lCode = new bitset<CENCUS_BIT>[ wid * hei ];
	bitset<CENCUS_BIT>* rCode = new bitset<CENCUS_BIT>[ wid * hei ];
	bitset<CENCUS_BIT>* pLCode = lCode;
	bitset<CENCUS_BIT>* pRCode = rCode;
	for( int y = 0; y < hei; y ++ ) {
		uchar* pLData = ( uchar* ) ( lGray.ptr<uchar>( y ) );
		uchar* pRData = ( uchar* ) ( rGray.ptr<uchar>( y ) );
		for( int x = 0; x < wid; x ++ ) {
			int bitCnt = 0;
			for( int wy = - H_WD; wy <= H_WD; wy ++ ) {
				int qy = ( y + wy + hei ) % hei;
				uchar* qLData = ( uchar* ) ( lGray.ptr<uchar>( qy ) );
				uchar* qRData = ( uchar* ) ( rGray.ptr<uchar>( qy ) );
				for( int wx = - H_WD; wx <= H_WD; wx ++ ) {
					if( wy != 0 || wx != 0 ) {
						int qx = ( x + wx + wid ) % wid;
						( *pLCode )[ bitCnt ] = ( pLData[ x ] > qLData[ qx ] );
						( *pRCode )[ bitCnt ] = ( pRData[ x ] > qRData[ qx ] );
						bitCnt ++;
					}

				}
			}
			pLCode ++;
			pRCode ++;
		}
	}
	// build cost volume
	bitset<CENCUS_BIT> lB;
	bitset<CENCUS_BIT> rB;
	pRCode = rCode;
	for( int y = 0; y < hei; y ++ ) {
		int index = y * wid;
		for( int x = 0; x < wid; x ++ ) {
			rB = *pRCode;
			for( int d = 0; d < maxDis; d ++ ) {
				double* cost   = ( double* ) rCostVol[ d ].ptr<double>( y );
				cost[ x ] = CENCUS_BIT;
				if( x + d < wid ) {
					lB = lCode[ index + x + d ];
					cost[ x ] = ( rB ^ lB ).count();
				}

			}
			pRCode ++;
		}
	}
	delete [] lCode;
	delete [] rCode;

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
					//fout1 << cost[x] <<endl;
					//fout2 << Dp_sum <<endl;
					cost[x] = cost[x] * W1 + Dp_sum * (1 - W1);
					//cost[x] = (1 - exp(-cost[x]/W1))  + (1 - exp(-Dp_sum/W2));
				}
			}
		}
	}
	//fout1.close();
	//fout2.close();
}
#endif
