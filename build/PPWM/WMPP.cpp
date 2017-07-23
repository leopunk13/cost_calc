#include "WMPP.h"

namespace WMPP_FUNC{
	void lrCheck( Mat& lDis, Mat& rDis, int* lValid, int* rValid, const int disSc, const int maxDis)
	{
		int hei = lDis.rows;
		int wid = lDis.cols;
		int imgSize = hei * wid;
		memset( lValid, 0, imgSize * sizeof( int ) );
		memset( rValid, 0, imgSize * sizeof( int ) );
		int* pLValid = lValid;
		int* pRValid = rValid;
		for( int y = 0; y < hei; y ++ ) {
			double* lDisData = ( double* ) lDis.ptr<double>( y );
			double* rDisData = ( double* ) rDis.ptr<double>( y );
			for( int x = 0; x < wid; x ++ ) {
				// check left image
				int lDep = lDisData[ x ] / disSc;
				//assert( ( x - lDep ) >= 0 && ( x - lDep ) < wid );
				int rLoc = ( x - lDep + wid ) % wid;
				//int rLoc = x - max(lDep, 0);
				int rDep = rDisData[ rLoc ] / disSc;
				// disparity should not be zero
				if( lDep == rDep && lDep >= 5) {
					*pLValid = 1;
				}
				// check right image
				rDep = rDisData[ x ] / disSc;
				//assert( ( x + rDep ) >= 0 && ( x + rDep ) < wid );
				int lLoc = ( x + rDep + wid ) % wid;
				lDep = lDisData[ lLoc ] / disSc;
				// disparity should not be zero
				if( rDep == lDep && rDep >= 5 ) {
					*pRValid = 1;
				}
				pLValid ++;
				pRValid ++;
			}
		}
	}
	void fillInv( Mat& lDis, Mat& rDis, int* lValid, int* rValid)
	{
		int hei = lDis.rows;
		int wid = lDis.cols;
		// fill left dep
		int* pLValid = lValid;
		for( int y = 0; y < hei; y ++ ) {
			int* yLValid = lValid + y * wid;
			double* lDisData = ( double* ) lDis.ptr<double>( y );
			for( int x = 0; x < wid; x ++ ) {
				if( *pLValid == 0) {
					// find left first valid pixel
					int lFirst = x;
					int lFind = 0;
					while( lFirst >= 0 ) {
						if( yLValid[ lFirst ] ) {
							lFind = 1;
							break;
						}
						lFirst --;
					}
					int rFind = 0;
					// find right first valid pixel
					int rFirst = x;
					while( rFirst < wid ) {
						if( yLValid[ rFirst ] ) {
							rFind = 1;
							break;
						}
						rFirst ++;
					}
					// set x's depth to the lowest one
					if( lFind && rFind ) {
						if( lDisData[ lFirst ] <= lDisData[ rFirst ] ) {
							lDisData[ x ] = lDisData[ lFirst ];
						} else {
							lDisData[ x ] = lDisData[ rFirst ];
						}
					} else if( lFind ) {
						lDisData[ x ] = lDisData[ lFirst ];
					} else if ( rFind ) {
						lDisData[ x ] = lDisData[ rFirst ];
					}

				}
				pLValid ++;
			}
		}
		// fill right dep
		int* pRValid = rValid;
		for( int y = 0; y < hei; y ++ ) {
			int* yRValid = rValid + y * wid;
			double* rDisData = ( double* ) ( rDis.ptr<double>( y ) );
			for( int x = 0; x < wid; x ++ ) {
				if( *pRValid == 0) {
					// find left first valid pixel
					int lFirst = x;
					int lFind = 0;
					while( lFirst >= 0 ) {
						if( yRValid[ lFirst ] ) {
							lFind = 1;
							break;
						}
						lFirst --;
					}
					// find right first valid pixel
					int rFirst = x;
					int rFind = 0;
					while( rFirst < wid ) {
						if( yRValid[ rFirst ] ) {
							rFind = 1;
							break;
						}
						rFirst ++;
					}
					if( lFind && rFind ) {
						// set x's depth to the lowest one
						if( rDisData[ lFirst ] <= rDisData[ rFirst ] ) {
							rDisData[ x ] = rDisData[ lFirst ];
						} else {
							rDisData[ x ] = rDisData[ rFirst ];
						}
					} else if( lFind ) {
						rDisData[ x ] = rDisData[ lFirst ];
					} else if ( rFind )  {
						rDisData[ x ] = rDisData[ rFirst ];
					}

				}
				pRValid ++;
			}
		}
	}

	void ItRegionVoting(Mat lImg, Mat rImg, Mat lDis, Mat rDis, int* lValid, int* rValid, const int maxDis)
	{
		int hei = lDis.rows;
		int wid = rDis.cols;
		Mat ltmp = lImg.clone();
		Mat LeftLab, RightLab;
		ltmp.convertTo(ltmp, CV_32F, 255);
		cvtColor(ltmp, LeftLab, CV_RGB2Lab);
		Mat rtmp = rImg.clone();
		rtmp.convertTo(rtmp, CV_32F, 255);
		cvtColor(rtmp, RightLab, CV_RGB2Lab);
		int* pLValid = lValid;
		//for (int it = 0; it < 5; it++)
		//{
			for (int y = 0; y < hei; y++)
			{
				for (int x = 0; x < wid; x++)
				{
					int temp1 = 0;
					int temp2 = 0;
					double hist[maxDis];
					memset(hist, 0, sizeof(hist));
					for (int i = 0; i < maxDis; i++)
					{
						hist[i] = 0.0f;
					}
					if (*pLValid == 0 && x < maxDis)
					{
						int left = x - 1;
						int right = x + 1;
						int up = y + 1;
						int down = y - 1;
						int* upLValid = lValid + up * wid;
						int* downLValid = lValid + down * wid;
						while (up < hei && sqrt(pow(LeftLab.at<Vec3b>(y, x)[0] - LeftLab.at<Vec3b>(up, x)[0], 2) + pow(LeftLab.at<Vec3b>(y, x)[1] - LeftLab.at<Vec3b>(up, x)[1], 2) + pow(LeftLab.at<Vec3b>(y, x)[2] - LeftLab.at<Vec3b>(up, x)[2], 2)) < Tc && abs(up - y) < Ts)
						{
							while (left >= 0 && sqrt(pow(LeftLab.at<Vec3b>(up, x)[0] - LeftLab.at<Vec3b>(up, left)[0], 2) + pow(LeftLab.at<Vec3b>(up, x)[1] - LeftLab.at<Vec3b>(up, left)[1], 2) + pow(LeftLab.at<Vec3b>(up, x)[2] - LeftLab.at<Vec3b>(up, left)[2], 2)) < Tc && abs(left - x) < Ts)
							{
								if (upLValid[left])
								{
									temp1 = lDis.ptr<double>(up)[left];
									temp2 = int(temp1);//������ת��  
									hist[temp2]++;
								}
								left--;
							}
							while (right < wid && sqrt(pow(LeftLab.at<Vec3b>(up, x)[0] - LeftLab.at<Vec3b>(up, right)[0], 2) + pow(LeftLab.at<Vec3b>(up, x)[1] - LeftLab.at<Vec3b>(up, right)[1], 2) + pow(LeftLab.at<Vec3b>(up, x)[2] - LeftLab.at<Vec3b>(up, right)[2], 2)) < Tc && abs(right - x) < Ts)
							{
								if (upLValid[right])
								{
									temp1 = lDis.ptr<double>(up)[right];
									temp2 = int(temp1);//������ת��  
									hist[temp2]++;
								}
								right++;
							}
							up++;
						}

						while (down >= 0 && sqrt(pow(LeftLab.at<Vec3b>(y, x)[0] - LeftLab.at<Vec3b>(down, x)[0], 2) + pow(LeftLab.at<Vec3b>(y, x)[1] - LeftLab.at<Vec3b>(down, x)[1], 2) + pow(LeftLab.at<Vec3b>(y, x)[2] - LeftLab.at<Vec3b>(down, x)[2], 2)) < Tc && abs(down - y) < Ts)
						{
							while (left >= 0 && sqrt(pow(LeftLab.at<Vec3b>(down, x)[0] - LeftLab.at<Vec3b>(down, left)[0], 2) + pow(LeftLab.at<Vec3b>(down, x)[1] - LeftLab.at<Vec3b>(down, left)[1], 2) + pow(LeftLab.at<Vec3b>(down, x)[2] - LeftLab.at<Vec3b>(down, left)[2], 2)) < Tc && abs(left - x) < Ts)
							{
								if (downLValid[left])
								{
									temp1 = lDis.ptr<double>(down)[left];
									temp2 = int(temp1);//������ת��  
									hist[temp2]++;
								}
								left--;
							}
							while (right < wid && sqrt(pow(LeftLab.at<Vec3b>(down, x)[0] - LeftLab.at<Vec3b>(down, right)[0], 2) + pow(LeftLab.at<Vec3b>(down, x)[1] - LeftLab.at<Vec3b>(down, right)[1], 2) + pow(LeftLab.at<Vec3b>(down, x)[2] - LeftLab.at<Vec3b>(down, right)[2], 2)) < Tc && abs(right - x) < Ts)
							{
								if (downLValid[right])
								{
									temp1 = lDis.ptr<double>(down)[right];
									temp2 = int(temp1);//������ת��  
									hist[temp2]++;
								}
								right++;
							}
							down--;
						}

						double max_hist = 0.0f;
						int max_value = 0;
						double Sp = 0.0f;
						for (int i = 0; i < maxDis; i++)
						{
							if (hist[i] != 0)
							{
								if (max_hist <= hist[i])
								{
									max_hist = hist[i];
									max_value = i;
								}
								Sp += hist[i];
							}
						}
						if (Sp > 5 && max_hist / Sp > 0.33)
							lDis.ptr<double>(y)[x] = max_value;
					}
					pLValid ++;
				}
			}
		//}
		int* pRValid = rValid;
		//for (int it = 0; it < 5; it++)
		//{
			for (int y = 0; y < hei; y++)
			{
				for (int x = 0; x < wid; x++)
				{
					int temp1 = 0;
					int temp2 = 0;
					double hist[maxDis];
					memset(hist, 0, sizeof(hist));
					for (int i = 0; i < maxDis; i++)
					{
						hist[i] = 0.0f;
					}
					if (*pRValid == 0 && x < wid - maxDis)
					{
						int left = x - 1;
						int right = x + 1;
						int up = y + 1;
						int down = y - 1;
						int* upRValid = rValid + up * wid;
						int* downRValid = rValid + down * wid;
						while (up < hei && sqrt(pow(RightLab.at<Vec3b>(y, x)[0] - RightLab.at<Vec3b>(up, x)[0], 2) + pow(RightLab.at<Vec3b>(y, x)[1] - RightLab.at<Vec3b>(up, x)[1], 2) + pow(RightLab.at<Vec3b>(y, x)[2] - RightLab.at<Vec3b>(up, x)[2], 2)) < Tc && abs(up - y) < Ts)
						{
							while (left >= 0 && sqrt(pow(RightLab.at<Vec3b>(up, x)[0] - RightLab.at<Vec3b>(up, left)[0], 2) + pow(RightLab.at<Vec3b>(up, x)[1] - RightLab.at<Vec3b>(up, left)[1], 2) + pow(RightLab.at<Vec3b>(up, x)[2] - RightLab.at<Vec3b>(up, left)[2], 2)) < Tc && abs(left - x) < Ts)
							{
								if (upRValid[left])
								{
									temp1 = rDis.ptr<double>(up)[left];
									temp2 = int(temp1);//������ת��  
									hist[temp2]++;
								}
								left--;
							}
							while (right < wid && sqrt(pow(RightLab.at<Vec3b>(up, x)[0] - RightLab.at<Vec3b>(up, right)[0], 2) + pow(RightLab.at<Vec3b>(up, x)[1] - RightLab.at<Vec3b>(up, right)[1], 2) + pow(RightLab.at<Vec3b>(up, x)[2] - RightLab.at<Vec3b>(up, right)[2], 2)) < Tc && abs(right - x) < Ts)
							{
								if (upRValid[right])
								{
									temp1 = rDis.ptr<double>(up)[right];
									temp2 = int(temp1);//������ת��  
									hist[temp2]++;
								}
								right++;
							}
							up++;
						}

						while (down >= 0 && sqrt(pow(RightLab.at<Vec3b>(y, x)[0] - RightLab.at<Vec3b>(down, x)[0], 2) + pow(RightLab.at<Vec3b>(y, x)[1] - RightLab.at<Vec3b>(down, x)[1], 2) + pow(RightLab.at<Vec3b>(y, x)[2] - RightLab.at<Vec3b>(down, x)[2], 2)) < Tc && abs(down - y) < Ts)
						{
							while (left >= 0 && sqrt(pow(RightLab.at<Vec3b>(down, x)[0] - RightLab.at<Vec3b>(down, left)[0], 2) + pow(RightLab.at<Vec3b>(down, x)[1] - RightLab.at<Vec3b>(down, left)[1], 2) + pow(RightLab.at<Vec3b>(down, x)[2] - RightLab.at<Vec3b>(down, left)[2], 2)) < Tc && abs(left - x) < Ts)
							{
								if (downRValid[left])
								{
									temp1 = rDis.ptr<double>(down)[left];
									temp2 = int(temp1);//������ת��  
									hist[temp2]++;
								}
								left--;
							}
							while (right < wid && sqrt(pow(RightLab.at<Vec3b>(down, x)[0] - RightLab.at<Vec3b>(down, right)[0], 2) + pow(RightLab.at<Vec3b>(down, x)[1] - RightLab.at<Vec3b>(down, right)[1], 2) + pow(RightLab.at<Vec3b>(down, x)[2] - RightLab.at<Vec3b>(down, right)[2], 2)) < Tc && abs(right - x) < Ts)
							{
								if (downRValid[right])
								{
									temp1 = rDis.ptr<double>(down)[right];
									temp2 = int(temp1);//������ת��  
									hist[temp2]++;
								}
								right++;
							}
							down--;
						}

						double max_hist = 0.0f;
						int max_value = 0;
						double Sp = 0.0f;
						for (int i = 0; i < maxDis; i++)
						{
							if (hist[i] != 0)
							{
								if (max_hist <= hist[i])
								{
									max_hist = hist[i];
									max_value = i;
								}
								Sp += hist[i];
							}
						}
						if (Sp > 5 && max_hist / Sp > 0.33)
							rDis.ptr<double>(y)[x] = max_value;
					}
					pRValid ++;
				}
			}
		//}
	}

	Mat SmallHoleFilling(Mat leftDisp)
	{
		int nc = leftDisp.cols;
		int nr = leftDisp.rows;
		for (int it = 0; it < 5; it++)
		{
			for (int y = 0; y < nr; y++)
			{
				for (int x = 0; x < nc; x++)
				{
					int left = x - 1;
					int right = x + 1;
					int up = y + 1;
					int down = y - 1;
					double dl = 255.0f;
					double dr = 255.0f;
					double du = 255.0f;
					double dd = 255.0f;
					if (leftDisp.ptr<double>(y)[x] == 0)
					{
						while (left >= 0 && (leftDisp.ptr<double>(y)[left] == 0 || leftDisp.ptr<double>(y)[left] == 255))
						{

							left--;
						}
						if (left >= 0)
							dl = leftDisp.ptr<double>(y)[left];
						while (right < nc && (leftDisp.ptr<double>(y)[right] == 0 || leftDisp.ptr<double>(y)[right] == 255))
						{

							right++;
						}
						if (right < nc)
							dr = leftDisp.ptr<double>(y)[right];
						while (up < nr && (leftDisp.ptr<double>(up)[x] == 0 || leftDisp.ptr<double>(up)[x] == 255))
						{

							up++;
						}
						if (up < nr)
							du = leftDisp.ptr<double>(up)[x];
						while (down >= 0 && (leftDisp.ptr<double>(down)[x] == 0 || leftDisp.ptr<double>(down)[x] == 255))
						{
							down--;
						}
						if (down >= 0)
							dd = leftDisp.ptr<double>(down)[x];
						if (min(dl, dr) == 255.0 && min(du, dd) == 255.0)
							leftDisp.ptr<double>(y)[x] = 0;
						else if (min(dl, dr) == 255.0 && min(du, dd) != 255.0)
							leftDisp.ptr<double>(y)[x] = min(du, dd);
						else if (min(du, dd) == 255 && min(dl, dr) != 255)
							leftDisp.ptr<double>(y)[x] = min(dl, dr);
						else if (abs(min(dl, dr) - min(du, dd) <= 2))
							leftDisp.ptr<double>(y)[x] = (min(dl, dr) + min(du, dd)) / 2;
						else
							leftDisp.ptr<double>(y)[x] = 0;

					}
					if (leftDisp.ptr<double>(y)[x] == 255)
					{
						while (left >= 0 && (leftDisp.ptr<double>(y)[left] == 0 || leftDisp.ptr<double>(y)[left] == 255))
						{
							left--;
						}
						if (left >= 0)
							dl = leftDisp.ptr<double>(y)[left];
						while (right < nc && (leftDisp.ptr<double>(y)[right] == 0 || leftDisp.ptr<double>(y)[right] == 255))
						{
							right++;
						}
						if (right < nc)
							dr = leftDisp.ptr<double>(y)[right];
						while (up < nr && (leftDisp.ptr<double>(up)[x] == 0 || leftDisp.ptr<double>(up)[x] == 255))
						{
							up++;
						}
						if (up < nr)
							du = leftDisp.ptr<double>(up)[x];
						while (down >= 0 && (leftDisp.ptr<double>(down)[x] == 0 || leftDisp.ptr<double>(down)[x] == 255))
						{
							down--;
						}
						if (down >= 0)
							dd = leftDisp.ptr<double>(down)[x];
						if (min(dl, dr) == 255 && min(du, dd) == 255)
							leftDisp.ptr<double>(y)[x] = 255;
						else if (min(dl, dr) == 255 && min(du, dd) != 255)
							leftDisp.ptr<double>(y)[x] = min(du, dd);
						else if (min(du, dd) == 255 && min(dl, dr) != 255)
							leftDisp.ptr<double>(y)[x] = min(dl, dr);
						else
							leftDisp.ptr<double>(y)[x] = min(min(dl, dr), min(du, dd));
					}

				}
			}
		}

		return leftDisp;
	}

	void wgtMedian( const Mat& lImg, const Mat& rImg, Mat& lDis, Mat& rDis, int* lValid, int* rValid, const int maxDis, const int disSc )
	{
		int hei = lDis.rows;
		int wid = lDis.cols;
		int wndR = MED_SZ / 2;
		double* disHist = new double[ maxDis ];

		// filter left
		int* pLValid = lValid;
		for( int y = 0; y < hei; y ++  ) {
			double* lDisData = ( double* ) lDis.ptr<double>( y );
			float* pL = ( float* ) lImg.ptr<float>( y );
			for( int x = 0; x < wid; x ++ ) {
				if( *pLValid == 0 ) {
					// just filter invalid pixels
					memset( disHist, 0, sizeof( double ) * maxDis );
					double sumWgt = 0.0f;
					// set disparity histogram by bilateral weight
					for( int wy = - wndR; wy <= wndR; wy ++ ) {
						int qy = ( y + wy + hei ) % hei;
						// int* qLValid = lValid + qy * wid;
						float* qL = ( float* ) lImg.ptr<float>( qy );
						double* qDisData = ( double* ) lDis.ptr<double>( qy );
						for( int wx = - wndR; wx <= wndR; wx ++ ) {
							int qx = ( x + wx + wid ) % wid;
							// invalid pixel also used
							// if( qLValid[ qx ] && wx != 0 && wy != 0 ) {
							int qDep = qDisData[ qx ] / disSc;
							if( qDep != 0 ) {

								double disWgt = wx * wx + wy * wy;
								// disWgt = sqrt( disWgt );
								double clrWgt = ( pL[ 3 * x ] - qL[ 3 * qx ] ) * ( pL[ 3 * x ] - qL[ 3 * qx ] ) +
									( pL[ 3 * x + 1 ] - qL[ 3 * qx + 1 ] ) * ( pL[ 3 * x + 1 ] - qL[ 3 * qx + 1 ] ) +
									( pL[ 3 * x + 2 ] - qL[ 3 * qx + 2 ] ) * ( pL[ 3 * x + 2 ] - qL[ 3 * qx + 2 ] );
								// clrWgt = sqrt( clrWgt );
								double biWgt = exp( - disWgt / ( SIG_DIS * SIG_DIS ) - clrWgt / ( SIG_CLR * SIG_CLR ) );
								disHist[ qDep ] += biWgt;
								sumWgt += biWgt;
							}
							// }
						}
					}
					double halfWgt = sumWgt / 2.0f;
					sumWgt = 0.0f;
					int filterDep = 0;
					for( int d = 0; d < maxDis; d ++ ) {
						sumWgt += disHist[ d ];
						if( sumWgt >= halfWgt ) {
							filterDep = d;
							break;
						}
					}
					// set new disparity
					lDisData[ x ] = filterDep * disSc;
				}
				pLValid ++;
			}
		}
		// filter right depth
		int* pRValid = rValid;
		for( int y = 0; y < hei; y ++  ) {
			double* rDisData = ( double* ) rDis.ptr<double>( y );
			float* pR = ( float* ) rImg.ptr<float>( y );
			for( int x = 0; x < wid; x ++ ) {
				if( *pRValid == 0 ) {
					// just filter invalid pixels
					memset( disHist, 0, sizeof( double ) * maxDis );
					double sumWgt = 0.0f;
					// set disparity histogram by bilateral weight
					for( int wy = - wndR; wy <= wndR; wy ++ ) {
						int qy = ( y + wy + hei ) % hei;
						// int* qRValid = rValid + qy * wid;
						float* qR = ( float* ) rImg.ptr<float>( qy );
						double* qDisData = ( double* ) rDis.ptr<double>( qy );
						for( int wx = - wndR; wx <= wndR; wx ++ ) {
							int qx = ( x + wx + wid ) % wid;
							// if( qRValid[ qx ] && wx != 0 && wy != 0 ) {
							int qDep = qDisData[ qx ] / disSc;
							if( qDep != 0 ) {

								double disWgt = wx * wx + wy * wy;
								disWgt = sqrt( disWgt );
								double clrWgt =
									( pR[ 3 * x ] - qR[ 3 * qx ] ) * ( pR[ 3 * x ] - qR[ 3 * qx ] ) +
									( pR[ 3 * x + 1 ] - qR[ 3 * qx + 1 ] ) * ( pR[ 3 * x + 1 ] - qR[ 3 * qx + 1 ] ) +
									( pR[ 3 * x + 2 ] - qR[ 3 * qx + 2 ] ) * ( pR[ 3 * x + 2 ] - qR[ 3 * qx + 2 ] );
								clrWgt = sqrt( clrWgt );
								double biWgt = exp( - disWgt / ( SIG_DIS * SIG_DIS ) - clrWgt / ( SIG_CLR * SIG_CLR ) );
								disHist[ qDep ] += biWgt;
								sumWgt += biWgt;
							}
							// }
						}
					}
					double halfWgt = sumWgt / 2.0f;
					sumWgt = 0.0f;
					int filterDep = 0;
					for( int d = 0; d < maxDis; d ++ ) {
						sumWgt += disHist[ d ];
						if( sumWgt >= halfWgt ) {
							filterDep = d;
							break;
						}
					}
					// set new disparity
					rDisData[ x ] = filterDep * disSc;
				}
				pRValid ++;
			}
		}
		Mat l_chk, r_chk;
		lDis.convertTo(l_chk,CV_16UC1);
		rDis.convertTo(r_chk,CV_16UC1);
		imwrite("disp0_pp.png", l_chk);
		imwrite("disp1_pp.png", r_chk);
		delete [] disHist;
	}
	void saveChk( const int hei, const int wid,  int* lValid, int* rValid, const int maxDis )
	{
		Mat lChk = Mat::zeros( hei, wid, CV_16UC1 );
		Mat rChk = Mat::zeros( hei, wid, CV_16UC1 );
		int* pLV = lValid;
		int* pRV = rValid;
		for( int y = 0; y < hei; y ++ ) {
			double* lChkData = ( double* )( lChk.ptr<double>( y ) );
			double* rChkData = ( double* )( rChk.ptr<double>( y ) );
			for( int x = 0; x < wid; x ++ ) {
				if( *pLV ) {
					lChkData[ x ] = 0;
				} else{
					lChkData[ x ] = maxDis;
				}

				if( *pRV ) {
					rChkData[ x ] = 0;
				} else{
					rChkData[ x ] = maxDis;
				}
				pLV ++;
				pRV ++;
			}
		}
		Mat l_chk, r_chk;
		lChk.convertTo(l_chk,CV_16UC1);
		rChk.convertTo(r_chk,CV_16UC1);
		imwrite( "l_chk.png", l_chk );
		imwrite( "r_chk.png", r_chk );
	}
}

void WMPP::postProcess( const Mat& lImg, const Mat& rImg, const int maxDis, const int disSc, Mat& lDis, Mat& rDis,
	Mat& lSeg, Mat& lChk)
{
	// color image should be 3x3 median filtered
	// according to weightedMedianMatlab.m from CVPR11
	Mat lTmp, rTmp;
	lImg.convertTo( lTmp, CV_32F );
	rImg.convertTo( rTmp, CV_32F );
	int hei = lDis.rows;
	int wid = lDis.cols;
	int imgSize = hei * wid;
	int* lValid = new int[ imgSize ];
	int* rValid = new int[ imgSize ];
	
	//WMPP_FUNC::lrCheck( lDis, rDis, lValid, rValid, disSc );
	
	//WMPP_FUNC::wgtMedian( lTmp, rTmp, lDis, rDis, lValid, rValid, maxDis, disSc );
	// iter 3 times
    for( int i = 0; i < 3; i ++ ) {
		// save check results
	    //WMPP_FUNC::wgtMedian( lTmp, rTmp, lDis, rDis, lValid, rValid, maxDis, disSc );
		WMPP_FUNC::lrCheck( lDis, rDis, lValid, rValid, disSc, maxDis);
		WMPP_FUNC::fillInv(lDis, rDis, lValid, rValid);
		WMPP_FUNC::wgtMedian( lTmp, rTmp, lDis, rDis, lValid, rValid, maxDis, disSc );
		/*WMPP_FUNC::lrCheck( lDis, rDis, lValid, rValid, disSc, maxDis );		
		WMPP_FUNC::ItRegionVoting(lImg, rImg, lDis, rDis, lValid, rValid, maxDis);
		WMPP_FUNC::wgtMedian( lTmp, rTmp, lDis, rDis, lValid, rValid, maxDis, disSc );*/
	}
	//WMPP_FUNC:: fillInv(lDis, rDis, lValid, rValid );
	// see last check results
	WMPP_FUNC::lrCheck( lDis, rDis, lValid, rValid, disSc, maxDis );
	WMPP_FUNC::saveChk( hei, wid, lValid, rValid ,maxDis);
	delete [] lValid;
	delete [] rValid;
}
