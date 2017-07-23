Mat BoxFilter( const Mat& imSrc, const int r )
{
	int H = imSrc.rows;
	int W = imSrc.cols;
	// image size must large than filter size
	CV_Assert( W >= r && H >= r );
	Mat imDst = Mat::ones( H, W, imSrc.type() );

	// cumulative sum over Y axis
	Mat imCum = CumSum( imSrc, 1 );

	// difference along Y ( [ 0, r ], [r + 1, H - r - 1], [ H - r, H ] )
	for( int y = 0; y < r + 1; y ++ ) {
		double* dstData = ( double* ) imDst.ptr<double>( y );
		double* plusData = ( double* ) imCum.ptr<double>( y + r );
		for( int x = 0; x < W; x ++ ) {
			dstData[ x ] = plusData[ x ];
		}

	}
	for( int y = r + 1; y < H - r; y ++ ) {
		double* dstData = ( double* ) imDst.ptr<double>( y );
		double* minusData = ( double*  ) imCum.ptr<double>( y - r - 1);
		double* plusData = ( double* ) imCum.ptr<double>( y + r );
		for( int x = 0; x < W; x ++ ) {
			dstData[ x ] = plusData[ x ] - minusData[ x ];
		}
	}
	for( int y = H - r; y < H; y ++ ) {
		double* dstData = ( double* ) imDst.ptr<double>( y );
		double* minusData = ( double*  ) imCum.ptr<double>( y - r - 1);
		double* plusData = ( double* ) imCum.ptr<double>( H - 1 );
		for( int x = 0; x < W; x ++ ) {
			dstData[ x ] = plusData[ x ] - minusData[ x ];
		}
	}

	// cumulative sum over X axis
	imCum = CumSum( imDst, 2 );

	for( int y = 0; y < H; y ++ ) {
		double* dstData = ( double* ) imDst.ptr<double>( y );
		double* cumData = ( double* ) imCum.ptr<double>( y );
		for( int x = 0; x < r + 1; x ++ ) {
			dstData[ x ] = cumData[ x + r ];
		}
		for( int x = r + 1; x < W - r; x ++ ) {
			dstData[ x ] = cumData[ x + r ] - cumData[ x - r - 1 ];

		}
		for( int x = W - r; x < W; x ++ ) {
			dstData[ x ] = cumData[ W - 1 ] - cumData[ x - r - 1 ];
		}
	}

	return imDst;
}
