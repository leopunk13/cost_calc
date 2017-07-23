///////////////////////////////////////////////////////
// File: GetMethod
// Desc: Get all the method interface by name
//
// Author: Zhang Kang
// Date: 2013/09/06
///////////////////////////////////////////////////////
#pragma  once
#include "CC/GrdCC.h"
#include "CC/CenCC.h"
#include "CC/CGCC.h"
#include "CC/AECensusCC.h"
#include "CC/ADDCensusCC.h"
#include "CAFilter/GFCA.h"
#include "CAFilter/BFCA.h"
#include "CAFilter/BoxCA.h"
#include "CAST/STCA.h"
#include "PPWM/WMPP.h"
#include "PPSG/SGPP.h"

// get cost compuation method name
CCMethod* getCCType( const string name );

// get cost aggregation method name
CAMethod* getCAType( const string name );

// get cost compuation method name
PPMethod* getPPType( const string name );
