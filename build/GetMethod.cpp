#include "GetMethod.h"


// get cost compuation method name
CCMethod* getCCType( const string name ) 
{
	if( name == "GRD" ) {
		return new GrdCC();
	} else if( name == "CEN" ) {
		return new CenCC();
	} else if( name == "BSM" ) {
		return NULL;
	} else if ( name == "CG" ) {
		return new CGCC();
	}else if (name == "AECEN") {
		return new AECencusCC();
	}else if (name == "ADDCEN"){
		return new ADDCensusCC();
	}
}

// get cost aggregation method name
CAMethod* getCAType( const string name )
{
	if( name == "GF" ) {
		return new GFCA();
	} else if( name == "BF" ) {
		return new BFCA();
	} else if( name == "BOX" ) {
		return new BoxCA();
	} else if ( name == "ST" ) {
		return new STCA();
	}
}

// get cost compuation method name
PPMethod* getPPType( const string name ) 
{
    if( name == "WM" ) {
		return new WMPP();
	} else if( name == "NP" ) {
		return NULL;
	} else if( name == "SG" ) {
		return new SGPP();
	}
}
