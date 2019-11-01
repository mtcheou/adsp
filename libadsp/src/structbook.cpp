// structbook.cpp

#include "structbook.h"

// ============================
//  CStructBook Class
// ============================
CStructBook::CStructBook()
{
	m_iNumElement =0;
	m_norm = 0.0;
}

CStructBook::~CStructBook()
{
	m_iNumElement =0;
	m_norm = 0.0;
}

void CStructBook::setNumElement ( int numElement )
{
	m_iNumElement = numElement;
}

void CStructBook::setNorm( double norm )
{
	m_norm = norm;
}

int	CStructBook::getNumElement()
{
	return m_iNumElement;
}

double	CStructBook::getNorm()
{
	return m_norm;
}



// ============================
//  CStructBookExp Class
// ============================
CStructBookExp::CStructBookExp()
{
	m_structBook = NULL;
	m_opCurve = NULL;
	m_numElemOpCurve = 0;
    m_minamp =0.0;
    m_maxamp =0.0;
    m_minrho =0.0;
    m_maxrho =0.0;
    m_minphase =0.0;
    m_maxphase =0.0;
}

CStructBookExp::~CStructBookExp()
{
	if ( m_structBook!=NULL ) 
	{
	    
	    delete [] m_structBook;
	    m_structBook = NULL;
	    m_iNumElement = 0;
	    m_norm = 0.0;
	}
	if (m_opCurve != NULL)
	{
	    delete [] m_opCurve;
	    m_opCurve = NULL;
	}
}


// Intrinsic

//===============================================================
// Function: addElement
// Goal: Destructor
// Return:
//===============================================================

void CStructBookExp::addElement ( CParameter* parm )
{
	strtContinuousExp* aux;
	strtContinuousExp cExpStructure;

	if ( m_structBook==NULL )
	{
		m_structBook = new strtContinuousExp[1];
	}

	aux = new strtContinuousExp[m_iNumElement + 1];

	memcpy ( aux,m_structBook,m_iNumElement*sizeof ( strtContinuousExp ) );

	cExpStructure.innerProduct = ( ( CExpParm* ) parm )->innerProd;
	cExpStructure.rho = ( ( CExpParm* ) parm )->rho;
	cExpStructure.xi = ( ( CExpParm* ) parm )->xi;
	cExpStructure.phase = ( ( CExpParm* ) parm )->phase;
	cExpStructure.a = ( ( CExpParm* ) parm )->a;
	cExpStructure.b = ( ( CExpParm* ) parm )->b;

	aux[m_iNumElement] = cExpStructure;

	if ( m_structBook!=NULL )
	{
        delete [] m_structBook;
	}

	m_structBook = new strtContinuousExp[m_iNumElement + 1];

	memcpy ( m_structBook,aux, ( m_iNumElement+1 ) *sizeof ( strtContinuousExp ) );

	m_iNumElement++;

	delete [] aux;
}

//===============================================================
// Function: addElement
// Goal: Destructor
// Return:
//===============================================================

void CStructBookExp::addElement ( strtContinuousExp cExpStructure )
{
    strtContinuousExp* aux;

    if ( m_structBook==NULL )
    {
        m_structBook = new strtContinuousExp;
    }

    aux = new strtContinuousExp[m_iNumElement + 1];

    memcpy ( aux,m_structBook,m_iNumElement*sizeof ( strtContinuousExp ) );

    aux[m_iNumElement] = cExpStructure;

    if ( m_structBook!=NULL )
    {
        delete [] m_structBook;
    }

    m_structBook = new strtContinuousExp[m_iNumElement + 1];

    memcpy ( m_structBook,aux, ( m_iNumElement+1 ) *sizeof ( strtContinuousExp ) );

    m_iNumElement++;

    delete [] aux;
}

//===============================================================
// Function: addElement
// Goal: Destructor
// Return:
//===============================================================

void CStructBookExp::addElement ( strtContinuousExp* cExpStructure, int numEl )
{
    strtContinuousExp* aux;

    if ( m_structBook==NULL )
    {
        m_structBook = new strtContinuousExp;
    }

    aux = new strtContinuousExp[m_iNumElement + numEl];

    memcpy ( aux,m_structBook,m_iNumElement*sizeof ( strtContinuousExp ) );

    memcpy ( &aux[m_iNumElement],cExpStructure, ( numEl ) *sizeof ( strtContinuousExp ) );

    if ( m_structBook!=NULL )
    {
        delete [] m_structBook;
    }

    m_structBook = new strtContinuousExp[m_iNumElement + numEl];

    memcpy ( m_structBook,aux, ( m_iNumElement+numEl ) *sizeof ( strtContinuousExp ) );

    m_iNumElement += numEl;

    delete [] aux;
}

//===============================================================
// Function: setNextAtomIndex
// Goal: 
// Return:
//===============================================================
void CStructBookExp::setNextAtomIndex ( int atomIndex,
                                        int nextAtom)
{
	m_structBook[atomIndex].nextAtom = nextAtom;
}

//===============================================================
// Function: setPrevAtomIndex
// Goal: 
// Return:
//===============================================================
void CStructBookExp::setPrevAtomIndex ( int atomIndex,
                                        int prevAtom)
{
	m_structBook[atomIndex].prevAtom = prevAtom;
}

//===============================================================
// Function: setOrigAtomIndex
// Goal: 
// Return:
//===============================================================
void CStructBookExp::setOrigAtomIndex ( int atomIndex,
                                        int origAtom)
{
	m_structBook[atomIndex].origAtomIndex = origAtom;
}

//===============================================================
// Function: removeElement
// Goal: 
// Return:
//===============================================================

void CStructBookExp::removeElement ( int atomIndex )
{
    if ( m_structBook==NULL )
    {
        cout << "There is nothing to remove!!!" << endl;
        exit;
        
    }
    if (atomIndex>(m_iNumElement-1))
    {
        cout << "Atom index greater than SB number of elements!!!" << endl;
        exit;
    }
    if (m_iNumElement==1)
    {
        delete [] m_structBook;
        m_structBook=NULL;
        m_iNumElement=0;
        return;
    }
        
    //cout << atomIndex<< endl;
    //cout << m_iNumElement << endl;
    
    strtContinuousExp* aux;
    aux = new strtContinuousExp[m_iNumElement - 1];

    if (atomIndex!=0)
    {
        memcpy ( aux,m_structBook,(atomIndex)*sizeof ( strtContinuousExp ) );
        //cout << (atomIndex+1) << endl;
        //cout << "PAssou 1"<<endl;
    }
    if ( (m_iNumElement-1-atomIndex)!=0)
    {
        memcpy ( &aux[atomIndex],&m_structBook[atomIndex+1],(m_iNumElement-1-atomIndex)*sizeof ( strtContinuousExp ) );
        //cout << atomIndex << " " << (atomIndex+1) << " " << (m_iNumElement-1-atomIndex) << endl;
        //cout << "PAssou 2"<<endl;
    }
    if (m_structBook!=NULL) delete [] m_structBook;
    m_structBook = new strtContinuousExp[m_iNumElement - 1];
    memcpy ( m_structBook,aux, ( m_iNumElement-1 ) *sizeof ( strtContinuousExp ) );
    //cout << "PAssou 3"<<endl;

    m_iNumElement--;

    delete [] aux;
}

//===============================================================
// Function: saveElementASCII
// Goal:
// Return:
//===============================================================

void CStructBookExp::saveElementASCII ( FILE* stream,
                                        double meanApproxRatio,
                                        double approxRatio,
                                        double befSupInnerP,
                                        double aftSupInnerP,
                                        double normRatio)
{
    double snr = 20*(log(1.0/(normRatio))/log(10.0));
	fprintf ( stream,"%5i %15.8f %15.8f %15.8f %15.8f %5i %5i    %5i %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f\n",
	          m_iNumElement,
	          m_structBook[m_iNumElement-1].innerProduct,
	          m_structBook[m_iNumElement-1].rho,
	          m_structBook[m_iNumElement-1].xi,
	          m_structBook[m_iNumElement-1].phase,
	          m_structBook[m_iNumElement-1].a,
	          m_structBook[m_iNumElement-1].b,
              m_structBook[m_iNumElement-1].prevAtom+1,
	          approxRatio,
	          meanApproxRatio,
	          befSupInnerP,
	          aftSupInnerP,
              normRatio,
              snr);
}

void CStructBookExp::saveElementBin( FILE* stream)
{
    int numwritten;
    numwritten = fwrite( &m_iNumElement, sizeof(int), 1, stream );
    numwritten = fwrite( &m_structBook[m_iNumElement-1], sizeof(strtContinuousExp), 1, stream );
}

void CStructBookExp::saveElementBin( FILE* stream, int iNumElement, strtContinuousExp sb)
{
    int numwritten;
    numwritten = fwrite( &iNumElement, sizeof(int), 1, stream );
    numwritten = fwrite( &sb, sizeof(strtContinuousExp), 1, stream );
}

//===============================================================
// Function: loadElementASCII
// Goal:
// Return:
//===============================================================

void CStructBookExp::loadElementASCII ( FILE* stream )
{
//	char str[_MAX_PATH];
	strtContinuousExp* aux;
	strtContinuousExp cExpStructure;

	float innerProduct,rho,xi,phase;
	int a,b;

	//fgets(str,84,stream);
	fscanf ( stream,"%i %f %f %f %f %i %i %*f %*f",	&m_iNumElement,
	         &innerProduct,
	         &rho,
	         &xi,
	         &phase,
	         &a,
	         &b );

	cExpStructure.innerProduct = innerProduct;
	cExpStructure.rho = rho;
	cExpStructure.xi =xi;
	cExpStructure.phase=phase;
	cExpStructure.a = a;
	cExpStructure.b = b;

	if ( m_structBook==NULL )
	{
		m_structBook = new strtContinuousExp[m_iNumElement];
		m_structBook[m_iNumElement-1] = cExpStructure;
	}
	else
	{
		aux = new strtContinuousExp[m_iNumElement-1];
		memcpy ( aux,m_structBook, ( m_iNumElement-1 ) *sizeof ( strtContinuousExp ) );
		delete [] m_structBook;
		m_structBook = new strtContinuousExp[m_iNumElement];
		memcpy ( m_structBook,aux, ( m_iNumElement-1 ) *sizeof ( strtContinuousExp ) );
		m_structBook[m_iNumElement-1] = cExpStructure;
		delete [] aux;
	}


}

//===============================================================
// Function: loadElementASCII
// Goal:
// Return:
//===============================================================

void CStructBookExp::loadElementASCII ( char* str )
{
	strtContinuousExp* aux;
	strtContinuousExp cExpStructure;

	float innerProduct,rho,xi,phase;
	int a,b;

	sscanf ( str,"%i %f %f %f %f %i %i %*f %*f",	&m_iNumElement,
	         &innerProduct,
	         &rho,
	         &xi,
	         &phase,
	         &a,
	         &b );

	cExpStructure.innerProduct = innerProduct;
	cExpStructure.rho = rho;
	cExpStructure.xi =xi;
	cExpStructure.phase=phase;
	cExpStructure.a = a;
	cExpStructure.b = b;

	if ( m_structBook==NULL )
	{
		m_structBook = new strtContinuousExp[m_iNumElement];
		m_structBook[m_iNumElement-1] = cExpStructure;
	}
	else
	{
		aux = new strtContinuousExp[m_iNumElement-1];
		memcpy ( aux,m_structBook, ( m_iNumElement-1 ) *sizeof ( strtContinuousExp ) );
		delete [] m_structBook;
		m_structBook = new strtContinuousExp[m_iNumElement];
		memcpy ( m_structBook,aux, ( m_iNumElement-1 ) *sizeof ( strtContinuousExp ) );
		m_structBook[m_iNumElement-1] = cExpStructure;
		delete [] aux;
	}


}

//===============================================================
// Function: adjustStructBook
// Goal:
// Return:
//===============================================================

void CStructBookExp::adjustStructBook ( int L )
{
	strtContinuousExp* aux;

	m_iNumElement = m_iNumElement-L;

	aux = new strtContinuousExp[m_iNumElement];

	memcpy ( aux,m_structBook,m_iNumElement*sizeof ( strtContinuousExp ) );

	if ( m_structBook!=NULL )
	{
		delete [] m_structBook;
	}

	m_structBook = new strtContinuousExp[m_iNumElement];

	memcpy ( m_structBook,aux,m_iNumElement*sizeof ( strtContinuousExp ) );

	delete [] aux;
}

//===============================================================
// Function: saveStructBook
// Goal:
// Return:
//===============================================================

void CStructBookExp::saveStructBook ( char* file )
{
	FILE* stream;
	int numwritten;

	if ( ( stream = fopen ( file, "w+t" ) ) != NULL )
	{
		numwritten = fwrite (	m_structBook,
		                      sizeof ( strtContinuousExp ),
		                      m_iNumElement,
		                      stream );
		fclose ( stream );
	}
	else
		cout << "Problem opening the file\n" << endl;
}

//===============================================================
// Function: saveStructBook
// Goal:
// Return:
//===============================================================

void CStructBookExp::saveStructBook ( FILE* stream )
{
	int numwritten;

	if ( m_structBook != NULL )
	{
		fwrite ( &m_iNumElement, sizeof ( int ), 1, stream );
		numwritten = fwrite ( m_structBook,
		                      sizeof ( strtContinuousExp ),
		                      m_iNumElement,
		                      stream );
		//fprintf(stream,"\n");
	}
	else
	{
		fwrite ( &m_iNumElement, sizeof ( int ), 1, stream );
		//cout << "Unable to save structure book !!!" << endl;
		cout << "   * Null samples !!!" << endl;
	}
}

//===============================================================
// Function: saveStructBookASCII
// Goal:
// Return:
//===============================================================

void CStructBookExp::saveStructBookASCII ( FILE* stream )
{

	if ( m_structBook != NULL )
	{
		fprintf ( stream,"N. Struct: %5i\n",m_iNumElement );
		for ( int i =0;i<m_iNumElement;i++ )
		{
			fprintf ( stream,"%14.8f %14.8f %14.8f %14.8f %5i %5i\n",	m_structBook[i].innerProduct,
			          m_structBook[i].rho,
			          m_structBook[i].xi,
			          m_structBook[i].phase,
			          m_structBook[i].a,
			          m_structBook[i].b );
		}
	}
	else
	{
		fprintf ( stream,"###### No Atoms ######\n" );
	}
}


//===============================================================
// Function: loadStructBook
// Goal:
// Return:
//===============================================================

void CStructBookExp::loadStructBook ( FILE* stream )
{
	int numread;
//	float temp_thres,temp_norm;
	//fscanf(stream,"%i\n",&m_iNumElement);
	//fscanf(stream,"%f\n",&temp_thres);
	//m_threshold = (double)temp_thres;
	//fscanf(stream,"%f\n",&temp_norm);
	//m_norm = (double)temp_norm;

	fread ( &m_iNumElement, sizeof ( int ), 1, stream );
	if ( m_iNumElement!=0 )
	{

		if ( m_structBook == NULL )
		{
			m_structBook = new strtContinuousExp [m_iNumElement];
		}

		numread = fread ( m_structBook,
		                  sizeof ( strtContinuousExp ),
		                  m_iNumElement,
		                  stream );
	}
	//fscanf(stream,"\n");

}

//===============================================================
// Function: printToScreen
// Goal:
// Return:
//===============================================================

void CStructBookExp::printToScreen()
{
	printf ( "Structures\n| Step  | Amp \t\t| freq \t\t| rho \t\t| phase \t| a \t| b \t|\n" ); //u \t\t|

	for ( int i=0; i<m_iNumElement; i++ )
	{
		printf ( "| %i \t| %f \t| %f \t| %f \t| %f \t| %i \t| %i \t|\n", //| %f \t
		         ( i+1 ),
		         m_structBook[i].innerProduct,
		         m_structBook[i].xi,//(float)(structure_cont[n-1].xi*f_amost/(2*pi)),
		         m_structBook[i].rho,
		         m_structBook[i].phase,
		         //m_structBook[i].u,
		         m_structBook[i].a,
		         m_structBook[i].b );
	}
}

//===============================================================
// Function: printElementToScreen
// Goal:
// Return:
//===============================================================

void CStructBookExp::printElementToScreen(int index)
{
    printf("Coef: %f | rho: %f | xi: %f | phase: %f |a: %d| b: %d|na: %d|pa: %d|origatom: %d|\n",
            m_structBook[index].innerProduct,
            m_structBook[index].rho,
            m_structBook[index].xi,
            m_structBook[index].phase,
            m_structBook[index].a,
            m_structBook[index].b,
            m_structBook[index].nextAtom,
            m_structBook[index].prevAtom,
            m_structBook[index].origAtomIndex);

}

// ================================================================

void CStructBookExp::convertCoefNegativeToPositive()
{
    int i;
    for(i=0; i< m_iNumElement; i++)
    {
        if( m_structBook[i].innerProduct < 0.0)
        {
            m_structBook[i].phase += pi;
            m_structBook[i].innerProduct = - m_structBook[i].innerProduct;
        } 
    
        if (m_structBook[i].phase >= (2*pi) )
        {
            m_structBook[i].phase -= 2*pi;
        }
    
        if (m_structBook[i].phase < 0 )
        {
            m_structBook[i].phase += 2*pi;
        }
    }
}

void CStructBookExp::findParmRange(double& min_coef,
                                double& max_coef,
                                double& min_rho,
                                double& max_rho,
                                double& min_phase,
                                double& max_phase)
{
    int i,j;
    double mincoef_aux, minrho_aux, minphase_aux;
    double maxcoef_aux, maxrho_aux, maxphase_aux;
    mincoef_aux     = 10000;
    minrho_aux      = 10000;
    minphase_aux    = 10000;
    maxcoef_aux     = -10000;
    maxrho_aux      = -10000;
    maxphase_aux    = -10000;
    for (i=0; i < m_iNumElement; i++)
    {
        if (mincoef_aux > m_structBook[i].innerProduct)
            mincoef_aux =  m_structBook[i].innerProduct;
        if (maxcoef_aux < m_structBook[i].innerProduct)
            maxcoef_aux =  m_structBook[i].innerProduct;
        if (minphase_aux >  m_structBook[i].phase)
            minphase_aux =  m_structBook[i].phase;
        if (maxphase_aux <  m_structBook[i].phase)
            maxphase_aux =  m_structBook[i].phase;
        if (minrho_aux >  fabs(m_structBook[i].rho))
            minrho_aux =  fabs(m_structBook[i].rho);
        if (maxrho_aux <  fabs(m_structBook[i].rho))
            maxrho_aux =  fabs(m_structBook[i].rho);
    }
    //min_coef = mincoef_aux;
    min_coef = 0.0;
    max_coef = maxcoef_aux;
    //min_rho  = minrho_aux;
    min_rho  = 0.0;
    max_rho  = maxrho_aux;
    //min_phase = minphase_aux;
    //max_phase = maxphase_aux;
    min_phase= 0.0;
    max_phase= 2*pi;

    //
    m_minamp = min_coef;
    m_maxamp = max_coef;
    m_minrho = min_rho;
    m_maxrho = max_rho;
    m_minphase = min_phase;
    m_maxphase = max_phase;
}

int CStructBookExp::findDeltaSupMax()
{
    int i;
    int deltaSupMax =0;
    for (i=0; i < m_iNumElement; i++)
    {
        if (deltaSupMax < (m_structBook[i].b - m_structBook[i].a) )
        {
            deltaSupMax =  (m_structBook[i].b - m_structBook[i].a);
        }
    }
    return deltaSupMax;
}


void CStructBookExp::sepByAmp(	strtContinuousExp* pSB,
                                int numElement,
                                double lowerAmpRangeLimit,
                                double upperAmpRangeLimit)
{
    double coef;
    int i,j;
    for(i=0; i<numElement; i++ )
    {
        coef = pSB[i].innerProduct;
        if ( (coef>lowerAmpRangeLimit) && (coef<=upperAmpRangeLimit) )
        {
            addElement(pSB[i]);  
        }
    }
}

void CStructBookExp::sepBySubBlock(	strtContinuousExp* pSB,
                                    int numElement,
                                    int lowerSubBlockLimit,
                                    int upperSubBlockLimit)
{
    int a;
    int i,j;
    for(i=0; i<numElement; i++ )
    {
        a = pSB[i].a;
        if ( (a>=lowerSubBlockLimit) && (a<upperSubBlockLimit) )
        {
            addElement(pSB[i]);  
        }
    }
}


//=========================================
// Quantization function
//=========================================

//===============================================================
// Function:
// Goal:
// Return:
//===============================================================

void CStructBookExp::setQuantConfig(    strtContinuousExp* pSB,
                                        int sbNumElement,
                                        double norm,
                                        int nbits_amp,
						                int nbits_rho,
						                int nbits_phase,
						                double freq1, // (Ffund/ Finit) (electric/audio)
						                double freq2, // (Fs/ Fend) (electric/audio)
						                int signalSize,
                                        int sigType)
{
    sbQHeadExp.norm = norm;
    
    sbQHeadExp.nbits_amp = nbits_amp;
    sbQHeadExp.nbits_rho = nbits_rho;
    sbQHeadExp.nbits_phase = nbits_phase;
    if (sigType == 1)
    {
	    double rf;
	    rf = RF_COEF * ( (freq2/2) / freq1);
	    sbQHeadExp.nbits_xi = log(rf) / log(2.0);
        sbQHeadExp.min_xi = freq1;
        sbQHeadExp.max_xi = freq2;
    }
    if (sigType == 2)
    {
	    sbQHeadExp.min_xi = freq1;
        sbQHeadExp.max_xi = freq2;
        int Nfreq = (int)(24 * ceil( log10(freq2/freq1)/log10(2.0) ) )+1;
        sbQHeadExp.nbits_xi = (int)ceil(log(Nfreq) / log(2.0));
    }
    sbQHeadExp.nbits_a = (int)ceil( log(signalSize) / log(2.0) );
    sbQHeadExp.nbits_b = sbQHeadExp.nbits_a;
	
	sbQHeadExp.max_amp =-10000;
	// Find max_amp 
    int i;
	for( i=0;i < sbNumElement; i++)
	{
		if (pSB[i].innerProduct > sbQHeadExp.max_amp) 
			sbQHeadExp.max_amp = pSB[i].innerProduct;
	}	
	sbQHeadExp.min_amp =  0;
	sbQHeadExp.max_rho = -100000;
	sbQHeadExp.min_rho =  100000;
	sbQHeadExp.max_phase = -2*pi; 
	sbQHeadExp.min_phase =  2*pi;

    double	nlevel_amp;
	nlevel_amp = pow(2.0, (double)sbQHeadExp.nbits_amp) - 1;
	double step_amp;
	step_amp = fabs( (sbQHeadExp.max_amp - sbQHeadExp.min_amp) / nlevel_amp );

	int k = 0;

	double innerProductQ = 0;

	// Find max_rho, min_rho, max_phase, min_phase and ajust xi
	for(i=0;i< sbNumElement;i++)
	{
		innerProductQ = (int)( (pSB[i].innerProduct + step_amp/2) / step_amp);
		
		if(innerProductQ != 0)
		{
			// rho
			if (fabs(pSB[i].rho) > sbQHeadExp.max_rho)
				sbQHeadExp.max_rho = fabs(pSB[i].rho);
			if (fabs(pSB[i].rho) < sbQHeadExp.min_rho)
				sbQHeadExp.min_rho = fabs(pSB[i].rho);
			// phase
			if (pSB[i].phase > sbQHeadExp.max_phase)
				sbQHeadExp.max_phase = pSB[i].phase;
			if (pSB[i].phase < sbQHeadExp.min_phase)
				sbQHeadExp.min_phase = pSB[i].phase;

			// xi ajustment
			if (pSB[i].xi < 0)	
				pSB[i].xi = - pSB[i].xi;
			if (pSB[i].xi > 2*pi)
				pSB[i].xi = pSB[i].xi - 2*pi;

			k++;
		}
	}

	m_iNumElement = k;
    sbQHeadExp.numberStruct = m_iNumElement;
	
	if (m_structBook == NULL)
	{
		m_structBook = new strtContinuousExp[m_iNumElement];
	}

}

//===============================================================
// Function: quantStructBook
// Goal:
// Return:
//===============================================================

void CStructBookExp::quantStructBook(   strtContinuousExp* pSB,
                                        int sbNumElement)
{
    // =========================================
	// Quantization steps definition
    double	nlevel_amp = pow(2.0, (double)sbQHeadExp.nbits_amp) - 1;
	double step_amp = fabs( (sbQHeadExp.max_amp - sbQHeadExp.min_amp) / nlevel_amp );

    double nlevel_rho = pow(2.0, (double)sbQHeadExp.nbits_rho) - 1;
	double step_rho = fabs( (sbQHeadExp.max_rho - sbQHeadExp.min_rho) / nlevel_rho );

    double nlevel_phase = pow(2.0, (double)sbQHeadExp.nbits_phase) - 1;
	double step_phase = fabs( (sbQHeadExp.max_phase - sbQHeadExp.min_phase) / nlevel_phase );
	

	//===================================================
	// Quantize structure book

    double innerProductQ;
    int i;
    int k =0;
    
	for(i=0; i < sbNumElement;i++)
	{
		
		innerProductQ = (int)( (pSB[i].innerProduct + step_amp/2) / step_amp) * step_amp;
		
		if(innerProductQ != 0)
		{
			// innerProduct
			m_structBook[k].innerProduct = innerProductQ;
			
			// rho
			m_structBook[k].rho = (int)((pSB[i].rho + step_rho/2)/step_rho) * step_rho;

            // xi
            m_structBook[k].xi = pSB[i].xi;
			
			// phase
			m_structBook[k].phase = (int)((pSB[i].phase + step_phase/2)/step_phase) * step_phase;

			// a
			m_structBook[k].a = pSB[i].a;

			// b
			m_structBook[k].b = pSB[i].b;

            k++;
		
		}
	}

}

//===============================================================
// Function: computeNumBits
// Goal:
// Return:
//===============================================================

int CStructBookExp::computeNumBits()
{
    int numBits = sbQHeadExp.numberStruct * (  sbQHeadExp.nbits_amp +
                                               sbQHeadExp.nbits_rho +
                                               sbQHeadExp.nbits_xi +
                                               sbQHeadExp.nbits_phase +
                                               sbQHeadExp.nbits_a +
                                               sbQHeadExp.nbits_b);
    return numBits;
}


//===============================================================
// Function: computeRate
// Goal:
// Return:
//===============================================================

double CStructBookExp::computeRate(int numSamples)
{
    double rate= (double)computeNumBits() / (double) numSamples;
    return rate;
}

//===============================================================
// Function: setOpCurve
// Goal:
// Return:
//=============================================================== 

void CStructBookExp::setOpCurve(strtOpCurveExp* opCurve,int numElement)
{
    if (m_opCurve==NULL) m_opCurve = new strtOpCurveExp[numElement];
    
    memcpy(m_opCurve,opCurve,sizeof(strtOpCurveExp)*numElement);
    
    m_numElemOpCurve = numElement;
}

void CStructBookExp::setMinAmp(double minamp)
{
    m_minamp = minamp;
} 

void  CStructBookExp::setMaxAmp(double maxamp)
{
    m_maxamp = maxamp;
}

void  CStructBookExp::setMinRho(double minrho)
{
    m_minrho = minrho;
}

void  CStructBookExp::setMaxRho(double maxrho)
{
    m_maxrho = maxrho;
}

void  CStructBookExp::setMinPhase(double minphase)
{
    m_minphase = minphase;
}

void  CStructBookExp::setMaxPhase(double maxphase)
{
    m_maxphase = maxphase;
}


//===============================================================
// Function: getStructBook()
// Goal:
// Return:
//===============================================================

strtContinuousExp*	CStructBookExp::getStructBook() const
{
	return m_structBook;
}

//===============================================================
// Function: getOpCurve()
// Goal:
// Return:
//===============================================================

strtOpCurveExp*	CStructBookExp::getOpCurve() const
{
	return m_opCurve;
}


//===============================================================
// Function: getNextAtomIndex
// Goal: 
// Return:
//===============================================================
int CStructBookExp::getNextAtomIndex ( int atomIndex)
{
	return m_structBook[atomIndex].nextAtom;
}

//===============================================================
// Function: getPrevAtomIndex
// Goal: 
// Return:
//===============================================================
int CStructBookExp::getPrevAtomIndex ( int atomIndex)
{
	return m_structBook[atomIndex].prevAtom;
}

//===============================================================
// Function: getOrigAtomIndex
// Goal: 
// Return:
//===============================================================
int CStructBookExp::getOrigAtomIndex ( int atomIndex)
{
	return m_structBook[atomIndex].origAtomIndex;
}

int CStructBookExp::getNumElemOpCurve()
{
    return m_numElemOpCurve;
}

double CStructBookExp::getMinAmp()
{
    return m_minamp;
}

double CStructBookExp::getMaxAmp()
{
    return m_maxamp;
}

double CStructBookExp::getMinRho()
{
    return m_minrho;
}

double CStructBookExp::getMaxRho()
{
    return m_maxrho;
}

double CStructBookExp::getMinPhase()
{
    return m_minphase;
}

double CStructBookExp::getMaxPhase()
{
    return m_maxphase;
}




/*
// ********************************************************
//  CStructBookQExp Class
// ********************************************************


//===============================================================
// Function: CStructBookQExp
// Goal: Constructor
// Return:
//===============================================================

CStructBookQExp::CStructBookQExp()
{
	m_structBookQ = NULL;
}

//===============================================================
// Function: ~CStructBookQExp
// Goal: Destructor
// Return:
//===============================================================

CStructBookQExp::~CStructBookQExp()
{
	if (m_structBookQ != NULL)
		delete [] m_structBookQ;
}

//===============================================================
// Function: setStructBookQ
// Goal: 
// Return:
//===============================================================

void CStructBookQExp::setStructBookQ(	strtStructBookQHeader sBookQHeader,
                                    strtStructBookQ* structBookQ,
                                    int iNumElementQ)
{
	structBookQuantizedHeader = sBookQHeader;

	m_iNumElementQ = iNumElementQ;

	if (m_structBookQ == NULL)
	{
		m_structBookQ = new strtStructBookQ [m_iNumElementQ];
	}
	memcpy( m_structBookQ,
            structBookQ,
            m_iNumElementQ * sizeof(strtStructBookQ));

}

//===============================================================
// Function: setStructBookQ
// Goal: 
// Return:
//===============================================================

void CStructBookQExp::setStructBookQ()
{
	m_iNumElementQ = structBookQuantizedHeader.numberStruct;
	if (m_structBookQ==NULL)
	{
		m_structBookQ = new strtStructBookQ[m_iNumElementQ];
	}
}

//===============================================================
// Function: saveStructBookQHeader
// Goal: 
// Return:
//===============================================================

void CStructBookQExp::saveStructBookQHeader(FILE* stream)
{
	int numwritten;

	numwritten = fwrite(	&structBookQuantizedHeader, 
				sizeof( strtStructBookQHeader ), 
				1, 
				stream );

	setNumElementQ(structBookQuantizedHeader.numberStruct);

	if (m_structBookQ == NULL)
	{
		m_structBookQ = new strtStructBookQ [m_iNumElementQ];
	}
}

//===============================================================
// Function: saveStructBookQHeader
// Goal: 
// Return:
//===============================================================

void CStructBookQExp::saveStructBookQHeader(CBitStream& bitStream)
{	
	strtFpBinForm aux;

	if (structBookQuantizedHeader.norm!=0)
	{	
		// chnull flag
		bitStream.writeBuffer(0,1);

		// norm
		aux = bitStream.double2bin(	structBookQuantizedHeader.norm,
						NORM_INT_NBITS, NORM_DEC_NBITS);
		bitStream.writeBuffer(	aux.integer,
					aux.int_nbits);
		bitStream.writeBuffer(	aux.decimal,
					aux.dec_nbits);
		
		// numberStruct			
		bitStream.writeBuffer(	structBookQuantizedHeader.numberStruct,
					NUM_STRUCT_NBITS);
		// max_amp	
		aux = bitStream.double2bin(	structBookQuantizedHeader.max_amp,
						AMP_INT_NBITS, AMP_DEC_NBITS);
		bitStream.writeBuffer(	aux.integer,
					aux.int_nbits);
		bitStream.writeBuffer(	aux.decimal,
					aux.dec_nbits);
		// nbits_amp
		bitStream.writeBuffer(	structBookQuantizedHeader.nbits_amp,
					AMP_2NBITS);
		// max_rho	
		aux = bitStream.double2bin(	structBookQuantizedHeader.max_rho,
						RHO_INT_NBITS, RHO_DEC_NBITS);
		bitStream.writeBuffer(	aux.integer,
					aux.int_nbits);
		bitStream.writeBuffer(	aux.decimal,
					aux.dec_nbits);
		// min_rho
		aux = bitStream.double2bin(	structBookQuantizedHeader.min_rho,
						RHO_INT_NBITS, RHO_DEC_NBITS);
		bitStream.writeBuffer(	aux.integer,
					aux.int_nbits);
		bitStream.writeBuffer(	aux.decimal,
					aux.dec_nbits);
		// nbits_rho
		bitStream.writeBuffer(	structBookQuantizedHeader.nbits_rho,
					RHO_2NBITS);
		// max_phase			
		aux = bitStream.double2bin(	structBookQuantizedHeader.max_phase,
						PHASE_INT_NBITS, PHASE_DEC_NBITS);
		bitStream.writeBuffer(	aux.integer,
					aux.int_nbits);
		bitStream.writeBuffer(	aux.decimal,
					aux.dec_nbits);
		// min_phase			
		aux = bitStream.double2bin(	structBookQuantizedHeader.min_phase,
						PHASE_INT_NBITS, PHASE_DEC_NBITS);
		bitStream.writeBuffer(	aux.integer,
					aux.int_nbits);
		bitStream.writeBuffer(	aux.decimal,
					aux.dec_nbits);
		// nbits_phase
		bitStream.writeBuffer(	structBookQuantizedHeader.nbits_phase,
					PHASE_2NBITS);
	}
	else
	{
		// chnull flag
		bitStream.writeBuffer(1,1);
	}
	
		
}


//===============================================================
// Function: saveStructBookQ
// Goal: 
// Return:
//===============================================================

void CStructBookQExp::saveStructBookQ(	CBitStream& bitStream,
					double Ffund,
					double Fs,
					int signalSize)
{

	double rf;
	//rf = RF_COEF * (Fs / Ffund);
	rf = RF_COEF * ( (Fs/2) / Ffund);
	double aux;
	aux = log(rf) / log(2.0);
	int nbits_xi = (int)(aux+1);
	int nbits_u_a_b = (int)ceil(((log10((double)(signalSize)))/(log10((double)(2)))));

	for (int i=0; i<m_iNumElementQ; i++)
	{
		// amp
		bitStream.writeBuffer(	m_structBookQ[i].innerProduct,
					structBookQuantizedHeader.nbits_amp);
		// rho signal
		bitStream.writeBuffer(	m_structBookQ[i].rhoSignal,
					1);
		// rho
		bitStream.writeBuffer(	m_structBookQ[i].rho,
					structBookQuantizedHeader.nbits_rho);
		// xi
		bitStream.writeBuffer(	m_structBookQ[i].xi,
					nbits_xi);
		// phase
		bitStream.writeBuffer(	m_structBookQ[i].phase,
					structBookQuantizedHeader.nbits_phase);
		// u
		//bitStream.writeBuffer(m_structBookQ[i].u,
		//			nbits_u_a_b);
		//
		// a
		bitStream.writeBuffer(	m_structBookQ[i].a,
					nbits_u_a_b);
		// b
		bitStream.writeBuffer(	m_structBookQ[i].b,
					nbits_u_a_b);
	}
}

//===============================================================
// Function: loadStructBookQHeader
// Goal: 
// Return:
//===============================================================

void CStructBookQExp::loadStructBookQHeader(FILE* stream)
{
	int numread;
	
	numread = fread(&structBookQuantizedHeader, 
			sizeof( strtStructBookQHeader ), 
			1, 
			stream );
	m_iNumElementQ = structBookQuantizedHeader.numberStruct;
	if (m_structBookQ == NULL)
	{
		m_structBookQ = new strtStructBookQ [m_iNumElementQ];
	}
}

//===============================================================
// Function: loadStructBookQHeader
// Goal: 
// Return:
//===============================================================

void CStructBookQExp::loadStructBookQHeader(CBitStream& bitStream)
{	
	strtFpBinForm aux;
	int chnull;
	
	// chnull flag
	bitStream.readBuffer(chnull,1);

	if ( chnull != 1)
	{
		// norm
		aux.int_nbits = NORM_INT_NBITS;
		aux.dec_nbits = NORM_DEC_NBITS;
		bitStream.readBuffer(	aux.integer,
					aux.int_nbits);
		bitStream.readBuffer(	aux.decimal,
					aux.dec_nbits);
		structBookQuantizedHeader.norm = bitStream.bin2double(aux);
		
		// numberStruct
		bitStream.readBuffer(	structBookQuantizedHeader.numberStruct,
					NUM_STRUCT_NBITS);
		
		// max_amp
		aux.int_nbits = AMP_INT_NBITS;
		aux.dec_nbits = AMP_DEC_NBITS;
		bitStream.readBuffer(	aux.integer,
					aux.int_nbits);
		bitStream.readBuffer(	aux.decimal,
					aux.dec_nbits);
		structBookQuantizedHeader.max_amp = bitStream.bin2double(aux);
		
		// nbits_amp
		bitStream.readBuffer(	structBookQuantizedHeader.nbits_amp,
					AMP_2NBITS);
		
		// max_rho
		aux.int_nbits = RHO_INT_NBITS;
		aux.dec_nbits = RHO_DEC_NBITS;
		bitStream.readBuffer(	aux.integer,
					aux.int_nbits);
		bitStream.readBuffer(	aux.decimal,
					aux.dec_nbits);
		structBookQuantizedHeader.max_rho = bitStream.bin2double(aux);
		
		// min_rho
		bitStream.readBuffer(	aux.integer,
					aux.int_nbits);
		bitStream.readBuffer(	aux.decimal,
					aux.dec_nbits);
		structBookQuantizedHeader.min_rho = bitStream.bin2double(aux);
		
		// nbits_rho
		bitStream.readBuffer(	structBookQuantizedHeader.nbits_rho,
					RHO_2NBITS);
		
		// max_phase	
		aux.int_nbits = PHASE_INT_NBITS;
		aux.dec_nbits = PHASE_DEC_NBITS;
		bitStream.readBuffer(	aux.integer,
					aux.int_nbits);
		bitStream.readBuffer(	aux.decimal,
					aux.dec_nbits);
		structBookQuantizedHeader.max_phase = bitStream.bin2double(aux);
		
		// min_phase
		bitStream.readBuffer(	aux.integer,
					aux.int_nbits);
		bitStream.readBuffer(	aux.decimal,
					aux.dec_nbits);
		structBookQuantizedHeader.min_phase = bitStream.bin2double(aux);
		
		// nbits_phase
		bitStream.readBuffer(	structBookQuantizedHeader.nbits_phase,
					PHASE_2NBITS);
	}
	else
	{
		structBookQuantizedHeader.numberStruct=0;
	}
	
//	cout << "Header reconstruido do canal:" << endl;
//	cout << structBookQuantizedHeader.norm << endl;
//	cout << structBookQuantizedHeader.numberStruct << endl;
//	cout << structBookQuantizedHeader.max_amp << endl;
//	cout << structBookQuantizedHeader.nbits_amp << endl;
//	cout << structBookQuantizedHeader.max_rho << endl;
//	cout << structBookQuantizedHeader.min_rho << endl;
//	cout << structBookQuantizedHeader.nbits_rho << endl;
//	cout << structBookQuantizedHeader.max_phase << endl;
//	cout << structBookQuantizedHeader.min_phase << endl;
//	cout << structBookQuantizedHeader.nbits_phase << endl;	
}


//===============================================================
// Function: loadStructBookQ
// Goal: 
// Return:
//===============================================================

void CStructBookQExp::loadStructBookQ(	CBitStream& bitStream,
					double Ffund,
					double Fs,
					int signalSize)
{
	if (structBookQuantizedHeader.numberStruct != 0)
	{

		double rf;
		rf = RF_COEF * ( (Fs/2) / Ffund);
		double aux;
		aux = log(rf) / log(2.0);
		int nbits_xi = (int)(aux+1);
		int nbits_u_a_b = (int)ceil(((log10((double)(signalSize)))/(log10((double)(2)))));

		for (int i=0; i<m_iNumElementQ; i++)
		{
			// amp
			bitStream.readBuffer(	m_structBookQ[i].innerProduct,
						structBookQuantizedHeader.nbits_amp);
			// rho signal
			bitStream.readBuffer(	m_structBookQ[i].rhoSignal,
						1);
			// rho
			bitStream.readBuffer(	m_structBookQ[i].rho,
						structBookQuantizedHeader.nbits_rho);
			// xi
			bitStream.readBuffer(	m_structBookQ[i].xi,
						nbits_xi);
			// phase
			bitStream.readBuffer(	m_structBookQ[i].phase,
						structBookQuantizedHeader.nbits_phase);
			
			// a
			bitStream.readBuffer(	m_structBookQ[i].a,
						nbits_u_a_b);
			// b
			bitStream.readBuffer(	m_structBookQ[i].b,
						nbits_u_a_b);
		}
	}
}

//===============================================================
// Function:  quantizeFloat
// Goal: 
// Return:
//===============================================================

float CStructBookQExp::quantizeFloat(float input, int int_nbits, int dec_nbits)
{
	float output;
	CBitStream bitStream;
	strtFpBinForm aux;
	
	if (input>=0)
	{
		aux = bitStream.float2bin(input,int_nbits,dec_nbits);
		output = bitStream.bin2float(aux);
	}
	else
	{
		aux = bitStream.float2bin(-input,int_nbits,dec_nbits);
		output = -bitStream.bin2float(aux);
	}

	return output;
}

//===============================================================
// Function:  quantizeDouble
// Goal: 
// Return:
//===============================================================

double CStructBookQExp::quantizeDouble(double input, int int_nbits, int dec_nbits)
{
	double output;
	CBitStream bitStream;
	strtFpBinForm aux;
	
	if (input>=0)
	{
		aux = bitStream.double2bin(input,int_nbits,dec_nbits);
		output = bitStream.bin2double(aux);
	}
	else
	{
		aux = bitStream.double2bin(-input,int_nbits,dec_nbits);
		output = -bitStream.bin2double(aux);
	}
	
	return output;
}

//===============================================================
// Function:  quantizeStructBookQHeader
// Goal:
// Return:
//===============================================================

void CStructBookQExp::quantizeStructBookQHeader()
{
	structBookQuantizedHeader.norm = quantizeDouble(structBookQuantizedHeader.norm,
							NORM_INT_NBITS,
							NORM_DEC_NBITS);
	structBookQuantizedHeader.max_amp = quantizeDouble(structBookQuantizedHeader.max_amp,
							AMP_INT_NBITS,
							AMP_DEC_NBITS);
	structBookQuantizedHeader.max_rho = quantizeDouble(structBookQuantizedHeader.max_rho,
							RHO_INT_NBITS,
							RHO_DEC_NBITS);
	structBookQuantizedHeader.min_rho = quantizeDouble(structBookQuantizedHeader.min_rho,
							RHO_INT_NBITS,
							RHO_DEC_NBITS);	
	structBookQuantizedHeader.max_phase = quantizeDouble(structBookQuantizedHeader.max_phase,
							PHASE_INT_NBITS,
							PHASE_DEC_NBITS);
	structBookQuantizedHeader.min_phase = quantizeDouble(structBookQuantizedHeader.min_phase,
							PHASE_INT_NBITS,
							PHASE_DEC_NBITS);
}



//===============================================================
// Function:  quantizeStructureBook
// Goal: squared amp
// Return:
//===============================================================

void CStructBookQExp::quantizeStructureBook(	CStructBook* structBook ,
						int nbits_amp,
						int nbits_rho,
						int nbits_phase,
						double Ffund,
						double Fs,
						int signalSize)
{
	//structBook.printToScreen();

	int i=0;
	double rf;
	rf = RF_COEF * (Fs / Ffund);
	double aux;
	aux = log(rf) / log(2.0);
	int nbits_xi = (int)(aux+1);
	//int nbits_u  = (int)(((log10((double)(signalSize)))/(log10((double)(2))))+0.5);
	int nbits_a  = (int)(((log10((double)(signalSize)))/(log10((double)(2))))+0.5);
	int nbits_b  = (int)(((log10((double)(signalSize)))/(log10((double)(2))))+0.5);

	//cout << "Number of bits:" << endl;
	//cout << "amp: " << nbits_amp << " rho: " << nbits_rho << " ph: " << nbits_phase << endl;
	//cout << "xi: " << nbits_xi << " a: " << nbits_a << " b: " << nbits_b << endl;

	double	nlevel_amp, nlevel_rho,//	nlevel_u,
		nlevel_xi, nlevel_phase;// nlevel_a,nlevel_b;

	nlevel_amp = pow(2.0, (double)nbits_amp) - 1;
	nlevel_rho = pow(2.0, (double)nbits_rho) - 1;
	//nlevel_u = pow(2.0, (double)nbits_u) - 1;
	nlevel_xi = pow(2.0, (double)nbits_xi) - 1;
	nlevel_phase = pow(2.0, (double)nbits_phase) - 1;
	//nlevel_a = pow(2.0, (double)nbits_a) - 1;
	//nlevel_b = pow(2.0, (double)nbits_b) - 1;

	//cout << "Number of levels:" << endl;
	//cout << "amp: " << nlevel_amp << " rho: " << nlevel_rho << " ph: " << nlevel_phase << endl;
	//cout << "xi: " << nlevel_xi << endl;
	
	// Maximum and minimum values definition
	double	max_amp, min_amp, max_rho, min_rho,
		max_xi, min_xi, max_phase, min_phase;
	
	max_amp = -100000; 
	//min_amp =  100000; 
	min_amp = 0;
	max_rho = -100000;
	min_rho =  100000;
	max_xi  =  2*pi;
	min_xi  =  0;
	max_phase = -2*pi; 
	min_phase =  2*pi;

	strtContinuousExp* pSBook;
	pSBook = new strtContinuousExp[structBook->getNumElement()];
	memcpy(	pSBook,
		structBook->getStructBook(),
		structBook->getNumElement() * sizeof(strtContinuousExp) );

	// Find max_amp and min_amp
	for( i=0;i < structBook->getNumElement(); i++)
    	{
		if (pow (pSBook[i].innerProduct, 2.0) > max_amp) 
			max_amp = pow (pSBook[i].innerProduct, 2.0);
		//if (pow (pSBook[i].innerProduct, 2.0) < min_amp) 
		//	min_amp = pow (pSBook[i].innerProduct, 2.0);
	} // end for
	

	double step_amp;
	step_amp = fabs( (max_amp - min_amp) / nlevel_amp );

	int k = 0;

	double innerProductQSqr = 0;

	// Find max_rho, min_rho, max_phase, min_phase, and ajust xi
	for(i=0;i< structBook->getNumElement();i++)
    	{
		innerProductQSqr = (int)( (pow (pSBook[i].innerProduct, 2.0) + step_amp/2) / step_amp);
		if(innerProductQSqr != 0)
		{
			// rho
			if ( fabs(pSBook[i].rho) > max_rho)
				max_rho = fabs(pSBook[i].rho);
			if ( fabs(pSBook[i].rho) < min_rho)
				min_rho = fabs(pSBook[i].rho);
			// phase
			if (pSBook[i].phase > max_phase)
				max_phase = pSBook[i].phase;
			if (pSBook[i].phase < min_phase)
				min_phase = pSBook[i].phase;

			// xi ajustment
			if (pSBook[i].xi < 0)	
				pSBook[i].xi = - pSBook[i].xi;
			if (pSBook[i].xi > 2*pi)
				pSBook[i].xi = pSBook[i].xi - 2*pi;

			k++;
		}
    	} // end for

	m_iNumElementQ = k;
	
	if (m_structBookQ == NULL)
	{
		m_structBookQ = new strtStructBookQ[m_iNumElementQ];
	}

	//cout << "Range: " << endl;
	//cout << "min amp: " << min_amp << " min rho: " << min_rho << " min ph: " << min_phase << endl;
	//cout << "max amp: " << max_amp << " max rho: " << max_rho << " max ph: " << max_phase << endl;

	// =========================================
	// Quantization steps definition
	// step_amp already defined above

	double step_rho;
	step_rho = fabs( (max_rho - min_rho) / nlevel_rho );
	//double step_u;
	//step_u = (double)signalSize / nlevel_u;	
	double step_xi;
	step_xi = (2*pi) / rf;
	double step_phase;
	step_phase = fabs( (max_phase - min_phase) / nlevel_phase );
	//double step_a;
	//step_a = (double)signalSize / nlevel_a;	
	//double step_b;
	//step_b = (double)signalSize / nlevel_b;	

	//cout << "Step: " << endl;
	//cout << "step amp: " << step_amp << " step rho: " << step_rho << " step ph: " << step_phase << endl;
	//cout << "step xi: " << step_xi << endl;

	//===================================================
	// Quantize structure book

	int min_rhoQ, min_phaseQ; // min_ampQ,

	// Obtaining offset
	if(min_rho >= 0)
		min_rhoQ = (int)(((min_rho) + step_rho/2)/step_rho);
	if(min_rho < 0)
		min_rhoQ = -(int)(-((min_rho) + step_rho/2)/step_rho);
	if(min_phase >= 0)
		min_phaseQ = (int)(((min_phase) + step_phase/2)/step_phase);
	if(min_phase < 0)
		min_phaseQ = -(int)(-((min_phase) + step_phase/2)/step_phase);
	
	k=0;

	for(i=0; i<structBook.getNumElement();i++)
    	{
		//if(pow (pSBook[i].innerProduct, 2.0) >= 0)
		//{
		innerProductQSqr = (int)( (pow (pSBook[i].innerProduct, 2.0) + step_amp/2) / step_amp);
		//}
		//if(pow (pSBook[i].innerProduct, 2.0)<0)
		//{
		//	innerProductQSqr = -(int)((-(pow (pSBook[i].innerProduct, 2.0))-step_amp/2)/step_amp);
		//}
		if(innerProductQSqr != 0)
		{
			// innerProduct - amp
			m_structBookQ[k].innerProduct = (int)innerProductQSqr;
			// rho signal
			if(pSBook[i].rho >= 0)
			{
				m_structBookQ[k].rhoSignal = 1;
			}
			else
			{
				m_structBookQ[k].rhoSignal = 0;
			}
			// rho
			m_structBookQ[k].rho = (int)((fabs(pSBook[i].rho) + step_rho/2)/step_rho);
			// u
			//m_structBookQ[k].u = (int)((pSBook[i].u + step_u/2)/step_u);
			// xi
			m_structBookQ[k].xi = (int)((pSBook[i].xi + step_xi/2)/step_xi);
			// phase
			//if(pSBook[i].phase >= 0)
			m_structBookQ[k].phase = (int)((pSBook[i].phase + step_phase/2)/step_phase);
			//if(pSBook[i].phase < 0)
			//	m_structBookQ[k].phase = -(int)((-pSBook[i].phase + step_phase/2)/step_phase);
			// a
			m_structBookQ[k].a = pSBook[i].a;
			// b
			m_structBookQ[k].b = pSBook[i].b;

			k++;
		}
    	} // end for

	// Indexing -> offset
	for(i=0; i<m_iNumElementQ;i++)
	{
		//m_structBookQ[i].innerProduct = m_structBookQ[i].innerProduct - min_ampQ;
		m_structBookQ[i].rho = m_structBookQ[i].rho - min_rhoQ;
		m_structBookQ[i].phase = m_structBookQ[i].phase - min_phaseQ;
	}

	// Fill header

	//structBookQuantizedHeader.signalSize = signalSize;
	//structBookQuantizedHeader.threshold = (float)structBook.getThreshold();
	structBookQuantizedHeader.norm = (float)structBook.getNorm();
	structBookQuantizedHeader.numberStruct = m_iNumElementQ;
	structBookQuantizedHeader.max_amp = (float)max_amp;
	//structBookQuantizedHeader.min_amp = min_amp;
	//structBookQuantizedHeader.step_amp = step_amp;
	structBookQuantizedHeader.nbits_amp = nbits_amp;
	structBookQuantizedHeader.max_rho = (float)max_rho;
	structBookQuantizedHeader.min_rho = (float)min_rho;
	//structBookQuantizedHeader.step_rho = step_rho;
	structBookQuantizedHeader.nbits_rho = nbits_rho;
	//structBookQuantizedHeader.step_xi = (float)step_xi;
	//structBookQuantizedHeader.nbits_xi = nbits_xi;
	structBookQuantizedHeader.max_phase = (float)max_phase;
	structBookQuantizedHeader.min_phase = (float)min_phase;
	//structBookQuantizedHeader.step_phase = step_phase;
	structBookQuantizedHeader.nbits_phase = nbits_phase;

	calcNumBits(Ffund,Fs,signalSize);

	delete [] pSBook;
}

//===============================================================
// Function: quantizeStructureBook
// Goal: linear amp 
// Return:
//===============================================================

void CStructBookQExp::quantizeStructureBook(	CStructBook& structBook,
						int nbits_amp,
						int nbits_rho,
						int nbits_phase,
						double Ffund,
						double Fs,
						int signalSize,
						int dummy)
{
	//cout << "Original Structure Book" << endl;
	//structBook.printToScreen();

	int i=0;
	double rf;
	rf = RF_COEF * ( (Fs/2) / Ffund);
	double aux;
	aux = log(rf) / log(2.0);
	//int nbits_xi = (int)(aux+1);
	//int nbits_u  = (int)(((log10((double)(signalSize)))/(log10((double)(2))))+0.5);
	//int nbits_a  = (int)(((log10((double)(signalSize)))/(log10((double)(2))))+0.5);
	//int nbits_b  = (int)(((log10((double)(signalSize)))/(log10((double)(2))))+0.5);


	double	nlevel_amp, nlevel_rho, nlevel_phase; //nlevel_xi,
	
	nlevel_amp = pow(2.0, (double)nbits_amp) - 1;
	nlevel_rho = pow(2.0, (double)nbits_rho) - 1;
	//nlevel_xi = pow(2.0, (double)nbits_xi) - 1;
	nlevel_phase = pow(2.0, (double)nbits_phase) - 1;
	
	// Maximum and minimum values definition
	double	max_amp, min_amp, max_rho, min_rho, max_phase, min_phase;

	strtContinuousExp* pSBook;
	pSBook = new strtContinuousExp[structBook.getNumElement()];
	memcpy(	pSBook,
		structBook.getStructBook(),
		structBook.getNumElement() * sizeof(strtContinuousExp) );
	
	max_amp =-10000;
	// Find max_amp 
	for( i=0;i < structBook.getNumElement(); i++)
	{
		if (pSBook[i].innerProduct > max_amp) 
			max_amp = pSBook[i].innerProduct;
	} // end for
	
	min_amp =  0;
	max_rho = -100000;
	min_rho =  100000;
	max_phase = -2*pi; 
	min_phase =  2*pi;
	
	double max_amp_quantized = quantizeDouble(max_amp, AMP_INT_NBITS,AMP_DEC_NBITS);

	double step_amp;
	step_amp = fabs( (max_amp_quantized - min_amp) / nlevel_amp );

	int k = 0;

	double innerProductQ = 0;

	// Find max_rho, min_rho, max_phase, min_phase and ajust xi
	for(i=0;i< structBook.getNumElement();i++)
	{
		if (pSBook[i].innerProduct==max_amp)
		{
			innerProductQ = (int)( (max_amp_quantized + step_amp/2) / step_amp);
		}
		else
		{
			innerProductQ = (int)( (pSBook[i].innerProduct + step_amp/2) / step_amp);
		}
		
		if(innerProductQ != 0)
		{
			// rho
			if (fabs(pSBook[i].rho) > max_rho)
				max_rho = fabs(pSBook[i].rho);
			if (fabs(pSBook[i].rho) < min_rho)
				min_rho = fabs(pSBook[i].rho);
			// phase
			if (pSBook[i].phase > max_phase)
				max_phase = pSBook[i].phase;
			if (pSBook[i].phase < min_phase)
				min_phase = pSBook[i].phase;

			// xi ajustment
			if (pSBook[i].xi < 0)	
				pSBook[i].xi = - pSBook[i].xi;
			if (pSBook[i].xi > 2*pi)
				pSBook[i].xi = pSBook[i].xi - 2*pi;

			k++;
		}
	} // end for

	m_iNumElementQ = k;
	
	if (m_structBookQ == NULL)
	{
		m_structBookQ = new strtStructBookQ[m_iNumElementQ];
	}

	double max_rho_quantized = quantizeDouble(max_rho, RHO_INT_NBITS,RHO_DEC_NBITS);
	double min_rho_quantized = quantizeDouble(min_rho, RHO_INT_NBITS,RHO_DEC_NBITS);
	double max_phase_quantized = quantizeDouble(max_phase, PHASE_INT_NBITS,PHASE_DEC_NBITS);
	double min_phase_quantized = quantizeDouble(min_phase, PHASE_INT_NBITS,PHASE_DEC_NBITS);
	
	// =========================================
	// Quantization steps definition
	// step_amp already defined above

	double step_rho;
	step_rho = fabs( (max_rho_quantized - min_rho_quantized) / nlevel_rho );
	double step_xi;
	step_xi = (pi) / rf;
	double step_phase;
	step_phase = fabs( (max_phase_quantized - min_phase_quantized) / nlevel_phase );
	

	//===================================================
	// Quantize structure book

	k=0;

	for(i=0; i<structBook.getNumElement();i++)
	{
		if (pSBook[i].innerProduct==max_amp)
		{
			innerProductQ = (int)( (max_amp_quantized + step_amp/2) / step_amp);
		}
		else
		{
			innerProductQ = (int)( (pSBook[i].innerProduct + step_amp/2) / step_amp);
		}
		
		if(innerProductQ != 0)
		{
			// innerProduct
			m_structBookQ[k].innerProduct = (int)innerProductQ;
			// rho signal
			if(pSBook[i].rho >= 0)
			{
				m_structBookQ[k].rhoSignal = 1;
			}
			else
			{
				m_structBookQ[k].rhoSignal = 0;
			}
			// rho
			if (pSBook[i].rho==max_rho)
			{
				m_structBookQ[k].rho = (int)((fabs(max_rho_quantized) + step_rho/2)/step_rho);
			}
			else if (pSBook[i].rho==min_rho)
			{
				m_structBookQ[k].rho = (int)((fabs(min_rho_quantized) + step_rho/2)/step_rho);
			}
			else
			{
				m_structBookQ[k].rho = (int)((fabs(pSBook[i].rho) + step_rho/2)/step_rho);
			}
			// xi
			m_structBookQ[k].xi = (int)((pSBook[i].xi + step_xi/2)/step_xi);
			// phase
			if (pSBook[i].phase==max_phase)
			{
				m_structBookQ[k].phase = (int)((fabs(max_phase_quantized) + step_phase/2)/step_phase);
			}
			else if (pSBook[i].phase==min_phase)
			{
				m_structBookQ[k].phase = (int)((fabs(min_phase_quantized) + step_phase/2)/step_phase);
			}
			else
			{
				m_structBookQ[k].phase = (int)((pSBook[i].phase + step_phase/2)/step_phase);
			}
			// a
			m_structBookQ[k].a = pSBook[i].a;
			// b
			m_structBookQ[k].b = pSBook[i].b;

			k++;
		}
	} // end for

	int min_rhoQ, min_phaseQ; // min_ampQ,

	if(min_rho >= 0)
		min_rhoQ = (int)(((min_rho) + step_rho/2)/step_rho);
	if(min_rho < 0)
		min_rhoQ = -(int)(-((min_rho) + step_rho/2)/step_rho);
	if(min_phase >= 0)
		min_phaseQ = (int)(((min_phase) + step_phase/2)/step_phase);
	if(min_phase < 0)
		min_phaseQ = -(int)(-((min_phase) + step_phase/2)/step_phase);
	
	// Indexing -> offset
	for(i=0; i<m_iNumElementQ;i++)
	{
		m_structBookQ[i].rho = m_structBookQ[i].rho - min_rhoQ;
		m_structBookQ[i].phase = m_structBookQ[i].phase - min_phaseQ;
	}

	// Fill header

	
	structBookQuantizedHeader.norm = structBook.getNorm();
	structBookQuantizedHeader.numberStruct = m_iNumElementQ;
	structBookQuantizedHeader.max_amp = max_amp;
	structBookQuantizedHeader.nbits_amp = nbits_amp;
	structBookQuantizedHeader.max_rho = max_rho;
	structBookQuantizedHeader.min_rho = min_rho;
	structBookQuantizedHeader.nbits_rho = nbits_rho;
	structBookQuantizedHeader.max_phase = max_phase;
	structBookQuantizedHeader.min_phase = min_phase;
	structBookQuantizedHeader.nbits_phase = nbits_phase;
	
//	cout << "Header do canal:" << endl;
//	cout << structBookQuantizedHeader.norm << endl;
//	cout << structBookQuantizedHeader.numberStruct << endl;
//	cout << structBookQuantizedHeader.max_amp << endl;
//	cout << structBookQuantizedHeader.nbits_amp << endl;
//	cout << structBookQuantizedHeader.max_rho << endl;
//	cout << structBookQuantizedHeader.min_rho << endl;
//	cout << structBookQuantizedHeader.nbits_rho << endl;
//	cout << structBookQuantizedHeader.max_phase << endl;
//	cout << structBookQuantizedHeader.min_phase << endl;
//	cout << structBookQuantizedHeader.nbits_phase << endl;

	calcNumBits(Ffund,Fs,signalSize);

	delete [] pSBook;
}



//===============================================================
// Function: dequantizeStructureBook
// Goal: squared amp
// Return:
//===============================================================
void CStructBookQExp::dequantizeStructureBook(	CStructBook& structBook, 
						double Ffund,
						double Fs,
						int signalSize)
{	
	
	structBook.setNorm( (double)structBookQuantizedHeader.norm);
	
	double	nlevel_amp, nlevel_rho, nlevel_phase;

	nlevel_amp = pow(2.0, (double)structBookQuantizedHeader.nbits_amp) -1;
	nlevel_rho = pow(2.0, (double)structBookQuantizedHeader.nbits_rho) -1;
	nlevel_phase = pow(2.0, (double)structBookQuantizedHeader.nbits_phase) -1;
	

	double step_amp;
	step_amp = fabs( (double)(structBookQuantizedHeader.max_amp - 0) / nlevel_amp );
	double step_rho;
	step_rho = fabs( (double)(structBookQuantizedHeader.max_rho - structBookQuantizedHeader.min_rho) / nlevel_rho );
	double rf;
	rf = RF_COEF * ( (Fs/2) / Ffund);
	double step_xi;
	//step_xi = (2*pi) / rf;
	step_xi = (pi) / rf;
	double step_phase;
	step_phase = fabs( (double)(structBookQuantizedHeader.max_phase - structBookQuantizedHeader.min_phase) / nlevel_phase );

	strtContinuousExp* pSBook;
	pSBook = new strtContinuousExp [structBookQuantizedHeader.numberStruct]; 

	for (int i=0; i<structBookQuantizedHeader.numberStruct; i++)
	{
		pSBook[i].innerProduct = sqrt( (double)(m_structBookQ[i].innerProduct * step_amp) );
		pSBook[i].rho = (double)(m_structBookQ[i].rho * step_rho) + structBookQuantizedHeader.min_rho;
		if (m_structBookQ[i].rhoSignal==0)
		{
			pSBook[i].rho = - pSBook[i].rho;
		}
		pSBook[i].xi = (double)(m_structBookQ[i].xi * step_xi);
		pSBook[i].phase = (double)(m_structBookQ[i].phase * step_phase) + structBookQuantizedHeader.min_phase;
		pSBook[i].a = m_structBookQ[i].a;
		pSBook[i].b = m_structBookQ[i].b;

		structBook.addElement(pSBook[i]);
	}
	//structBook.printToScreen();
	

	delete [] pSBook;
	
}

//===============================================================
// Function: dequantizeStructureBook
// Goal: linear amp
// Return:
//===============================================================
void CStructBookQExp::dequantizeStructureBook(	CStructBook& structBook,
						double Ffund,
						double Fs,
						int signalSize,
						int dummy)
{	
	structBook.setNorm( (double)structBookQuantizedHeader.norm);

	double	nlevel_amp, nlevel_rho, nlevel_phase;

	nlevel_amp = pow(2.0, (double)structBookQuantizedHeader.nbits_amp) -1;
	nlevel_rho = pow(2.0, (double)structBookQuantizedHeader.nbits_rho) -1;
	nlevel_phase = pow(2.0, (double)structBookQuantizedHeader.nbits_phase) -1;
	
	double step_amp;
	step_amp = fabs( (double)(structBookQuantizedHeader.max_amp - 0) / nlevel_amp );
	double step_rho;
	step_rho = fabs( (double)(structBookQuantizedHeader.max_rho - structBookQuantizedHeader.min_rho) / nlevel_rho );
	double rf;
	rf = RF_COEF * ( (Fs/2) / Ffund);
	double step_xi;
	step_xi = pi / rf;
	double step_phase;
	step_phase = fabs( (double)(structBookQuantizedHeader.max_phase - structBookQuantizedHeader.min_phase) / nlevel_phase );

	strtContinuousExp* pSBook;
	pSBook = new strtContinuousExp [structBookQuantizedHeader.numberStruct]; 

	for (int i=0; i<structBookQuantizedHeader.numberStruct; i++)
	{
		pSBook[i].innerProduct =  (double)(m_structBookQ[i].innerProduct * step_amp);
		pSBook[i].rho = (double)(m_structBookQ[i].rho * step_rho) + structBookQuantizedHeader.min_rho;
		if (m_structBookQ[i].rhoSignal==0)
		{
			pSBook[i].rho = - pSBook[i].rho;
		}
		pSBook[i].xi = (double)(m_structBookQ[i].xi * step_xi);
		pSBook[i].phase = (double)(m_structBookQ[i].phase * step_phase) + structBookQuantizedHeader.min_phase;
		pSBook[i].a = m_structBookQ[i].a;
		pSBook[i].b = m_structBookQ[i].b;

		structBook.addElement(pSBook[i]);
	}
	//cout << "Quantized Structure Book" << endl;
	//structBook.printToScreen();
	
	delete [] pSBook;
}




//===============================================================
// Function: getSBQHeader()
// Goal: 
// Return:
//===============================================================

strtStructBookQHeader CStructBookQExp::getSBQHeader() const
{
	return structBookQuantizedHeader;
}

//===============================================================
// Function: getStructBookQ()
// Goal: 
// Return:
//===============================================================

strtStructBookQ* CStructBookQExp::getStructBookQ() const
{
	return m_structBookQ;
}

//===============================================================
// Function: getNumElementQ()
// Goal: 
// Return:
//===============================================================

int	CStructBookQExp::getNumElementQ() const
{
	return m_iNumElementQ;
}

//===============================================================
// Function: getNumBits()
// Goal: 
// Return:
//===============================================================

int	CStructBookQExp::getNumBits() const
{
	return m_iNumBits;
}

//===============================================================
// Function: getNorm()
// Goal: 
// Return:
//===============================================================

double	CStructBookQExp::getNorm() const
{
	return structBookQuantizedHeader.norm;
}


//===============================================================
// Function: setNumElementQ()
// Goal: 
// Return:
//===============================================================

void CStructBookQExp::setNumElementQ(int num)
{

	m_iNumElementQ = num;
}


//===============================================================
// Function: calcNumBits()
// Goal: 
// Return:
//===============================================================

void CStructBookQExp::calcNumBits( double Ffund,
                                double Fs,
                                int signalSize)
{
	double rf;
	rf = RF_COEF * ( (Fs/2) / Ffund);
	double aux;
	aux = log(rf) / log(2.0);
	
	//rf = RF_COEF * (Fs / Ffund);
	//double aux;
	//aux = log(rf) / log(2.0);
	
	int nbits_xi = (int)(aux+1);
	int nbits_a_b = (int)ceil(((log10((double)(signalSize)))/(log10(2.0))));
	int nbits_rho_sig = 1 ;
	int nbits_header = 	NORM_INT_NBITS  + NORM_DEC_NBITS + NUM_STRUCT_NBITS + AMP_INT_NBITS + 
				AMP_DEC_NBITS   + 2*(RHO_INT_NBITS  + RHO_DEC_NBITS) + 
				2*(PHASE_INT_NBITS + PHASE_DEC_NBITS) + AMP_2NBITS + RHO_2NBITS + 
				PHASE_2NBITS;
	int nbits_chnull = 1;


	if (structBookQuantizedHeader.numberStruct!=0)
	{
		m_iNumBits =	nbits_chnull + 
						nbits_header +
						structBookQuantizedHeader.numberStruct * (
						structBookQuantizedHeader.nbits_amp + 
						nbits_rho_sig +
						structBookQuantizedHeader.nbits_rho + 
						structBookQuantizedHeader.nbits_phase +
						nbits_xi + 
						2 * nbits_a_b );
	}
	else
	{
		m_iNumBits =	nbits_chnull; 
	}
					
}

//===============================================================
// Function: operator= 
// Goal: 
// Return:
//===============================================================

const CStructBookQExp& CStructBookQ::operator=(const CStructBookQ& structBookQ)
{
	if (structBookQ.getStructBookQ() != NULL)
	{
		this->structBookQuantizedHeader = structBookQ.getSBQHeader();

		this->m_iNumElementQ = structBookQ.getNumElementQ();
		if (this->m_structBookQ == NULL)
		{
			this->m_structBookQ = new strtStructBookQ[this->m_iNumElementQ];
		}

		memcpy(	this->m_structBookQ,
				structBookQ.getStructBookQ(),
				(this->m_iNumElementQ) * sizeof(strtStructBookQ));	
		
		m_iNumBits = structBookQ.getNumBits();
		
	}
	else
	{
		this->m_iNumElementQ = 0;
		this->m_structBookQ = NULL;
	}

	return *this;
}
*/

// -------------------------------------------------------------------------------------------

void separateByAmp(CStructBook* sbIn,CStructBook* sbOut, int NAmpRange, double* ampRangeLimit )
{
    strtContinuousExp* pSBIn;
    pSBIn = ((CStructBookExp*)sbIn)->getStructBook();
    int numElementIn = ((CStructBookExp*)sbIn)->getNumElement();
    double coef;
    int i,j;
    for(i=0; i<numElementIn; i++ )
    {
        coef = pSBIn[i].innerProduct;
        for(j=0; j<NAmpRange; j++ )
        {
            if ( (coef>ampRangeLimit[j]) && (coef<=ampRangeLimit[j+1]) )
            {
                ((CStructBookExp*)sbOut)[j].addElement(pSBIn[i]);  
            }
        }
    }
}

