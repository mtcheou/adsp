//dictionary.cpp

#include "dictionary.h"


// ========================================
//  CDictionary Class
// ========================================

//===============================================================
// Function: Contructors
// Goal:
// Return:
//===============================================================
CDictionary::CDictionary()
{
    m_signalSize = 0;
    m_complexAtom = NULL;
    m_realAtom = NULL;
    m_rAtomNorm = 0.0;
}


//===============================================================
// Function: Destructor
// Goal:
// Return:
//===============================================================

CDictionary::~CDictionary()
{
    if (m_complexAtom) delete [] m_complexAtom;
    if (m_realAtom) delete [] m_realAtom;
}

double CDictionary::getRAtomNorm()
{
    return m_rAtomNorm;
}

// ========================================
//  CExpDictionary Class
// ========================================


// Interface
void CExpDictionary::setSignalSize(int signalSize)
{
    m_signalSize = signalSize;
    if (m_complexAtom) delete [] m_complexAtom;
    if (m_realAtom) delete [] m_realAtom;
    m_complexAtom = new CComplex[signalSize];
    m_realAtom = new double[signalSize];
}

void CExpDictionary::setComplexAtom(CParameter* parm)
{
    int i;
    int u;
    double real_peak=0;
    double imag_peak=0;

    double rho = ((CExpParm*)parm)->rho;
    double xi = ((CExpParm*)parm)->xi;
    int a = ((CExpParm*)parm)->a;
    int b = ((CExpParm*)parm)->b;

    if (rho>=0)
    {
        u = a;
        for( i=0;i<m_signalSize; i++)
        {
            m_complexAtom[i] = CComplex(0,0,1);
            if(i>= u)
            {
                if(xi!=0)
                {
                    m_complexAtom[i] = CComplex( exp( (-rho*(double)( (double)i-u )) ) * ( cos( (xi*(double)i) ) ), // sqrt(2*rho),//(double)signalSize,
                                        exp( (-rho*(double)( (double)i-u )) ) * ( sin( (xi*(double)i) ) ), // sqrt(2*rho),//(double)signalSize,
                                        1);
                }
                else
                {
                    m_complexAtom[i] = CComplex( exp( (-rho*(double)( (double)i-u )) ), // sqrt(2*rho), // (double)signalSize,
                                        0.0,//exp( (-rho*(double)( (double)i-u )) ),// sqrt(2*rho), // (double)signalSize,
                                        1);
                }
            }

            if ( fabs(real_peak) < fabs(m_complexAtom[i].Real()) ) real_peak = m_complexAtom[i].Real();
            if ( fabs(imag_peak) < fabs(m_complexAtom[i].Imag()) ) imag_peak = m_complexAtom[i].Imag();
        }
    }
    else
    {

        u = m_signalSize - 1 - b;
        for( i=0;i<m_signalSize; i++)
        {
            m_complexAtom[m_signalSize-1-i] = CComplex(0,0,1);
            if (i>=u)
            {
                if(xi!=0)
                {
                    m_complexAtom[m_signalSize-1-i] =   CComplex( exp( (rho*(double)( (double)i-u )) ) * ( cos( (xi*(double)(m_signalSize-1-i)) ) ), // sqrt(2*rho),//(double)signalSize,
                                            exp( (rho*(double)( (double)i-u )) ) * ( sin( (xi*(double)(m_signalSize-1-i)) ) ), // sqrt(2*rho),//(double)signalSize,
                                            1);
                }
                else
                {
                    m_complexAtom[m_signalSize-1-i] =   CComplex( exp( (rho*(double)( (double)i-u )) ), // sqrt(2*rho), // (double)signalSize,
                                            0.0,//exp( (rho*(double)( (double)i-u )) ),// sqrt(2*rho), // (double)signalSize,
                                            1);
                }
            }
            if ( fabs(real_peak) < fabs(m_complexAtom[m_signalSize-1-i].Real()) ) real_peak = m_complexAtom[m_signalSize-1-i].Real();
            if ( fabs(imag_peak) < fabs(m_complexAtom[m_signalSize-1-i].Imag()) ) imag_peak = m_complexAtom[m_signalSize-1-i].Imag();
        }
    }

    for(i = 0; i < a; i++)
    {
        m_complexAtom[i] = CComplex(0,0,1);
    }

    for(i = m_signalSize - 1; i > b; i--)
    {
        m_complexAtom[i] = CComplex(0,0,1);
    }

/*  if ( fabs(real_peak) < 1e-10 )
    {
        for( i=0; i < signalSize; i++)
        {
            atom[i].setReal(0);
        }
    }
    if ( fabs(imag_peak) < 1e-10 )
    {
        for( i=0; i < signalSize; i++)
        {
            atom[i].setImag(0);
        }
    }
*/

}
void CExpDictionary::setRealAtom(CParameter* parm)
{
    int i;
    double norm;
    int u;
    double peak=0;

    double rho = ((CExpParm*)parm)->rho;
    double xi = ((CExpParm*)parm)->xi;
    double phi = ((CExpParm*)parm)->phase;
    int a = ((CExpParm*)parm)->a;
    int b = ((CExpParm*)parm)->b;

    if (rho>=0)
    {
        u = a;
        for(i=0; i < m_signalSize; i++)
        {
            m_realAtom[i]=0;
            if(i >= u)
            {
                if(xi!=0)
                    m_realAtom[i] = exp( -rho * (double)(i-u) ) * cos ( (xi*(double)i) + phi );
                else
                    m_realAtom[i] = exp( -rho * (double)(i-u) ) * cos ( phi );
            }

            if ( fabs(peak) < fabs(m_realAtom[i]) ) peak = m_realAtom[i];
        }
    }
    else
    {
        u = m_signalSize - 1 - b;
        for(i=0; i < m_signalSize; i++)
        {
            m_realAtom[m_signalSize-1-i]=0;
            if( i >= u )
            {
                if(xi!=0)
                    m_realAtom[m_signalSize-1-i] = exp( rho * (double)(i-u) ) * cos ( (xi*(double)(m_signalSize-1-i)) + phi );
                else
                    m_realAtom[m_signalSize-1-i] = exp( rho * (double)(i-u) ) * cos ( phi );
            }

            if ( fabs(peak) < fabs(m_realAtom[m_signalSize-1-i]) ) peak = m_realAtom[m_signalSize-1-i];
        }
    }

    if (fabs(peak) < 1e-10)
    {
        for(i = 0; i < m_signalSize ; i++)
        {
            m_realAtom[i] = 0;
        }
    }
    else
    {
        //========================
        for(i = 0; i < a; i++)
        {
            m_realAtom[i] = 0;
        }

        for(i = m_signalSize - 1; i > b; i--)
        {
            m_realAtom[i] = 0;
        }
        // ======================
        // Calculate norm
        norm = 0;
        for(i=0; i < m_signalSize; i++)
        {
            norm += m_realAtom[i] * m_realAtom[i];
        }
        norm = sqrt(norm);

        m_rAtomNorm = norm;

        // ======================
        if ( norm == 0 ) norm = 1;

        // ======================

        for(i=0; i<m_signalSize; i++)
        {
            m_realAtom[i] = m_realAtom[i] / norm;
        }

    }
}

void CExpDictionary::setRealAtom(strtContinuousExp sb)
{
    int i;
    double norm;
    int u;
    double peak=0;

    double rho = sb.rho;
    double xi = sb.xi;
    double phi = sb.phase;
    int a = sb.a;
    int b = sb.b;

    if (rho>=0)
    {
        u = a;
        for(i=0; i < m_signalSize; i++)
        {
            m_realAtom[i]=0;
            if(i >= u)
            {
                if(xi!=0)
                    m_realAtom[i] = exp( -rho * (double)(i-u) ) * cos ( (xi*(double)i) + phi );
                else
                    m_realAtom[i] = exp( -rho * (double)(i-u) ) * cos ( phi );
            }

            if ( fabs(peak) < fabs(m_realAtom[i]) ) peak = m_realAtom[i];
        }
    }
    else
    {
        u = m_signalSize - 1 - b;
        for(i=0; i < m_signalSize; i++)
        {
            m_realAtom[m_signalSize-1-i]=0;
            if( i >= u )
            {
                if(xi!=0)
                    m_realAtom[m_signalSize-1-i] = exp( rho * (double)(i-u) ) * cos ( (xi*(double)(m_signalSize-1-i)) + phi );
                else
                    m_realAtom[m_signalSize-1-i] = exp( rho * (double)(i-u) ) * cos ( phi );
            }

            if ( fabs(peak) < fabs(m_realAtom[m_signalSize-1-i]) ) peak = m_realAtom[m_signalSize-1-i];
        }
    }

    if (fabs(peak) < 1e-10)
    {
        for(i = 0; i < m_signalSize ; i++)
        {
            m_realAtom[i] = 0;
        }
    }
    else
    {
        //========================
        for(i = 0; i < a; i++)
        {
            m_realAtom[i] = 0;
        }

        for(i = m_signalSize - 1; i > b; i--)
        {
            m_realAtom[i] = 0;
        }
        // ======================
        // Calculate norm
        norm = 0;
        for(i=0; i < m_signalSize; i++)
        {
            norm += m_realAtom[i] * m_realAtom[i];
        }
        norm = sqrt(norm);

        m_rAtomNorm = norm;

        // ======================
        if ( norm == 0 ) norm = 1;

        // ======================

        for(i=0; i<m_signalSize; i++)
        {
            m_realAtom[i] = m_realAtom[i] / norm;
        }

    }
}

void CExpDictionary::setRealAtomOnSupport(CParameter* parm)
{
    int i;
    double norm;
    int u;
    double peak=0;

    double rho = ((CExpParm*)parm)->rho;
    double xi = ((CExpParm*)parm)->xi;
    double phi = ((CExpParm*)parm)->phase;
    int a = ((CExpParm*)parm)->a;
    int b = ((CExpParm*)parm)->b;

    if (rho>=0)
    {
        u = a;
        for(i=a; i <=b; i++)
        {
            if(xi!=0)
                m_realAtom[i] = exp( -rho * (double)(i-u) ) * cos ( (xi*(double)i) + phi );
            else
                m_realAtom[i] = exp( -rho * (double)(i-u) ) * cos ( phi );

            if ( fabs(peak) < fabs(m_realAtom[i]) ) peak = m_realAtom[i];
        }
    }
    else
    {
        u = m_signalSize - 1 - b;
        for(i=a; i <= b; i++)
        {
            if(xi!=0)
                m_realAtom[i] = exp( rho * (double)(m_signalSize-1-i-u) ) * cos ( (xi*(double)(i)) + phi );
            else
                m_realAtom[i] = exp( rho * (double)(m_signalSize-1-i-u) ) * cos ( phi );

            if ( fabs(peak) < fabs(m_realAtom[i]) ) peak = m_realAtom[i];
        }
    }

    if (fabs(peak) < 1e-10)
    {
        for(i = a; i <= b ; i++)
        {
            m_realAtom[i] = 0;
        }
    }
    else
    {
        // ======================
        // Calculate norm
        norm = 0;
        for(i=a; i <= b; i++)
        {
            norm += m_realAtom[i] * m_realAtom[i];
        }
        norm = sqrt(norm);

        m_rAtomNorm = norm;

        // ======================
        if ( norm == 0 ) norm = 1;

        // ======================

        for(i=a; i<=b; i++)
        {
            m_realAtom[i] = m_realAtom[i] / norm;
        }

    }
}

void CExpDictionary::setRealAtomOnSupport(strtContinuousExp sb)
{
    int i;
    double norm;
    int u;
    double peak=0;

    double rho = sb.rho;
    double xi = sb.xi;
    double phi = sb.phase;
    int a = sb.a;
    int b = sb.b;

    if (rho>=0)
    {
        u = a;
        for(i=a; i <=b; i++)
        {
            if(xi!=0)
                m_realAtom[i] = exp( -rho * (double)(i-u) ) * cos ( (xi*(double)i) + phi );
            else
                m_realAtom[i] = exp( -rho * (double)(i-u) ) * cos ( phi );

            if ( fabs(peak) < fabs(m_realAtom[i]) ) peak = m_realAtom[i];
        }
    }
    else
    {
        u = m_signalSize - 1 - b;
        for(i=a; i <= b; i++)
        {
            if(xi!=0)
                m_realAtom[i] = exp( rho * (double)(m_signalSize-1-i-u) ) * cos ( (xi*(double)(i)) + phi );
            else
                m_realAtom[i] = exp( rho * (double)(m_signalSize-1-i-u) ) * cos ( phi );

            if ( fabs(peak) < fabs(m_realAtom[i]) ) peak = m_realAtom[i];
        }
    }

    if (fabs(peak) < 1e-10)
    {
        for(i = a; i <= b ; i++)
        {
            m_realAtom[i] = 0;
        }
    }
    else
    {
        // ======================
        // Calculate norm
        norm = 0;
        for(i=a; i <= b; i++)
        {
            norm += m_realAtom[i] * m_realAtom[i];
        }
        norm = sqrt(norm);

        m_rAtomNorm = norm;

        // ======================
        if ( norm == 0 ) norm = 1;

        // ======================

        for(i=a; i<=b; i++)
        {
            m_realAtom[i] = m_realAtom[i] / norm;
        }

    }
}

int CExpDictionary::getSignalSize()
{
    return m_signalSize;
}
CComplex* CExpDictionary::getComplexAtom()
{
    return m_complexAtom;
}
double* CExpDictionary::getRealAtom()
{
    return m_realAtom;
}
double CExpDictionary::computeUnormRAtomSample(CParameter* parm,int t)
{

    int u;
    int i;
    double sample=0.0;

    double rho = ((CExpParm*)parm)->rho;
    double xi = ((CExpParm*)parm)->xi;
    double phi = ((CExpParm*)parm)->phase;
    int a = ((CExpParm*)parm)->a;
    int b = ((CExpParm*)parm)->b;

    i=t;
    if (rho>=0)
    {
        u = a;
        sample=0;
        if(i >= u)
        {
            if(xi!=0)
                sample = exp( -rho * (double)(i-u) ) * cos ( (xi*(double)i) + phi );
            else
                sample = exp( -rho * (double)(i-u) ) * cos ( phi );
        }
    }
    else
    {
        u = m_signalSize - 1 - b;
        sample=0;
//         if( i >= u )
//         {
            if(xi!=0)
                sample = exp( rho * (double)(m_signalSize - 1 - i - u) ) * cos ( (xi*(double)(i)) + phi );
            else
                sample = exp( rho * (double)(m_signalSize - 1 - i - u) ) * cos ( phi );
//         }
    }
    return sample;
}

double CExpDictionary::computeUnormRAtomSampleSigSize(CParameter* parm,int t, int signalSize)
{

    int u;
    int i;
    double sample=0.0;

    double rho = ((CExpParm*)parm)->rho;
    double xi = ((CExpParm*)parm)->xi;
    double phi = ((CExpParm*)parm)->phase;
    int a = ((CExpParm*)parm)->a;
    int b = ((CExpParm*)parm)->b;

    i=t;
    if (rho>=0)
    {
        u = a;
        sample=0;
        if(i >= u)
        {
            if(xi!=0)
                sample = exp( -rho * (double)(i-u) ) * cos ( (xi*(double)i) + phi );
            else
                sample = exp( -rho * (double)(i-u) ) * cos ( phi );
        }
    }
    else
    {
        u = signalSize - 1 - b;
        sample=0;
//         if( i >= u )
//         {
            if(xi!=0)
                sample = exp( rho * (double)(signalSize - 1 - i - u) ) * cos ( (xi*(double)(i)) + phi );
            else
                sample = exp( rho * (double)(signalSize - 1 - i - u) ) * cos ( phi );
//         }
    }
    return sample;
}


///////////////////////////////////
// Intrinsic
///////////////////////////////////

//===============================================================
// Function: executeDecomp
// Goal:
// Return:
//===============================================================

void CExpDictionary::executeDecomp( CDataSignal* dataSignal,
                                    CStructBook** structBook,
                                    CFileGenData* genData,
                                    CFileDictionary* dicData)
{
    //CFileDictionary* dicData;
    //dicData = new CFileDictionary;
    //dicData->setFileName("dictionary.dat");
    //dicData->loadData();
    //dicData->printToScreen();

    int initBlock = genData->getInitBlock();
    int finalBlock = genData->getEndBlock();
    int nMaxStep = genData->getNumMaxStep();


    cgMatrix<double> residue(1,m_signalSize,0.0);
    cgMatrix<double> cgRealAtom(1,m_signalSize,0.0);
    cgMatrix<double> cgRealAtomAux(1,m_signalSize,0.0);
    double** pSignal = dataSignal->getSignal();

    // ====================================

    int step = 0;
    double norm =0;
    CExpParm expParm;

    int numBlock = (int)ceil((double)dataSignal->getSignalSize()/(double)dataSignal->getBlockHop());
    if (finalBlock==9999)
    {
        finalBlock = numBlock;
    }
    // Allocate memory for candidate atoms with tume continuity
    int numSignal = dataSignal->getNumSignal();
    CStructBook** sbContinuity;
    sbContinuity = new CStructBook* [numSignal];
    int k;
    for (k=0; k<numSignal; k++)
    {
        sbContinuity[k] = new CStructBookExp[numBlock];
    }

    int L = (int)ceil(((log10((double)(m_signalSize)))/(log10((double)(2)))));
    double* approxRatio;
    approxRatio = new double[L];

    for (k=0;k<L;k++)
    {
        approxRatio[k] = 0.0;
    }
    double meanApproxRatio;

    //double tolAppRatio = getApproxRatio(m_signalSize);
    double tolAppRatio =  genData->getApproxRatioTarget();
    double snrTarget = genData->getSNRTarget();

    double befSupInnerP,aftSupInnerP;

    char fileName[_MAX_PATH];
    strcpy(fileName, dataSignal->getFileName());
    char* pos;
    pos = strrchr( fileName, '.');
    char aux[_MAX_PATH];
    sprintf(aux,"_b%d-%d.sba",initBlock,finalBlock);
    strcpy( &pos[0], aux);

    // Writing the Main Header
    FILE* iosba;
    iosba = fopen(fileName,"w");
    fprintf(iosba,"Sign. Type :          %5i\n", dataSignal->getType());
    fprintf(iosba,"Dict. Type :          %5i\n", genData->getDicType());
    fprintf(iosba,"No. Signals:          %5i\n", dataSignal->getNumSignal());
    fprintf(iosba,"Signal Size:       %8i\n", dataSignal->getSignalSize());
    fprintf(iosba,"Block Hop:            %5i\n", dataSignal->getBlockHop());
    fprintf(iosba,"Block Size :          %5i\n", dataSignal->getBlockSize());
    if (genData->getSigType()==1)
    {
        fprintf(iosba,"Samp. Freq :     %10.2f\n", ((CComtradeSignal*)dataSignal)->getSamplingRate(1));
    }
    if (genData->getSigType()==2)
    {
        fprintf(iosba,"Samp. Freq :     %10.2f\n", ((CAudioSignal*)dataSignal)->getSamplingRate());
    }
    fflush(iosba);
    fclose(iosba);

    //FILE* stream2;
    //stream2 = fopen("eval_prevblocksb.out","w");
    //fclose(stream2);

    char sbbFName[_MAX_PATH];
    strcpy(sbbFName, dataSignal->getFileName());
    pos = strrchr( sbbFName, '.');
    sprintf(aux,"_b%d-%d.sbb",initBlock,finalBlock);
    strcpy( &pos[0], aux);

    FILE* iosbb;
    iosbb = fopen("header.sbb","wb");
    int dummyint;
    double dummydouble;
    dummyint = dataSignal->getType();
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = genData->getDicType();
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = dataSignal->getNumSignal();
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = dataSignal->getSignalSize();
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = dataSignal->getBlockHop();
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = dataSignal->getBlockSize();
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    if (genData->getSigType()==1)
        dummydouble = ((CComtradeSignal*)dataSignal)->getSamplingRate(1);
    if (genData->getSigType()==2)
        dummydouble = ((CAudioSignal*)dataSignal)->getSamplingRate();
    fwrite(&dummydouble, sizeof(double), 1, iosbb);
    fclose(iosbb);

    iosbb = fopen(sbbFName,"wb");
    fclose(iosbb);

    fstream file_stage;
    if (genData->getPrintDecompStage()==1)
    {
        //file pointers
        file_stage.open("decomp_stages.out",ios::out);
        // file header
        file_stage  <<  setw (10) << setfill(' ') << "Signal" << " "
                    <<  setw (10) << setfill(' ') << "Block" << " "
                    <<  setw (10) << setfill(' ') << "No." << " "
                    <<  setw (10) << setfill(' ') << "Stage" << " "
                    <<  setw (20) << setfill(' ') << "Coef." << " "
                    <<  setw (20) << setfill(' ') << "Decay" << " "
                    <<  setw (20) << setfill(' ') << "Freq" << " "
                    <<  setw (20) << setfill(' ') << "Phase"<< " "
                    <<  setw (10) << setfill(' ') << "Ti"<< " "
                    <<  setw (10) << setfill(' ') << "Tf"<< " "
                    << endl;
    }

    int iSignal;
    double sigNorm;
    int iBlock;
    int i;
    for (i=0; i < dataSignal->getNumSignal(); i++ )
    {
        // Writing the Signal Header
        iosba = fopen(fileName,"a");
        fprintf(iosba,"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
        fprintf(iosba,"Signal:               %5i\n",i+1);
        fprintf(iosba,"Norm:            %10.5f\n",dataSignal->getNorm(i));
        fflush(iosba);
        fclose(iosba);

        iosbb = fopen(sbbFName,"ab");
        iSignal = i+1;
        fwrite(&iSignal, sizeof(int), 1, iosbb);
        sigNorm = dataSignal->getNorm(i);
        fwrite(&sigNorm, sizeof(double), 1, iosbb);
        fclose(iosbb);

        if (dataSignal->getNorm(i)!=0.0)
        {
            for (int j=initBlock-1; j < finalBlock; j++ )
            {

                residue.zeros();
                // Loading signal into the vector
                if ((int)(j*(double)dataSignal->getBlockHop())+m_signalSize < dataSignal->getSignalSize())
                {
                    residue.fillVector((pSignal[i])+(int)(j*(double)dataSignal->getBlockHop()),
                                        m_signalSize);
                }
                else
                {
                    residue.fillVector((pSignal[i])+(int)(j*(double)dataSignal->getBlockHop()),
                                        dataSignal->getSignalSize() - (int)(j*(double)dataSignal->getBlockHop()));
                }
                // Normalizing initial residue (with signal norm)
                residue /= dataSignal->getNorm(i);


                norm = residue.norm();
                ((CStructBookExp*)structBook[i])[j].setNorm(norm);
                double initBlockNorm = norm;

                // Writing the Block Header
                iosba = fopen(fileName,"a");
                fprintf(iosba,"--------------------------------------------------------------\n");
                fprintf(iosba,"Block:                %5i\n",j+1);
                fprintf(iosba,"Norm:            %10.5f\n",initBlockNorm);
                fprintf(iosba,"No.    Coef.           Decaying        Freq            Phase           Ti   Tf    PrevAtom AppRatio   meanAppRat befSup     aftSup     normRatio  SNR(dB)     \n");
                fflush(iosba);
                fclose(iosba);

                iosbb = fopen(sbbFName,"ab");
                iBlock = j+1;
                fwrite(&iBlock, sizeof(int), 1, iosbb);
                fwrite(&initBlockNorm, sizeof(double), 1, iosbb);
                fclose(iosbb);

                //stream2 = fopen("eval_prevblocksb.out","a");
                //fprintf(stream2,"Signal %i; Block %i\n",i+1,j+1);
                //fprintf(stream2,"current     prev\n");
                //fclose(stream2);
                int decomp_stage;


                // Beginning of the decomposition
                if (norm!=0)
                {
                    if ((j!=0) &&
                        (genData->getFlagEvalAtomCont()==1))
                    {
                        evalAtomContinuity( residue,
                                            sbContinuity[i],
                                            structBook[i],
                                            i,
                                            j,
                                            j-1,
                                            step,
                                            0,
                                            fileName,
                                            sbbFName,
                                            initBlockNorm,
                                            file_stage,
                                            genData->getPrintDecompStage());
                    }
                    int nAtomCont = step;
                    cout << "- nAtomCont: " << nAtomCont << endl;
                    do
                    {
                        if (step>=nMaxStep) break;
#ifdef DBG_WATCH_DCMP_STEP
                        residue.PrintToFile("residue.dat");
#endif
                        cout << "##########################################################################" << endl;
                        cout << "->Decomposing Signal: " << i+1 << "; Block: "<< j+1 << "; Atom: "<< step+1 << endl;
                        cout << "##########################################################################" << endl;
                        // Signal projection over a discrete dictionary
                        cout << "Signal projection over a discrete dictionary" << endl;
                        //expParm = fastMPKolasa(residue);
                        fastMPKolasaModified(residue,dataSignal,dicData,&expParm);
                        expParm.printParm2Screen();
                        if (genData->getPrintDecompStage()==1)
                        {
                            decomp_stage=1;
                            file_stage  <<  setw (10) << setfill(' ') << i+1 << " "
                                        <<  setw (10) << setfill(' ') << j+1 << " "
                                        <<  setw (10) << setfill(' ') << step+1 << " "
                                        <<  setw (10) << setfill(' ') << decomp_stage << " "
                                        <<  setw (20) << setfill(' ') << expParm.innerProd << " "
                                        <<  setw (20) << setfill(' ') << expParm.rho << " "
                                        <<  setw (20) << setfill(' ') << expParm.xi << " "
                                        <<  setw (20) << setfill(' ') << expParm.phase << " "
                                        <<  setw (10) << setfill(' ') << expParm.a << " "
                                        <<  setw (10) << setfill(' ') << expParm.b << " "
                                        << endl;
                        }

                        // Calculate approximation ratio referred to this step
                        //cout<< "->Calculating approximation ratio..." << endl;
                        approxRatio[step%L] = fabs(expParm.innerProd)/norm;
                        meanApproxRatio = 0.0;
                        for (k=0;k<L;k++)
                        {
                            meanApproxRatio += approxRatio[k]/(double)L;
                        }
                        cout << "Mean Approx. Ratio: " << meanApproxRatio << endl;

                        // Maximize approximation by finding optimum parameters
                        //optimizeContinuousParms(residue, &expParm);

                        // Heuristics
                        befSupInnerP = 0.0;
                        aftSupInnerP = 0.0;

                        if (expParm.a != expParm.b)  // damp and pure cases
                        {
                            if (genData->getFlagFindSupport()==1)
                            {
                                befSupInnerP = expParm.innerProd;
                                cout<< "->Finding the best time support ..." << endl;
                                if (expParm.rho==0.0)
                                {
                                    // Two-way search
                                    findFastBestTimeSupport(residue,&expParm,1,
                                                            genData->getCoefTempSup());
                                }
                                else
                                {
                                    // One-way search
                                    findFastBestTimeSupport(residue,&expParm,0,
                                                            genData->getCoefTempSup());
                                }
                                expParm.printParm2Screen();
                                if (genData->getPrintDecompStage()==1)
                                {
                                    decomp_stage=2;
                                    file_stage  <<  setw (10) << setfill(' ') << i+1 << " "
                                                <<  setw (10) << setfill(' ') << j+1 << " "
                                                <<  setw (10) << setfill(' ') << step+1 << " "
                                                <<  setw (10) << setfill(' ') << decomp_stage << " "
                                                <<  setw (20) << setfill(' ') << expParm.innerProd << " "
                                                <<  setw (20) << setfill(' ') << expParm.rho << " "
                                                <<  setw (20) << setfill(' ') << expParm.xi << " "
                                                <<  setw (20) << setfill(' ') << expParm.phase << " "
                                                <<  setw (10) << setfill(' ') << expParm.a << " "
                                                <<  setw (10) << setfill(' ') << expParm.b << " "
                                                << endl;
                                }
                                aftSupInnerP = expParm.innerProd;
                            }
                            if (genData->getFlagOptDecay()==1)
                            {
                                optimizeDecaying(residue,&expParm);
                                expParm.printParm2Screen();
                                if (genData->getPrintDecompStage()==1)
                                {
                                    decomp_stage=3;
                                    file_stage  <<  setw (10) << setfill(' ') << i+1 << " "
                                                <<  setw (10) << setfill(' ') << j+1 << " "
                                                <<  setw (10) << setfill(' ') << step+1 << " "
                                                <<  setw (10) << setfill(' ') << decomp_stage << " "
                                                <<  setw (20) << setfill(' ') << expParm.innerProd << " "
                                                <<  setw (20) << setfill(' ') << expParm.rho << " "
                                                <<  setw (20) << setfill(' ') << expParm.xi << " "
                                                <<  setw (20) << setfill(' ') << expParm.phase << " "
                                                <<  setw (10) << setfill(' ') << expParm.a << " "
                                                <<  setw (10) << setfill(' ') << expParm.b << " "
                                                << endl;
                                }
                            }
                        }
                        //  proceedHeuristics(residue,&expParm,dataSignal);
                        if ((j!=0) && (genData->getFlagSeqBlock()==1))
                        {
                            cout<< "->Search in the Previous Block Structure Book ..." << endl;
                            searchSBPreviousBlock(  residue,
                                                    &expParm,
                                                    sbContinuity[i],//structBook[i],
                                                    j-1,
                                                    0,
                                                    genData->getCoefSeqBlock());
                            expParm.printParm2Screen();
                        }



                        // Adjusting parameters
                        adjustParameters(residue,&expParm);

                        // Print to screen the chosen atom parameter
                        cout << "Parameters adjusted!!" << endl;
                        cout << "Chosen atom parameters: " << endl;
                        expParm.printParm2Screen();

                        // Updating residue
                        cout<< "->Updating residue..." << endl << endl;
                        setRealAtom((CParameter*)&expParm);
                        cgRealAtom.fillVector(m_realAtom);
                        cgRealAtomAux = cgRealAtom*expParm.innerProd;
                        //cgRealAtomAux.PrintToFile("scaled_atom.dat");
                        residue = residue - cgRealAtom*expParm.innerProd;
                        norm = residue.norm();

                        // Adding the hop to time support
                        //expParm.a = expParm.a + (int)(j*(double)dataSignal->getBlockHop());
                        //expParm.b = expParm.b + (int)(j*(double)dataSignal->getBlockHop());

                        // Add element to structure book
                        ((CStructBookExp*)structBook[i])[j].addElement(&expParm);
                        int indorig = ((CStructBookExp*)structBook[i])[j].getNumElement();
                        ((CStructBookExp*)structBook[i])[j].setNextAtomIndex(indorig-1,-1);
                        ((CStructBookExp*)structBook[i])[j].setPrevAtomIndex(indorig-1,-1);
                        ((CStructBookExp*)structBook[i])[j].setOrigAtomIndex(indorig-1,indorig-1);

                        // Add element to structure book with candidate atom with continuity
                        // for the next block
                        if (expParm.b==m_signalSize-1)
                        {
                            ((CStructBookExp*)sbContinuity[i])[j].addElement(&expParm);
                            int ind = ((CStructBookExp*)sbContinuity[i])[j].getNumElement();
                            ((CStructBookExp*)sbContinuity[i])[j].setNextAtomIndex(ind-1,-1);
                            ((CStructBookExp*)sbContinuity[i])[j].setPrevAtomIndex(ind-1,-1);
                            ((CStructBookExp*)sbContinuity[i])[j].setOrigAtomIndex(ind-1,indorig-1);
                        }


                        iosba = fopen(fileName,"a");
                        ((CStructBookExp*)structBook[i])[j].saveElementASCII(   iosba,
                                                                                meanApproxRatio,
                                                                                approxRatio[step%L],
                                                                                befSupInnerP,
                                                                                aftSupInnerP,
                                                                                (norm/initBlockNorm));
                        fflush(iosba);
                        fclose(iosba);

                        cout << "SNR: " << 20*log10(initBlockNorm/norm) << " (dB)"<< endl;

                        iosbb = fopen(sbbFName,"ab");
                        ((CStructBookExp*)structBook[i])[j].saveElementBin(iosbb);
                        fclose(iosbb);
                        ((CStructBookExp*)structBook[i])[j].printElementToScreen(indorig-1);

                        step++;
                    }
                    //while( ( fabs(expParm.innerProd) > pow(2.0,-15) ) && (step<NUM_MAX_STEP)  );
                    while(      (   (meanApproxRatio > tolAppRatio) ||
                                    (step<(L+nAtomCont)) ||
                                    (20*log10(initBlockNorm/norm)<snrTarget)     )
                                //((norm/initBlockNorm)>1e-8)
                                && (step<nMaxStep)
                                //&& ( fabs(expParm.innerProd) > 1e-8 )
                         ); //( fabs(expParm.innerProd) > pow(2.0,-15) )   );

                }
                else
                {
                    cout << "  ### Block "<< j+1 <<" with null samples ### " << endl;
                    iosba = fopen(fileName,"a");
                    fprintf(iosba,"###### Block with null samples ######\n");
                    fflush(iosba);
                    fclose(iosba);
                }

                iosba = fopen(fileName,"a");
                fprintf(iosba,"99999\n");
                fflush(iosba);
                fclose(iosba);

                iosbb = fopen(sbbFName,"ab");
                dummyint = 99999;
                fwrite( &dummyint, sizeof( int ), 1, iosbb );
                fclose(iosbb);

                step = 0;
            }
            iosba = fopen(fileName,"a");
            fprintf(iosba,"88888\n");
            fflush(iosba);
            fclose(iosba);

            iosbb = fopen(sbbFName,"ab");
            dummyint = 88888;
            fwrite( &dummyint, sizeof( int ), 1, iosbb );
            fclose(iosbb);
        }
        else
        {
                cout << "  ### Signal "<< i+1 <<" with null samples ### " << endl;
                iosba = fopen(fileName,"a");
                fprintf(iosba,"###### Signal with null samples ######\n");
                fflush(iosba);
                fclose(iosba);
        }
        iosba = fopen(fileName,"a");
        fprintf(iosba,"77777\n");
        fflush(iosba);
        fclose(iosba);

        iosbb = fopen(sbbFName,"ab");
        dummyint = 77777;
        fwrite( &dummyint, sizeof( int ), 1, iosbb );
        fclose(iosbb);

        //delete [] sbPrevProj;
    }

    delete [] approxRatio;
    for (i=0; i<dataSignal->getNumSignal(); i++)
    {
        delete [] ((CStructBookExp*)sbContinuity[i]);
    }
    delete [] sbContinuity;
    //delete dicData;

    if (genData->getPrintDecompStage()==1)
    {
        file_stage.close();
    }

}

//===============================================================
// Function: executeDecompElectric
// Goal:
// Return:
//===============================================================

void CExpDictionary::executeDecompElectric( CDataSignal* dataSignal,
                                            CStructBook** structBook,
                                            CFileDecomp* genData,
                                            CFileDictionary* dicData,
                                            CFileDecompBlockRange* blockRange)
{
    //CFileDictionary* dicData;
    //dicData = new CFileDictionary;
    //dicData->setFileName("dictionary.dat");
    //dicData->loadData();
    //dicData->printToScreen();

    int initBlock = blockRange->getInitBlock();
    int finalBlock = blockRange->getEndBlock();
    int nMaxStep = genData->getNumMaxStep();


    cgMatrix<double> residue(1,m_signalSize,0.0);
    cgMatrix<double> cgRealAtom(1,m_signalSize,0.0);
    cgMatrix<double> cgRealAtomAux(1,m_signalSize,0.0);
    double** pSignal = dataSignal->getSignal();

    // ====================================

    int step = 0;
    double norm =0;
    CExpParm expParm;

    int numBlock = (int)ceil((double)dataSignal->getSignalSize()/(double)dataSignal->getBlockHop());
    if (finalBlock==9999)
    {
        finalBlock = numBlock;
    }
    // Allocate memory for candidate atoms with tume continuity
    int numSignal = dataSignal->getNumSignal();
    CStructBook** sbContinuity;
    sbContinuity = new CStructBook* [numSignal];
    int k;
    for (k=0; k<numSignal; k++)
    {
        sbContinuity[k] = new CStructBookExp[numBlock];
    }

    int L = (int)ceil(((log10((double)(m_signalSize)))/(log10((double)(2)))));
    double* approxRatio;
    approxRatio = new double[L];

    for (k=0;k<L;k++)
    {
        approxRatio[k] = 0.0;
    }
    double meanApproxRatio;

    //double tolAppRatio = getApproxRatio(m_signalSize);
    double tolAppRatio =  genData->getApproxRatioTarget();
    double snrTarget = genData->getSNRTarget();

    double befSupInnerP,aftSupInnerP;

    char fileName[_MAX_PATH];
    strcpy(fileName, dataSignal->getFileName());
    char* pos;
    pos = strrchr( fileName, '.');
    char aux[_MAX_PATH];
    sprintf(aux,"_b%d-%d.sba",initBlock,finalBlock);
    strcpy( &pos[0], aux);

    // Writing the Main Header
    FILE* iosba;
    iosba = fopen(fileName,"w");
    fprintf(iosba,"Sign. Type :          %5i\n", dataSignal->getType());
    fprintf(iosba,"Dict. Type :          %5i\n", genData->getDicType());
    fprintf(iosba,"No. Signals:          %5i\n", dataSignal->getNumSignal());
    fprintf(iosba,"Signal Size:       %8i\n", dataSignal->getSignalSize());
    fprintf(iosba,"Block Hop:            %5i\n", dataSignal->getBlockHop());
    fprintf(iosba,"Block Size :          %5i\n", dataSignal->getBlockSize());
    if (genData->getSigType()==1)
    {
        fprintf(iosba,"Samp. Freq :     %10.2f\n", ((CComtradeSignal*)dataSignal)->getSamplingRate(1));
    }
    if (genData->getSigType()==2)
    {
        fprintf(iosba,"Samp. Freq :     %10.2f\n", ((CAudioSignal*)dataSignal)->getSamplingRate());
    }
    fflush(iosba);
    fclose(iosba);

    //FILE* stream2;
    //stream2 = fopen("eval_prevblocksb.out","w");
    //fclose(stream2);

    char sbbFName[_MAX_PATH];
    strcpy(sbbFName, dataSignal->getFileName());
    pos = strrchr( sbbFName, '.');
    sprintf(aux,"_b%d-%d.sbb",initBlock,finalBlock);
    strcpy( &pos[0], aux);

    FILE* iosbb;
    iosbb = fopen("header.sbb","wb");
    int dummyint;
    double dummydouble;
    dummyint = dataSignal->getType();
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = genData->getDicType();
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = dataSignal->getNumSignal();
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = dataSignal->getSignalSize();
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = dataSignal->getBlockHop();
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = dataSignal->getBlockSize();
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    if (genData->getSigType()==1)
        dummydouble = ((CComtradeSignal*)dataSignal)->getSamplingRate(1);
    if (genData->getSigType()==2)
        dummydouble = ((CAudioSignal*)dataSignal)->getSamplingRate();
    fwrite(&dummydouble, sizeof(double), 1, iosbb);
    fclose(iosbb);

    iosbb = fopen(sbbFName,"wb");
    fclose(iosbb);

    fstream file_stage;
    if (genData->getPrintDecompStage()==1)
    {
        //file pointers
        file_stage.open("decomp_stages.out",ios::out);
        // file header
        file_stage  <<  setw (10) << setfill(' ') << "Signal" << " "
                    <<  setw (10) << setfill(' ') << "Block" << " "
                    <<  setw (10) << setfill(' ') << "No." << " "
                    <<  setw (10) << setfill(' ') << "Stage" << " "
                    <<  setw (20) << setfill(' ') << "Coef." << " "
                    <<  setw (20) << setfill(' ') << "Decay" << " "
                    <<  setw (20) << setfill(' ') << "Freq" << " "
                    <<  setw (20) << setfill(' ') << "Phase"<< " "
                    <<  setw (10) << setfill(' ') << "Ti"<< " "
                    <<  setw (10) << setfill(' ') << "Tf"<< " "
                    << endl;
    }

    int iSignal;
    double sigNorm;
    int iBlock;
    int i;
    for (i=0; i < dataSignal->getNumSignal(); i++ )
    {
        // Writing the Signal Header
        iosba = fopen(fileName,"a");
        fprintf(iosba,"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
        fprintf(iosba,"Signal:               %5i\n",i+1);
        fprintf(iosba,"Norm:            %10.5f\n",dataSignal->getNorm(i));
        fflush(iosba);
        fclose(iosba);

        iosbb = fopen(sbbFName,"ab");
        iSignal = i+1;
        fwrite(&iSignal, sizeof(int), 1, iosbb);
        sigNorm = dataSignal->getNorm(i);
        fwrite(&sigNorm, sizeof(double), 1, iosbb);
        fclose(iosbb);

        if (dataSignal->getNorm(i)!=0.0)
        {
            for (int j=initBlock-1; j < finalBlock; j++ )
            {

                residue.zeros();
                // Loading signal into the vector
                if ((int)(j*(double)dataSignal->getBlockHop())+m_signalSize < dataSignal->getSignalSize())
                {
                    residue.fillVector((pSignal[i])+(int)(j*(double)dataSignal->getBlockHop()),
                                        m_signalSize);
                }
                else
                {
                    residue.fillVector((pSignal[i])+(int)(j*(double)dataSignal->getBlockHop()),
                                        dataSignal->getSignalSize() - (int)(j*(double)dataSignal->getBlockHop()));
                }
                // Normalizing initial residue (with signal norm)
                residue /= dataSignal->getNorm(i);


                norm = residue.norm();
                ((CStructBookExp*)structBook[i])[j].setNorm(norm);
                double initBlockNorm = norm;

                // Writing the Block Header
                iosba = fopen(fileName,"a");
                fprintf(iosba,"--------------------------------------------------------------\n");
                fprintf(iosba,"Block:                %5i\n",j+1);
                fprintf(iosba,"Norm:            %10.5f\n",initBlockNorm);
                fprintf(iosba,"No.    Coef.           Decaying        Freq            Phase           Ti   Tf    PrevAtom AppRatio   meanAppRat befSup     aftSup     normRatio  SNR(dB)     \n");
                fflush(iosba);
                fclose(iosba);

                iosbb = fopen(sbbFName,"ab");
                iBlock = j+1;
                fwrite(&iBlock, sizeof(int), 1, iosbb);
                fwrite(&initBlockNorm, sizeof(double), 1, iosbb);
                fclose(iosbb);

                //stream2 = fopen("eval_prevblocksb.out","a");
                //fprintf(stream2,"Signal %i; Block %i\n",i+1,j+1);
                //fprintf(stream2,"current     prev\n");
                //fclose(stream2);
                int decomp_stage;


                // Beginning of the decomposition
                if (norm!=0)
                {
                    if ((j!=0) &&
                        (genData->getFlagEvalAtomCont()==1))
                    {
                        evalAtomContinuity( residue,
                                            sbContinuity[i],
                                            structBook[i],
                                            i,
                                            j,
                                            j-1,
                                            step,
                                            0,
                                            fileName,
                                            sbbFName,
                                            initBlockNorm,
                                            file_stage,
                                            genData->getPrintDecompStage());
                    }
                    int nAtomCont = step;
                    cout << "- nAtomCont: " << nAtomCont << endl;
                    do
                    {
                        if (step>=nMaxStep) break;
#ifdef DBG_WATCH_DCMP_STEP
                        residue.PrintToFile("residue.dat");
#endif
                        cout << "##########################################################################" << endl;
                        cout << "->Decomposing Signal: " << i+1 << "; Block: "<< j+1 << "; Atom: "<< step+1 << endl;
                        cout << "##########################################################################" << endl;
                        // Signal projection over a discrete dictionary
                        cout << "Signal projection over a discrete dictionary" << endl;
                        //expParm = fastMPKolasa(residue);
                        fastMPKolasaModified(residue,dataSignal,dicData,&expParm);
                        expParm.printParm2Screen();
                        if (genData->getPrintDecompStage()==1)
                        {
                            decomp_stage=1;
                            file_stage  <<  setw (10) << setfill(' ') << i+1 << " "
                                        <<  setw (10) << setfill(' ') << j+1 << " "
                                        <<  setw (10) << setfill(' ') << step+1 << " "
                                        <<  setw (10) << setfill(' ') << decomp_stage << " "
                                        <<  setw (20) << setfill(' ') << expParm.innerProd << " "
                                        <<  setw (20) << setfill(' ') << expParm.rho << " "
                                        <<  setw (20) << setfill(' ') << expParm.xi << " "
                                        <<  setw (20) << setfill(' ') << expParm.phase << " "
                                        <<  setw (10) << setfill(' ') << expParm.a << " "
                                        <<  setw (10) << setfill(' ') << expParm.b << " "
                                        << endl;
                        }

                        // Calculate approximation ratio referred to this step
                        //cout<< "->Calculating approximation ratio..." << endl;
                        approxRatio[step%L] = fabs(expParm.innerProd)/norm;
                        meanApproxRatio = 0.0;
                        for (k=0;k<L;k++)
                        {
                            meanApproxRatio += approxRatio[k]/(double)L;
                        }
                        cout << "Mean Approx. Ratio: " << meanApproxRatio << endl;

                        // Maximize approximation by finding optimum parameters
                        //optimizeContinuousParms(residue, &expParm);

                        // Heuristics
                        befSupInnerP = 0.0;
                        aftSupInnerP = 0.0;

                        if (expParm.a != expParm.b)  // damp and pure cases
                        {
                            if (genData->getFlagFindSupport()==1)
                            {
                                befSupInnerP = expParm.innerProd;
                                cout<< "->Finding the best time support ..." << endl;
                                if (expParm.rho==0.0)
                                {
                                    // Two-way search
                                    findFastBestTimeSupport(residue,&expParm,1,
                                                            genData->getCoefTempSup());
                                }
                                else
                                {
                                    // One-way search
                                    findFastBestTimeSupport(residue,&expParm,0,
                                                            genData->getCoefTempSup());
                                }
                                expParm.printParm2Screen();
                                if (genData->getPrintDecompStage()==1)
                                {
                                    decomp_stage=2;
                                    file_stage  <<  setw (10) << setfill(' ') << i+1 << " "
                                                <<  setw (10) << setfill(' ') << j+1 << " "
                                                <<  setw (10) << setfill(' ') << step+1 << " "
                                                <<  setw (10) << setfill(' ') << decomp_stage << " "
                                                <<  setw (20) << setfill(' ') << expParm.innerProd << " "
                                                <<  setw (20) << setfill(' ') << expParm.rho << " "
                                                <<  setw (20) << setfill(' ') << expParm.xi << " "
                                                <<  setw (20) << setfill(' ') << expParm.phase << " "
                                                <<  setw (10) << setfill(' ') << expParm.a << " "
                                                <<  setw (10) << setfill(' ') << expParm.b << " "
                                                << endl;
                                }
                                aftSupInnerP = expParm.innerProd;
                            }
                            if (genData->getFlagOptDecay()==1)
                            {
                                optimizeDecaying(residue,&expParm);
                                expParm.printParm2Screen();
                                if (genData->getPrintDecompStage()==1)
                                {
                                    decomp_stage=3;
                                    file_stage  <<  setw (10) << setfill(' ') << i+1 << " "
                                                <<  setw (10) << setfill(' ') << j+1 << " "
                                                <<  setw (10) << setfill(' ') << step+1 << " "
                                                <<  setw (10) << setfill(' ') << decomp_stage << " "
                                                <<  setw (20) << setfill(' ') << expParm.innerProd << " "
                                                <<  setw (20) << setfill(' ') << expParm.rho << " "
                                                <<  setw (20) << setfill(' ') << expParm.xi << " "
                                                <<  setw (20) << setfill(' ') << expParm.phase << " "
                                                <<  setw (10) << setfill(' ') << expParm.a << " "
                                                <<  setw (10) << setfill(' ') << expParm.b << " "
                                                << endl;
                                }
                            }

                        }
                        if (expParm.rho != 0.0) // damp case
                        {
                            // Discrimine Sine
                            cout<< "->Discrimining sine..." << endl;
                            discrimineSine(residue,&expParm,dataSignal);
                            expParm.printParm2Screen();
                            if (genData->getPrintDecompStage()==1)
                            {
                                decomp_stage=4;
                                file_stage  <<  setw (10) << setfill(' ') << i+1 << " "
                                            <<  setw (10) << setfill(' ') << j+1 << " "
                                            <<  setw (10) << setfill(' ') << step+1 << " "
                                            <<  setw (10) << setfill(' ') << decomp_stage << " "
                                            <<  setw (20) << setfill(' ') << expParm.innerProd << " "
                                            <<  setw (20) << setfill(' ') << expParm.rho << " "
                                            <<  setw (20) << setfill(' ') << expParm.xi << " "
                                            <<  setw (20) << setfill(' ') << expParm.phase << " "
                                            <<  setw (10) << setfill(' ') << expParm.a << " "
                                            <<  setw (10) << setfill(' ') << expParm.b << " "
                                            << endl;
                            }
                        }


                        if ((j!=0) && (genData->getFlagSeqBlock()==1))
                        {
                            cout<< "->Search in the Previous Block Structure Book ..." << endl;
                            searchSBPreviousBlock(  residue,
                                                    &expParm,
                                                    sbContinuity[i],
                                                    j-1,
                                                    0,
                                                    genData->getCoefSeqBlock());
                            expParm.printParm2Screen();
                        }



                        // Adjusting parameters
                        adjustParameters(residue,&expParm);

                        // Print to screen the chosen atom parameter
                        cout << "Parameters adjusted!!" << endl;
                        cout << "Chosen atom parameters: " << endl;
                        expParm.printParm2Screen();

                        // Updating residue
                        cout<< "->Updating residue..." << endl << endl;
                        setRealAtom((CParameter*)&expParm);
                        cgRealAtom.fillVector(m_realAtom);
                        cgRealAtomAux = cgRealAtom*expParm.innerProd;
                        //cgRealAtomAux.PrintToFile("scaled_atom.dat");
                        residue = residue - cgRealAtom*expParm.innerProd;
                        norm = residue.norm();

                        // Adding the hop to time support
                        //expParm.a = expParm.a + (int)(j*(double)dataSignal->getBlockHop());
                        //expParm.b = expParm.b + (int)(j*(double)dataSignal->getBlockHop());

                        // Add element to structure book
                        ((CStructBookExp*)structBook[i])[j].addElement(&expParm);
                        int indorig = ((CStructBookExp*)structBook[i])[j].getNumElement();
                        ((CStructBookExp*)structBook[i])[j].setNextAtomIndex(indorig-1,-1);
                        ((CStructBookExp*)structBook[i])[j].setPrevAtomIndex(indorig-1,-1);
                        ((CStructBookExp*)structBook[i])[j].setOrigAtomIndex(indorig-1,indorig-1);

                        // Add element to structure book with candidate atom with continuity
                        // for the next block
                        if (expParm.b==m_signalSize-1)
                        {
                            ((CStructBookExp*)sbContinuity[i])[j].addElement(&expParm);
                            int ind = ((CStructBookExp*)sbContinuity[i])[j].getNumElement();
                            ((CStructBookExp*)sbContinuity[i])[j].setNextAtomIndex(ind-1,-1);
                            ((CStructBookExp*)sbContinuity[i])[j].setPrevAtomIndex(ind-1,-1);
                            ((CStructBookExp*)sbContinuity[i])[j].setOrigAtomIndex(ind-1,indorig-1);
                        }


                        iosba = fopen(fileName,"a");
                        ((CStructBookExp*)structBook[i])[j].saveElementASCII(   iosba,
                                                                                meanApproxRatio,
                                                                                approxRatio[step%L],
                                                                                befSupInnerP,
                                                                                aftSupInnerP,
                                                                                (norm/initBlockNorm));
                        fflush(iosba);
                        fclose(iosba);

                        cout << "SNR: " << 20*log10(initBlockNorm/norm) << " (dB)"<< endl;

                        iosbb = fopen(sbbFName,"ab");
                        ((CStructBookExp*)structBook[i])[j].saveElementBin(iosbb);
                        fclose(iosbb);
                        ((CStructBookExp*)structBook[i])[j].printElementToScreen(indorig-1);

                        step++;
                    }
                    //while( ( fabs(expParm.innerProd) > pow(2.0,-15) ) && (step<NUM_MAX_STEP)  );
                    while(      (   (meanApproxRatio > tolAppRatio) ||
                                    (step<(L+nAtomCont)) ||
                                    (20*log10(initBlockNorm/norm)<snrTarget)     )
                                //((norm/initBlockNorm)>1e-8)
                                && (step<nMaxStep)
                                //&& ( fabs(expParm.innerProd) > 1e-8 )
                         ); //( fabs(expParm.innerProd) > pow(2.0,-15) )   );

                }
                else
                {
                    cout << "  ### Block "<< j+1 <<" with null samples ### " << endl;
                    iosba = fopen(fileName,"a");
                    fprintf(iosba,"###### Block with null samples ######\n");
                    fflush(iosba);
                    fclose(iosba);
                }

                iosba = fopen(fileName,"a");
                fprintf(iosba,"99999\n");
                fflush(iosba);
                fclose(iosba);

                iosbb = fopen(sbbFName,"ab");
                dummyint = 99999;
                fwrite( &dummyint, sizeof( int ), 1, iosbb );
                fclose(iosbb);

                step = 0;
            }
            iosba = fopen(fileName,"a");
            fprintf(iosba,"88888\n");
            fflush(iosba);
            fclose(iosba);

            iosbb = fopen(sbbFName,"ab");
            dummyint = 88888;
            fwrite( &dummyint, sizeof( int ), 1, iosbb );
            fclose(iosbb);
        }
        else
        {
                cout << "  ### Signal "<< i+1 <<" with null samples ### " << endl;
                iosba = fopen(fileName,"a");
                fprintf(iosba,"###### Signal with null samples ######\n");
                fflush(iosba);
                fclose(iosba);
        }
        iosba = fopen(fileName,"a");
        fprintf(iosba,"77777\n");
        fflush(iosba);
        fclose(iosba);

        iosbb = fopen(sbbFName,"ab");
        dummyint = 77777;
        fwrite( &dummyint, sizeof( int ), 1, iosbb );
        fclose(iosbb);

        //delete [] sbPrevProj;
    }

    delete [] approxRatio;
    for (i=0; i<dataSignal->getNumSignal(); i++)
    {
        delete [] ((CStructBookExp*)sbContinuity[i]);
    }
    delete [] sbContinuity;
    //delete dicData;

    if (genData->getPrintDecompStage()==1)
    {
        file_stage.close();
    }

}


//===============================================================
// Function: findMPKolasa
// Goal:
// Return:
//===============================================================

CExpParm CExpDictionary::fastMPKolasa(cgMatrix<double>& residue)
{
    CParameter* parm = new CExpParm;
    CExpParm expParm;
    int i,k;
    double opt_phase=0;
    int NC;
    int N;

    double *in1,*in2;
    fftw_complex *out1, *out2;
    fftw_plan plan_forward1, plan_forward2;

    FILE* stream;
//  ---------------------------------------------------
    stream = fopen( "results_kolasa.out", "w" );


    N = m_signalSize;
    NC = ( N/2 ) +1;

    in1 = (double*) fftw_malloc ( sizeof ( double ) * N);
    out1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * NC );
    in2 = (double*) fftw_malloc ( sizeof ( double ) * N );
    out2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * NC );

    plan_forward1 = fftw_plan_dft_r2c_1d ( N, in1, out1, FFTW_ESTIMATE );
    plan_forward2 = fftw_plan_dft_r2c_1d ( N, in2, out2, FFTW_ESTIMATE );

    double s, tau, xi;
    int k_pi;
    double C;
    double innerProd_xp;
    double innerProd_xq;
    double innerProd_pp;
    double innerProd_qq;
    double innerProd_pq;
    double a1,b1;
    double innerProd;
    double maxInnerProd = 0;


    // s=1 Impulse
    for(i=0;i<N;i++)
    {
        innerProd = residue[0][i];
        if (fabs(innerProd)>fabs(maxInnerProd))
        {
            maxInnerProd = innerProd;
            expParm.innerProd = maxInnerProd;
            expParm.rho = 1.0;
            expParm.xi = 0.0;
            expParm.phase = 0.0 ;
            expParm.a = i;
            expParm.b = i;

        }
        fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
                                  innerProd, innerProd, 0.0, 1.0, 0.0,//innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
                                  0.0, 1.0, (double)i, 0.0,0.0);
    }

    // Damped sinusoids
    s=2;
    while (s<N)
    {
        int delta_tau = s;
        // Decreasing damped sinusoids
        for (tau=0; tau<(double)N; tau=tau+delta_tau)
        {
            ((CExpParm*)parm)->rho = (double)(1/s);
            ((CExpParm*)parm)->xi = 0.0;
            ((CExpParm*)parm)->phase = 0.0;
            ((CExpParm*)parm)->a = tau;
            ((CExpParm*)parm)->b = N-1;
            setRealAtom(parm);

            for (i=0;i<N;i++)
            {
                in1[i] = residue[0][i] * m_realAtom[i];
                in2[i] = m_realAtom[i] * m_realAtom[i];
            }

            fftw_execute ( plan_forward1 );

            fftw_execute ( plan_forward2 );

            C = out2[0][0];
            k_pi=NC-1;
            int delta_k = N/s;
            for (k=0;k<NC;k=k+delta_k)
            {
                innerProd_xp = out1[k][0];  // /sqrt((double)N);
                innerProd_xq = -out1[k][1]; // /sqrt((double)N);

                if ( (2*k)<NC )
                {
                    innerProd_pp = 0.5 * (C + out2[2*k][0]);// /sqrt((double)N) );
                    innerProd_qq = 0.5 * (C - out2[2*k][0]);// /sqrt((double)N) );
                    innerProd_pq = -0.5 * (out2[2*k][1]);// /sqrt((double)N) );
                }
                else
                {
                    innerProd_pp = 0.5 * (C + out2[N-(2*k)][0]);// /sqrt((double)N) );
                    innerProd_qq = 0.5 * (C - out2[N-(2*k)][0]);// /sqrt((double)N) );
                    innerProd_pq = -0.5 * ( - out2[N-(2*k)][1]);// /sqrt((double)N) );
                }

                a1 = innerProd_xp * innerProd_qq - innerProd_xq * innerProd_pq;
                b1 = innerProd_xq * innerProd_pp - innerProd_xp * innerProd_pq;

                if ( (k == 0) || (k == k_pi) )
                {
                    opt_phase = 0;
                    innerProd = innerProd_xp  / sqrt(innerProd_pp);
                }
                else if (a1 == 0)
                {
                    opt_phase = (double)(pi/2);
                    innerProd = -innerProd_xq / sqrt(innerProd_qq);
                }
                else if ( (a1!=0) && (k != 0) )
                {
                    opt_phase = atan( -(b1/a1) );
                    innerProd = (a1/fabs(a1))*(innerProd_xp*a1 + innerProd_xq*b1) /
                                sqrt(a1*a1*innerProd_pp +
                                     b1*b1*innerProd_qq +
                                     2*a1*b1*innerProd_pq);
                }

                xi = (k*2*pi)/N;

                if (fabs(innerProd)>fabs(maxInnerProd))
                {
                    maxInnerProd = innerProd;
                    expParm.innerProd = maxInnerProd;
                    expParm.rho = (double)1/s;
                    expParm.xi = xi;
                    expParm.phase = opt_phase ;
                    expParm.a = tau;
                    expParm.b = N-1;
                }

                // For debug
                //fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
                //                innerProd, innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
                //                innerProd_pq, s, (double)tau, xi,opt_phase);
            }
        }

        // Increasing
        for (tau=1; tau<(double)N; tau=tau+delta_tau)
        {
            ((CExpParm*)parm)->rho = -(double)(1/s);
            ((CExpParm*)parm)->xi = 0.0;
            ((CExpParm*)parm)->phase = 0.0;
            ((CExpParm*)parm)->a = 0.0;
            ((CExpParm*)parm)->b = tau;
            setRealAtom(parm);

            for (i=0;i<N;i++)
            {
                in1[i] = residue[0][i] * m_realAtom[i];
                in2[i] = m_realAtom[i] * m_realAtom[i];
            }

            fftw_execute ( plan_forward1 );

            fftw_execute ( plan_forward2 );

            C = out2[0][0];
            k_pi=NC-1;
            int delta_k = N/s;
            for (k=0;k<NC;k=k+delta_k)
            {
                innerProd_xp = out1[k][0];  // /sqrt((double)N);
                innerProd_xq = -out1[k][1]; // /sqrt((double)N);

                if ( (2*k)<NC )
                {
                    innerProd_pp = 0.5 * (C + out2[2*k][0]);// /sqrt((double)N) );
                    innerProd_qq = 0.5 * (C - out2[2*k][0]);// /sqrt((double)N) );
                    innerProd_pq = -0.5 * (out2[2*k][1]);// /sqrt((double)N) );
                }
                else
                {
                    innerProd_pp = 0.5 * (C + out2[N-(2*k)][0]);// /sqrt((double)N) );
                    innerProd_qq = 0.5 * (C - out2[N-(2*k)][0]);// /sqrt((double)N) );
                    innerProd_pq = -0.5 * ( - out2[N-(2*k)][1]);// /sqrt((double)N) );
                }

                a1 = innerProd_xp * innerProd_qq - innerProd_xq * innerProd_pq;
                b1 = innerProd_xq * innerProd_pp - innerProd_xp * innerProd_pq;

                if ( (k == 0) || (k == k_pi) )
                {
                    opt_phase = 0;
                    innerProd = innerProd_xp  / sqrt(innerProd_pp);
                }
                else if (a1 == 0)
                {
                    opt_phase = (double)(pi/2);
                    innerProd = -innerProd_xq / sqrt(innerProd_qq);
                }
                else if ( (a1!=0) && (k != 0) )
                {
                    opt_phase = atan( -(b1/a1) );
                    innerProd = (a1/fabs(a1))*(innerProd_xp*a1 + innerProd_xq*b1) /
                                sqrt(a1*a1*innerProd_pp +
                                     b1*b1*innerProd_qq +
                                     2*a1*b1*innerProd_pq);
                }

                xi = (k*2*pi)/N;

                if (fabs(innerProd)>fabs(maxInnerProd))
                {
                    maxInnerProd = innerProd;
                    expParm.innerProd = maxInnerProd;
                    expParm.rho = -(double)1/s;
                    expParm.xi = xi;
                    expParm.phase = opt_phase ;
                    expParm.a = 0.0;
                    expParm.b = tau;
                }

                // For debug
                fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
                                  innerProd, innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
                                  innerProd_pq, s, (double)tau, xi,opt_phase);
            }
        }

        s=s*2;
    }

    // s = N Aprox. as Pure Sinusoid
    tau =0;

    ((CExpParm*)parm)->rho = 0.0;
    ((CExpParm*)parm)->xi = 0.0;
    ((CExpParm*)parm)->phase = 0.0;
    ((CExpParm*)parm)->a = 0;
    ((CExpParm*)parm)->b = N-1;
    setRealAtom(parm);

    for (i=0;i<N;i++)
    {
        in1[i] = residue[0][i] * m_realAtom[i];
        in2[i] = m_realAtom[i] * m_realAtom[i];
    }

    fftw_execute ( plan_forward1 );

    fftw_execute ( plan_forward2 );

    C = out2[0][0];
    k_pi=NC-1;
    int delta_k = N/s;
    for (k=0;k<NC;k=k+delta_k)
    {
        innerProd_xp = out1[k][0];  // /sqrt((double)N);
        innerProd_xq = -out1[k][1]; // /sqrt((double)N);

        if ( (2*k)<NC )
        {
            innerProd_pp = 0.5 * (C + out2[2*k][0]);// /sqrt((double)N) );
            innerProd_qq = 0.5 * (C - out2[2*k][0]);// /sqrt((double)N) );
            innerProd_pq = -0.5 * (out2[2*k][1]);// /sqrt((double)N) );
        }
        else
        {
            innerProd_pp = 0.5 * (C + out2[N-(2*k)][0]);// /sqrt((double)N) );
            innerProd_qq = 0.5 * (C - out2[N-(2*k)][0]);// /sqrt((double)N) );
            innerProd_pq = -0.5 * ( - out2[N-(2*k)][1]);// /sqrt((double)N) );
        }

        a1 = innerProd_xp * innerProd_qq - innerProd_xq * innerProd_pq;
        b1 = innerProd_xq * innerProd_pp - innerProd_xp * innerProd_pq;

        if ( (k == 0) || (k == k_pi) )
        {
            opt_phase = 0;
            innerProd = innerProd_xp / sqrt(innerProd_pp);
        }
        else if (a1 == 0)
        {
            opt_phase = (double)(pi/2);
            innerProd = -innerProd_xq / sqrt(innerProd_qq);
        }
        else if ( (a1!=0) && (k != 0) )
        {
            opt_phase = atan( -(b1/a1) );
            innerProd = (a1/fabs(a1))*(innerProd_xp*a1 + innerProd_xq*b1) /
                        sqrt(a1*a1*innerProd_pp +
                             b1*b1*innerProd_qq +
                             2*a1*b1*innerProd_pq);
        }

        xi = (k*2*pi)/N;

        if (fabs(innerProd)>fabs(maxInnerProd))
        {
            maxInnerProd = innerProd;
            expParm.innerProd = maxInnerProd;
            expParm.rho = 0.0;
            expParm.xi = xi;
            expParm.phase = opt_phase ;
            expParm.a = tau;
            expParm.b = N-1;
        }
        // For debug
        fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
                          innerProd, innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
                          innerProd_pq, s, (double)tau, xi,opt_phase);
    }

    fclose(stream);

    fftw_destroy_plan ( plan_forward1 );
    fftw_destroy_plan ( plan_forward2 );

    fftw_free ( in1 );
    fftw_free ( out1 );
    fftw_free ( in2 );
    fftw_free ( out2 );

    delete parm;

    return expParm;
}


//===============================================================
// Function: fastMPKolasaModified
// Goal:
// Return:
//==============================================================

void CExpDictionary::fastMPKolasaModified(  cgMatrix<double>& residue,
                                        	CDataSignal* dataSignal,
                                            CFileDictionary* dicData,
                                            CParameter* chosenParm)
{

    CParameter* parm;
    parm = new CExpParm;
    CExpParm expParm;
    expParm.innerProd = 0;
    expParm.rho = 0;
    expParm.xi = 0;
    expParm.phase = 0;
    expParm.a = 0;
    expParm.b = 0;
    //double            *in_vec1, *in_vec2, *out_vec3;
    fftw_complex    *z1, *z2;
    fftw_complex    *w1, *w2;
    fftw_complex    *z1_fft, *z2_fft, *w1_fft, *w2_fft ;
    fftw_complex    *conv_zw1, *conv_zw2, *conv_zxi0_w2, *conv_zxi0_w2_expinc;
    fftw_complex    *prod_cvec_zw1, *prod_cvec_zw2;

    fftw_plan       plan_forward1a, plan_forward1b, plan_backward1,
                    plan_forward2a, plan_forward2b, plan_backward2;

    int N = m_signalSize;

    //FILE* stream;
    //  ---------------------------------------------------
    //stream = fopen( "results_modified_kolasa.out", "w" );

    // Memory allocation
    z1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    z2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    //w1 = (double*) fftw_malloc ( sizeof ( double ) * (2*N));
    //w2 = (double*) fftw_malloc ( sizeof ( double ) * (2*N));
    w1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    w2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    z1_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    z2_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    w1_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    w2_fft = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    conv_zw1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    conv_zw2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    conv_zxi0_w2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    conv_zxi0_w2_expinc = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));

    prod_cvec_zw1 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));
    prod_cvec_zw2 = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * (2*N));

    // PLAN FFT
    plan_forward1a = fftw_plan_dft_1d ( 2*N, z1, z1_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_forward1b = fftw_plan_dft_1d ( 2*N, w1, w1_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    //plan_forward1b = fftw_plan_dft_r2c_1d ( 2*N, w1,  w1_fft, FFTW_ESTIMATE);
    plan_backward1 = fftw_plan_dft_1d ( 2*N, prod_cvec_zw1, conv_zw1, FFTW_BACKWARD, FFTW_ESTIMATE);

    plan_forward2a = fftw_plan_dft_1d ( 2*N, z2, z2_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_forward2b = fftw_plan_dft_1d ( 2*N, w2, w2_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    //plan_forward2b = fftw_plan_dft_r2c_1d ( 2*N, w2, w2_fft,FFTW_ESTIMATE );
    plan_backward2 = fftw_plan_dft_1d ( 2*N, prod_cvec_zw2, conv_zw2, FFTW_BACKWARD, FFTW_ESTIMATE);

    double innerProd_xp;
    double innerProd_xq;
    double innerProd_pp;
    double innerProd_qq;
    double innerProd_pq;
    double a1, b1;
    double innerProd = 0;
    double maxInnerProd = 0;
    double opt_phase;

    int i,j;
    int k;
    int s;
    double xi;
    int delta_tau;
    int tau;
    double Fs;
    double Ffund;
    double delta_f;
    int Nfreq;
    double* xi_vec;
    double freqi,freqf;
    int fdiscrtype;

    //dicData->printToScreen();

    for(j=0;j<dicData->getNumDicBlock();j++)
    {
        // Using the dictionary data
        s = (dicData->getDicBlock())[j].scale;
        delta_tau = (dicData->getDicBlock())[j].delta_tau;
        fdiscrtype = (dicData->getDicBlock())[j].fdiscrtype;
        freqi = (dicData->getDicBlock())[j].freqi;
        freqf = (dicData->getDicBlock())[j].freqf;

        //printf(" %d %d %d %f %f \n",s,delta_tau,fdiscrtype,freqi,freqf);




        if ((s > m_signalSize)&&(s!=88888))
        {
            cout << "Scale greater than block size!!!" << endl;
            printf("Scale %d - block size %d\n",s,m_signalSize);
            exit(1);
        }
        //////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////

        if ( (freqf==9999999999) && (freqi==9999999999) )
        {
            if (dataSignal->getType() == 1)
            {
                Fs = ((CComtradeSignal*)dataSignal)->getSamplingRate(1);
                freqi = ((CComtradeSignal*)dataSignal)->getFundamentalFrequency();
                freqf = Fs;
            }
            else if (dataSignal->getType() == 2) // Audio
            {
                Fs = ((CAudioSignal*)dataSignal)->getSamplingRate();
                freqi =  27.5; // A0 = 27.5 Hz
                freqf = Fs;
                //double A0 = 27.5; // Hz
                //double C8 = 4186.01;
                //double freqi= A0 * pow(2.0,-23/12);
                //double freqf= C8 * pow(2.0,23/12); // B9
            }
            else if (dataSignal->getType() == 3) // Noise for Audio
            {
                Fs = ((CNoiseSignal*)dataSignal)->getSamplingRate();
                freqi =  Fs/100;
                freqf = Fs;
            }

        }
        else if (freqf==9999999999)
        {
            if (dataSignal->getType() == 1)
            {
                Fs = ((CComtradeSignal*)dataSignal)->getSamplingRate(1);
                freqf = Fs;
            }
            else if (dataSignal->getType() == 2) // Audio
            {
                Fs = ((CAudioSignal*)dataSignal)->getSamplingRate();
                freqf = Fs;
            }
            else if (dataSignal->getType() == 3) // Noise for Audio
            {
                Fs = ((CNoiseSignal*)dataSignal)->getSamplingRate();
                freqf = Fs;
            }

        }

        if (freqf>Fs)
        {
            cout << "Final frequency greater than sampling frequency !!!" << endl;
            printf("final freq. %f - Fs %f\n",freqf,Fs);
            exit(1);
        }

        if (fdiscrtype==1) // linear
        {
            delta_f = (2*pi/freqf)* freqi;
            Nfreq = (int)(freqf/(2*freqi));
            //Nfreq = (int)ceil(freqf/freqi);
            xi_vec = new double[Nfreq];
            for (i=0;i<Nfreq;i++)
            {
                xi_vec[i] = (2*pi/Fs) * (freqi * i );
            }
        }
        else if (fdiscrtype==2) // geometric with quarter-tone discretization
        {


            Nfreq = (int)ceil(24 * ( log10(freqf/freqi)/log10(2.0) ) )+1;

            xi_vec = new double[Nfreq];
            xi_vec[0] = 0.0;
            for (i=1;i<Nfreq;i++)
            {
                xi_vec[i] = (2*pi/Fs) * (freqi * pow (2.0, (double)(i-1)/24) );
            }
        }

//        // Set linear frequencies for COMTRADE files
//        if (dataSignal->getType() == 1) // Comtrade
//        {
//            if ( (freqf==9999999999) && (freqi==9999999999) )
//            {
//                Ffund = ((CComtradeSignal*)dataSignal)->getFundamentalFrequency();
//                Fs = ((CComtradeSignal*)dataSignal)->getSamplingRate(1);
//            }
//            else if (freqf==9999999999)
//            {
//                Ffund = freqi;
//                Fs = ((CComtradeSignal*)dataSignal)->getSamplingRate(1);
//            }
//            delta_f = (2*pi/Fs)* Ffund;
//            Nfreq = (int)(Fs/(2*Ffund));
//            xi_vec = new double[Nfreq];
//            for (i=0;i<Nfreq;i++)
//            {
//                xi_vec[i] = (double)i * delta_f;
//            }
//        }
//        // Set linear/geometric frequencies space for audio or noise
//        if ( (dataSignal->getType() == 2) || // Audio
//             (dataSignal->getType() == 3) )   // Noise for Audio
//        {
//            if (dataSignal->getType() == 2) Fs = ((CAudioSignal*)dataSignal)->getSamplingRate();
//            if (dataSignal->getType() == 3) Fs = ((CNoiseSignal*)dataSignal)->getSamplingRate();
//            //double A0 = 27.5; // Hz
//            //double C8 = 4186.01;
//            //double freqi= A0 * pow(2.0,-23/12);
//            //double freqf= C8 * pow(2.0,23/12); // B9
//
//            if (fdiscrtype==1) // linear
//            {
//                Nfreq = (int)ceil(freqf/freqi);
//                xi_vec = new double[Nfreq];
//                for (i=0;i<Nfreq;i++)
//                {
//                    xi_vec[i] = (2*pi/Fs) * (freqi * i );
//                }
//            }
//            else if (fdiscrtype==2) // geometric with quarter-tone discretization
//            {
//
//
//                Nfreq = (int)ceil(24 * ( log10(freqf/freqi)/log10(2.0) ) )+1;
//
//                xi_vec = new double[Nfreq];
//                xi_vec[0] = 0.0;
//                for (i=1;i<Nfreq;i++)
//                {
//                    xi_vec[i] = (2*pi/Fs) * (freqi * pow (2.0, (double)(i-1)/24) );
//                }
//            }
//        }
        //////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////
        if (s==1)
        {

            // s=1 Impulse
            for(i=0;i<N;i+=delta_tau)
            {
                innerProd = residue[0][i];
                if (fabs(innerProd)>fabs(maxInnerProd))
                {
                    maxInnerProd = innerProd;
                    expParm.innerProd = maxInnerProd;
                    expParm.rho = 0.0;
                    expParm.xi = 0.0;
                    expParm.phase = 0.0 ;
                    expParm.a = i;
                    expParm.b = i;

                }
                //fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
                //                  innerProd, innerProd, 0.0, 1.0, 0.0,//innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
                //                  0.0, 1.0, (double)i, 0.0,0.0);
            }
        }
        else if (s==88888)
        {
            // s = N Aprox. as Pure Sinusoid
            cgMatrix<double> innerProduct;
            //k_pi = (int) s/2.0;
            //for(k=0;k<(int)s;k++)
            for(k=0;k<Nfreq;k++)
            {
                //xi = k * ((2*pi)/s);
                xi = xi_vec[k];

                ((CExpParm*)parm)->rho      = 0.0;
                ((CExpParm*)parm)->xi       = xi;
                ((CExpParm*)parm)->phase    = 0.0;
                ((CExpParm*)parm)->a        = 0;
                ((CExpParm*)parm)->b        = N-1;
                setComplexAtom(parm);

                // Computing optimum phase
                opt_phase = computeOptimumPhase(    residue,
                                                    xi,
                                                    innerProd);
        //      opt_phase= computeOptimumPhase( residue,
        //                                      xi);

        //      ((CExpParm*)parm)->rho      = 0.0;
        //      ((CExpParm*)parm)->xi       = xi;
        //      ((CExpParm*)parm)->phase    = opt_phase;
        //      ((CExpParm*)parm)->a        = 0;
        //      ((CExpParm*)parm)->b        = N-1;
        //      setRealAtom(parm);

        //      innerProduct = residue * m_realAtom;
            //  innerProd = innerProduct[0][0];
                //printf(" %12.10f %12.10f %12.10f %5d\n",innerProd,0.0,xi,0);
                if (fabs(innerProd)>fabs(maxInnerProd))
                {
                    maxInnerProd = innerProd;
                    expParm.innerProd = maxInnerProd;
                    expParm.rho = 0.0;
                    expParm.xi = xi;
                    expParm.phase = opt_phase ;
                    expParm.a = 0;
                    expParm.b = N-1;
                }

            }

        }
        else // s =[2,4,...,N/2] Damped Sinusoids
        {
            for (i=0;i<2*N;i++)
            {
                z1[i][0] = 0.0;
                z1[i][1] = 0.0;
                z2[i][0] = 0.0;
                z2[i][1] = 0.0;
                w1[i][0] = 0.0;
                w1[i][1] = 0.0;
                w2[i][0] = 0.0;
                w2[i][1] = 0.0;
            }
            double tol = 0;//1e-10;
            //  int k_pi;


            ((CExpParm*)parm)->rho      = 1.0/(double)s;
            ((CExpParm*)parm)->xi       = 0.0;
            ((CExpParm*)parm)->phase    = 0.0;
            ((CExpParm*)parm)->a        = 0;
            ((CExpParm*)parm)->b        = N-1;
            setRealAtom(parm);

            for (i=0;i<N;i++)
            {
                w1[i][0] = m_realAtom[N-1-i];
                w2[i][0] = m_realAtom[N-1-i] * m_realAtom[N-1-i];
            }
            //k_pi = (int) s/2.0;
            //for(k=0;k<(int)s;k++)
            for(k=0;k<Nfreq;k++)
            {
                //xi = k * ((2*pi)/s);
                xi = xi_vec[k];

                if ( (xi==0.0) || (xi>=((2*pi)/s) ))
                {

                    // z1
                    ((CExpParm*)parm)->rho      = 0.0;
                    ((CExpParm*)parm)->xi       = xi;
                    ((CExpParm*)parm)->phase    = 0.0;
                    ((CExpParm*)parm)->a        = 0;
                    ((CExpParm*)parm)->b        = N-1;
                    setComplexAtom(parm);

                    for (i=0;i<N;i++)
                    {
                        z1[i][0] = m_complexAtom[i].Real() * residue[0][i];
                        z1[i][1] = m_complexAtom[i].Imag() * residue[0][i];
                    }

                    // z2
                    ((CExpParm*)parm)->rho      = 0.0;
                    ((CExpParm*)parm)->xi       = 2*xi;
                    ((CExpParm*)parm)->phase    = 0.0;
                    ((CExpParm*)parm)->a        = 0;
                    ((CExpParm*)parm)->b        = N-1;
                    setComplexAtom(parm);

                    for (i=0;i<N;i++)
                    {
                        z2[i][0] = m_complexAtom[i].Real();
                        z2[i][1] = m_complexAtom[i].Imag();
                    }

                    // zw1
                    fftw_execute (plan_forward1a);
                    fftw_execute (plan_forward1b);

                    // product between complex vectors
                    for (i=0;i<2*N;i++)
                    {
                        prod_cvec_zw1[i][0]= (z1_fft[i][0]*w1_fft[i][0] - z1_fft[i][1]*(w1_fft[i][1]));
                        prod_cvec_zw1[i][1]= (z1_fft[i][1]*w1_fft[i][0] + z1_fft[i][0]*(w1_fft[i][1]));
                        //cout << prod_cvec_zw1[i][0]<< " + j* " << prod_cvec_zw1[i][1] << endl;
                    }

                    fftw_execute ( plan_backward1);

                    // zw2
                    fftw_execute (plan_forward2a);
                    fftw_execute (plan_forward2b);

                    // product between complex vectors
                    for (i=0;i<2*N;i++)
                    {
                        prod_cvec_zw2[i][0]= (z2_fft[i][0]*w2_fft[i][0] - z2_fft[i][1]*(w2_fft[i][1]));
                        prod_cvec_zw2[i][1]= (z2_fft[i][1]*w2_fft[i][0] + z2_fft[i][0]*(w2_fft[i][1]));
                    }

                    fftw_execute ( plan_backward2);

                    if (k==0) // (xi==0)
                    {
                        memcpy(conv_zxi0_w2, conv_zw2, sizeof(fftw_complex) * (2*N));
                    }

                    // Compute inner product and optimum phase
                    //for (int tau=-N+1; tau<N+1; tau+=delta_tau)
                    for (tau=0; tau<(double)N; tau=tau+delta_tau)
                    {
                        innerProd_xp = conv_zw1[tau+N-1][0]/(double)(2*N);
                        innerProd_xq = conv_zw1[tau+N-1][1]/(double)(2*N);
                        innerProd_pp = 0.5*(conv_zxi0_w2[tau+N-1][0] + conv_zw2[tau+N-1][0])/(double)(2*N);
                        innerProd_qq = 0.5*(conv_zxi0_w2[tau+N-1][0] - conv_zw2[tau+N-1][0])/(double)(2*N);
                        innerProd_pq = 0.5*(conv_zw2[tau+N-1][1]/(double)(2*N));

                        a1 = innerProd_xp * innerProd_qq - innerProd_xq * innerProd_pq;
                        b1 = innerProd_xq * innerProd_pp - innerProd_xp * innerProd_pq;

                        if (( k == 0 )||( (int)xi*10000 == (int)pi*10000))
                        {
                            opt_phase = 0;
                            if (fabs(innerProd_pp) > tol)
                            {
                                innerProd = innerProd_xp / sqrt(innerProd_pp);
                            }
                            else
                            {
                                innerProd = 0.0;
                                cout << " 1 - " << innerProd_pp << endl;
                            }
                        }
                        else if (a1 == 0)
                        {
                            opt_phase = (double)(pi/2);
                            if (fabs(innerProd_qq) > tol)
                            {
                                innerProd = -innerProd_xq / sqrt(innerProd_qq);
                            }
                            else
                            {
                                innerProd = 0.0;
                                cout << " 2 - " << innerProd_qq  << endl;
                            }
                        }
                        else //if ( (a1 != 0) && (xi != 0) )
                        {
                            opt_phase = atan( -(b1/a1) );
                            if ( fabs(a1*a1*innerProd_pp + b1*b1*innerProd_qq + 2*a1*b1*innerProd_pq) > tol)
                            {
                                innerProd =  (a1/fabs(a1))*(innerProd_xp*a1 + innerProd_xq*b1) /
                                             sqrt(a1*a1*innerProd_pp +
                                             b1*b1*innerProd_qq +
                                             2*a1*b1*innerProd_pq);
                            }
                            else
                            {
                                innerProd = 0.0;
                                cout << " 3 - " << (a1*a1*innerProd_pp + b1*b1*innerProd_qq + 2*a1*b1*innerProd_pq) << endl;
                            }

                        }

                        if (fabs(innerProd)>fabs(maxInnerProd))
                        {
                            //printf("Decreasing\n");
                            //printf(" %12.10f %12.10f %12.10f %5d\n",innerProd,1.0/(double)s,xi,tau);
                            maxInnerProd = innerProd;
                            expParm.innerProd = maxInnerProd;
                            expParm.rho = 1.0/(double)s;
                            expParm.xi = xi;
                            expParm.phase = opt_phase ;
                            expParm.a = tau;
                            expParm.b = N-1;
                        }
                        // For debug
                        //"%+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E\n"
                        //fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
                        //                  innerProd, innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
                        //                  innerProd_pq, s, (double)tau, xi,opt_phase);//, conv_zw1[tau+N-1][0]/(double)(2*N),
                                          //conv_zw1[tau+N-1][1]/(double)(2*N), conv_zw2[tau+N-1][0]/(double)(2*N), conv_zw2[tau+N-1][1]/(double)(2*N), conv_zxi0_w2[tau+N-1][0]/(double)(2*N),  z1[tau+N-1][0],
                                          //z1[tau+N-1][1], w1[tau+N-1], z2[tau+N-1][0], z2[tau+N-1][1], w2[tau+N-1]);
                    }

                    // Increasing exponential

                    // product between complex vectors (z1_fft e  conjugate(w1_fft))
                    for (i=0;i<2*N;i++)
                    {
                        prod_cvec_zw1[i][0]= (z1_fft[i][0]*w1_fft[i][0] - z1_fft[i][1]*(-w1_fft[i][1]));
                        prod_cvec_zw1[i][1]= (z1_fft[i][1]*w1_fft[i][0] + z1_fft[i][0]*(-w1_fft[i][1]));
                    }

                    fftw_execute ( plan_backward1);

                    // product between complex vectors (z2_fft e  conjugate(w2_fft))
                    for (i=0;i<2*N;i++)
                    {
                        prod_cvec_zw2[i][0]= (z2_fft[i][0]*w2_fft[i][0] - z2_fft[i][1]*(-w2_fft[i][1]));
                        prod_cvec_zw2[i][1]= (z2_fft[i][1]*w2_fft[i][0] + z2_fft[i][0]*(-w2_fft[i][1]));
                    }

                    fftw_execute ( plan_backward2);

                    if (k==0) // (xi==0)
                    {
                        memcpy(conv_zxi0_w2_expinc, conv_zw2, sizeof(fftw_complex) * (2*N));
                    }

                    // Compute inner product and optimum phase
                    for (tau=N-1; tau>0; tau=tau-delta_tau)
                    {
                        innerProd_xp = conv_zw1[(tau+N+1)%(2*N)][0]/(double)(2*N);
                        innerProd_xq = conv_zw1[(tau+N+1)%(2*N)][1]/(double)(2*N);
                        innerProd_pp = 0.5*(conv_zxi0_w2_expinc[(tau+N+1)%(2*N)][0] + conv_zw2[(tau+N+1)%(2*N)][0])/(double)(2*N);
                        innerProd_qq = 0.5*(conv_zxi0_w2_expinc[(tau+N+1)%(2*N)][0] - conv_zw2[(tau+N+1)%(2*N)][0])/(double)(2*N);
                        innerProd_pq = 0.5*(conv_zw2[(tau+N+1)%(2*N)][1]/(double)(2*N));

                        a1 = innerProd_xp * innerProd_qq - innerProd_xq * innerProd_pq;
                        b1 = innerProd_xq * innerProd_pp - innerProd_xp * innerProd_pq;

                        if (( k == 0 )||( (int)xi*10000 == (int)pi*10000))
                        {
                            opt_phase = 0;
                            if (fabs(innerProd_pp) > tol)
                            {
                                innerProd = innerProd_xp / sqrt(innerProd_pp);
                            }
                            else
                            {
                                innerProd = 0.0;
                                cout << " 1 - " << innerProd_pp << endl;
                            }
                        }
                        else if (a1 == 0)
                        {
                            opt_phase = (double)(pi/2);
                            if (fabs(innerProd_qq) > tol)
                            {
                                innerProd = -innerProd_xq / sqrt(innerProd_qq);
                            }
                            else
                            {
                                innerProd = 0.0;
                                cout << " 2 - " << innerProd_qq  << endl;
                            }
                        }
                        else //if ( (a1 != 0) && (xi != 0) )
                        {
                            opt_phase = atan( -(b1/a1) );
                            if ( fabs(a1*a1*innerProd_pp + b1*b1*innerProd_qq + 2*a1*b1*innerProd_pq) > tol)
                            {
                                innerProd =  (a1/fabs(a1))*(innerProd_xp*a1 + innerProd_xq*b1) /
                                             sqrt(a1*a1*innerProd_pp +
                                             b1*b1*innerProd_qq +
                                             2*a1*b1*innerProd_pq);
                            }
                            else
                            {
                                innerProd = 0.0;
                                cout << " 3 - " << (a1*a1*innerProd_pp + b1*b1*innerProd_qq + 2*a1*b1*innerProd_pq) << endl;
                            }

                        }
                        if (fabs(innerProd)>fabs(maxInnerProd))
                        {
                            //printf("Increasing\n");
                            //printf(" %12.10f %12.10f %12.10f %5d\n",innerProd,1.0/(double)s,xi,tau);
                            maxInnerProd = innerProd;
                            expParm.innerProd = maxInnerProd;
                            expParm.rho = -1.0/(double)s;
                            expParm.xi = xi;
                            expParm.phase = opt_phase ;
                            expParm.a = 0;
                            expParm.b = tau;
                        }
                        // For debug
                        //"%+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E %+10.5E\n"
                        //fprintf( stream, "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
                        //                 innerProd, innerProd_xp, innerProd_xq, innerProd_pp, innerProd_qq,
                        //                  innerProd_pq, s, (double)tau, xi,opt_phase);//, conv_zw1[tau+N-1][0]/(double)(2*N),
                                          //conv_zw1[tau+N-1][1]/(double)(2*N), conv_zw2[tau+N-1][0]/(double)(2*N), conv_zw2[tau+N-1][1]/(double)(2*N), conv_zxi0_w2[tau+N-1][0]/(double)(2*N),  z1[tau+N-1][0],
                                          //z1[tau+N-1][1], w1[tau+N-1], z2[tau+N-1][0], z2[tau+N-1][1], w2[tau+N-1]);
                    }
                }
            }
            //s=s*2;
        }
        if (xi_vec!=NULL)
        {
            delete [] xi_vec;
        }
    }


    ((CExpParm*)chosenParm)->innerProd = expParm.innerProd;
	((CExpParm*)chosenParm)->rho = expParm.rho;
	((CExpParm*)chosenParm)->xi = expParm.xi;
	((CExpParm*)chosenParm)->phase = expParm.phase;
	((CExpParm*)chosenParm)->a = expParm.a;
	((CExpParm*)chosenParm)->b = expParm.b;

    //fclose(stream);


    fftw_destroy_plan ( plan_forward1a );
    fftw_destroy_plan ( plan_forward1b);
    fftw_destroy_plan ( plan_backward1 );
    fftw_destroy_plan ( plan_forward2a );
    fftw_destroy_plan ( plan_forward2b);
    fftw_destroy_plan ( plan_backward2 );

    fftw_free(z1);
    fftw_free(z2);
    fftw_free(w1);
    fftw_free(w2);
    fftw_free(z1_fft);
    fftw_free(z2_fft);
    fftw_free(w1_fft);
    fftw_free(w2_fft);
    fftw_free(conv_zw1);
    fftw_free(conv_zw2);
    fftw_free(conv_zxi0_w2);
    fftw_free(conv_zxi0_w2_expinc);
    fftw_free(prod_cvec_zw1);
    fftw_free(prod_cvec_zw2);

    delete parm;
}


//===============================================================
// Function: computeOptimumPhase (With continuous parameters)
// Goal:
// Return:
//===============================================================

double CExpDictionary::computeOptimumPhase( cgMatrix<double>& residue,
                                            double xi)
{
    double opt_phase = 0;
    int i;
    cgMatrix<double> realPartComplexDic(1,m_signalSize,0.0);
    cgMatrix<double> imagPartComplexDic(1,m_signalSize,0.0);

    for ( i=0;i<m_signalSize;i++)
    {
        realPartComplexDic[0][i] = m_complexAtom[i].Real();
        imagPartComplexDic[0][i] = m_complexAtom[i].Imag();
    }

    double innerProdReal=0;
    double innerProdImag=0;
    double innerProdRealImag=0;
    for (i=0;i< m_signalSize;i++)
    {
        innerProdReal += residue.getData(0,i) * m_complexAtom[i].Real();
        innerProdImag += residue.getData(0,i) * m_complexAtom[i].Imag();
        innerProdRealImag += m_complexAtom[i].Real() * m_complexAtom[i].Imag();
    }

    cgMatrix<double> innerProductReal(1,1,innerProdReal);
    cgMatrix<double> innerProductImag(1,1,innerProdImag);
    //innerProductReal = residue * (realPartComplexDic.transpose());
    //innerProductImag = residue * (imagPartComplexDic.transpose());

    double p,q ;
    p = realPartComplexDic.norm();
    q = imagPartComplexDic.norm();

    cgMatrix<double> innerProductRealImag(1,1,innerProdRealImag);
    //innerProductRealImag = realPartComplexDic * (imagPartComplexDic.transpose());

    cgMatrix<double> a1;
    cgMatrix<double> b1;
    a1 = innerProductReal*(q*q) - innerProductImag * innerProductRealImag;
    b1 = innerProductImag*(p*p) - innerProductReal * innerProductRealImag;

    double innerProd;
    if ( (xi == 0) ||
        ((int)(10000*xi) == (int)(10000*pi))) //caso n�o haja sen�ide
    {
        opt_phase = 0;
        innerProd = innerProductReal[0][0]/p;
    }
    else if (a1.getData(0,0) == 0)
    {
        opt_phase = (double)(pi/2);
        innerProd = -innerProductImag[0][0]/q;
    }
    else if (   (a1.getData(0,0)!=0) && (xi!=0) )
    {
        opt_phase = atan( -(b1.getData(0,0)/a1.getData(0,0)) );
        innerProd = (a1[0][0]/fabs(a1[0][0]))*(innerProductReal[0][0]*a1[0][0] + innerProductImag[0][0]*b1[0][0])/
            sqrt(a1[0][0]*a1[0][0]*p*p+b1[0][0]*b1[0][0]*q*q+2*a1[0][0]*b1[0][0]*innerProductRealImag[0][0]);
    }

//  FILE* stream;
//  ---------------------------------------------------
/*  stream = fopen( "results_modified_kolasa.out", "a" );
    fprintf(stream,"%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
                                                    innerProd,
                                                    innerProductReal[0][0],
                                                    innerProductImag[0][0],
                                                    p*p,
                                                    q*q,
                                                    innerProductRealImag[0][0],
                                                    (double)m_signalSize,//pow(2.0, (double)(expDic.getGammaExp())[atomIndex].j),
                                                    0.0,//((expDic.getGammaExp())[atomIndex].p)*pow(2.0, (double)(expDic.getGammaExp())[atomIndex].j-1),
                                                    xi,//(expDic.getGammaExp())[atomIndex].k*2*pi*pow(2.0, -(double)(expDic.getGammaExp())[atomIndex].j));
                                                    opt_phase);
    fclose(stream);
*/

    return opt_phase;
}

//===============================================================
// Function: computeOptimumPhase (With continuous parameters)
// Goal:
// Return: opt_phase and innerProd
//===============================================================

double CExpDictionary::computeOptimumPhase( cgMatrix<double>& residue,
                                            double xi,
                                            double& innerProd)
{
    double opt_phase = 0;
    int i;
    cgMatrix<double> realPartComplexDic(1,m_signalSize,0.0);
    cgMatrix<double> imagPartComplexDic(1,m_signalSize,0.0);

    for ( i=0;i<m_signalSize;i++)
    {
        realPartComplexDic[0][i] = m_complexAtom[i].Real();
        imagPartComplexDic[0][i] = m_complexAtom[i].Imag();
    }

    double innerProdReal=0;
    double innerProdImag=0;
    double innerProdRealImag=0;
    for (i=0;i< m_signalSize;i++)
    {
        innerProdReal += residue.getData(0,i) * m_complexAtom[i].Real();
        innerProdImag += residue.getData(0,i) * m_complexAtom[i].Imag();
        innerProdRealImag += m_complexAtom[i].Real() * m_complexAtom[i].Imag();
    }

    cgMatrix<double> innerProductReal(1,1,innerProdReal);
    cgMatrix<double> innerProductImag(1,1,innerProdImag);
    //innerProductReal = residue * (realPartComplexDic.transpose());
    //innerProductImag = residue * (imagPartComplexDic.transpose());

    double p,q ;
    p = realPartComplexDic.norm();
    q = imagPartComplexDic.norm();

    cgMatrix<double> innerProductRealImag(1,1,innerProdRealImag);
    //innerProductRealImag = realPartComplexDic * (imagPartComplexDic.transpose());

    cgMatrix<double> a1;
    cgMatrix<double> b1;
    a1 = innerProductReal*(q*q) - innerProductImag * innerProductRealImag;
    b1 = innerProductImag*(p*p) - innerProductReal * innerProductRealImag;

    if ( (xi == 0) ||
        ((int)(10000*xi) == (int)(10000*pi))) //caso n�o haja sen�ide
    {
        opt_phase = 0;
        innerProd = innerProductReal[0][0]/p;
    }
    else if (a1.getData(0,0) == 0)
    {
        opt_phase = (double)(pi/2);
        innerProd = -innerProductImag[0][0]/q;
    }
    else if (   (a1.getData(0,0)!=0) && (xi!=0) )
    {
        opt_phase = atan( -(b1.getData(0,0)/a1.getData(0,0)) );
        innerProd = (a1[0][0]/fabs(a1[0][0]))*(innerProductReal[0][0]*a1[0][0] + innerProductImag[0][0]*b1[0][0])/
            sqrt(a1[0][0]*a1[0][0]*p*p+b1[0][0]*b1[0][0]*q*q+2*a1[0][0]*b1[0][0]*innerProductRealImag[0][0]);
    }

//  FILE* stream;
//  ---------------------------------------------------
/*  stream = fopen( "results_modified_kolasa.out", "a" );
    fprintf(stream,"%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
                                                    innerProd,
                                                    innerProductReal[0][0],
                                                    innerProductImag[0][0],
                                                    p*p,
                                                    q*q,
                                                    innerProductRealImag[0][0],
                                                    (double)m_signalSize,//pow(2.0, (double)(expDic.getGammaExp())[atomIndex].j),
                                                    0.0,//((expDic.getGammaExp())[atomIndex].p)*pow(2.0, (double)(expDic.getGammaExp())[atomIndex].j-1),
                                                    xi,//(expDic.getGammaExp())[atomIndex].k*2*pi*pow(2.0, -(double)(expDic.getGammaExp())[atomIndex].j));
                                                    opt_phase);
    fclose(stream);
*/

    return opt_phase;
}


//===============================================================
// Function: optimizeContinuousParms
// Goal:
// Return:
//===============================================================
void CExpDictionary::optimizeContinuousParms(   cgMatrix<double>& residue,
                                                CParameter* parm)
{
    CParameter* parm_aux;
    parm_aux = new CExpParm;
    if (((CExpParm*)parm)->a != ((CExpParm*)parm)->b) // NOT impulse case; if impulse case, do nothing
    {

        double rho, u, xi, phase;

        // decaying factor
        rho = ((CExpParm*)parm)->rho;
        // time shift
        if (rho>=0.0)
            u = ((CExpParm*)parm)->a;
        else
            u = ((CExpParm*)parm)->b;
        // frequency
        xi = ((CExpParm*)parm)->xi;

        phase = ((CExpParm*)parm)->phase;

        // Define inner product
        double highestProduct;
        highestProduct = ((CExpParm*)parm)->innerProd;

        // =======================================
        double increment_rho, rho_aux, rho_fix;

        increment_rho =  rho / 2.0;  //((double)expDic.getSignalSize()/2.0);
        rho_aux = increment_rho + rho;
        rho_fix = rho;

        double increment_u, u_aux, u_fix;

        if (rho!=0.0)
            increment_u = 1.0 / (2.0 * rho);
        else
            increment_u = (double)m_signalSize / 2.0;
        u_aux = increment_u + u;
        u_fix = u;

        if(u_aux>(m_signalSize-1))
            u_aux = m_signalSize-1;
        if(u_aux<0) u_aux =0;

        double increment_xi, xi_aux, xi_fix;

        increment_xi = rho * (pi) ;
        xi_aux = increment_xi + xi;
        xi_fix = xi;

        // =========================
        // Pseudo-Newton
        cout << "Pseudo-Newton" << endl;

        int time = 0;

        double ph_fix = phase;

        double u_lim  = 1.0 /  (double) ( 2 * ( m_signalSize ) * ( m_signalSize ) );
        double xi_lim = 1.0 /  (double) ( 2 * ( m_signalSize ) * ( m_signalSize ) );
        double rho_lim  = 1.0 / (double) ( 2 * ( m_signalSize ) * ( m_signalSize ) );

        //cout << "rho_lim: " << rho_lim << "u_lim: " << u_lim << "xi_lim: " << xi_lim << endl;

        /*ofstream outOptFile("epss_opt_atom_cont_param.dat",ios::out);

        outOptFile << step << ' ' << s_fix << ' ' << u_fix  << ' '
            << xi_fix << ' ' << increment_s << ' ' << increment_u << ' '
            << increment_xi << '\n';
        */

        int incr_count=0;
        cgMatrix<double> innerProduct;
        double opt_phase = phase;

        while(  (fabs(increment_u)  > u_lim) ||
                (fabs(increment_xi) > xi_lim) ||
                (fabs(increment_rho)  > rho_lim) )
        {
            if (time==0) // xi_fix, rho_aux, u_fix
            {
                ((CExpParm*)parm_aux)->rho      = rho_aux;
                ((CExpParm*)parm_aux)->xi       = xi_fix;
                ((CExpParm*)parm_aux)->phase    = 0.0;
                if (rho>=0.0)
                {
                    ((CExpParm*)parm_aux)->a        = (int)u_fix;
                    ((CExpParm*)parm_aux)->b        = m_signalSize-1;
                }
                else
                {
                    ((CExpParm*)parm_aux)->a        = 0;
                    ((CExpParm*)parm_aux)->b        = (int)u_fix;
                }
                setComplexAtom(parm_aux);

                // Computing optimum phase
                opt_phase = computeOptimumPhase(residue,
                                                xi);

                ((CExpParm*)parm_aux)->phase    = opt_phase;
                setRealAtom(parm_aux);

            }
            if (time==1) // xi_aux, rho_fix, u_fix
            {
                ((CExpParm*)parm_aux)->rho      = rho_fix;
                ((CExpParm*)parm_aux)->xi       = xi_aux;
                ((CExpParm*)parm_aux)->phase    = 0.0;
                if (rho>=0.0)
                {
                    ((CExpParm*)parm_aux)->a        = (int)u_fix;
                    ((CExpParm*)parm_aux)->b        = m_signalSize-1;
                }
                else
                {
                    ((CExpParm*)parm_aux)->a        = 0;
                    ((CExpParm*)parm_aux)->b        = (int)u_fix;
                }
                setComplexAtom(parm_aux);

                // Computing optimum phase
                opt_phase = computeOptimumPhase(residue,
                                                xi);

                ((CExpParm*)parm_aux)->phase    = opt_phase;
                setRealAtom(parm_aux);

            }
            if (time==2) // xi_fix, rho_fix, u_aux
            {
                ((CExpParm*)parm_aux)->rho      = rho_fix;
                ((CExpParm*)parm_aux)->xi       = xi_fix;
                ((CExpParm*)parm_aux)->phase    = 0.0;
                if (rho>=0.0)
                {
                    ((CExpParm*)parm_aux)->a        = (int)u_aux;
                    ((CExpParm*)parm_aux)->b        = m_signalSize-1;
                }
                else
                {
                    ((CExpParm*)parm_aux)->a        = 0;
                    ((CExpParm*)parm_aux)->b        = (int)u_aux;
                }
                setComplexAtom(parm_aux);

                // Computing optimum phase
                opt_phase = computeOptimumPhase(residue,
                                                xi);

                ((CExpParm*)parm_aux)->phase    = opt_phase;
                setRealAtom(parm_aux);

            }

            //realAtomVector.fillVector(realAtom);

            //realAtomVector.PrintToFile("atom.dat");
            //residue.PrintToFile("residue.dat");

            innerProduct = residue * m_realAtom;//( realAtomVector.transpose() );

            //cout << "rho_fix: " << rho_fix << "u_fix: " << u_fix << "xi_fix: " << xi_fix << endl;
            //cout << "rho_aux: " << rho_aux << "u_aux: " << u_aux << "xi_aux: " << xi_aux << endl;
            //cout << "PP: " << innerProduct.getData(0,0) << endl;

            if( fabs(innerProduct.getData(0,0)) > fabs(highestProduct) )
            {
                // Greater Product
                highestProduct = innerProduct.getData(0,0);

                if(time==0)
                {
                    rho_fix = rho_aux;
                    rho_aux = rho_fix + increment_rho;
                    ph_fix = opt_phase;

                    increment_rho = pow(4.0,(double)incr_count) * increment_rho ;
                }
                if(time==1)
                {
                    xi_fix = xi_aux;
                    xi_aux = xi_fix + increment_xi;
                    ph_fix = opt_phase;

                    increment_xi = /*pow(2.0,(double)incr_count)*/2.0 * increment_xi ;
                }
                if(time==2)
                {
                    u_fix = u_aux;
                    u_aux = u_fix + increment_u;
                    ph_fix = opt_phase;

                    if (u_aux > (m_signalSize - 1) )
                        u_aux = m_signalSize - 1;
                    if (u_aux<0) u_aux=0;

                    increment_u = pow(2.0,(double)incr_count) * increment_u ;
                }

                incr_count++;
            }
            else
            {
                // Smaller product
                if(time==0)
                {
                    increment_rho = increment_rho/pow(4.0,(double)incr_count);
                    increment_rho= - (increment_rho/2.0);
                    rho_aux = rho_fix + increment_rho;
                }
                if(time==1)
                {
                    increment_xi = - (increment_xi/2.0);
                    xi_aux = xi_fix + increment_xi;
                }
                if(time==2)
                {
                    increment_u = - (increment_u/2.0);
                    u_aux = u_fix + increment_u - 0.5;

                    if (u_aux>(m_signalSize - 1))
                        u_aux = m_signalSize - 1;

                    if (u_aux<0) u_aux=0;
                }

                time++;
                time = time % 3;

                incr_count = 0;
            }

            if (fabs(rho_fix) < 1e-20) break;

            /*outOptFile << step << ' ' << s_fix << ' ' << u_fix  << ' '
            << xi_fix << ' ' << increment_s << ' ' << increment_u << ' '
            << increment_xi << '\n';
            */
        }

        // Ajust parameters
        if (xi_fix < 0) xi_fix = xi_fix + 2*pi;
        if (xi_fix > 2*pi) xi_fix = xi_fix - 2*pi;

        if (fabs(rho_fix) < 1e-20) // pure sinusoid
        {
            ((CExpParm*)parm_aux)->rho      = 0.0;
            ((CExpParm*)parm_aux)->xi       = xi_fix;
            ((CExpParm*)parm_aux)->phase    = 0.0;
            if (rho>=0.0)
            {
                ((CExpParm*)parm_aux)->a        = (int)u_aux;
                ((CExpParm*)parm_aux)->b        = m_signalSize-1;
            }
            else
            {
                ((CExpParm*)parm_aux)->a        = 0;
                ((CExpParm*)parm_aux)->b        = (int)u_aux;
            }
            setComplexAtom(parm_aux);

            // Computing optimum phase
            opt_phase = computeOptimumPhase(residue,
                                            xi);

            ((CExpParm*)parm_aux)->phase    = opt_phase;
            setRealAtom(parm_aux);

            //realAtomVector.fillVector(realAtom);

            innerProduct = residue * m_realAtom;


        }
        ((CExpParm*)parm)->innerProd = innerProduct.getData(0,0);
        ((CExpParm*)parm)->rho = ((CExpParm*)parm_aux)->rho;
        ((CExpParm*)parm)->xi = ((CExpParm*)parm_aux)->xi;
        ((CExpParm*)parm)->phase = ((CExpParm*)parm_aux)->phase;
        ((CExpParm*)parm)->a = ((CExpParm*)parm_aux)->a;
        ((CExpParm*)parm)->b = ((CExpParm*)parm_aux)->b;
    }
    delete parm_aux;
}


//===============================================================
// Function: proceedHeuristics
// Goal:
// Return:
//===============================================================

void CExpDictionary::proceedHeuristics( cgMatrix<double>& residue,
                                        CParameter* parm,
                                        CDataSignal* dataSignal)
{

    if (((CExpParm*)parm)->a != ((CExpParm*)parm)->b)   // damp and pure cases
    {

        // Find the best time support, and recalculating rho and xi
        //cout<< "->Finding the best time support, and recalculating rho and xi..." << endl;
        //findExpAtomBestTimeSupport(residue,cExpStructure);
        //findFastBestTimeSupport(residue,parm);

        // Quantize frequency
        //cout<< "->Quantizing the frequency..." << endl;
        //quantizeFrequency(residue,parm,dataSignal);

    /*  if ( (((CExpParm*)parm)->a!=0) || (((CExpParm*)parm)->b!=m_signalSize-1) )
        {
            // Find the best time support
            cout<< "->Finding the best time support, and recalculating rho ..." << endl;
            //findExpAtomBestTimeSupport(residue,cExpStructure,1);
            findFastBestTimeSupport(residue,parm,1);
            ((CExpParm*)parm)->printParm2Screen();
        }*/

        cout<< "->Finding the best time support ..." << endl;
        findFastBestTimeSupport(residue,parm,2,1);
        ((CExpParm*)parm)->printParm2Screen();

        optimizeDecaying(residue,parm);
        ((CExpParm*)parm)->printParm2Screen();

        if ( (((CExpParm*)parm)->rho != 0.0) && (dataSignal->getType()==1)) // damp case
        {
            // Discrimine Sine
            cout<< "->Discrimining sine..." << endl;
            discrimineSine(residue,parm,dataSignal);
            ((CExpParm*)parm)->printParm2Screen();
        }
    }
}

//===============================================================
// Function: findFastBestTimeSupport
// Goal: with rho and xi optimization
// Return:
//===============================================================
void CExpDictionary::findFastBestTimeSupport(   cgMatrix<double>& residue,
                                                CParameter* parm)
{
    CParameter* parm_aux;
    parm_aux = new CExpParm;

    double highestProduct = ((CExpParm*)parm)->innerProd;
    double rho = ((CExpParm*)parm)->rho;
    double xi = ((CExpParm*)parm)->xi;
    if (xi > pi) xi = (2*pi) - xi; // xi adjustment
    int a = ((CExpParm*)parm)->a;
    int b = ((CExpParm*)parm)->b;

    double opt_phase = 0;
    double phase = ((CExpParm*)parm)->phase;
    int new_a = a;
    int new_b = b;

    cgMatrix<double> innerProduct;
    cgMatrix<double> innerProductAux;

    // Decreasing Exponential
    if(rho>=0)
    {
        // From right to left
        double innerProd=0;
        double innerProdReal=0;
        double innerProdImag=0;
        double innerProdRealImag=0;
        double sqrNormReal =0;
        double sqrNormImag =0;
        double a1=0;
        double b1=0;

        ((CExpParm*)parm_aux)->rho      = rho;
        ((CExpParm*)parm_aux)->xi       = xi;
        ((CExpParm*)parm_aux)->phase    = 0.0;
        ((CExpParm*)parm_aux)->a        = a;
        ((CExpParm*)parm_aux)->b        = b;
        setComplexAtom(parm_aux);


        // Initialize inner product variables
        int i;
        for (i=a;i<=b;i++)
        {
            innerProdReal       += residue.getData(0,i) * m_complexAtom[i].Real();
            innerProdImag       += residue.getData(0,i) * m_complexAtom[i].Imag();
            innerProdRealImag   += m_complexAtom[i].Real() * m_complexAtom[i].Imag();
            sqrNormReal         += m_complexAtom[i].Real() * m_complexAtom[i].Real();
            sqrNormImag         += m_complexAtom[i].Imag() * m_complexAtom[i].Imag();
        }

        for(int b_aux=b; b_aux>a; b_aux--)
        {

            a1 = innerProdReal*sqrNormImag - innerProdImag * innerProdRealImag;
            b1 = innerProdImag*sqrNormReal - innerProdReal * innerProdRealImag;

            if ( (xi == 0) ||
                ((int)(10000*xi) == (int)(10000*pi))) //caso n�o haja sen�ide
            {
                opt_phase = 0;
            }
            else if (a1 == 0)
            {
                opt_phase = (double)(pi/2);
            }
            else if ( (a1!=0) && (xi!=0) )
            {
                opt_phase = atan( -(b1/a1) );
            }

            //cout << opt_phase << endl;

            innerProd = ( (innerProdReal * cos(opt_phase)) - (innerProdImag * sin(opt_phase)) )/
                        sqrt( (sqrNormReal*cos(opt_phase)*cos(opt_phase)) -
                        (sin(2*opt_phase)* innerProdRealImag) +
                        (sqrNormImag*sin(opt_phase)*sin(opt_phase)) );

            //update inner products
            innerProdReal       -= residue.getData(0,b_aux) * m_complexAtom[b_aux].Real();
            innerProdImag       -= residue.getData(0,b_aux) * m_complexAtom[b_aux].Imag();
            innerProdRealImag   -= m_complexAtom[b_aux].Real() * m_complexAtom[b_aux].Imag();
            sqrNormReal         -= m_complexAtom[b_aux].Real() * m_complexAtom[b_aux].Real();
            sqrNormImag         -= m_complexAtom[b_aux].Imag() * m_complexAtom[b_aux].Imag();


            if ( fabs(highestProduct) <= fabs(innerProd) )
            {
                highestProduct = innerProd;
                new_b = b_aux;
                b = new_b;
                //u = a;
                phase = opt_phase;
            }
        }

        innerProdReal = innerProdImag = innerProdRealImag = sqrNormReal = sqrNormImag =0;
        for (i=a;i<=b;i++)
        {
            innerProdReal       += residue.getData(0,i) * m_complexAtom[i].Real();
            innerProdImag       += residue.getData(0,i) * m_complexAtom[i].Imag();
            innerProdRealImag   += m_complexAtom[i].Real() * m_complexAtom[i].Imag();
            sqrNormReal         += m_complexAtom[i].Real() * m_complexAtom[i].Real();
            sqrNormImag         += m_complexAtom[i].Imag() * m_complexAtom[i].Imag();
        }

        // From left to right
        for(int a_aux=a; a_aux<b; a_aux++)
        {
            a1 = innerProdReal*sqrNormImag - innerProdImag * innerProdRealImag;
            b1 = innerProdImag*sqrNormReal - innerProdReal * innerProdRealImag;

            if ( (xi == 0) ||
                ((int)(10000*xi) == (int)(10000*pi))) //caso n�o haja sen�ide
            {
                opt_phase = 0;
            }
            else if (a1 == 0)
            {
                opt_phase = (double)(pi/2);
            }
            else if ( (a1!=0) && (xi!=0) )
            {
                opt_phase = atan( -(b1/a1) );
            }

            //cout << opt_phase << endl;

            innerProd = ( (innerProdReal * cos(opt_phase)) - (innerProdImag * sin(opt_phase)) )/
                        sqrt( (sqrNormReal*cos(opt_phase)*cos(opt_phase)) -
                        (sin(2*opt_phase)* innerProdRealImag) +
                        (sqrNormImag*sin(opt_phase)*sin(opt_phase)) );

            //update inner products
            innerProdReal       -= residue.getData(0,a_aux) * m_complexAtom[a_aux].Real();
            innerProdImag       -= residue.getData(0,a_aux) * m_complexAtom[a_aux].Imag();
            innerProdRealImag   -= m_complexAtom[a_aux].Real() * m_complexAtom[a_aux].Imag();
            sqrNormReal         -= m_complexAtom[a_aux].Real() * m_complexAtom[a_aux].Real();
            sqrNormImag         -= m_complexAtom[a_aux].Imag() * m_complexAtom[a_aux].Imag();

            if ( fabs(highestProduct) <= fabs(innerProd) )
            {
                highestProduct = innerProd;
                new_a = a_aux;
                //u = new_a;
                phase = opt_phase;
            }
        }
    }

    // Increasing exponential case
    if(rho<0)
    {
        // From left to right
        double innerProd=0;
        double innerProdReal=0;
        double innerProdImag=0;
        double innerProdRealImag=0;
        double sqrNormReal =0;
        double sqrNormImag =0;
        double a1=0;
        double b1=0;

        ((CExpParm*)parm_aux)->rho      = rho;
        ((CExpParm*)parm_aux)->xi       = xi;
        ((CExpParm*)parm_aux)->phase    = 0.0;
        ((CExpParm*)parm_aux)->a        = a;
        ((CExpParm*)parm_aux)->b        = b;
        setComplexAtom(parm_aux);

        // Initialize inner product variables
        int i;
        for (i=a;i<=b;i++)
        {
            innerProdReal       += residue.getData(0,i) * m_complexAtom[i].Real();
            innerProdImag       += residue.getData(0,i) * m_complexAtom[i].Imag();
            innerProdRealImag   += m_complexAtom[i].Real() * m_complexAtom[i].Imag();
            sqrNormReal         += m_complexAtom[i].Real() * m_complexAtom[i].Real();
            sqrNormImag         += m_complexAtom[i].Imag() * m_complexAtom[i].Imag();
        }

        for(int a_aux=a; a_aux<b; a_aux++)
        {
            a1 = innerProdReal*sqrNormImag - innerProdImag * innerProdRealImag;
            b1 = innerProdImag*sqrNormReal - innerProdReal * innerProdRealImag;

            if ( (xi == 0) ||
                ((int)(10000*xi) == (int)(10000*pi))) //caso n�o haja sen�ide
            {
                opt_phase = 0;
            }
            else if (a1 == 0)
            {
                opt_phase = (double)(pi/2);
            }
            else if ( (a1!=0) && (xi!=0) )
            {
                opt_phase = atan( -(b1/a1) );
            }

            //cout << opt_phase << endl;

            innerProd = ( (innerProdReal * cos(opt_phase)) - (innerProdImag * sin(opt_phase)) )/
                        sqrt( (sqrNormReal*cos(opt_phase)*cos(opt_phase)) -
                        (sin(2*opt_phase)* innerProdRealImag) +
                        (sqrNormImag*sin(opt_phase)*sin(opt_phase)) );

            //update inner products
            innerProdReal       -= residue.getData(0,a_aux) * m_complexAtom[a_aux].Real();
            innerProdImag       -= residue.getData(0,a_aux) * m_complexAtom[a_aux].Imag();
            innerProdRealImag   -= m_complexAtom[a_aux].Real() * m_complexAtom[a_aux].Imag();
            sqrNormReal         -= m_complexAtom[a_aux].Real() * m_complexAtom[a_aux].Real();
            sqrNormImag         -= m_complexAtom[a_aux].Imag() * m_complexAtom[a_aux].Imag();

            if ( fabs(highestProduct) <= fabs(innerProd) )
            {
                highestProduct = innerProd;
                new_a = a_aux;
                a = a_aux;
                //u = (double)(expDic.getSignalSize()-1-b);
                phase = opt_phase;
            }
        }

        innerProdReal = innerProdImag = innerProdRealImag = sqrNormReal = sqrNormImag =0;
        for (i=a;i<=b;i++)
        {
            innerProdReal       += residue.getData(0,i) * m_complexAtom[i].Real();
            innerProdImag       += residue.getData(0,i) * m_complexAtom[i].Imag();
            innerProdRealImag   += m_complexAtom[i].Real() * m_complexAtom[i].Imag();
            sqrNormReal         += m_complexAtom[i].Real() * m_complexAtom[i].Real();
            sqrNormImag         += m_complexAtom[i].Imag() * m_complexAtom[i].Imag();
        }

        // From right to left
        for(int b_aux=b; b_aux>a; b_aux--)
        {
            a1 = innerProdReal*sqrNormImag - innerProdImag * innerProdRealImag;
            b1 = innerProdImag*sqrNormReal - innerProdReal * innerProdRealImag;

            if ( (xi == 0) ||
                ((int)(10000*xi) == (int)(10000*pi))) //caso n�o haja sen�ide
            {
                opt_phase = 0;
            }
            else if (a1 == 0)
            {
                opt_phase = (double)(pi/2);
            }
            else if ( (a1!=0) && (xi!=0) )
            {
                opt_phase = atan( -(b1/a1) );
            }

            //cout << opt_phase << endl;

            innerProd = ( (innerProdReal * cos(opt_phase)) - (innerProdImag * sin(opt_phase)) )/
                        sqrt( (sqrNormReal*cos(opt_phase)*cos(opt_phase)) -
                        (sin(2*opt_phase)* innerProdRealImag) +
                        (sqrNormImag*sin(opt_phase)*sin(opt_phase)) );

            //update inner products
            innerProdReal       -= residue.getData(0,b_aux) * m_complexAtom[b_aux].Real();
            innerProdImag       -= residue.getData(0,b_aux) * m_complexAtom[b_aux].Imag();
            innerProdRealImag   -= m_complexAtom[b_aux].Real() * m_complexAtom[b_aux].Imag();
            sqrNormReal         -= m_complexAtom[b_aux].Real() * m_complexAtom[b_aux].Real();
            sqrNormImag         -= m_complexAtom[b_aux].Imag() * m_complexAtom[b_aux].Imag();


            if ( fabs(highestProduct) <= fabs(innerProd) )
            {
                highestProduct = innerProd;
                new_b = b_aux;
                //u = (double)(expDic.getSignalSize()-1-b_aux);
                phase = opt_phase;
            }
        }
    }


    // Modified down to this point: u(maybe), new_a, new_b.

    // Recalculate rho and xi for the new time support
    // by a pseudo-newton method

    double rho_fix, rho_aux, increment_rho;
    rho_fix = rho;
    increment_rho = rho * 0.5;//-rho / 10;
    rho_aux = rho_fix + increment_rho;

    double xi_fix, xi_aux, increment_xi;
    xi_fix = xi;
    increment_xi = xi * 0.5;
    xi_aux = xi_fix + increment_xi;

    // =======================
    // Pseudo-Newton method

    cout << "Pseudo-Newton" << endl;
    double ph_fix = phase;
    int time =0;
    int incr_count = 0;

    double xi_lim = 1.0 / (double) ( 20 * ( m_signalSize ) * ( m_signalSize ) );
    double rho_lim  = 1.0 / (double) ( 100 * ( m_signalSize ) * ( m_signalSize ) * ( m_signalSize ) );

    while(  (fabs(increment_xi) > xi_lim) ||
            (fabs(increment_rho) > rho_lim) )
    {
        if(time==0) // xi_aux
        {
            ((CExpParm*)parm_aux)->rho      = rho_fix;
            ((CExpParm*)parm_aux)->xi       = xi_aux;
            ((CExpParm*)parm_aux)->phase    = 0.0;
            ((CExpParm*)parm_aux)->a        = new_a;
            ((CExpParm*)parm_aux)->b        = new_b;

            setComplexAtom(parm_aux);

            // Computing optimum phase
            opt_phase = computeOptimumPhase(residue,
                                            xi);

            ((CExpParm*)parm_aux)->phase    = opt_phase;
            setRealAtom(parm_aux);

        }
        if(time==1) // rho_aux
        {
            ((CExpParm*)parm_aux)->rho      = rho_aux;
            ((CExpParm*)parm_aux)->xi       = xi_fix;
            ((CExpParm*)parm_aux)->phase    = 0.0;
            ((CExpParm*)parm_aux)->a        = new_a;
            ((CExpParm*)parm_aux)->b        = new_b;

            setComplexAtom(parm_aux);

            // Computing optimum phase
            opt_phase = computeOptimumPhase(residue,
                                            xi);

            ((CExpParm*)parm_aux)->phase    = opt_phase;
            setRealAtom(parm_aux);
        }

        //realAtomVector.fillVector(realAtom);

        innerProduct = residue * m_realAtom;//( realAtomVector.transpose() );

        //cout << rho_fix << " " << xi_fix << endl;

        if( fabs(highestProduct) < fabs(innerProduct.getData(0,0)) )
        {
            highestProduct = innerProduct.getData(0,0);

            if(time==0)
            {
                /*if(fabs(xi_aux)<0.001)
                {
                    xi_aux=0;
                    xi_fix=0;
                    ph_fix=0;
                }
                else
                {*/
                xi_fix = xi_aux;
                xi_aux = xi_fix + increment_xi;
                ph_fix = opt_phase;
                //}

                increment_xi = 2.0 * increment_xi ;

            }
            if(time==1)
            {
                rho_fix = rho_aux;
                rho_aux = rho_fix + increment_rho;
                ph_fix = opt_phase;

                increment_rho = pow(4.0,(double)incr_count) * increment_rho ;
            }
            incr_count++;
        }
        else
        {
            if(time==0)
            {
                /*if(fabs(xi_aux)<0.001)
                {
                    xi_fix = xi_aux; //= 0;
                    increment_xi = - increment_xi / 2.0;
                    xi_aux = xi_fix + increment_xi;
                }
                else
                {*/
                increment_xi = - increment_xi / 2.0;
                xi_aux = xi_fix + increment_xi;
                //}
            }
            if(time==1)
            {
                increment_rho =  increment_rho / pow(4.0,(double)incr_count);
                increment_rho = - increment_rho / 2.0;
                rho_aux = rho_fix + increment_rho;
            }

            time++;
            time = time % 2;

            incr_count = 0;
        }
    }

    if (xi_fix > pi) xi_fix = (2*pi) - xi_fix; // xi adjustment

    ((CExpParm*)parm_aux)->rho      = rho_fix;
    ((CExpParm*)parm_aux)->xi       = xi_fix;
    ((CExpParm*)parm_aux)->phase    = 0.0;
    ((CExpParm*)parm_aux)->a        = new_a;
    ((CExpParm*)parm_aux)->b        = new_b;

    setComplexAtom(parm_aux);

    // Computing optimum phase
    opt_phase = computeOptimumPhase(residue,
                                    xi);

    ((CExpParm*)parm_aux)->phase    = opt_phase;
    setRealAtom(parm_aux);

    innerProduct = residue * m_realAtom;

    ((CExpParm*)parm)->innerProd =  innerProduct.getData(0,0);
    ((CExpParm*)parm)->rho   =  ((CExpParm*)parm_aux)->rho;
    ((CExpParm*)parm)->xi    =  ((CExpParm*)parm_aux)->xi;
    ((CExpParm*)parm)->phase =  ((CExpParm*)parm_aux)->phase;
    ((CExpParm*)parm)->a     =  ((CExpParm*)parm_aux)->a;
    ((CExpParm*)parm)->b     =  ((CExpParm*)parm_aux)->b;


    delete parm_aux;

}

//===============================================================
// Function: findFastBestTimeSupport
// Goal: with rho optimization
// Return:
// Obs: If dummy = 0 -> One-way search - decreasing exp: rigth to left; increasing exp: left to right
//      If dummy = 1 -> Two-way search
//      If dummy = 2 -> runs only from left to the right
//===============================================================

void CExpDictionary::findFastBestTimeSupport(   cgMatrix<double>& residue,
                                                CParameter* parm,
                                                int dummy,
                                                double coefHeur)
{
    CParameter* parm_aux;
    parm_aux = new CExpParm;

    double highestProduct = ((CExpParm*)parm)->innerProd;
    double rho = ((CExpParm*)parm)->rho;
    double xi = ((CExpParm*)parm)->xi;
    if (xi > pi) xi = (2*pi) - xi; // xi adjustment
    int a = ((CExpParm*)parm)->a;
    int b = ((CExpParm*)parm)->b;

    double phase = ((CExpParm*)parm)->phase;
    int new_a = a;
    int new_b = b;

    cgMatrix<double> innerProduct;
//  double innerProd;
    double innerProdAux;

    ((CExpParm*)parm_aux)->rho      = rho;
    ((CExpParm*)parm_aux)->xi       = xi;
    ((CExpParm*)parm_aux)->phase    = phase;
    ((CExpParm*)parm_aux)->a        = a;
    ((CExpParm*)parm_aux)->b        = b;
    setRealAtom(parm_aux);

    // Decreasing Exponential
    if(rho>=0)
    {

        if ( (dummy==0) || (dummy==1) || (dummy==2) )
        {
          innerProdAux = highestProduct;
          // From right to left
          for(int b_aux=b; b_aux>a; b_aux--)
          {

             innerProdAux = (innerProdAux - (residue[0][b_aux]* m_realAtom[b_aux]) ) /
                             sqrt(1-(m_realAtom[b_aux]*m_realAtom[b_aux]));


             /*((CExpParm*)parm_aux)->b = b_aux-1;

             //setComplexAtom(parm_aux);

             // Computing optimum phase
             //opt_phase = computeOptimumPhase(residue,
             //                             ((CExpParm*)parm_aux)->xi,
             //                             innerProd);

             //((CExpParm*)parm_aux)->phase = opt_phase;
             setRealAtom(parm_aux);

             innerProduct = residue * m_realAtom;
             */



             //if ( fabs(highestProduct) <= fabs(innerProd) )
             //if ( fabs(highestProduct) <= fabs(innerProduct[0][0]) )
             if ( fabs(highestProduct)*coefHeur <= fabs(innerProdAux) )
             {
                    //highestProduct = innerProd;
                    //highestProduct = innerProduct[0][0];
                    highestProduct = innerProdAux;
                    new_b = b_aux-1;
                    b = new_b;
                    //u = a;
             }
          }
        }
        ((CExpParm*)parm_aux)->b = b;
        setRealAtom(parm_aux);

        if (dummy==1)
        {
            innerProdAux = highestProduct;
            // From left to right
            for(int a_aux=a; a_aux<b; a_aux++)
            {

                innerProdAux = (innerProdAux - (residue[0][a_aux]*m_realAtom[a_aux]) )/
                                sqrt(1 - (m_realAtom[a_aux]*m_realAtom[a_aux]));

                /*((CExpParm*)parm_aux)->a = a_aux+1;

                //setComplexAtom(parm_aux);

                // Computing optimum phase
                //opt_phase = computeOptimumPhase(residue,
                //                              ((CExpParm*)parm_aux)->xi,
                //                              innerProd);

                //((CExpParm*)parm_aux)->phase  = opt_phase;
                setRealAtom(parm_aux);

                innerProduct = residue * m_realAtom;
                */

                if ( fabs(highestProduct)*coefHeur <= fabs(innerProdAux) )
                {
                    //highestProduct = innerProd;
                    highestProduct = innerProdAux;
                    new_a = a_aux+1;
                    //u = new_a;
                }
            }
        }
    }


    // Increasing exponential case
    if(rho<0)
    {
        if ( (dummy==0) || (dummy==1) )
        {
            innerProdAux = highestProduct;
            // From left to right
            for(int a_aux=a; a_aux<b; a_aux++)
            {
                innerProdAux = (innerProdAux - (residue[0][a_aux]*m_realAtom[a_aux]) )/
                                sqrt(1 - (m_realAtom[a_aux]*m_realAtom[a_aux]));

        /*      ((CExpParm*)parm_aux)->a = a_aux;

                setComplexAtom(parm_aux);

                // Computing optimum phase
                opt_phase = computeOptimumPhase(residue,
                                                ((CExpParm*)parm_aux)->xi,
                                                innerProd);

                ((CExpParm*)parm_aux)->phase    = opt_phase;
                setRealAtom(parm_aux);

                //innerProduct = residue * m_realAtom;
        */
                if ( fabs(highestProduct)*coefHeur <= fabs(innerProdAux) )
                {
                    highestProduct = innerProdAux;
                    new_a = a_aux+1;
                    a = new_a;
                    //u = new_a;
                }
            }
            ((CExpParm*)parm_aux)->a = a;
            setRealAtom(parm_aux);
        }
        if ( (dummy==1) || (dummy==2))
        {
          innerProdAux = highestProduct;
          // From right to left
          for(int b_aux=b; b_aux>a; b_aux--)
          {

             innerProdAux = (innerProdAux - (residue[0][b_aux]* m_realAtom[b_aux]) ) /
                             sqrt(1-(m_realAtom[b_aux]*m_realAtom[b_aux]));

          /*    ((CExpParm*)parm_aux)->b = b_aux;

             setComplexAtom(parm_aux);

             // Computing optimum phase
             opt_phase = computeOptimumPhase(residue,
                                             ((CExpParm*)parm_aux)->xi,
                                             innerProd);

             ((CExpParm*)parm_aux)->phase   = opt_phase;
             setRealAtom(parm_aux);

             //innerProduct = residue * m_realAtom;
    */
             if ( fabs(highestProduct)*coefHeur <= fabs(innerProdAux) )
             {
                    highestProduct = innerProdAux;
                    new_b = b_aux-1;
                    //u = a;
             }
          }
        }

    }



    ((CExpParm*)parm_aux)->innerProd = highestProduct;
    //((CExpParm*)parm_aux)->phase     = phase;
    ((CExpParm*)parm_aux)->a         = new_a;
    ((CExpParm*)parm_aux)->b         = new_b;


    ((CExpParm*)parm)->innerProd = ((CExpParm*)parm_aux)->innerProd;
    ((CExpParm*)parm)->rho   =  ((CExpParm*)parm_aux)->rho;
    ((CExpParm*)parm)->xi    =  ((CExpParm*)parm_aux)->xi;
    ((CExpParm*)parm)->phase =  ((CExpParm*)parm_aux)->phase;
    ((CExpParm*)parm)->a     =  ((CExpParm*)parm_aux)->a;
    ((CExpParm*)parm)->b     =  ((CExpParm*)parm_aux)->b;


    delete parm_aux;

}

/* void CExpDictionary::findFastBestTimeSupport(    cgMatrix<double>& residue,
                                                CParameter* parm,
                                                int dummy)
{
    CParameter* parm_aux;
    parm_aux = new CExpParm;

    double highestProduct = ((CExpParm*)parm)->innerProd;
    double rho = ((CExpParm*)parm)->rho;
    double xi = ((CExpParm*)parm)->xi;
    if (xi > pi) xi = (2*pi) - xi; // xi adjustment
    int a = ((CExpParm*)parm)->a;
    int b = ((CExpParm*)parm)->b;

    double opt_phase = 0;
    double phase = ((CExpParm*)parm)->phase;
    int new_a = a;
    int new_b = b;

    cgMatrix<double> innerProduct;
    cgMatrix<double> innerProductAux;

    // Decreasing Exponential
    if(rho>=0)
    {
        // From right to left
        double innerProd=0;
        double innerProdReal=0;
        double innerProdImag=0;
        double innerProdRealImag=0;
        double sqrNormReal =0;
        double sqrNormImag =0;
        double a1=0;
        double b1=0;

        ((CExpParm*)parm_aux)->rho      = rho;
        ((CExpParm*)parm_aux)->xi       = xi;
        ((CExpParm*)parm_aux)->phase    = 0.0;
        ((CExpParm*)parm_aux)->a        = a;
        ((CExpParm*)parm_aux)->b        = b;
        setComplexAtom(parm_aux);


        // Initialize inner product variables
        int i;

        for (i=a;i<=b;i++)
        {
            innerProdReal       += residue.getData(0,i) * m_complexAtom[i].Real();
            innerProdImag       += residue.getData(0,i) * m_complexAtom[i].Imag();
            innerProdRealImag   += m_complexAtom[i].Real() * m_complexAtom[i].Imag();
            sqrNormReal         += m_complexAtom[i].Real() * m_complexAtom[i].Real();
            sqrNormImag         += m_complexAtom[i].Imag() * m_complexAtom[i].Imag();
        }

        for(int b_aux=b; b_aux>a; b_aux--)
        {

            a1 = innerProdReal*sqrNormImag - innerProdImag * innerProdRealImag;
            b1 = innerProdImag*sqrNormReal - innerProdReal * innerProdRealImag;

            if ( (xi == 0) ||
                ((int)(10000*xi) == (int)(10000*pi))) //caso n�o haja sen�ide
            {
                opt_phase = 0;
            }
            else if (a1 == 0)
            {
                opt_phase = (double)(pi/2);
            }
            else if ( (a1!=0) && (xi!=0) )
            {
                opt_phase = atan( -(b1/a1) );
            }

            //cout << opt_phase << endl;

            innerProd = ( (innerProdReal * cos(opt_phase)) - (innerProdImag * sin(opt_phase)) )/
                        sqrt( (sqrNormReal*cos(opt_phase)*cos(opt_phase)) -
                        (sin(2*opt_phase)* innerProdRealImag) +
                        (sqrNormImag*sin(opt_phase)*sin(opt_phase)) );

            //update inner products
            innerProdReal       -= residue.getData(0,b_aux) * m_complexAtom[b_aux].Real();
            innerProdImag       -= residue.getData(0,b_aux) * m_complexAtom[b_aux].Imag();
            innerProdRealImag   -= m_complexAtom[b_aux].Real() * m_complexAtom[b_aux].Imag();
            sqrNormReal         -= m_complexAtom[b_aux].Real() * m_complexAtom[b_aux].Real();
            sqrNormImag         -= m_complexAtom[b_aux].Imag() * m_complexAtom[b_aux].Imag();


            if ( fabs(highestProduct) <= fabs(innerProd) )
            {
                highestProduct = innerProd;
                new_b = b_aux;
                b = new_b;
                //u = a;
                phase = opt_phase;
            }
        }

        innerProdReal = innerProdImag = innerProdRealImag = sqrNormReal = sqrNormImag =0;
        for (i=a;i<=b;i++)
        {
            innerProdReal       += residue.getData(0,i) * m_complexAtom[i].Real();
            innerProdImag       += residue.getData(0,i) * m_complexAtom[i].Imag();
            innerProdRealImag   += m_complexAtom[i].Real() * m_complexAtom[i].Imag();
            sqrNormReal         += m_complexAtom[i].Real() * m_complexAtom[i].Real();
            sqrNormImag         += m_complexAtom[i].Imag() * m_complexAtom[i].Imag();
        }

        if (dummy!=2)
        {
            // From left to right
            for(int a_aux=a; a_aux<b; a_aux++)
            {
                a1 = innerProdReal*sqrNormImag - innerProdImag * innerProdRealImag;
                b1 = innerProdImag*sqrNormReal - innerProdReal * innerProdRealImag;

                if ( (xi == 0) ||
                    ((int)(10000*xi) == (int)(10000*pi))) //caso n�o haja sen�ide
                {
                    opt_phase = 0;
                }
                else if (a1 == 0)
                {
                    opt_phase = (double)(pi/2);
                }
                else if ( (a1!=0) && (xi!=0) )
                {
                    opt_phase = atan( -(b1/a1) );
                }

                //cout << opt_phase << endl;

                innerProd = ( (innerProdReal * cos(opt_phase)) - (innerProdImag * sin(opt_phase)) )/
                            sqrt( (sqrNormReal*cos(opt_phase)*cos(opt_phase)) -
                            (sin(2*opt_phase)* innerProdRealImag) +
                            (sqrNormImag*sin(opt_phase)*sin(opt_phase)) );

                //update inner products
                innerProdReal       -= residue.getData(0,a_aux) * m_complexAtom[a_aux].Real();
                innerProdImag       -= residue.getData(0,a_aux) * m_complexAtom[a_aux].Imag();
                innerProdRealImag   -= m_complexAtom[a_aux].Real() * m_complexAtom[a_aux].Imag();
                sqrNormReal         -= m_complexAtom[a_aux].Real() * m_complexAtom[a_aux].Real();
                sqrNormImag         -= m_complexAtom[a_aux].Imag() * m_complexAtom[a_aux].Imag();

                if ( fabs(highestProduct) <= fabs(innerProd) )
                {
                    highestProduct = innerProd;
                    new_a = a_aux;
                    //u = new_a;
                    phase = opt_phase;
                }
            }
        }
    }

    // Increasing exponential case
    if(rho<0)
    {

        // From left to right
        double innerProd=0;
        double innerProdReal=0;
        double innerProdImag=0;
        double innerProdRealImag=0;
        double sqrNormReal =0;
        double sqrNormImag =0;
        double a1=0;
        double b1=0;

        ((CExpParm*)parm_aux)->rho      = rho;
        ((CExpParm*)parm_aux)->xi       = xi;
        ((CExpParm*)parm_aux)->phase    = 0.0;
        ((CExpParm*)parm_aux)->a        = a;
        ((CExpParm*)parm_aux)->b        = b;
        setComplexAtom(parm_aux);

        // Initialize inner product variables
        int i;
        if (dummy!=2)
        {
            for (i=a;i<=b;i++)
            {
                innerProdReal       += residue.getData(0,i) * m_complexAtom[i].Real();
                innerProdImag       += residue.getData(0,i) * m_complexAtom[i].Imag();
                innerProdRealImag   += m_complexAtom[i].Real() * m_complexAtom[i].Imag();
                sqrNormReal         += m_complexAtom[i].Real() * m_complexAtom[i].Real();
                sqrNormImag         += m_complexAtom[i].Imag() * m_complexAtom[i].Imag();
            }

            for(int a_aux=a; a_aux<b; a_aux++)
            {
                a1 = innerProdReal*sqrNormImag - innerProdImag * innerProdRealImag;
                b1 = innerProdImag*sqrNormReal - innerProdReal * innerProdRealImag;

                if ( (xi == 0) ||
                    ((int)(10000*xi) == (int)(10000*pi))) //caso n�o haja sen�ide
                {
                    opt_phase = 0;
                }
                else if (a1 == 0)
                {
                    opt_phase = (double)(pi/2);
                }
                else if ( (a1!=0) && (xi!=0) )
                {
                    opt_phase = atan( -(b1/a1) );
                }

                //cout << opt_phase << endl;

                innerProd = ( (innerProdReal * cos(opt_phase)) - (innerProdImag * sin(opt_phase)) )/
                            sqrt( (sqrNormReal*cos(opt_phase)*cos(opt_phase)) -
                            (sin(2*opt_phase)* innerProdRealImag) +
                            (sqrNormImag*sin(opt_phase)*sin(opt_phase)) );

                //update inner products
                innerProdReal       -= residue.getData(0,a_aux) * m_complexAtom[a_aux].Real();
                innerProdImag       -= residue.getData(0,a_aux) * m_complexAtom[a_aux].Imag();
                innerProdRealImag   -= m_complexAtom[a_aux].Real() * m_complexAtom[a_aux].Imag();
                sqrNormReal         -= m_complexAtom[a_aux].Real() * m_complexAtom[a_aux].Real();
                sqrNormImag         -= m_complexAtom[a_aux].Imag() * m_complexAtom[a_aux].Imag();

                if ( fabs(highestProduct) <= fabs(innerProd) )
                {
                    highestProduct = innerProd;
                    new_a = a_aux;
                    a = a_aux;
                    //u = (double)(expDic.getSignalSize()-1-b);
                    phase = opt_phase;
                }
            }
        }

        innerProdReal = innerProdImag = innerProdRealImag = sqrNormReal = sqrNormImag =0;
        for (i=a;i<=b;i++)
        {
            innerProdReal       += residue.getData(0,i) * m_complexAtom[i].Real();
            innerProdImag       += residue.getData(0,i) * m_complexAtom[i].Imag();
            innerProdRealImag   += m_complexAtom[i].Real() * m_complexAtom[i].Imag();
            sqrNormReal         += m_complexAtom[i].Real() * m_complexAtom[i].Real();
            sqrNormImag         += m_complexAtom[i].Imag() * m_complexAtom[i].Imag();
        }

        // From right to left
        for(int b_aux=b; b_aux>a; b_aux--)
        {
            a1 = innerProdReal*sqrNormImag - innerProdImag * innerProdRealImag;
            b1 = innerProdImag*sqrNormReal - innerProdReal * innerProdRealImag;

            if ( (xi == 0) ||
                ((int)(10000*xi) == (int)(10000*pi))) //caso n�o haja sen�ide
            {
                opt_phase = 0;
            }
            else if (a1 == 0)
            {
                opt_phase = (double)(pi/2);
            }
            else if ( (a1!=0) && (xi!=0) )
            {
                opt_phase = atan( -(b1/a1) );
            }

            //cout << opt_phase << endl;

            innerProd = ( (innerProdReal * cos(opt_phase)) - (innerProdImag * sin(opt_phase)) )/
                        sqrt( (sqrNormReal*cos(opt_phase)*cos(opt_phase)) -
                        (sin(2*opt_phase)* innerProdRealImag) +
                        (sqrNormImag*sin(opt_phase)*sin(opt_phase)) );

            //update inner products
            innerProdReal       -= residue.getData(0,b_aux) * m_complexAtom[b_aux].Real();
            innerProdImag       -= residue.getData(0,b_aux) * m_complexAtom[b_aux].Imag();
            innerProdRealImag   -= m_complexAtom[b_aux].Real() * m_complexAtom[b_aux].Imag();
            sqrNormReal         -= m_complexAtom[b_aux].Real() * m_complexAtom[b_aux].Real();
            sqrNormImag         -= m_complexAtom[b_aux].Imag() * m_complexAtom[b_aux].Imag();


            if ( fabs(highestProduct) <= fabs(innerProd) )
            {
                highestProduct = innerProd;
                new_b = b_aux;
                //u = (double)(expDic.getSignalSize()-1-b_aux);
                phase = opt_phase;
            }
        }
    }


    // Modified down to this point: u(maybe), new_a, new_b.

    ((CExpParm*)parm_aux)->innerProd = highestProduct;
    ((CExpParm*)parm_aux)->phase     = phase;
    ((CExpParm*)parm_aux)->a         = new_a;
    ((CExpParm*)parm_aux)->b         = new_b;


    ((CExpParm*)parm)->innerProd = ((CExpParm*)parm_aux)->innerProd;
    ((CExpParm*)parm)->rho   =  ((CExpParm*)parm_aux)->rho;
    ((CExpParm*)parm)->xi    =  ((CExpParm*)parm_aux)->xi;
    ((CExpParm*)parm)->phase =  ((CExpParm*)parm_aux)->phase;
    ((CExpParm*)parm)->a     =  ((CExpParm*)parm_aux)->a;
    ((CExpParm*)parm)->b     =  ((CExpParm*)parm_aux)->b;


    delete parm_aux;

}
*/
//===============================================================
// Function: findBestTimeSupport
// Goal: with rho optimization
// Return:
// Obs: If dummy = 2 -> runs only from left to the right
//===============================================================
void CExpDictionary::findBestTimeSupport(   cgMatrix<double>& residue,
                                            CParameter* parm,
                                            int dummy)
{
    CParameter* parm_aux;
    parm_aux = new CExpParm;

    double highestProduct = ((CExpParm*)parm)->innerProd;
    double rho = ((CExpParm*)parm)->rho;
    double xi = ((CExpParm*)parm)->xi;
    if (xi > pi) xi = (2*pi) - xi; // xi adjustment
    int a = ((CExpParm*)parm)->a;
    int b = ((CExpParm*)parm)->b;

    double opt_phase = 0;
    double phase = ((CExpParm*)parm)->phase;
    int new_a = a;
    int new_b = b;

    cgMatrix<double> innerProduct;
    double innerProd;

    ((CExpParm*)parm_aux)->rho      = rho;
    ((CExpParm*)parm_aux)->xi       = xi;
    ((CExpParm*)parm_aux)->phase    = 0.0;
    ((CExpParm*)parm_aux)->a        = a;
    ((CExpParm*)parm_aux)->b        = b;

    // Decreasing Exponential
    if(rho>=0)
    {
        // From right to left
        for(int b_aux=b; b_aux>a; b_aux--)
        {

            ((CExpParm*)parm_aux)->b = b_aux;

            setComplexAtom(parm_aux);

            // Computing optimum phase
            opt_phase = computeOptimumPhase(residue,
                                            ((CExpParm*)parm_aux)->xi,
                                            innerProd);

            ((CExpParm*)parm_aux)->phase    = opt_phase;
            setRealAtom(parm_aux);

            //innerProduct = residue * m_realAtom;

            if ( fabs(highestProduct) <= fabs(innerProd) )
            {
                highestProduct = innerProd;
                new_b = b_aux;
                b = new_b;
                //u = a;
                phase = opt_phase;
            }
        }

        if (dummy!=2)
        {
            // From left to right
            for(int a_aux=a; a_aux<b; a_aux++)
            {

                ((CExpParm*)parm_aux)->a = a_aux;

                setComplexAtom(parm_aux);

                // Computing optimum phase
                opt_phase = computeOptimumPhase(residue,
                                                ((CExpParm*)parm_aux)->xi,
                                                innerProd);

                ((CExpParm*)parm_aux)->phase    = opt_phase;
                setRealAtom(parm_aux);

                //innerProduct = residue * m_realAtom;

                if ( fabs(highestProduct) <= fabs(innerProd) )
                {
                    highestProduct = innerProd;
                    new_a = a_aux;
                    //u = new_a;
                    phase = opt_phase;
                }
            }
        }
    }

    // Increasing exponential case
    if(rho<0)
    {
        if (dummy!=2)
        {
            // From left to right
            for(int a_aux=a; a_aux<b; a_aux++)
            {

                ((CExpParm*)parm_aux)->a = a_aux;

                setComplexAtom(parm_aux);

                // Computing optimum phase
                opt_phase = computeOptimumPhase(residue,
                                                ((CExpParm*)parm_aux)->xi,
                                                innerProd);

                ((CExpParm*)parm_aux)->phase    = opt_phase;
                setRealAtom(parm_aux);

                //innerProduct = residue * m_realAtom;

                if ( fabs(highestProduct) <= fabs(innerProd) )
                {
                    highestProduct = innerProd;
                    new_a = a_aux;
                    //u = new_a;
                    phase = opt_phase;
                }
            }
        }

        // From right to left
        for(int b_aux=b; b_aux>a; b_aux--)
        {

            ((CExpParm*)parm_aux)->b = b_aux;

            setComplexAtom(parm_aux);

            // Computing optimum phase
            opt_phase = computeOptimumPhase(residue,
                                            ((CExpParm*)parm_aux)->xi,
                                            innerProd);

            ((CExpParm*)parm_aux)->phase    = opt_phase;
            setRealAtom(parm_aux);

            //innerProduct = residue * m_realAtom;

            if ( fabs(highestProduct) <= fabs(innerProd) )
            {
                highestProduct = innerProd;
                new_b = b_aux;
                b = new_b;
                //u = a;
                phase = opt_phase;
            }
        }

    }


    // Modified down to this point: u(maybe), new_a, new_b.

    ((CExpParm*)parm_aux)->innerProd = highestProduct;
    ((CExpParm*)parm_aux)->phase     = phase;
    ((CExpParm*)parm_aux)->a         = new_a;
    ((CExpParm*)parm_aux)->b         = new_b;


    ((CExpParm*)parm)->innerProd = ((CExpParm*)parm_aux)->innerProd;
    ((CExpParm*)parm)->rho   =  ((CExpParm*)parm_aux)->rho;
    ((CExpParm*)parm)->xi    =  ((CExpParm*)parm_aux)->xi;
    ((CExpParm*)parm)->phase =  ((CExpParm*)parm_aux)->phase;
    ((CExpParm*)parm)->a     =  ((CExpParm*)parm_aux)->a;
    ((CExpParm*)parm)->b     =  ((CExpParm*)parm_aux)->b;


    delete parm_aux;

}


void CExpDictionary::optimizeDecaying(  cgMatrix<double>& residue,
                                        CParameter* parm_aux)
{

    // Recalculate rho for the new time support
    // by a pseudo-newton method
    cgMatrix<double> innerProduct;
    double highestProduct = ((CExpParm*)parm_aux)->innerProd;

    double rho_fix, rho_aux, increment_rho;
    if (((CExpParm*)parm_aux)->rho == 0)
    {
        rho_fix = 1e-7;
    }
    else
    {
        rho_fix = ((CExpParm*)parm_aux)->rho;
    }
    increment_rho = rho_fix * 0.5;//-rho / 10;
    rho_aux = rho_fix + increment_rho;

    // =======================
    // Pseudo-Newton method

    cout << "Optimizing the Decaying with Pseudo-Newton algorithm" << endl;
    double ph_fix = ((CExpParm*)parm_aux)->phase;
    int time =0;
    int incr_count = 0;

    double rho_lim  = 1e-8;
    double opt_phase;
    while (fabs(increment_rho) > rho_lim)
    {

        ((CExpParm*)parm_aux)->rho      = rho_aux;

        setComplexAtom(parm_aux);

        // Computing optimum phase
        opt_phase = computeOptimumPhase(residue,
                                        ((CExpParm*)parm_aux)->xi);

        ((CExpParm*)parm_aux)->phase    = opt_phase;
        setRealAtom(parm_aux);


        //realAtomVector.fillVector(realAtom);

        innerProduct = residue * m_realAtom;//( realAtomVector.transpose() );

        //cout << rho_fix << " " << xi_fix << endl;

        if( fabs(highestProduct) < fabs(innerProduct.getData(0,0)) )
        {
            highestProduct = innerProduct.getData(0,0);

            rho_fix = rho_aux;
            rho_aux = rho_fix + increment_rho;
            ph_fix = opt_phase;

            increment_rho = pow(4.0,(double)incr_count) * increment_rho;

            incr_count++;
        }
        else
        {
            increment_rho =  increment_rho / pow(4.0,(double)incr_count);
            increment_rho = - increment_rho / 2.0;
            rho_aux = rho_fix + increment_rho;
            incr_count = 0;
        }
    }

    if (fabs(rho_fix)<1e-7)
    {
        ((CExpParm*)parm_aux)->rho      = 0.0;
        setComplexAtom(parm_aux);

        // Computing optimum phase
        opt_phase = computeOptimumPhase(residue,
                                        ((CExpParm*)parm_aux)->xi);

        ((CExpParm*)parm_aux)->phase    = opt_phase;
        setRealAtom(parm_aux);

        innerProduct = residue * m_realAtom;
        ((CExpParm*)parm_aux)->innerProd = innerProduct[0][0];
    }
    else
    {
        ((CExpParm*)parm_aux)->innerProd = highestProduct;
        ((CExpParm*)parm_aux)->rho = rho_fix;
        ((CExpParm*)parm_aux)->phase = ph_fix;
    }

}

void CExpDictionary::optimizeDecayingErrorNorm(  cgMatrix<double>& residue,
                                                 CParameter* parm_aux,
                                                 double& minResNorm,
                                                 double coef_xi0,
                                                 double ratomsample_xi0)
{

    // =======================
    // Pseudo-Newton method

    cout << " -- Optimizing the Decaying with Pseudo-Newton algorithm regarding error norm" << endl;
    cout << "    INIT: ";
    ((CExpParm*)parm_aux)->printParm2Screen();
    // Recalculate rho for the new time support
    // by a pseudo-newton method
    CParameter* parm_aux2;
    parm_aux2 = new CExpParm;
    CParameter* parm_aux3;
    parm_aux3 = new CExpParm;
    double res_norm;
    cgMatrix<double> cgRealAtom(1,m_signalSize,0.0);
    cgMatrix<double> cgRealAtomAux(1,m_signalSize,0.0);
    cgMatrix<double> res(1,m_signalSize,0.0);

    double rho_fix, rho_aux, increment_rho;
    if (((CExpParm*)parm_aux)->rho == 0)
    {
        rho_fix = 1e-7;
    }
    else
    {
        rho_fix = ((CExpParm*)parm_aux)->rho;
    }
    increment_rho = rho_fix * 0.5;//-rho / 10;
    rho_aux = rho_fix + increment_rho;


 //   double innerProd = ((CExpParm*)parm_aux)->innerProd;
//    cout << "Olhai Primeiro!!!" << endl;
//    ((CExpParm*)parm_aux)->printParm2Screen();

    int time = 0;
    int incr_count = 0;



    double norm2, norm3;


    double rho_lim  = 1e-8; // 1/(2^16) 16 bit representation
    int nmaxiter = 5000;
    int imaxiter = 0;
    cout << "rho_fix" << " " << "rho_aux" << " " << "res_norm" << " " << "minResNorm" << " " << "coef_xi0" << " " << "ratomsample_xi0" << endl;
    while ( (fabs(increment_rho) > rho_lim) )
    {
        //cout << fabs(increment_rho) << " " << rho_lim <<endl;
        ((CExpParm*)parm_aux2)->rho = rho_aux;
        ((CExpParm*)parm_aux2)->xi = 0.0;
        ((CExpParm*)parm_aux2)->phase = 0.0;
        ((CExpParm*)parm_aux2)->a = 0;
        ((CExpParm*)parm_aux2)->b = m_signalSize-1;
        setRealAtom(parm_aux2);
        norm2 = m_rAtomNorm;
        ((CExpParm*)parm_aux2)->innerProd = (coef_xi0*ratomsample_xi0)/ m_realAtom[0];

        ((CExpParm*)parm_aux3)->rho = ((CExpParm*)parm_aux2)->rho;
        ((CExpParm*)parm_aux3)->xi = ((CExpParm*)parm_aux)->xi;
        ((CExpParm*)parm_aux3)->phase = ((CExpParm*)parm_aux)->phase;
        ((CExpParm*)parm_aux3)->a = 0;
        ((CExpParm*)parm_aux3)->b = m_signalSize-1;
        setRealAtom(parm_aux3);
        norm3 = m_rAtomNorm;
        ((CExpParm*)parm_aux3)->innerProd = ((CExpParm*)parm_aux2)->innerProd*(norm3/norm2);

        cgRealAtom.fillVector(m_realAtom);
        cgRealAtomAux = cgRealAtom*((CExpParm*)parm_aux3)->innerProd;
        res = residue - cgRealAtomAux;
        res_norm = res.norm();
        cout << rho_fix << " "<< rho_aux << " " << res_norm << " " << minResNorm << " " << coef_xi0 << " " << ratomsample_xi0 << endl;
        if( fabs(res_norm) < fabs(minResNorm) )
        {
            minResNorm = res_norm;
            ((CExpParm*)parm_aux)->innerProd = ((CExpParm*)parm_aux3)->innerProd;
            rho_fix = rho_aux;
            ((CExpParm*)parm_aux)->rho = rho_fix;
            rho_aux = rho_fix + increment_rho;
            increment_rho = pow(4.0,(double)incr_count) * increment_rho;
            incr_count++;
//            cout << "Olhai!!!" << endl;
//            ((CExpParm*)parm_aux)->printParm2Screen();
        }
        else
        {
            increment_rho =  increment_rho / pow(4.0,(double)incr_count);
            increment_rho = - increment_rho / 2.0;
            rho_aux = rho_fix + increment_rho;
            incr_count = 0;
        }
        imaxiter++;
        if (imaxiter==nmaxiter) break;
    }

//    cout << "innerProd depois da otim." << innerProd << endl;

/*    if (fabs(rho_fix) < 1e-7)
    {
        cout << "= Adjust to pure sine" << endl;
        ((CExpParm*)parm_aux)->rho      = 0.0;
        setRealAtom(parm_aux);
        ((CExpParm*)parm_aux)->innerProd =ratom_sample/m_realAtom[0];
    }
    else if(rho_fix>4)
    {
        cout << "= Adjust to impulse from decreasing" << endl;
        ((CExpParm*)parm_aux)->b = ((CExpParm*)parm_aux)->a;
        ((CExpParm*)parm_aux)->xi = 0.0;
        ((CExpParm*)parm_aux)->rho = 0.0;
        ((CExpParm*)parm_aux)->phase = 0.0;
        setRealAtom(parm_aux);
        ((CExpParm*)parm_aux)->innerProd =ratom_sample/m_realAtom[0];
    }
    else if(rho_fix<-4)
    {
        cout << "= Adjust to impulse from increasing" << endl;
        ((CExpParm*)parm_aux)->a = ((CExpParm*)parm_aux)->b;
        ((CExpParm*)parm_aux)->xi = 0.0;
        ((CExpParm*)parm_aux)->rho = 0.0;
        ((CExpParm*)parm_aux)->phase = 0.0;
        ((CExpParm*)parm_aux)->innerProd = residue[0][((CExpParm*)parm_aux)->b];
    }
*/
//    else
//    {
//        cout << "IF 4" << endl;
//        ((CExpParm*)parm_aux)->rho = rho_fix;
//        setRealAtom(parm_aux);
//        ((CExpParm*)parm_aux)->innerProd =ratom_sample/m_realAtom[0];
//    }

//    cout << "innerProd depois do ajuste." << ((CExpParm*)parm_aux)->innerProd << endl;
//    cout << "Olhai de novo!!" << endl;
//    ((CExpParm*)parm_aux)->printParm2Screen();

    cout << "    END:  ";
    ((CExpParm*)parm_aux)->printParm2Screen();

    delete parm_aux2;
    delete parm_aux3;
}


//===============================================================
// Function: quantizeFrequency
// Goal:
// Return:
//===============================================================

void CExpDictionary::quantizeFrequency( cgMatrix<double>& residue,
                                        CParameter* parm,
                                        CDataSignal* dataSignal)
{
    CParameter* parm_aux;
    parm_aux = new CExpParm;

    double Ffund = ((CComtradeSignal*)dataSignal)->getFundamentalFrequency();
    double Fs = ((CComtradeSignal*)dataSignal)->getSamplingRate(1);

    double rf;
    rf = RF_COEF * (Fs / Ffund);

    double aux;
    aux = log(rf) / log(2.0);

    int num_bits_xi;
    num_bits_xi = (int)(aux);

    double step_xi;
    step_xi = (2*pi) / rf;

    double xi = ((CExpParm*)parm)->xi;

    // Quantization
    xi =  floor( (xi + step_xi/2) / step_xi) * step_xi;



    cgMatrix<double> innerProduct;

    double opt_phase = 0;

    ((CExpParm*)parm)->xi = xi;
    setComplexAtom(parm);

    opt_phase = computeOptimumPhase(residue,
                                    xi);
    cout << opt_phase << endl;

    ((CExpParm*)parm)->phase = opt_phase;
    setRealAtom(parm);

    //realAtomVector.fillVector(realAtom);

    innerProduct = residue * m_realAtom;//( realAtomVector.transpose() );

    double highestProduct = innerProduct.getData(0,0);

    // Recalculating rho by Pseudo-Newton method
    double rho = ((CExpParm*)parm)->rho;
    double rho_fix, rho_aux, increment_rho;

    int cont=0;

    /*if( fabs(rho)< (.1 / (double)expDic.getSignalSize() ) )
    {
        // rho close to zero
        if(rho>0)
            increment_rho = -.05 / (double)expDic.getSignalSize();
        else
            increment_rho =  .05 / (double)expDic.getSignalSize();
        cont = 0;
    }
    else
    {
        cont = -1000;
        increment_rho = -rho / 2;
    }*/

    rho_fix = rho;
    increment_rho = rho * 0.5; // 2;
    rho_aux = rho_fix + increment_rho;

    // ===================
    // Pseudo-Newton
    cout<< "Pseudo-Newton" << endl;

    double ph_fix = opt_phase;

    int signalSize = m_signalSize;

    double rho_lim  = rho_aux / (double) ( /*10000*/2 * pow(signalSize,3.0));

    int incr_count=0;

    while(  //(cont < 1000) &&
            (fabs(increment_rho)  > fabs(rho_lim)) )
    {
        cont++;

        ((CExpParm*)parm_aux)->rho      = rho_aux;
        ((CExpParm*)parm_aux)->xi       = xi;
        ((CExpParm*)parm_aux)->phase    = 0.0;
        ((CExpParm*)parm_aux)->a        = ((CExpParm*)parm)->a;
        ((CExpParm*)parm_aux)->b        = ((CExpParm*)parm)->b;

        setComplexAtom(parm_aux);

        // Computing optimum phase
        opt_phase = computeOptimumPhase(residue,
                                            xi);
        ((CExpParm*)parm_aux)->phase    = opt_phase;
        setRealAtom(parm_aux);

        //realAtomVector.fillVector(realAtom);

        innerProduct = residue * m_realAtom; //( realAtomVector.transpose() );

        //cout << "rho_fix:" << rho_fix << " increment rho: " << increment_rho<< endl;

        if ( fabs(innerProduct.getData(0,0)) > fabs(highestProduct) )
        {
            highestProduct = innerProduct.getData(0,0);
            rho_fix = rho_aux;
            rho_aux = rho_fix + increment_rho;
            ph_fix = opt_phase;
            increment_rho = pow(4.0,(double)incr_count) * increment_rho;
            incr_count++;
        }
        else
        {
            increment_rho /= pow(4.0,(double)incr_count);
            increment_rho = /*- .75*/- 0.5 * increment_rho;
            rho_aux = rho_fix + increment_rho;
            incr_count =0;
        }

    }

    // Fill the damped sinusoid structure
    ((CExpParm*)parm)->innerProd = highestProduct;
    ((CExpParm*)parm)->rho = rho_fix;
    ((CExpParm*)parm)->xi = xi;
    ((CExpParm*)parm)->phase = ph_fix;

    delete parm_aux;
}


//===============================================================
// Function: discrimineSine
// Goal:
// Return:
//===============================================================

void CExpDictionary::discrimineSine(    cgMatrix<double>& residue,
                                        CParameter* parm,
                                        CDataSignal* dataSignal)
{
    CParameter* parm_aux;
    parm_aux = new CExpParm;

    cgMatrix<double> realAtomVector(    1,
                                        m_signalSize,
                                        0.0);
    cgMatrix<double> residueAuxVector(  1,
                                        m_signalSize,
                                        0.0);

    double dampSin_opt_phase;

    cgMatrix<double> dampSinInnerProduct;

    setRealAtom(parm);

    realAtomVector.fillVector(m_realAtom);
    int i;
    // Calculating error per sample
    for(i = ((CExpParm*)parm)->a;
        i < (((CExpParm*)parm)->b + 1);
        i++)
    {
        residueAuxVector[0][i] = residue[0][i] -
        (((CExpParm*)parm)->innerProd * realAtomVector[0][i]);
    }

    double dampSinErrorPerSample;
    dampSinErrorPerSample = residueAuxVector.norm() /
                            (double) (((CExpParm*)parm)->b -
                            ((CExpParm*)parm)->a + 1);

    double previousErrorPerSample;
    previousErrorPerSample = dampSinErrorPerSample;

    // Generate pure sinusoid (rho = 0)
    double pureSin_opt_phase = 0;

    ((CExpParm*)parm_aux)->rho      = 0.0;
    ((CExpParm*)parm_aux)->xi       = ((CExpParm*)parm)->xi;
    ((CExpParm*)parm_aux)->phase    = 0.0;
    ((CExpParm*)parm_aux)->a        = ((CExpParm*)parm)->a;
    ((CExpParm*)parm_aux)->b        = ((CExpParm*)parm)->b;

    setComplexAtom(parm_aux);

    // Computing optimum phase
    pureSin_opt_phase = computeOptimumPhase(residue,
                                    ((CExpParm*)parm_aux)->xi);

    ((CExpParm*)parm_aux)->phase    = pureSin_opt_phase;
    setRealAtom(parm_aux);

    cgMatrix<double> pureSinInnerProduct;

    realAtomVector.fillVector(m_realAtom);

    pureSinInnerProduct = residue * m_realAtom;//( realAtomVector.transpose() );

    // Calculating error per sample
    residueAuxVector.zeros();
    for(i = ((CExpParm*)parm)->a;
        i < (((CExpParm*)parm)->b + 1);
        i++)
    {
        residueAuxVector[0][i] = residue[0][i] -
        (pureSinInnerProduct.getData(0,0) * realAtomVector[0][i]);
    }

    double pureSinErrorPerSample;
    pureSinErrorPerSample = residueAuxVector.norm() /
                            (double) (((CExpParm*)parm)->b -
                            ((CExpParm*)parm)->a + 1);

    //============================

    double currentProduct = ((CExpParm*)parm)->innerProd;
    double rho;
    double rho_fix = ((CExpParm*)parm)->rho;
    double xi = ((CExpParm*)parm)->xi;
    //double u = cExpStructure.u;
    //double u_fix = u;
    int new_a = ((CExpParm*)parm)->a;
    int new_b = ((CExpParm*)parm)->b;
    double phase = ((CExpParm*)parm)->phase;

    double increment_t;
    double Ffund = ((CComtradeSignal*)dataSignal)->getFundamentalFrequency();
    double Fs = ((CComtradeSignal*)dataSignal)->getSamplingRate(1);
    double rf = RF_COEF * (Fs / Ffund);
    increment_t = rf;

    int nSampleCycle = (int)floor((2*pi)/xi);

    cout <<"Pure Sin Inner Prod: "<< pureSinInnerProduct.getData(0,0) << endl;

    /*char file1[40];
    char file2[40];
    char file3[40];
    char step_str[4];

    sprintf(step_str,"%03d",step);

    strcpy(file1,"psin_prod_step");
    strcat(file1,step_str);
    strcpy(file2,"psin_eps_step");
    strcat(file2,step_str);
    strcpy(file3,"dsin_eps_step");
    strcat(file3,step_str);
    */

    // First pure sinusoid test
    if (fabs(pureSinInnerProduct.getData(0,0)) >
        fabs(0.9999* /*20 */  ((CExpParm*)parm)->innerProd) )
    {
        ((CExpParm*)parm)->innerProd = pureSinInnerProduct.getData(0,0);
        ((CExpParm*)parm)->rho = 0;
        ((CExpParm*)parm)->phase = pureSin_opt_phase;
        //cExpStructure.u = new_a;
        //cExpStructure.a = new_a;
        //cExpStructure.b = new_b;
    }
    // Second sinusoid test
    else if(fabs(pureSinInnerProduct.getData(0,0)) >
            fabs(0.75 * ((CExpParm*)parm)->innerProd) )
    {
        rho = ((CExpParm*)parm)->rho;
        currentProduct = fabs (0.75 * ((CExpParm*)parm)->innerProd);
        //currentProduct = 0.95*pureSinInnerProduct.getData(0,0);

        //increment_t = (2*pi) / xi;
        /*new_a = cExpStructure.a - (int)floor(increment_t/2.0);
        new_b = cExpStructure.b + (int)floor(increment_t/2.0);
        if ( new_a < 0 )
            new_a = 0;
        if ( new_b > (expDic.getSignalSize() - 1) )
            new_b = expDic.getSignalSize() - 1;
        cout << "new_a_init: " << new_a << endl;
        cout << "new_b_init: " << new_b << endl;
        */
        // For rho > 0 (Decreasing Sinusoid)
        if(((CExpParm*)parm)->rho>0)
        {
            //cout << "rho > 0!!!" << endl;
            new_b = ((CExpParm*)parm)->b + (int)floor(increment_t/2.0);
            if ( new_b > (m_signalSize - 1) ) new_b = m_signalSize - 1;
            new_a = ((CExpParm*)parm)->a;
            //cout << "new_a_init: " << new_a << endl;
            //cout << "new_b_init: " << new_b << endl;
            /*strcat(file1,"_decr.dat");
            ofstream outFile1(file1);
            strcat(file2,"_decr.dat");
            ofstream outFile2(file2);
            strcat(file3,"_decr.dat");
            ofstream outFile3(file3);
            */

            // From right to left
            for(int b_aux = new_b; b_aux >= new_a; b_aux--)
            {
                // Procedures for Pure Sinusoid
                ((CExpParm*)parm_aux)->rho      = 0.0;
                ((CExpParm*)parm_aux)->xi       = ((CExpParm*)parm)->xi;
                ((CExpParm*)parm_aux)->phase    = 0.0;
                ((CExpParm*)parm_aux)->a        = new_a;
                ((CExpParm*)parm_aux)->b        = b_aux;

                setComplexAtom(parm_aux);

                // Computing optimum phase
                pureSin_opt_phase = computeOptimumPhase(residue,
                                                ((CExpParm*)parm_aux)->xi);

                ((CExpParm*)parm_aux)->phase    = pureSin_opt_phase;
                setRealAtom(parm_aux);

                realAtomVector.fillVector(m_realAtom);

                //realAtomVector.PrintToFile("pure_sin.dat");

                pureSinInnerProduct = residue * m_realAtom;//( realAtomVector.transpose() );

                // Calculating error per sample
                residueAuxVector.zeros();
                for(i = new_a;
                    i < (b_aux + 1);
                    i++)
                {
                    residueAuxVector[0][i] = residue[0][i] -
                    (pureSinInnerProduct.getData(0,0) * realAtomVector[0][i]);
                }

                pureSinErrorPerSample = residueAuxVector.norm() /
                                        (double) (b_aux - new_a + 1);

                // Procedures for damping sinusoid
                ((CExpParm*)parm_aux)->rho      = ((CExpParm*)parm)->rho;
                ((CExpParm*)parm_aux)->xi       = ((CExpParm*)parm)->xi;
                ((CExpParm*)parm_aux)->phase    = 0.0;
                ((CExpParm*)parm_aux)->a        = new_a;
                ((CExpParm*)parm_aux)->b        = b_aux;

                setComplexAtom(parm_aux);

                // Computing optimum phase
                dampSin_opt_phase = computeOptimumPhase(residue,
                                                ((CExpParm*)parm_aux)->xi);

                ((CExpParm*)parm_aux)->phase    = dampSin_opt_phase;
                setRealAtom(parm_aux);

                realAtomVector.fillVector(m_realAtom);

                //realAtomVector.PrintToFile("damp_sin.dat");

                dampSinInnerProduct = residue * m_realAtom;//( realAtomVector.transpose() );

                // Calculating error per sample
                residueAuxVector.zeros();
                for ( i = new_a; i < (b_aux + 1); i++)
                {
                    residueAuxVector[0][i] = residue[0][i] -
                    (dampSinInnerProduct.getData(0,0) * realAtomVector[0][i]);
                }

                dampSinErrorPerSample = residueAuxVector.norm() /
                                        (double) (b_aux - new_a + 1);

                /*outFile1 << fabs(pureSinInnerProduct.getData(0,0)) << ' ';
                outFile2 << pureSinErrorPerSample << ' ';
                outFile3 << dampSinErrorPerSample << ' ';
                */

                if (( fabs(/*0.95*/currentProduct) <= fabs(pureSinInnerProduct.getData(0,0)) ) &&
                    ( fabs(pureSinErrorPerSample) <= fabs(dampSinErrorPerSample) )  &&
                    ( fabs(pureSinErrorPerSample) <= fabs(previousErrorPerSample/2.0) ) )
                {
                    currentProduct = pureSinInnerProduct.getData(0,0);
                    ((CExpParm*)parm)->innerProd = currentProduct;

                    ((CExpParm*)parm)->b = b_aux;
                    //u_fix = new_a;
                    //cExpStructure.u = new_a;
                    //cExpStructure.a = new_a;
                    //rho_fix = 0;
                    ((CExpParm*)parm)->rho = 0;
                    //phase = pureSin_opt_phase;
                    ((CExpParm*)parm)->phase = pureSin_opt_phase;
                    previousErrorPerSample = pureSinErrorPerSample;
                }
            }

            /*outFile1 << endl;
            outFile2 << endl;
            outFile3 << endl;
            */
            new_a = ((CExpParm*)parm)->a - (int)floor(increment_t/2.0);
            if ( new_a < 0 ) new_a = 0;
            new_b = ((CExpParm*)parm)->b;
            //cout << "new_a_intermediate:" << new_a << endl;
            //cout << "new_b_intermediate:" << new_b << endl;

            //=================================
            // From left to right
            for (int a_aux = new_a; a_aux <= new_b; a_aux++)
            {
                // Procedures for Pure Sinusoid
                ((CExpParm*)parm_aux)->rho      = 0.0;
                ((CExpParm*)parm_aux)->xi       = ((CExpParm*)parm)->xi;
                ((CExpParm*)parm_aux)->phase    = 0.0;
                ((CExpParm*)parm_aux)->a        = a_aux;
                ((CExpParm*)parm_aux)->b        = new_b;

                setComplexAtom(parm_aux);

                // Computing optimum phase
                pureSin_opt_phase = computeOptimumPhase(residue,
                                                ((CExpParm*)parm_aux)->xi);

                ((CExpParm*)parm_aux)->phase    = pureSin_opt_phase;
                setRealAtom(parm_aux);


                realAtomVector.fillVector(m_realAtom);

                pureSinInnerProduct = residue * m_realAtom;//( realAtomVector.transpose() );

                // Calculating error per sample
                residueAuxVector.zeros();
                for(i = a_aux;
                    i < (new_b + 1);
                    i++)
                {
                    residueAuxVector[0][i] = residue[0][i] -
                    (pureSinInnerProduct.getData(0,0) * realAtomVector[0][i]);
                }

                pureSinErrorPerSample = residueAuxVector.norm() /
                                        (double) (new_b - a_aux + 1);

                // Procedures for damping sinusoid
                ((CExpParm*)parm_aux)->rho      = ((CExpParm*)parm)->rho;
                ((CExpParm*)parm_aux)->xi       = ((CExpParm*)parm)->xi;
                ((CExpParm*)parm_aux)->phase    = 0.0;
                ((CExpParm*)parm_aux)->a        = a_aux;
                ((CExpParm*)parm_aux)->b        = new_b;

                setComplexAtom(parm_aux);

                // Computing optimum phase
                dampSin_opt_phase = computeOptimumPhase(residue,
                                                ((CExpParm*)parm_aux)->xi);

                ((CExpParm*)parm_aux)->phase    = dampSin_opt_phase;
                setRealAtom(parm_aux);

                realAtomVector.fillVector(m_realAtom);

                dampSinInnerProduct = residue * m_realAtom;// ( realAtomVector.transpose() );

                // Calculating error per sample
                residueAuxVector.zeros();
                for ( i = a_aux; i < (new_b + 1); i++)
                {
                    residueAuxVector[0][i] = residue[0][i] -
                    (dampSinInnerProduct.getData(0,0) * realAtomVector[0][i]);
                }

                dampSinErrorPerSample = residueAuxVector.norm() /
                                        (double) (new_b - a_aux + 1);

                /*outFile1 << pureSinInnerProduct.getData(0,0) << endl;
                outFile2 << pureSinErrorPerSample << ' ';
                outFile3 << dampSinErrorPerSample << ' ';
                */

                if (( fabs(/*0.95*/currentProduct) <= fabs(pureSinInnerProduct.getData(0,0)) ) &&
                    ( fabs(pureSinErrorPerSample) <= fabs(dampSinErrorPerSample) ) &&
                    ( fabs(pureSinErrorPerSample) <= fabs(previousErrorPerSample/2.0) ) )
                {
                    currentProduct = pureSinInnerProduct.getData(0,0);
                    ((CExpParm*)parm)->innerProd = currentProduct;
                    //new_a = a_aux;
                    ((CExpParm*)parm)->a = a_aux;
                    //u_fix = a_aux;
                    //cExpStructure.u = new_a;
                    //rho_fix = 0;
                    ((CExpParm*)parm)->rho = 0;
                    //phase = pureSin_opt_phase;
                    ((CExpParm*)parm)->phase = pureSin_opt_phase;
                    previousErrorPerSample = pureSinErrorPerSample;
                }

            } // end a_aux

            new_a = ((CExpParm*)parm)->a;
            //cout << "new_a_final:" << new_a << endl;
            new_b = ((CExpParm*)parm)->b;
            //cout << "new_b_final:" << new_b << endl;

        } // end if rho >=0
        else if (((CExpParm*)parm)->rho < 0) // rho <= 0 (Increasing sinusoid)
        {
            //cout << "rho < 0!!!" << endl;
            new_a = ((CExpParm*)parm)->a - (int)floor(increment_t/2.0);
            if ( new_a < 0 ) new_a = 0;
            new_b = ((CExpParm*)parm)->b;
            //cout << "new_a_init:" << new_a << endl;
            //cout << "new_b_init:" << new_b << endl;
            /*strcat(file1,"_incr.dat");
            ofstream outFile1(file1);
            strcat(file2,"_incr.dat");
            ofstream outFile2(file2);
            strcat(file3,"_incr.dat");
            ofstream outFile3(file3);
            */

            //=================================
            // From left to right
            for (int a_aux = new_a; a_aux <= new_b; a_aux++)
            {
                // Procedures for Pure Sinusoid
                ((CExpParm*)parm_aux)->rho      = 0.0;
                ((CExpParm*)parm_aux)->xi       = ((CExpParm*)parm)->xi;
                ((CExpParm*)parm_aux)->phase    = 0.0;
                ((CExpParm*)parm_aux)->a        = new_a;
                ((CExpParm*)parm_aux)->b        = a_aux;

                setComplexAtom(parm_aux);

                // Computing optimum phase
                pureSin_opt_phase = computeOptimumPhase(residue,
                                                ((CExpParm*)parm_aux)->xi);

                ((CExpParm*)parm_aux)->phase    = pureSin_opt_phase;
                setRealAtom(parm_aux);

                realAtomVector.fillVector(m_realAtom);

                //realAtomVector.PrintToFile("pure_sin.dat");

                pureSinInnerProduct = residue * m_realAtom;//( realAtomVector.transpose() );

                // Calculating error per sample
                residueAuxVector.zeros();
                for(i = a_aux;
                    i < (new_b + 1);
                    i++)
                {
                    residueAuxVector[0][i] = residue[0][i] -
                    (pureSinInnerProduct.getData(0,0) * realAtomVector[0][i]);
                }

                pureSinErrorPerSample = residueAuxVector.norm() /
                                        (double) (new_b - a_aux + 1);

                // Procedures for damping sinusoid
                ((CExpParm*)parm_aux)->rho      = ((CExpParm*)parm)->rho;
                ((CExpParm*)parm_aux)->xi       = ((CExpParm*)parm)->xi;
                ((CExpParm*)parm_aux)->phase    = 0.0;
                ((CExpParm*)parm_aux)->a        = new_a;
                ((CExpParm*)parm_aux)->b        = a_aux;

                setComplexAtom(parm_aux);

                // Computing optimum phase
                dampSin_opt_phase = computeOptimumPhase(residue,
                                                ((CExpParm*)parm_aux)->xi);

                ((CExpParm*)parm_aux)->phase    = dampSin_opt_phase;
                setRealAtom(parm_aux);

                realAtomVector.fillVector(m_realAtom);

                //realAtomVector.PrintToFile("damp_sin.dat");

                dampSinInnerProduct = residue * m_realAtom;//( realAtomVector.transpose() );

                // Calculating error per sample
                residueAuxVector.zeros();
                for ( i = a_aux; i < (new_b + 1); i++)
                {
                    residueAuxVector[0][i] = residue[0][i] -
                    (dampSinInnerProduct.getData(0,0) * realAtomVector[0][i]);
                }

                dampSinErrorPerSample = residueAuxVector.norm() /
                                        (double) (new_b - a_aux + 1);

                /*tFile1 << fabs(pureSinInnerProduct.getData(0,0)) << ' ';
                outFile2 << pureSinErrorPerSample << ' ';
                outFile3 << dampSinErrorPerSample << ' ';
                */
                if (( fabs(/*0.95*/currentProduct) <= fabs(pureSinInnerProduct.getData(0,0)) ) &&
                    ( fabs(pureSinErrorPerSample) <= fabs(dampSinErrorPerSample) )  &&
                    ( fabs(pureSinErrorPerSample) <= fabs(previousErrorPerSample/2.0) ) )
                {
                    currentProduct = pureSinInnerProduct.getData(0,0);
                    ((CExpParm*)parm)->innerProd = currentProduct;
                    //new_a = a_aux;
                    ((CExpParm*)parm)->a = a_aux;
                    //u_fix = (expDic.getSignalSize() - 1 - new_b);
                    //cExpStructure.u = (expDic.getSignalSize() - 1 - new_b);
                    //rho_fix = 0;
                    ((CExpParm*)parm)->rho = 0;
                    //phase = pureSin_opt_phase;
                    ((CExpParm*)parm)->phase = pureSin_opt_phase;
                    previousErrorPerSample = pureSinErrorPerSample;
                }

            } // end for a_aux

            new_b = ((CExpParm*)parm)->b + (int)floor(increment_t/2.0);
            if ( new_b > (m_signalSize - 1) ) new_b = m_signalSize - 1;
            new_a = ((CExpParm*)parm)->a;
            //cout << "new_a_intermediate: " << new_a << endl;
            //cout << "new_b_intermediate: " << new_b << endl;

            /*outFile1 << endl;
            outFile2 << endl;
            outFile3 << endl;*/

            // From right to left
            for(int b_aux = new_b; b_aux >= new_a; b_aux--)
            {
                // Procedures for Pure Sinusoid
                ((CExpParm*)parm_aux)->rho      = 0.0;
                ((CExpParm*)parm_aux)->xi       = ((CExpParm*)parm)->xi;
                ((CExpParm*)parm_aux)->phase    = 0.0;
                ((CExpParm*)parm_aux)->a        = new_a;
                ((CExpParm*)parm_aux)->b        = b_aux;

                setComplexAtom(parm_aux);

                // Computing optimum phase
                pureSin_opt_phase = computeOptimumPhase(residue,
                                                ((CExpParm*)parm_aux)->xi);

                ((CExpParm*)parm_aux)->phase    = pureSin_opt_phase;
                setRealAtom(parm_aux);

                realAtomVector.fillVector(m_realAtom);

                //realAtomVector.PrintToFile("pure_sin.dat");

                pureSinInnerProduct = residue * m_realAtom;//( realAtomVector.transpose() );

                // Calculating error per sample
                residueAuxVector.zeros();
                for(i = new_a;
                    i < (b_aux + 1);
                    i++)
                {
                    residueAuxVector[0][i] = residue[0][i] -
                    (pureSinInnerProduct.getData(0,0) * realAtomVector[0][i]);
                }

                pureSinErrorPerSample = residueAuxVector.norm() /
                                        (double) (b_aux - new_a + 1);

                // Procedures for damping sinusoid
                ((CExpParm*)parm_aux)->rho      = ((CExpParm*)parm)->rho;
                ((CExpParm*)parm_aux)->xi       = ((CExpParm*)parm)->xi;
                ((CExpParm*)parm_aux)->phase    = 0.0;
                ((CExpParm*)parm_aux)->a        = new_a;
                ((CExpParm*)parm_aux)->b        = b_aux;

                setComplexAtom(parm_aux);

                // Computing optimum phase
                dampSin_opt_phase = computeOptimumPhase(residue,
                                                ((CExpParm*)parm_aux)->xi);

                ((CExpParm*)parm_aux)->phase    = dampSin_opt_phase;
                setRealAtom(parm_aux);

                realAtomVector.fillVector(m_realAtom);

                //realAtomVector.PrintToFile("damp_sin.dat");

                dampSinInnerProduct = residue * m_realAtom;//( realAtomVector.transpose() );

                // Calculating error per sample
                residueAuxVector.zeros();
                for ( i = new_a; i < (b_aux + 1); i++)
                {
                    residueAuxVector[0][i] = residue[0][i] -
                    (dampSinInnerProduct.getData(0,0) * realAtomVector[0][i]);
                }

                dampSinErrorPerSample = residueAuxVector.norm() /
                                        (double) (b_aux - new_a + 1);


                if (( fabs(/*0.95*/currentProduct) <= fabs(pureSinInnerProduct.getData(0,0)) ) &&
                    ( fabs(pureSinErrorPerSample) <= fabs(dampSinErrorPerSample) )  &&
                    ( fabs(pureSinErrorPerSample) <= fabs(previousErrorPerSample/2.0) ) )
                {
                    currentProduct = pureSinInnerProduct.getData(0,0);
                    ((CExpParm*)parm)->innerProd = currentProduct;
                    ((CExpParm*)parm)->b = b_aux;
                    ((CExpParm*)parm)->rho = 0;
                    ((CExpParm*)parm)->phase = pureSin_opt_phase;
                    previousErrorPerSample = pureSinErrorPerSample;
                }
            } // end for b_aux

            new_a = ((CExpParm*)parm)->a;
            new_b = ((CExpParm*)parm)->b;
            //cout << "new_a_final: " << new_a << endl;
            //cout << "new_b_final: " << new_b << endl;
        } // end if rho <= 0

        //if    ((fabs(currentProduct / (double)(new_b - new_a +1)) > fabs(cExpStructure.innerProduct / (double)(cExpStructure.b - cExpStructure.a + 1) ) ) &&
        /*if    ( fabs(currentProduct) > fabs(0.75 * cExpStructure.innerProduct) ) //)
        {
            cExpStructure.innerProduct = currentProduct;
            cExpStructure.rho = rho_fix;
            cExpStructure.u = u_fix;
            cExpStructure.phase = phase;
            cExpStructure.a = new_a;
            cExpStructure.b = new_b;
        }*/
    } // end second test

    //cout << "cExpStructure.a :" << ((CExpParm*)parm)->a << endl;
    //cout << "cExpStructure.b :" << ((CExpParm*)parm)->b << endl;

    delete parm_aux;
}

//===============================================================
// Function: searchSBPreviousBlock
// Goal:
// Return:
//===============================================================
void CExpDictionary::searchSBPreviousBlock(		cgMatrix<double>& residue,
												CParameter* parm,
												CStructBook* sbPreviousBlock,
												int iPrevBlock,
												int flagFile,
												double coefHeur)
{
	CParameter* parm_aux;
	CParameter* parm_aux2;
	parm_aux = new CExpParm;
	parm_aux2 = new CExpParm;

	CStructBook* sbPreviousAux;


	if (flagFile == 0)
	{
		sbPreviousAux = sbPreviousBlock;
	}
	else if (flagFile == 1)
	{
		iPrevBlock=0;
		sbPreviousAux = new CStructBookExp;
		char fileName[40];
		strcpy(fileName,"previous_block_sb.dat");
		FILE* stream;
		stream = fopen(fileName,"r");
		while (!feof(stream))
		{
			((CStructBookExp*)sbPreviousAux)[iPrevBlock].loadElementASCII(stream);
		}
		if (stream!=NULL) fclose(stream);
	}


	double maxInnerProd=0;
	double opt_phase = 0;
	((CExpParm*)parm_aux2)->innerProd=0.0;
	((CExpParm*)parm_aux2)->rho = 0.0;
	((CExpParm*)parm_aux2)->xi = 0.0;
	((CExpParm*)parm_aux2)->phase = 0.0;
	((CExpParm*)parm_aux2)->a = 0;
	((CExpParm*)parm_aux2)->b = 0;

	cgMatrix<double> innerProduct;
	double innerProd;

	int sbNumElement = ((CStructBookExp*)sbPreviousAux)[iPrevBlock].getNumElement();
	int i,b;
	for (i=0; i< sbNumElement; i++)
	{
		b = (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].b;
		if (b==m_signalSize-1)
		{
			((CExpParm*)parm_aux)->rho		= (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].rho;
			((CExpParm*)parm_aux)->xi		= (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].xi;
			((CExpParm*)parm_aux)->phase	= 0.0;
			((CExpParm*)parm_aux)->a		= 0;
			((CExpParm*)parm_aux)->b		= m_signalSize -1;

			setComplexAtom(parm_aux);

			// Computing optimum phase
			opt_phase = computeOptimumPhase(residue,
											((CExpParm*)parm_aux)->xi,
											innerProd);

			((CExpParm*)parm_aux)->phase	= opt_phase;


			if (fabs(innerProd)>fabs(maxInnerProd))
			{
				maxInnerProd = innerProd;
				((CExpParm*)parm_aux2)->innerProd = innerProd;
				((CExpParm*)parm_aux2)->rho = ((CExpParm*)parm_aux)->rho;
				((CExpParm*)parm_aux2)->xi = ((CExpParm*)parm_aux)->xi;
				((CExpParm*)parm_aux2)->phase = ((CExpParm*)parm_aux)->phase ;
				((CExpParm*)parm_aux2)->a = ((CExpParm*)parm_aux)->a;
				((CExpParm*)parm_aux2)->b = ((CExpParm*)parm_aux)->b;
			}
		}
	}

	if (((CExpParm*)parm_aux2)->innerProd!=0.0)
	{
		findFastBestTimeSupport(residue,parm_aux2,2,1.02);
	}

	cout << "$$ PARM $$" << endl;
	((CExpParm*)parm)->printParm2Screen();
	cout << "$$ PARM_PREVIOUS $$" << endl;
	((CExpParm*)parm_aux2)->printParm2Screen();

	FILE* stream;
	stream = fopen("eval_prevblocksb.out","a");
	double innerPPrev = ((CExpParm*)parm_aux2)->innerProd;
	double innerPCurrent = ((CExpParm*)parm)->innerProd;



	if ( fabs( ((CExpParm*)parm_aux2)->innerProd ) >  coefHeur*fabs( ((CExpParm*)parm)->innerProd) )
	{
		((CExpParm*)parm)->innerProd = ((CExpParm*)parm_aux2)->innerProd;
		((CExpParm*)parm)->rho		 = ((CExpParm*)parm_aux2)->rho;
		((CExpParm*)parm)->xi        = ((CExpParm*)parm_aux2)->xi;
		((CExpParm*)parm)->phase     = ((CExpParm*)parm_aux2)->phase ;
		((CExpParm*)parm)->a         = ((CExpParm*)parm_aux2)->a;
		((CExpParm*)parm)->b         = ((CExpParm*)parm_aux2)->b;
		fprintf(stream,"%10.7f %10.7f %i\n",innerPCurrent ,innerPPrev,1);
	}
	else
	{
		fprintf(stream,"%10.7f %10.7f %i\n",innerPCurrent ,innerPPrev,0);
	}
	fclose(stream);

	delete parm_aux;
	delete parm_aux2;
	if (flagFile==1) delete sbPreviousAux;
}

//===============================================================
// Function: searchSBPreviousBlock
// Goal: search for minimum error norm with phase and amplitude continuity
// Return:
//===============================================================
/*void CExpDictionary::searchSBPreviousBlock(     cgMatrix<double>& residue,
                                                CParameter* parm,
                                                CStructBook* sbPreviousBlock,
                                                int iPrevBlock,
                                                int flagFile,
                                                double coefHeur)
{
    CParameter* parm_aux;
    CParameter* parm_aux2;
    parm_aux = new CExpParm;
    parm_aux2 = new CExpParm;

    CStructBook* sbPreviousAux;


    if (flagFile == 0)
    {
        sbPreviousAux = sbPreviousBlock;
    }
    else if (flagFile == 1)
    {
        iPrevBlock=0;
        sbPreviousAux = new CStructBookExp;
        char fileName[40];
        strcpy(fileName,"previous_block_sb.dat");
        FILE* stream;
        stream = fopen(fileName,"r");
        while (!feof(stream))
        {
            ((CStructBookExp*)sbPreviousAux)[iPrevBlock].loadElementASCII(stream);
        }
        if (stream!=NULL) fclose(stream);
    }


    double opt_phase = 0.0;
    double delta_phase=0.0;
    ((CExpParm*)parm_aux2)->innerProd=0.0;
    ((CExpParm*)parm_aux2)->rho = 0.0;
    ((CExpParm*)parm_aux2)->xi = 0.0;
    ((CExpParm*)parm_aux2)->phase = 0.0;
    ((CExpParm*)parm_aux2)->a = 0;
    ((CExpParm*)parm_aux2)->b = 0;

    cgMatrix<double> innerProduct;
    double innerProd, ratom_sample;
    double res_norm;
    double minResNorm = 10000000000.0;
    cgMatrix<double> cgRealAtom(1,m_signalSize,0.0);
    cgMatrix<double> cgRealAtomAux(1,m_signalSize,0.0);
    cgMatrix<double> res(1,m_signalSize,0.0);

    int sbNumElement = ((CStructBookExp*)sbPreviousAux)[iPrevBlock].getNumElement();
    if (sbNumElement==0)
    {
        cout << "No candidates atoms with continuity in previous block" << endl;
        return;
    }
    int i,j,chosenContAtom;
    // Search for the candidate atom with continuity that leads the smallest error regarding the residue
    for (i=0; i< sbNumElement; i++)
    {
        // Compute the sample m_signalSize of the previous block
        ((CExpParm*)parm_aux)->innerProd= (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].innerProduct;
        ((CExpParm*)parm_aux)->rho      = (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].rho;
        ((CExpParm*)parm_aux)->xi       = (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].xi;
        ((CExpParm*)parm_aux)->phase    = (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].phase;
        ((CExpParm*)parm_aux)->a        = (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].a;
        ((CExpParm*)parm_aux)->b        = (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].b;
        setRealAtom(parm_aux);
        //cout << (((CExpParm*)parm_aux)->innerProd)<< endl;
        //cout << computeUnormRAtomSample(parm_aux,m_signalSize) <<endl;
        //cout << m_rAtomNorm <<endl;

        ratom_sample = (((CExpParm*)parm_aux)->innerProd) *
                    (computeUnormRAtomSample(parm_aux,m_signalSize)/m_rAtomNorm);
        //cout << "ratom_sample: " << ratom_sample << endl;

        // Compute the phase shift and amplitude for continuity
        delta_phase = m_signalSize*((CExpParm*)parm_aux)->xi +
                      (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].phase;
        delta_phase = delta_phase - floor(delta_phase/(2*pi))*(2*pi);
        ((CExpParm*)parm_aux)->phase    = delta_phase;
        ((CExpParm*)parm_aux)->a        = 0;
        ((CExpParm*)parm_aux)->b        = m_signalSize -1;// (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].b;

        if (((CExpParm*)parm_aux)->rho<0.0)
        {
            setRealAtom(parm_aux);
            ((CExpParm*)parm_aux)->innerProd =ratom_sample/m_realAtom[0];
            cgRealAtom.fillVector(m_realAtom);
            cgRealAtomAux = cgRealAtom*((CExpParm*)parm_aux)->innerProd;
            // Reduce atom support from right to left if decreasing exponential atom
            int binit = ((CExpParm*)parm_aux)->b;
            for (j=binit;j>=0;j--)
            {
                if (j!=binit) cgRealAtomAux[0][j+1] =0.0;
                res = residue - cgRealAtomAux;
                res_norm = res.norm();
                ((CExpParm*)parm_aux)->b = j;
                setRealAtom(parm_aux);
                ((CExpParm*)parm_aux)->innerProd =ratom_sample/m_realAtom[0];
                if ( (fabs(res_norm)<fabs(minResNorm)) && (fabs(((CExpParm*)parm_aux)->innerProd)>1e-7))
                {
                    minResNorm = res_norm;
                    ((CExpParm*)parm_aux2)->rho = ((CExpParm*)parm_aux)->rho;
                    ((CExpParm*)parm_aux2)->xi = ((CExpParm*)parm_aux)->xi;
                    ((CExpParm*)parm_aux2)->phase = ((CExpParm*)parm_aux)->phase ;
                    ((CExpParm*)parm_aux2)->a = ((CExpParm*)parm_aux)->a;
                    ((CExpParm*)parm_aux2)->b = ((CExpParm*)parm_aux)->b;
                    ((CExpParm*)parm_aux2)->innerProd = ((CExpParm*)parm_aux)->innerProd;
                    chosenContAtom = i;
                }
    //     cout << i << res.norm() << endl;
            }
            optimizeDecayingErrorNorm(residue,parm_aux2,minResNorm,ratom_sample);
        }
        else
        {
            setRealAtom(parm_aux);
            ((CExpParm*)parm_aux)->innerProd = ratom_sample/m_realAtom[0];
            cgRealAtom.fillVector(m_realAtom);
            cgRealAtomAux = cgRealAtom*((CExpParm*)parm_aux)->innerProd;
            res = residue - cgRealAtomAux;
            res_norm = res.norm();
            optimizeDecayingErrorNorm(residue,parm_aux,res_norm,ratom_sample);
            if ( (fabs(res_norm)<fabs(minResNorm)) && (fabs(((CExpParm*)parm_aux)->innerProd)>1e-7))
            {
                minResNorm = res_norm;
                ((CExpParm*)parm_aux2)->innerProd = ((CExpParm*)parm_aux)->innerProd;
                ((CExpParm*)parm_aux2)->rho = ((CExpParm*)parm_aux)->rho;
                ((CExpParm*)parm_aux2)->xi = ((CExpParm*)parm_aux)->xi;
                ((CExpParm*)parm_aux2)->phase = ((CExpParm*)parm_aux)->phase ;
                ((CExpParm*)parm_aux2)->a = ((CExpParm*)parm_aux)->a;
                ((CExpParm*)parm_aux2)->b = ((CExpParm*)parm_aux)->b;
                chosenContAtom = i;
            }
        }
    }

    // Compute residue of the best matching in the block
    setRealAtom(parm);
    cgRealAtom.fillVector(m_realAtom);
    res = residue - cgRealAtom*((CExpParm*)parm)->innerProd;
    double res_norm_bestmatch=0;
    res_norm_bestmatch = res.norm();
    cgRealAtomAux= cgRealAtom*((CExpParm*)parm)->innerProd;
//    cgRealAtomAux.PrintToFile("bestmatch_atom.dat");

//    cout << "minResNorm" << minResNorm << endl;
//    cout << "res_norm_bestmatch" << res_norm_bestmatch << endl;

//    cout << "Tecle qq coisa!!!" << endl;
//    int qqcoisa;
//    cin >> qqcoisa;

    cout << "$$ PARM $$" << endl;
    ((CExpParm*)parm)->printParm2Screen();
    cout << "$$ PARM_PREVIOUS $$" << endl;
    ((CExpParm*)parm_aux2)->printParm2Screen();
    cout << "$$$$$$$$$$$$$$$$$$$" << endl;

    FILE* stream;
    stream = fopen("eval_prevblocksb.out","a");

    if ( minResNorm <  coefHeur*res_norm_bestmatch )
    {
        ((CExpParm*)parm)->innerProd = ((CExpParm*)parm_aux2)->innerProd;
        ((CExpParm*)parm)->rho       = ((CExpParm*)parm_aux2)->rho;
        ((CExpParm*)parm)->xi        = ((CExpParm*)parm_aux2)->xi;
        ((CExpParm*)parm)->phase     = ((CExpParm*)parm_aux2)->phase;
        ((CExpParm*)parm)->a         = ((CExpParm*)parm_aux2)->a;
        ((CExpParm*)parm)->b         = ((CExpParm*)parm_aux2)->b;
        fprintf(stream,"%10.7f %10.7f %i %10.7f\n",res_norm_bestmatch ,minResNorm,1,coefHeur*res_norm_bestmatch);
        ((CStructBookExp*)sbPreviousAux)[iPrevBlock].removeElement(chosenContAtom);
    }
    else
    {
        fprintf(stream,"%10.7f %10.7f %i %10.7f\n",res_norm_bestmatch ,minResNorm,0,coefHeur*res_norm_bestmatch);
    }
    fclose(stream);

    delete parm_aux;
    delete parm_aux2;
    if (flagFile==1) delete sbPreviousAux;
}
*/
//===============================================================
// Function: evalAtomContinuity
// Goal:
// Return:
//===============================================================
void CExpDictionary::evalAtomContinuity(cgMatrix<double>& residue,
                                        CStructBook* sbContinuity,
                                        CStructBook* structBook,
                                        int iSignal,
                                        int iCurrBlock,
                                        int iPrevBlock,
                                        int& step,
                                        int flagFile,
                                        char* fileName,
                                        char* sbbFName,
                                        double initBlockNorm,
                                        fstream& file_stage,
                                        int flag_stage)
{
    CParameter* parm_aux;
    parm_aux = new CExpParm;
    CParameter* parm_aux2;
    parm_aux2 = new CExpParm;
    CParameter* parm_aux3;
    parm_aux3 = new CExpParm;
    CParameter* parm_aux4;
    parm_aux4 = new CExpParm;

    CStructBook* sbPreviousAux;

    if (flagFile == 0)
    {
        sbPreviousAux = sbContinuity;
    }
    else if (flagFile == 1)
    {
        iPrevBlock=0;
        sbPreviousAux = new CStructBookExp;
        char fileName[40];
        strcpy(fileName,"previous_block_sb.dat");
        FILE* stream;
        stream = fopen(fileName,"r");
        while (!feof(stream))
        {
            ((CStructBookExp*)sbPreviousAux)[iPrevBlock].loadElementASCII(stream);
        }
        if (stream!=NULL) fclose(stream);
    }

    double delta_phase=0.0;
    ((CExpParm*)parm_aux)->innerProd=0.0;
    ((CExpParm*)parm_aux)->rho = 0.0;
    ((CExpParm*)parm_aux)->xi = 0.0;
    ((CExpParm*)parm_aux)->phase = 0.0;
    ((CExpParm*)parm_aux)->a = 0;
    ((CExpParm*)parm_aux)->b = 0;

    cgMatrix<double> innerProduct;
    double innerProd, ratom_sample;
    double res_norm;
    double atom_norm,atom_norm_xi0;
    cgMatrix<double> cgRealAtom(1,m_signalSize,0.0);
    cgMatrix<double> cgRealAtomAux(1,m_signalSize,0.0);
    cgMatrix<double> res(1,m_signalSize,0.0);

    cgMatrix<double> cgRealAtom1(1,m_signalSize,0.0);
    cgMatrix<double> cgRealAtom2(1,2*m_signalSize,0.0);
    cgMatrix<double> cgRealAtom3(1,m_signalSize,0.0);
    cgMatrix<double> cgRealAtom4(1,m_signalSize,0.0);
    double norm1=0;
    double norm2=0;
    double norm3=0;
    double norm4=0;

    cgMatrix<double> cgVectorAux1(1,m_signalSize,0.0);
    cgMatrix<double> cgVectorAux2(1,2*m_signalSize,0.0);


    int sbNumElement = ((CStructBookExp*)sbPreviousAux)[iPrevBlock].getNumElement();
    if (sbNumElement==0)
    {
        cout << "No atoms with continuity in previous block" << endl;
        return;
    }
    int i=0,j=0;
    FILE* iosba;
    FILE* iosbb;
    for (i=0; i< sbNumElement; i++)
    {
        cout << "########################################################################################" << endl;
        cout << "->Decomposing Signal: " << iSignal+1 << "; Block: "<< iCurrBlock+1 << "; Atom: "<< step+1 << endl;
        cout << "########################################################################################" << endl;

        ((CExpParm*)parm_aux)->innerProd= (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].innerProduct;
        ((CExpParm*)parm_aux)->rho      = (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].rho;
        ((CExpParm*)parm_aux)->xi       = (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].xi;
        ((CExpParm*)parm_aux)->phase    = (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].phase;
        ((CExpParm*)parm_aux)->a        = (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].a;
        ((CExpParm*)parm_aux)->b        = (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].b;
        // Calculate atom norm
        setRealAtom(parm_aux);
        norm1 = m_rAtomNorm;
        cgRealAtom1.fillVector(m_realAtom);

//         cgVectorAux1 = cgRealAtom1*((CExpParm*)parm_aux)->innerProd;
//         cgVectorAux1.PrintToFile("ratom1.out");

        cout << "Bloco anterior" << endl;
        ((CExpParm*)parm_aux)->printParm2Screen();

        if (flag_stage==1)
        {
            int decomp_stage=91;
            file_stage  <<  setw (10) << setfill(' ') << iSignal+1 << " "
                        <<  setw (10) << setfill(' ') << iCurrBlock+1 << " "
                        <<  setw (10) << setfill(' ') << step+1 << " "
                        <<  setw (10) << setfill(' ') << decomp_stage << " "
                        <<  setw (20) << setfill(' ') << ((CExpParm*)parm_aux)->innerProd << " "
                        <<  setw (20) << setfill(' ') << ((CExpParm*)parm_aux)->rho << " "
                        <<  setw (20) << setfill(' ') << ((CExpParm*)parm_aux)->xi << " "
                        <<  setw (20) << setfill(' ') << ((CExpParm*)parm_aux)->phase << " "
                        <<  setw (10) << setfill(' ') << ((CExpParm*)parm_aux)->a << " "
                        <<  setw (10) << setfill(' ') << ((CExpParm*)parm_aux)->b << " "
                        << endl;
        }

        // Calculate norm with xi and phase equal zero


        ((CExpParm*)parm_aux2)->rho = ((CExpParm*)parm_aux)->rho;
        ((CExpParm*)parm_aux2)->xi = 0.0;
        ((CExpParm*)parm_aux2)->phase = 0.0;
        ((CExpParm*)parm_aux2)->a =  ((CExpParm*)parm_aux)->a;
        ((CExpParm*)parm_aux2)->b = m_signalSize-1;
        for (j=0;j<2*m_signalSize;j++)
        {
            cgRealAtom2[0][j] = computeUnormRAtomSampleSigSize(parm_aux2,j,2*m_signalSize);
        }
        norm2 = cgRealAtom2.norm();
        cgRealAtom2 /= norm2;
        ((CExpParm*)parm_aux2)->innerProd = ((CExpParm*)parm_aux)->innerProd*(norm2/norm1);
//         cgVectorAux2 = cgRealAtom2*((CExpParm*)parm_aux2)->innerProd;
//         cgVectorAux2.PrintToFile("ratom2.out");

        ((CExpParm*)parm_aux3)->rho = ((CExpParm*)parm_aux2)->rho;
        ((CExpParm*)parm_aux3)->xi = ((CExpParm*)parm_aux2)->xi;
        ((CExpParm*)parm_aux3)->phase = ((CExpParm*)parm_aux2)->phase;
        ((CExpParm*)parm_aux3)->a = 0;
        ((CExpParm*)parm_aux3)->b = m_signalSize-1;
        setRealAtom(parm_aux3);
        norm3 = m_rAtomNorm;
        cgRealAtom3.fillVector(m_realAtom);
        ((CExpParm*)parm_aux3)->innerProd = cgRealAtom2.dprodInterval(m_realAtom,m_signalSize,2*m_signalSize-1) *
                                            ((CExpParm*)parm_aux2)->innerProd;
//         cgVectorAux1 = cgRealAtom3*((CExpParm*)parm_aux3)->innerProd;
//         cgVectorAux1.PrintToFile("ratom3.out");

        ((CExpParm*)parm_aux4)->rho = ((CExpParm*)parm_aux)->rho;
        ((CExpParm*)parm_aux4)->xi = ((CExpParm*)parm_aux)->xi;
        delta_phase = m_signalSize*((CExpParm*)parm_aux)->xi +
                      ((CExpParm*)parm_aux)->phase;
        delta_phase = delta_phase - floor(delta_phase/(2*pi))*(2*pi);
        ((CExpParm*)parm_aux4)->phase = delta_phase;
        ((CExpParm*)parm_aux4)->a = 0;
        ((CExpParm*)parm_aux4)->b = m_signalSize-1;
        setRealAtom(parm_aux4);
        norm4 = m_rAtomNorm;
        cgRealAtom4.fillVector(m_realAtom);
        ((CExpParm*)parm_aux4)->innerProd= ((CExpParm*)parm_aux3)->innerProd * (norm4/norm3);
//         cgVectorAux1 = cgRealAtom4*((CExpParm*)parm_aux4)->innerProd;
//         cgVectorAux1.PrintToFile("ratom4.out");

        // Rho optimization
        cgRealAtomAux = cgRealAtom4*((CExpParm*)parm_aux4)->innerProd;
        res = residue - cgRealAtomAux;
        res_norm = res.norm();
        optimizeDecayingErrorNorm(  residue, parm_aux4, res_norm,
                                    ((CExpParm*)parm_aux3)->innerProd,
                                    cgRealAtom3[0][0]);




        if (flag_stage==1)
        {
            int decomp_stage=92;
            file_stage  <<  setw (10) << setfill(' ') << iSignal+1 << " "
                        <<  setw (10) << setfill(' ') << iCurrBlock+1 << " "
                        <<  setw (10) << setfill(' ') << step+1 << " "
                        <<  setw (10) << setfill(' ') << decomp_stage << " "
                        <<  setw (20) << setfill(' ') << ((CExpParm*)parm_aux4)->innerProd << " "
                        <<  setw (20) << setfill(' ') << ((CExpParm*)parm_aux4)->rho << " "
                        <<  setw (20) << setfill(' ') << ((CExpParm*)parm_aux4)->xi << " "
                        <<  setw (20) << setfill(' ') << ((CExpParm*)parm_aux4)->phase << " "
                        <<  setw (10) << setfill(' ') << ((CExpParm*)parm_aux4)->a << " "
                        <<  setw (10) << setfill(' ') << ((CExpParm*)parm_aux4)->b << " "
                        << endl;
        }

        if ( (fabs(((CExpParm*)parm_aux4)->innerProd) > 1e-8) && (fabs(((CExpParm*)parm_aux4)->rho) < 4) )
        {
            double norm = residue.norm();
            // adjust parameters (coef/phase)
            if(((CExpParm*)parm_aux4)->innerProd < 0.0)
            {
                ((CExpParm*)parm_aux4)->phase += pi;
                ((CExpParm*)parm_aux4)->innerProd = - ((CExpParm*)parm_aux4)->innerProd;
            }

            if (((CExpParm*)parm_aux4)->phase >= (2*pi) )
            {
                ((CExpParm*)parm_aux4)->phase -= 2*pi;
            }

            if (((CExpParm*)parm_aux4)->phase < 0 )
            {
                ((CExpParm*)parm_aux4)->phase += 2*pi;
            }
            //////////////////////////////////////
            //update residue
            setRealAtom(parm_aux4);
            cgRealAtom.fillVector(m_realAtom);
            cgRealAtomAux = cgRealAtom*((CExpParm*)parm_aux4)->innerProd;
            residue = residue - cgRealAtomAux;

            //cout << "- Atom parameters after opt. decaying!!!" << endl;
            //((CExpParm*)parm_aux4)->printParm2Screen();

            // Add element to structure book
            ((CStructBookExp*)sbContinuity)[iCurrBlock].addElement(parm_aux4);
            int ind = ((CStructBookExp*)sbContinuity)[iCurrBlock].getNumElement();
			((CStructBookExp*)sbContinuity)[iPrevBlock].setNextAtomIndex(i,ind-1);
			((CStructBookExp*)sbContinuity)[iCurrBlock].setPrevAtomIndex(ind-1,i);
            ///////////////////
            int origAtomIndex = ((CStructBookExp*)sbContinuity)[iPrevBlock].getOrigAtomIndex(i);
            ((CStructBookExp*)structBook)[iCurrBlock].addElement(parm_aux4);
            int indorig = ((CStructBookExp*)structBook)[iCurrBlock].getNumElement();
			((CStructBookExp*)structBook)[iPrevBlock].setNextAtomIndex(origAtomIndex,indorig-1);
			((CStructBookExp*)structBook)[iCurrBlock].setPrevAtomIndex(indorig-1,origAtomIndex);
            ((CStructBookExp*)structBook)[iCurrBlock].setNextAtomIndex(indorig-1,-1);
            ((CStructBookExp*)structBook)[iCurrBlock].setOrigAtomIndex(indorig-1,indorig-1);
            ((CStructBookExp*)sbContinuity)[iCurrBlock].setOrigAtomIndex(ind-1,indorig-1);

            iosba = fopen(fileName,"a");
            ((CStructBookExp*)structBook)[iCurrBlock].saveElementASCII( iosba,
                                                                        0.0,//meanApproxRatio,
                                                                        fabs(((CExpParm*)parm_aux4)->innerProd)/norm,//approxRatio[step%L],
                                                                        0.0,//befSupInnerP,
                                                                        0.0,//aftSupInnerP,
                                                                        residue.norm()/initBlockNorm);//(norm/initBlockNorm));
            fflush(iosba);
            fclose(iosba);

            iosbb = fopen(sbbFName,"ab");
            ((CStructBookExp*)structBook)[iCurrBlock].saveElementBin(iosbb);
            fclose(iosbb);

            step++;
        }
        else
        {
            cout << "!!!!!! Valor absoluto do coeficiente menor que 1e-8 ou decaimento > 4 !!!!!!" << endl;
            cout << "- abs(coef): "<< fabs(((CExpParm*)parm_aux4)->innerProd) << endl;
            cout << "- abs(rho): "<< fabs(((CExpParm*)parm_aux4)->rho) << endl;
        }
    }

    delete parm_aux;
    delete parm_aux2;
    delete parm_aux3;
    delete parm_aux4;
    if (flagFile==1) delete sbPreviousAux;
}

//===============================================================
// Function: evalAtomContinuity
// Goal:
// Return:
//===============================================================
/*void CExpDictionary::evalAtomContinuity(        cgMatrix<double>& residue,
                                                CParameter* parm,
                                                CStructBook* sbContinuity,
                                                int iPrevBlock)
{
    CParameter* parm_aux;
    parm_aux = new CExpParm;

    FILE* iosba;
    FILE* iosbb;

    cgMatrix<double> cgRealAtom(1,m_signalSize,0.0);

    double opt_phase = 0.0;
    double delta_phase=0.0;

    cgMatrix<double> innerProduct;
    double innerProd;

    int sbNumElement = ((CStructBookExp*)sbPreviousAux)[iPrevBlock].getNumElement();
    int i,b;
    double rho,ratom_sample;
    for (i=0; i< sbNumElement; i++)
    {
        rho = (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].rho;
        b = (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].b;
        if ((rho>=0.0) && (b==m_signalSize-1))
        {
            // Compute the sample m_signalSize of the previous block
            ((CExpParm*)parm_aux)->rho      = (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].rho;
            ((CExpParm*)parm_aux)->xi       = (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].xi;
            ((CExpParm*)parm_aux)->phase    = (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].phase;
            ((CExpParm*)parm_aux)->a        = (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].a;
            ((CExpParm*)parm_aux)->b        = (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].b;
            setRealAtom(parm_aux);
            ((CExpParm*)parm_aux)->b = m_signalSize;
            ratom_sample = (((CExpParm*)parm_aux)->innerProd) *
                           (computeUnormRAtomSample(parm_aux,m_signalSize)/m_rAtomNorm);

            // Compute the phase shift and amplitude for continuity
            delta_phase = m_signalSize*((CExpParm*)parm_aux)->xi;
            delta_phase = delta_phase - floor(delta_phase/(2*pi))*(2*pi);
            ((CExpParm*)parm_aux)->phase    = delta_phase +
                                              (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].phase;
            ((CExpParm*)parm_aux)->a        = 0;
            ((CExpParm*)parm_aux)->b        = (((CStructBookExp*)sbPreviousAux)[iPrevBlock].getStructBook())[i].b;
            setRealAtom(parm_aux);
            ((CExpParm*)parm_aux)->innerProd = ratom_sample/m_realAtom[0];

        }
    }

    delete parm_aux;
    if (flagFile==1) delete sbPreviousAux;
}
*/
//===============================================================
// Function: adjustParameters
// Goal:
// Return:
//===============================================================

void CExpDictionary::adjustParameters(  cgMatrix<double>& residue,
                                        CParameter* parm)
{

    cgMatrix<double> innerProd;

    // Adjust impulse atom
    int pos=0;
    double highestValue=0;
    if (fabs(((CExpParm*)parm)->rho) > 4)
    {
        setRealAtom(parm);

        for(int i=((CExpParm*)parm)->a; i<=((CExpParm*)parm)->b; i++)
        {
            if ( fabs(m_realAtom[i])>fabs(highestValue) )
            {
                highestValue = m_realAtom[i];
                pos = i;
            }
        }
        ((CExpParm*)parm)->a = pos;
        ((CExpParm*)parm)->b = pos;
        ((CExpParm*)parm)->xi = 0.0;
        ((CExpParm*)parm)->rho = 0.0;
        ((CExpParm*)parm)->phase = 0.0;

        setRealAtom(parm);

        innerProd = residue * m_realAtom;// (realAtomVector.transpose());

        ((CExpParm*)parm)->innerProd = innerProd.getData(0,0);

    }

    if (fabs(((CExpParm*)parm)->rho) < 1e-7)
    {
        ((CExpParm*)parm)->rho=0.0;
        setRealAtom(parm);

        innerProd = residue * m_realAtom;// (realAtomVector.transpose());

        ((CExpParm*)parm)->innerProd = innerProd.getData(0,0);
     }

    // Adjust innerProd signal => phase

    if(((CExpParm*)parm)->innerProd < 0.0)
    {
        ((CExpParm*)parm)->phase += pi;
        ((CExpParm*)parm)->innerProd = - ((CExpParm*)parm)->innerProd;
    }

    if (((CExpParm*)parm)->phase >= (2*pi) )
    {
        ((CExpParm*)parm)->phase -= 2*pi;
    }

    if (((CExpParm*)parm)->phase < 0 )
    {
        ((CExpParm*)parm)->phase += 2*pi;
    }

    //((CExpParm*)parm)->phase = fabs(((CExpParm*)parm)->phase);



    /*vec_aux.fillVector(realAtom);

    vec_aux *= cExpStructure.innerProduct;

    vec_aux.PrintToFile("aft_adjust.dat");*/

}

double CExpDictionary::getApproxRatio(int signalSize)
{
    double tolAppRatio;
    //          Dimensions:   64      128     256     512   1024    2048    4096
    double lambda_med_ger[7]={0.4439, 0.3256, 0.2361, 0.17, 0.1203, 0.0881, 0.07046518};

    if (signalSize==64)
    {
        tolAppRatio = lambda_med_ger[0];
    }
    else if (signalSize==128)
    {
        tolAppRatio = lambda_med_ger[1];
    }
    else if (signalSize==256)
    {
        tolAppRatio = lambda_med_ger[2];
    }
    else if (signalSize==512)
    {
        tolAppRatio = lambda_med_ger[3];
    }
    else if (signalSize==1024)
    {
        tolAppRatio = lambda_med_ger[4];
    }
    else if (signalSize==2048)
    {
        tolAppRatio = lambda_med_ger[5];
    }
    else if (signalSize==4096)
    {
        tolAppRatio = lambda_med_ger[6];
    }
    else
    {
        tolAppRatio = 0.0;
    }
    return tolAppRatio;
}

void CExpDictionary::computeNextBlockCoefPhase(CParameter* parm)
{
    // Calculate atom norm
    setRealAtom(parm);
    double atom_norm = m_rAtomNorm;

    // Calculate norm with xi and phase equal zero
    double xi_orig = ((CExpParm*)parm)->xi;
    double phase_orig = ((CExpParm*)parm)->phase;
    ((CExpParm*)parm)->xi = 0.0;
    ((CExpParm*)parm)->phase = 0.0;

    // Compute the sample m_signalSize of the previous block
    double ratom_sample = (((CExpParm*)parm)->innerProd) *
                            (computeUnormRAtomSample(parm,m_signalSize)/atom_norm);

    // Compute the phase shift and amplitude for continuity
    ((CExpParm*)parm)->a        = 0;
    ((CExpParm*)parm)->b        = m_signalSize -1;
    setRealAtom(parm);
    double atom_norm_xi0 = m_rAtomNorm;
    // Compute the inner product
    ((CExpParm*)parm)->innerProd = ratom_sample/m_realAtom[0];
    // Return xi value and phase adjustment
    ((CExpParm*)parm)->xi = xi_orig;
    double delta_phase =    m_signalSize*((CExpParm*)parm)->xi +
                            phase_orig;
    delta_phase = delta_phase - floor(delta_phase/(2*pi))*(2*pi);
    ((CExpParm*)parm)->phase    = delta_phase;
    setRealAtom(parm);
    atom_norm = m_rAtomNorm;
    ((CExpParm*)parm)->innerProd = ((CExpParm*)parm)->innerProd*
                                    (atom_norm/atom_norm_xi0);
}


