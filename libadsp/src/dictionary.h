// dictionary.h
#ifndef DICTIONARY_H
#define DICTIONARY_H

#include <iostream>
#include <stdlib.h>
#include <math.h>

#include <fstream>
#include <iomanip> 
#include <istream>

#include "fftw3.h"
#include "cgfunc.h"
#include "cgmatrix.h"
#include "complex.h"
#include "defines.h"
#include "datasignal.h"
#include "parameter.h"
#include "filemgr.h"
#include "structbook.h"

using namespace std;

// ========================================
//  CDictionary Class
// ========================================

/**
 * \brief Dictionary class
 */

class CDictionary {

protected:
    int         m_signalSize;
    CComplex*   m_complexAtom;
    double*     m_realAtom;
    double      m_rAtomNorm; // Eh atribuido um valor ao se gerar o Atomo

public:
    // Constructors
    CDictionary ();

    // Destructor
    virtual ~CDictionary();

    // Set methods
    virtual void setSignalSize(int signalSize)=0;
    virtual void setComplexAtom(CParameter* parm)=0;
    virtual void setRealAtom(CParameter* parm)=0;
    virtual void setRealAtomOnSupport(CParameter* parm)=0;
    // Get methods
    virtual int         getSignalSize()=0;
    virtual CComplex*   getComplexAtom()=0;
    virtual double*     getRealAtom()=0;
    double              getRAtomNorm();
    virtual double      computeUnormRAtomSample(CParameter* parm,int t)=0;
    
    // Signal decomposition
    virtual void executeDecomp( CDataSignal* dataSignal, 
                                CStructBook** structBook,
                                CFileGenData* genData,                                
                                CFileDictionary* dicData)=0;

    virtual void executeDecompElectric( CDataSignal* dataSignal,
                                        CStructBook** structBook,
                                        CFileDecomp* genData,            
                                        CFileDictionary* dicData,
                                        CFileDecompBlockRange* blockRange)=0;

};

/**
 * \brief Exponential dictionary class
 */

class CExpDictionary: public CDictionary {

public:


    // Constructors
    CExpDictionary (){};

    // Destructor
    virtual ~CExpDictionary(){};

    // Interface
    virtual void setSignalSize(int signalSize);
    virtual void setComplexAtom(CParameter* parm);
    virtual void setRealAtom(CParameter* parm);
    
    void setRealAtom(strtContinuousExp sb);
    
    virtual void setRealAtomOnSupport(CParameter* parm);
    
    void setRealAtomOnSupport(strtContinuousExp sb);

    virtual int         getSignalSize();
    virtual CComplex*   getComplexAtom();
    virtual double*     getRealAtom();
    double              computeUnormRAtomSample(CParameter* parm,int t);
    double computeUnormRAtomSampleSigSize(CParameter* parm,int t, int signalSize);

	// Signal decomposition
    virtual void executeDecomp( CDataSignal* dataSignal,
                                CStructBook** structBook,
                                CFileGenData* genData,                                
                                CFileDictionary* dicData);
    
    virtual void executeDecompElectric( CDataSignal* dataSignal,
                                        CStructBook** structBook,
                                        CFileDecomp* genData,                                
                                        CFileDictionary* dicData,
                                        CFileDecompBlockRange* blockRange);
    // *** Intrinsic ***
    // ---- Projection over discrete dictionary
    CExpParm fastMPKolasa(cgMatrix<double>& residue);
    void fastMPKolasaModified(  cgMatrix<double>& residue, 
								CDataSignal* dataSignal,
								CFileDictionary* dicData,
								CParameter* chosenParm);
    // ---- Optimum phase computation
    double computeOptimumPhase( cgMatrix<double>& residue, 
                                double xi);
    double computeOptimumPhase( cgMatrix<double>& residue, 
                                double xi,
                                double& innerProd);
    // Maximize approximation by finding optimum parameters
    void optimizeContinuousParms(   cgMatrix<double>& residue,
                                    CParameter* parm);
    // ---- Comtrade heuristics
    void proceedHeuristics( cgMatrix<double>& residue,
                            CParameter* parm,
                            CDataSignal* dataSignal);
        // with rho and xi optimization
    void findFastBestTimeSupport(   cgMatrix<double>& residue,
                                    CParameter* parm);
        // with rho optimization
    void findFastBestTimeSupport(   cgMatrix<double>& residue,
                                    CParameter* parm,
                                    int dummy,
                                    double coefHeur);
        // with rho optimization
    void findBestTimeSupport(   cgMatrix<double>& residue,
                                CParameter* parm,
                                int dummy);

    void optimizeDecaying(  cgMatrix<double>& residue,
                            CParameter* parm_aux);

    void optimizeDecayingErrorNorm(  cgMatrix<double>& residue,
                                     CParameter* parm_aux,
                                     double& minResNorm,
									 double coef_xi0,
									 double ratomsample_xi0);

    void quantizeFrequency( cgMatrix<double>& residue,
                            CParameter* parm,
                            CDataSignal* dataSignal);

    void discrimineSine(cgMatrix<double>& residue,
                        CParameter* parm,
                        CDataSignal* dataSignal);

    void searchSBPreviousBlock(	cgMatrix<double>& residue,
                                CParameter* parm,
                                CStructBook* sbPreviousBlock,
                                int iPrevBlock,
                                int flagFile,
                                double coefHeur);
//    void searchSBPreviousBlock( cgMatrix<double>& residue,
//                                CParameter* parm,
//                                CStructBook* structbook,
//                                int iPrevBlock,
//                                int flagFile,
//                                double coefHeur);

    void evalAtomContinuity(cgMatrix<double>& residue,
                            CStructBook* sbPreviousBlock,
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
							int flag_stage);

    void adjustParameters(  cgMatrix<double>& residue,
                            CParameter* parm);

    double getApproxRatio(int signalSize);

    //=============================================================

    void computeNextBlockCoefPhase(CParameter* parm);
    
   
};

#endif
