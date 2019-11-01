
#include "decomp.h"

void decompElectric(CFileDecomp* genData,
		    		CFileDecompBlockRange* blockRange,
	    			CFileDictionary* dicData,
	    			char* InputFile)
{
	//////////////////////////////////////////
    // Configuring and loading signals
    CDataSignal* dataSignal;
	dataSignal = new CComtradeSignal;

    dataSignal->setFileName(InputFile);
    dataSignal->setBlockSize(blockRange->getBlockSize());
    dataSignal->setBlockHop(blockRange->getBlockHop());
    dataSignal->setSignal();
    dataSignal->setNorm();

    // -----------------------------------------------------------
    //  Decomposing signal
    //
    // Allocating memory for Structure Book Set
    int numBlock = (int)ceil((double)dataSignal->getSignalSize()/(double)dataSignal->getBlockHop());
    int numSignal = dataSignal->getNumSignal();
    CStructBook** structBook;
    structBook = new CStructBook* [numSignal];
    int i;
    for (i=0; i<numSignal; i++)
    {
        structBook[i] = new CStructBookExp[numBlock];
    }

    // Setting Dictionary
    CDictionary* dic;
    if (genData->getDicType()==1)
    {
	    dic = new CExpDictionary;
	}
    int nbits = (int)ceil( log10( (double)(dataSignal->getBlockSize()) )  /  log10( (double)(2) ) );
    int dicSize = (int) pow(2.0,(double)nbits);
    dic->setSignalSize(dicSize);

    // Decomposing Blocks of Signals

    int initBlock = blockRange->getInitBlock();
    int finalBlock = blockRange->getEndBlock();
    int nMaxStep = genData->getNumMaxStep();


    cgMatrix<double> residue(1,dicSize,0.0);
    cgMatrix<double> cgRealAtom(1,dicSize,0.0);
    cgMatrix<double> cgRealAtomAux(1,dicSize,0.0);
    double** pSignal = dataSignal->getSignal();

    // ====================================

    int step = 0;
    double norm =0;


    CParameter* chosenParm;

    if (genData->getDicType()==1)
    {
    	chosenParm = new CExpParm;
    }

    if (finalBlock==9999)
    {
        finalBlock = numBlock;
    }
    // Allocate memory for candidate atoms with tume continuity
    CStructBook** sbContinuity;
    sbContinuity = new CStructBook* [numSignal];
    int k;
    for (k=0; k<numSignal; k++)
    {
        sbContinuity[k] = new CStructBookExp[numBlock];
    }

    int L = (int)ceil(((log10((double)(dicSize)))/(log10((double)(2)))));
    double* approxRatio;
    approxRatio = new double[L];
    for (k=0;k<L;k++)
    {
        approxRatio[k] = 0.0;
    }
    double meanApproxRatio;

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
    fprintf(iosba,"Block Size:           %5i\n", dataSignal->getBlockSize());
    if (genData->getSigType()==1)
    {
        fprintf(iosba,"Samp. Freq :     %10.2f\n", ((CComtradeSignal*)dataSignal)->getSamplingRate(1));
    }
    if (genData->getSigType()==2)
    {
        fprintf(iosba,"Samp. Freq :     %10.2f\n", ((CAudioSignal*)dataSignal)->getSamplingRate());
    }
    fprintf(iosba,"Init. Block:          %5i\n", initBlock);
    fprintf(iosba,"Final Block:          %5i\n", finalBlock);
    fflush(iosba);
    fclose(iosba);

    char sbbFName[_MAX_PATH];
    strcpy(sbbFName, dataSignal->getFileName());
    pos = strrchr( sbbFName, '.');
    sprintf(aux,"_b%d-%d_header.sbb",initBlock,finalBlock);
    strcpy( &pos[0], aux);

    FILE* iosbb;
    iosbb = fopen(sbbFName,"wb");
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
    dummyint = initBlock;
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    dummyint = finalBlock;
    fwrite(&dummyint, sizeof(int), 1, iosbb);
    fclose(iosbb);

	strcpy(sbbFName, dataSignal->getFileName());
    pos = strrchr( sbbFName, '.');
    sprintf(aux,"_b%d-%d.sbb",initBlock,finalBlock);
    strcpy( &pos[0], aux);
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
                if ((int)(j*(double)dataSignal->getBlockHop())+dicSize < dataSignal->getSignalSize())
                {
                    residue.fillVector((pSignal[i])+(int)(j*(double)dataSignal->getBlockHop()),
                                        dicSize);
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


                int decomp_stage;

                // Beginning of the decomposition
                if (norm!=0)
                {
                	if ((j!=0) &&
                        (genData->getFlagEvalAtomCont()==1))
                    {
                        ((CExpDictionary*)dic)->evalAtomContinuity( residue,
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
                        ((CExpDictionary*)dic)->fastMPKolasaModified(residue,dataSignal,dicData,chosenParm);
                        ((CExpParm*)chosenParm)->printParm2Screen();
                        if (genData->getPrintDecompStage()==1)
                        {
                            decomp_stage=1;
                            file_stage  <<  setw (10) << setfill(' ') << i+1 << " "
                                        <<  setw (10) << setfill(' ') << j+1 << " "
                                        <<  setw (10) << setfill(' ') << step+1 << " "
                                        <<  setw (10) << setfill(' ') << decomp_stage << " "
                                        <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->innerProd << " "
                                        <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->rho << " "
                                        <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->xi << " "
                                        <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->phase << " "
                                        <<  setw (10) << setfill(' ') << ((CExpParm*)chosenParm)->a << " "
                                        <<  setw (10) << setfill(' ') << ((CExpParm*)chosenParm)->b << " "
                                        << endl;
                        }

                        // Calculate approximation ratio referred to this step
                        //cout<< "->Calculating approximation ratio..." << endl;
                        approxRatio[step%L] = fabs(((CExpParm*)chosenParm)->innerProd)/norm;
                        meanApproxRatio = 0.0;
                        for (k=0;k<L;k++)
                        {
                            meanApproxRatio += approxRatio[k]/(double)L;
                        }
                        cout << "Mean Approx. Ratio: " << meanApproxRatio << endl;
                        cout << "Tol. Approx. Ratio: " << tolAppRatio << endl;


                        // Maximize approximation by finding optimum parameters
                        //optimizeContinuousParms(residue, chosenParm);

                        // Heuristics
                        befSupInnerP = 0.0;
                        aftSupInnerP = 0.0;

                        if (((CExpParm*)chosenParm)->a != ((CExpParm*)chosenParm)->b)  // damp and pure cases
                        {
                            if (genData->getFlagFindSupport()==1)
                            {
                                befSupInnerP = ((CExpParm*)chosenParm)->innerProd;
                                cout<< "->Finding the best time support ..." << endl;
                                if (((CExpParm*)chosenParm)->rho==0.0)
                                {
                                    // Two-way search
                                    ((CExpDictionary*)dic)->findFastBestTimeSupport(residue,chosenParm,1,
                                                            						genData->getCoefTempSup());
                                }
                                else
                                {
                                    // One-way search
                                    ((CExpDictionary*)dic)->findFastBestTimeSupport(residue,chosenParm,0,
                                                            						genData->getCoefTempSup());
                                }
                                ((CExpParm*)chosenParm)->printParm2Screen();
                                if (genData->getPrintDecompStage()==1)
                                {
                                    decomp_stage=2;
                                    file_stage  <<  setw (10) << setfill(' ') << i+1 << " "
                                                <<  setw (10) << setfill(' ') << j+1 << " "
                                                <<  setw (10) << setfill(' ') << step+1 << " "
                                                <<  setw (10) << setfill(' ') << decomp_stage << " "
                                                <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->innerProd << " "
                                                <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->rho << " "
                                                <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->xi << " "
                                                <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->phase << " "
                                                <<  setw (10) << setfill(' ') << ((CExpParm*)chosenParm)->a << " "
                                                <<  setw (10) << setfill(' ') << ((CExpParm*)chosenParm)->b << " "
                                                << endl;
                                }
                                aftSupInnerP = ((CExpParm*)chosenParm)->innerProd;
                            }
                            if (genData->getFlagOptDecay()==1)
                            {
                                ((CExpDictionary*)dic)->optimizeDecaying(residue,chosenParm);
                                ((CExpParm*)chosenParm)->printParm2Screen();
                                if (genData->getPrintDecompStage()==1)
                                {
                                    decomp_stage=3;
                                    file_stage  <<  setw (10) << setfill(' ') << i+1 << " "
                                                <<  setw (10) << setfill(' ') << j+1 << " "
                                                <<  setw (10) << setfill(' ') << step+1 << " "
                                                <<  setw (10) << setfill(' ') << decomp_stage << " "
                                                <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->innerProd << " "
                                                <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->rho << " "
                                                <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->xi << " "
                                                <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->phase << " "
                                                <<  setw (10) << setfill(' ') << ((CExpParm*)chosenParm)->a << " "
                                                <<  setw (10) << setfill(' ') << ((CExpParm*)chosenParm)->b << " "
                                                << endl;
                                }
                            }

                        }
                        if (((CExpParm*)chosenParm)->rho != 0.0) // damp case
                        {
                            // Discrimine Sine
                            cout<< "->Discrimining sine..." << endl;
                            ((CExpDictionary*)dic)->discrimineSine(residue,chosenParm,dataSignal);
                            ((CExpParm*)chosenParm)->printParm2Screen();
                            if (genData->getPrintDecompStage()==1)
                            {
                                decomp_stage=4;
                                file_stage  <<  setw (10) << setfill(' ') << i+1 << " "
                                            <<  setw (10) << setfill(' ') << j+1 << " "
                                            <<  setw (10) << setfill(' ') << step+1 << " "
                                            <<  setw (10) << setfill(' ') << decomp_stage << " "
                                            <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->innerProd << " "
                                            <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->rho << " "
                                            <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->xi << " "
                                            <<  setw (20) << setfill(' ') << ((CExpParm*)chosenParm)->phase << " "
                                            <<  setw (10) << setfill(' ') << ((CExpParm*)chosenParm)->a << " "
                                            <<  setw (10) << setfill(' ') << ((CExpParm*)chosenParm)->b << " "
                                            << endl;
                            }
                        }


                        if ((j!=0) && (genData->getFlagSeqBlock()==1))
                        {
                            cout<< "->Search in the Previous Block Structure Book ..." << endl;
                            ((CExpDictionary*)dic)->searchSBPreviousBlock(  residue,
																			chosenParm,
																			sbContinuity[i],
																			j-1,
																			0,
																			genData->getCoefSeqBlock());
                            ((CExpParm*)chosenParm)->printParm2Screen();
                        }


                        // Adjusting parameters
                        ((CExpDictionary*)dic)->adjustParameters(residue,chosenParm);

                        // Print to screen the chosen atom parameter
                        cout << "Parameters adjusted!!" << endl;
                        cout << "Chosen atom parameters: " << endl;
                        ((CExpParm*)chosenParm)->printParm2Screen();

                        // Updating residue
                        cout<< "->Updating residue..." << endl << endl;

                        ((CExpDictionary*)dic)->setRealAtom(chosenParm);
                        cgRealAtom.fillVector(((CExpDictionary*)dic)->getRealAtom());
                        cgRealAtomAux = cgRealAtom*((CExpParm*)chosenParm)->innerProd;
                        //cgRealAtomAux.PrintToFile("scaled_atom.dat");
                        residue = residue - cgRealAtom*((CExpParm*)chosenParm)->innerProd;
                        norm = residue.norm();

                         // Add element to structure book
                        ((CStructBookExp*)structBook[i])[j].addElement(chosenParm);
                        int indorig = ((CStructBookExp*)structBook[i])[j].getNumElement();
                        ((CStructBookExp*)structBook[i])[j].setNextAtomIndex(indorig-1,-1);
                        ((CStructBookExp*)structBook[i])[j].setPrevAtomIndex(indorig-1,-1);
                        ((CStructBookExp*)structBook[i])[j].setOrigAtomIndex(indorig-1,indorig-1);

                        // Add element to structure book with candidate atom with continuity
                        // for the next block
                        if (((CExpParm*)chosenParm)->b==dicSize-1)
                        {
                            ((CStructBookExp*)sbContinuity[i])[j].addElement(chosenParm);
                            int ind = ((CStructBookExp*)sbContinuity[i])[j].getNumElement();
                            ((CStructBookExp*)sbContinuity[i])[j].setNextAtomIndex(ind-1,-1);
                            ((CStructBookExp*)sbContinuity[i])[j].setPrevAtomIndex(ind-1,-1);
                            ((CStructBookExp*)sbContinuity[i])[j].setOrigAtomIndex(ind-1,indorig-1);
                        }




                        iosba = fopen(fileName,"a");
                        ((CStructBookExp*)structBook[i])[j].saveElementASCII(   iosba,
                                                                                meanApproxRatio,
                                                                                approxRatio[step%L],
                                                                                0,
                                                                                0,
                                                                                (norm/initBlockNorm));
                        fflush(iosba);
                        fclose(iosba);

                        cout << "SNR: " << 20*log10(initBlockNorm/norm) << " (dB)"<< endl;
						cout << "SNR Target: " << snrTarget<< endl;

                        iosbb = fopen(sbbFName,"ab");
                        ((CStructBookExp*)structBook[i])[j].saveElementBin(iosbb);
                        fclose(iosbb);

                        step++;
                    }
                    while(      (   (meanApproxRatio > tolAppRatio) ||
                                    (step<(L+nAtomCont)) )
                                //((norm/initBlockNorm)>1e-8)
                                && (step<nMaxStep)
                                && (20*log10(initBlockNorm/norm)<snrTarget)
                                //&& ( fabs(expParm.innerProd) > 1e-8 )
                         );
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

    if (genData->getPrintDecompStage()==1)
    {
        file_stage.close();
    }


    // Deallocating objects
    delete [] approxRatio;

    for (i=0; i<dataSignal->getNumSignal(); i++)
    {
        delete [] ((CStructBookExp*)structBook[i]);
    }
    delete [] structBook;
    if (genData->getDicType()==1)
    {
	    delete (CExpDictionary*)dic;
	}

    delete (CComtradeSignal*)dataSignal;

    return;
}

void decompAudio(	CFileDecomp* genData,
		    		CFileDecompBlockRange* blockRange,
	    			CFileDictionary* dicData,
	    			char* InputFile)
{
	//////////////////////////////////////////
    // Configuring and loading signals
    CDataSignal* dataSignal;
    dataSignal = new CAudioSignal;

    dataSignal->setFileName(InputFile);
    dataSignal->setBlockSize(blockRange->getBlockSize());
    dataSignal->setBlockHop(blockRange->getBlockHop());
    dataSignal->setSignal();
    dataSignal->setNorm();

    // Allocating memory for Structure Book Set
    int numBlock = (int)ceil((double)dataSignal->getSignalSize()/(double)dataSignal->getBlockHop());
    int numSignal = dataSignal->getNumSignal();
    CStructBook** structBook;
    structBook = new CStructBook* [numSignal];
    int i;
    for (i=0; i<numSignal; i++)
    {
        structBook[i] = new CStructBookExp[numBlock];
    }

    // Setting Dictionary
    CExpDictionary* expDic = new CExpDictionary;
    int nbits = (int)ceil( log10( (double)(dataSignal->getBlockSize()) )  /  log10( (double)(2) ) );
    int sigSize = (int) pow(2.0,(double)nbits);
    expDic->setSignalSize(sigSize);

    // Decomp Stage



    // Deallocating objects
    delete expDic;
    for (i=0; i<dataSignal->getNumSignal(); i++)
    {
        delete [] ((CStructBookExp*)structBook[i]);
    }
    delete [] structBook;

    delete (CAudioSignal*)dataSignal;

    return;
}

void decompNoise(	CFileDecomp* genData,
		    		CFileDecompBlockRange* blockRange,
	    			CFileDictionary* dicData)
{
	//////////////////////////////////////////
    // Configuring and loading signals
    CDataSignal* dataSignal;
	dataSignal = new CNoiseSignal;

	dataSignal->setNumSignal(1000);
	dataSignal->setSignalSize(blockRange->getBlockSize());
	dataSignal->setSamplingRate(44100);
    dataSignal->setBlockSize(blockRange->getBlockSize());
    dataSignal->setBlockHop(blockRange->getBlockHop());
    dataSignal->setSignal();
    dataSignal->setNorm();

    // Allocating memory for Structure Book Set
    int numBlock = (int)ceil((double)dataSignal->getSignalSize()/(double)dataSignal->getBlockHop());
    int numSignal = dataSignal->getNumSignal();
    CStructBook** structBook;
    structBook = new CStructBook* [numSignal];
    int i;
    for (i=0; i<numSignal; i++)
    {
        structBook[i] = new CStructBookExp[numBlock];
    }

    // Setting Dictionary
    CExpDictionary* expDic = new CExpDictionary;
    int nbits = (int)ceil( log10( (double)(dataSignal->getBlockSize()) )  /  log10( (double)(2) ) );
    int sigSize = (int) pow(2.0,(double)nbits);
    expDic->setSignalSize(sigSize);

    // Decomp Stage



    // Deallocating objects
    delete expDic;
    for (i=0; i<dataSignal->getNumSignal(); i++)
    {
        delete [] ((CStructBookExp*)structBook[i]);
    }
    delete [] structBook;

    delete (CNoiseSignal*)dataSignal;

    return;
}
