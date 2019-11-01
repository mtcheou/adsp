#include "mpursuit.h"

//===============================================================
// Function: fastMPKolasaModified
// Goal: 
// Return:
//==============================================================

CParameter* fastMPKolasaModified(	cgMatrix<double>& residue,
									CDataSignal* dataSignal,
									CFileDictionary* dicData,
									CDictionary* dic)
{
    CParameter* parm;
    CParameter* chosenParm;
    if (dicData->getDicType()==1)
    {
    	parm = new CExpParm;
		chosenParm = new CExpParm;
		((CExpParm*)chosenParm)->innerProd = 0;
		((CExpParm*)chosenParm)->rho = 0;
		((CExpParm*)chosenParm)->xi = 0;
		((CExpParm*)chosenParm)->phase = 0;
		((CExpParm*)chosenParm)->a = 0;
		((CExpParm*)chosenParm)->b = 0;
	}
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
        
        if ((s > dic->getSignalSize())&&(s!=88888))
        {
            cout << "Scale greater than block size!!!" << endl;
            printf("Scale %d - block size %d\n",s,dic->getSignalSize());
            exit(1);
        }

        // Set linear frequencies for COMTRADE files
        if (dataSignal->getType() == 1) // Comtrade
        {
            if ( (freqf==9999999999) && (freqi==9999999999) )
            {
                Ffund = ((CComtradeSignal*)dataSignal)->getFundamentalFrequency();
                Fs = ((CComtradeSignal*)dataSignal)->getSamplingRate(1);
            } 
            else if (freqf==9999999999)
            {
                Ffund = freqi;
                Fs = ((CComtradeSignal*)dataSignal)->getSamplingRate(1);
            }
            delta_f = (2*pi/Fs)* Ffund;
            Nfreq = (int)(Fs/(2*Ffund));
            xi_vec = new double[Nfreq];
            for (i=0;i<Nfreq;i++)
            {
                xi_vec[i] = (double)i * delta_f;
            }
        }
        // Set linear/geometric frequencies space for audio or noise
        if ( (dataSignal->getType() == 2) || // Audio
             (dataSignal->getType() == 3) )   // Noise for Audio
        {
            if (dataSignal->getType() == 2) Fs = ((CAudioSignal*)dataSignal)->getSamplingRate();
            if (dataSignal->getType() == 3) Fs = ((CNoiseSignal*)dataSignal)->getSamplingRate();
            //double A0 = 27.5; // Hz
            //double C8 = 4186.01;
            //double freqi= A0 * pow(2.0,-23/12);
            //double freqf= C8 * pow(2.0,23/12); // B9

            if (fdiscrtype==1) // linear
            {
                Nfreq = (int)ceil(freqf/freqi);
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
        }

        if (s==1)
        {

            // s=1 Impulse
            for(i=0;i<N;i+=delta_tau)
            {
                innerProd = residue[0][i];
                if (fabs(innerProd)>fabs(maxInnerProd))
                {
                    maxInnerProd = innerProd;
                    ((CExpParm*)chosenParm)->innerProd = maxInnerProd;
                    ((CExpParm*)chosenParm)->rho = 0.0;
                    ((CExpParm*)chosenParm)->xi = 0.0;
                    ((CExpParm*)chosenParm)->phase = 0.0 ;
                    ((CExpParm*)chosenParm)->a = i;
                    ((CExpParm*)chosenParm)->b = i;

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

                if (dicData->getDicType()==1)
    			{
					((CExpParm*)parm)->rho      = 0.0;
					((CExpParm*)parm)->xi       = xi;
					((CExpParm*)parm)->phase    = 0.0;
					((CExpParm*)parm)->a        = 0;
					((CExpParm*)parm)->b        = N-1;
					((CExpDictionary*)dic)->setComplexAtom(parm);
				}                
                                    
                // Computing optimum phase
                opt_phase = computeOptimumPhase(    residue,
                                                    xi,
                                                    innerProd);

                if (fabs(innerProd)>fabs(maxInnerProd))
                {        
                    maxInnerProd = innerProd;
                    if (dicData->getDicType()==1)
    				{
						((CExpParm*)chosenParm)->innerProd = maxInnerProd;
						((CExpParm*)chosenParm)->rho = 0.0;
						((CExpParm*)chosenParm)->xi = xi;
						((CExpParm*)chosenParm)->phase = opt_phase ;
						((CExpParm*)chosenParm)->a = 0;
						((CExpParm*)chosenParm)->b = N-1;
					}
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

			if (dicData->getDicType()==1)
    		{
				((CExpParm*)parm)->rho      = 1.0/(double)s;
				((CExpParm*)parm)->xi       = 0.0;
				((CExpParm*)parm)->phase    = 0.0;
				((CExpParm*)parm)->a        = 0;
				((CExpParm*)parm)->b        = N-1;
				((CExpDictionary*)dic)->setRealAtom(parm);
			}	

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
                    if (dicData->getDicType()==1)
    				{
						((CExpParm*)parm)->rho      = 0.0;
						((CExpParm*)parm)->xi       = xi;
						((CExpParm*)parm)->phase    = 0.0;
						((CExpParm*)parm)->a        = 0;
						((CExpParm*)parm)->b        = N-1;
						((CExpDictionary*)dic)->setComplexAtom(parm);
					}
                    

                    for (i=0;i<N;i++)
                    {
                        z1[i][0] = m_complexAtom[i].Real() * residue[0][i];
                        z1[i][1] = m_complexAtom[i].Imag() * residue[0][i];
                    }

                    // z2
                    if (dicData->getDicType()==1)
    				{
						((CExpParm*)parm)->rho      = 0.0;
						((CExpParm*)parm)->xi       = 2*xi;
						((CExpParm*)parm)->phase    = 0.0;
						((CExpParm*)parm)->a        = 0;
						((CExpParm*)parm)->b        = N-1;
						((CExpDictionary*)dic)->setComplexAtom(parm);
					}

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
                            if (dicData->getDicType()==1)
    						{
								((CExpParm*)chosenParm)->innerProd = maxInnerProd;
								((CExpParm*)chosenParm)->rho = 1.0/(double)s;
								((CExpParm*)chosenParm)->xi = xi;
								((CExpParm*)chosenParm)->phase = opt_phase ;
								((CExpParm*)chosenParm)->a = tau;
								((CExpParm*)chosenParm)->b = N-1;
							}
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
                            if (dicData->getDicType()==1)
    						{
								((CExpParm*)chosenParm)->innerProd = maxInnerProd;
								((CExpParm*)chosenParm)->rho = -1.0/(double)s;
								((CExpParm*)chosenParm)->xi = xi;
								((CExpParm*)chosenParm)->phase = opt_phase ;
								((CExpParm*)chosenParm)->a = 0;
								((CExpParm*)chosenParm)->b = tau;
							}
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
    return chosenParm;
}