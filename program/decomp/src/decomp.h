#ifndef DECOMP_H
#define DECOMP_H


#include "filemgr.h"
#include "structbook.h"
#include "dictionary.h"
#include "datasignal.h"

void decompElectric(CFileDecomp* genData,
		    		CFileDecompBlockRange* blockRange,
	    			CFileDictionary* dicData,
	    			char* InputFile);

void decompAudio(	CFileDecomp* genData,
		    		CFileDecompBlockRange* blockRange,
	    			CFileDictionary* dicData,
	    			char* InputFile);
	    			
void decompNoise(	CFileDecomp* genData,
		    		CFileDecompBlockRange* blockRange,
	    			CFileDictionary* dicData);	    				    			

#endif