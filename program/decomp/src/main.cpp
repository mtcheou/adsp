
#include "decomp.h"

using namespace std;

int main(int argc, char *argv[])
{

    // Define the ".cfg" file path
    // _MAX_PATH defined in stdlib.h
    char InputFile[_MAX_PATH];
    if (argc > 1)
    {
        strcpy( InputFile, argv[1]);
    }
    else
    {
    	strcpy( InputFile, "x001.cfg");
    }


	/////////////////////////////
    // Loading General Input File (for decompositon)
    CFileDecomp* genData;
    genData = new CFileDecomp;
    genData->setFileName("panelDecomp.dat");
    genData->loadData();
    genData->echo();
    
    // Loading block range input file
    CFileDecompBlockRange* blockRange;
    blockRange = new CFileDecompBlockRange;
    blockRange->setFileName("panelBlockRange.dat");
    blockRange->loadData();
    
    // Loading the dictionary file
    CFileDictionary* dicData;
    dicData = new CFileDictionary;
    dicData->setFileName("panelDictionary.dat");
    dicData->loadData();
    dicData->printToScreen();
    
    
    CDataSignal* dataSignal;
    if (genData->getSigType()==1)
    {
    	decompElectric(	genData,
		    			blockRange,
	    				dicData,
	    				InputFile);
    
    }
    if (genData->getSigType()==2)
    {
    
    }
    if (genData->getSigType()==3)
    {
    
    }
    
	
	///////
    delete genData;
    delete blockRange;
    delete dicData;
  	//////
    return 0;
}