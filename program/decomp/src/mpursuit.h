#ifndef MPURSUIT_H
#define MPURSUIT_H

#include "filemgr.h"
#include "structbook.h"
#include "dictionary.h"
#include "datasignal.h"

CParameter* fastMPKolasaModified(	cgMatrix<double>& residue,
									CDataSignal* dataSignal,
									CFileDictionary* dicData);
									
#endif