/// \file TimeDelayAlgorithms.cpp
/// Implementation file for TimeDelayAlgorithms class

#include "stdafx.h"
#include "TimeDelayAlgorithms.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "nsmfunctions.h"

#include "resource.h"
#include ".\NtGraphDialog.h"
#include ".\PlotUtility.h"
#include ".\spline\spline.h"
#include ".\MATLAB\MatlabFunctions.h"

////////////////////////////////////////////////////////////////////////////////////
///
/// Default constructor for the TimeDelayAlgorithms class.
/// Initializes all pointers to NULL, and sets the length to 0
///
////////////////////////////////////////////////////////////////////////////////////
TimeDelayAlgorithms::TimeDelayAlgorithms(void)
{
	m_pdInputI = NULL;
	m_pdInputQ = NULL;
	m_pdOutputI = NULL;
	m_pdOutputQ = NULL;
	m_iLength = 0;
	m_dSampleFrequency = 0.0;
	m_iBlockSize = 0;
	m_iInterpRate = 0;
	m_iInterpOrder = 0;
	m_dTimeDelay = 0.0;
	m_bOpenedFiles = FALSE;
	
	// Just assign some default filenames.
	m_csInputFilename = "Input.txt";
	m_csOutputFilename = "Output.txt";
}


////////////////////////////////////////////////////////////////////////////////////
///
/// Default destructor for the TimeDelayAlgorithms class.
/// Deallocates memory if class has allocated memory.
///
////////////////////////////////////////////////////////////////////////////////////
TimeDelayAlgorithms::~TimeDelayAlgorithms(void)
{
	if(m_bOpenedFiles)
	{
		delete [] m_pdInputI;
		delete [] m_pdInputQ;
		delete [] m_pdOutputI;
		delete [] m_pdOutputQ;
	}
}


////////////////////////////////////////////////////////////////////////////////////
///
/// Computes the estimated time delay between input and output vectors.
/// 
/// Reference MATLAB files:
/// - <a href = "../../NonlinearSystemModeler/Delay/timedelayestimationFunc.m"> timedelayestimationFunc.m </a>
/// - <a href = "../../NonlinearSystemModeler/Delay/timedelayFunc.m"> timedelayFunc.m </a>
///
/// \return the time (in ns) of the estimated
///
////////////////////////////////////////////////////////////////////////////////////
double TimeDelayAlgorithms::TimeDelayEstimation(void)
{
	int iStartIndex;
	int iOriginalDataBlockx2;
	int iLagrangeLength;
	int iMaxLags;
	int ixcovLength;
	double* pdMagIn = NULL;
	double* pdMagOut = NULL;
	double* pdTimeIn = NULL;
	double* pdTimeInLagrange = NULL;
	double* pdTimeOutLagrange = NULL;
	double* pdMagInLagrange = NULL;
	double* pdMagOutLagrange = NULL;
	double* pdCxy = NULL;
	int* piLags = NULL;
	double* pdLags = NULL;
	int iMaxCxyIndex;
	int iMaxCxyLag;
	double dTimeDelay;
	int iIndex;
	CNtGraphDialog* pPlotMagIn;
	CNtGraphDialog* pPlotMagOut;
	CNtGraphDialog* pPlotLagrange;
	CNtGraphDialog* pPlotCorrelation;


	// If invalid parameters, exit
	if(m_dSampleFrequency == 0.0 || m_iBlockSize < 1 || m_iInterpOrder < 0 ||
		m_iInterpRate < 0)
	{
		AfxMessageBox(_T("Error: Invalid Parameters"));
		return 0.0;
	}

	// Started to run of memory, place here.
	iOriginalDataBlockx2 = 2 * (m_iBlockSize + m_iInterpOrder);
	iStartIndex = 1000;
	
	// If invalid data, exit
	if(m_iLength <= (iStartIndex + iOriginalDataBlockx2) || m_pdInputI == NULL ||
		m_pdInputQ == NULL || m_pdOutputI == NULL || m_pdOutputQ == NULL)
	{
		AfxMessageBox(_T("Error: Unset Input/Output Data"));
		return 0.0;
	}

	// Almost direct port from matlab code
	if(m_iLength > 46080)
	{
		/// \note in MATLAB, the function does something weird that should give an error.
	}


	// Discard the first 1000 data to remove effects of system

	pdMagIn = new double[iOriginalDataBlockx2];
	pdMagOut = new double[iOriginalDataBlockx2];
	
	// Convert the data into magnitdue
	for(iIndex = iStartIndex; iIndex < (iStartIndex + iOriginalDataBlockx2); iIndex++)
	{
		pdMagIn[iIndex - iStartIndex] = sqrt(m_pdInputI[iIndex] * m_pdInputI[iIndex] +
											 m_pdInputQ[iIndex] * m_pdInputQ[iIndex]);

		pdMagOut[iIndex - iStartIndex] = sqrt(m_pdOutputI[iIndex] * m_pdOutputI[iIndex] +
											  m_pdOutputQ[iIndex] * m_pdOutputQ[iIndex]);
	}

	pdTimeIn = new double[iOriginalDataBlockx2];
	for(iIndex = 0; iIndex < iOriginalDataBlockx2; iIndex++)
	{
		pdTimeIn[iIndex] = iIndex + 1;
	}

	// Make lagrange Interpolation
	iLagrangeLength = LagrangeInterpolation(pdMagIn, iOriginalDataBlockx2, m_iInterpRate, m_iInterpOrder, pdTimeOutLagrange, pdMagOutLagrange, 1);

	// ComputeLagrangeLength does not do it properly, matlab code requires this many space allocated
	// Fixed in ComputeLagrangeLength
	//iLagrangeLength = ((iLagrangeLength / 26)+1) * 26;

	pdMagInLagrange = new double[iLagrangeLength];
	pdMagOutLagrange = new double[iLagrangeLength];
	pdTimeInLagrange = new double[iLagrangeLength];
	pdTimeOutLagrange = new double[iLagrangeLength];

	LagrangeInterpolation(pdMagOut, iOriginalDataBlockx2, m_iInterpRate, m_iInterpOrder, pdTimeOutLagrange, pdMagOutLagrange);
	LagrangeInterpolation(pdMagIn, iOriginalDataBlockx2, m_iInterpRate, m_iInterpOrder, pdTimeInLagrange, pdMagInLagrange);

	
	// Compute Time
	for(iIndex = 0; iIndex < iOriginalDataBlockx2; iIndex++)
	{
		pdTimeIn[iIndex] = pdTimeIn[iIndex] * 1000 / m_dSampleFrequency;
	}

	for(iIndex = 0; iIndex < iLagrangeLength; iIndex++)
	{
		pdTimeInLagrange[iIndex] = pdTimeInLagrange[iIndex] * 1000 / m_dSampleFrequency;
		pdTimeOutLagrange[iIndex] = pdTimeOutLagrange[iIndex] * 1000 / m_dSampleFrequency;
	}


	// Plot the Interpolated Results
	CREATE_NTDIALOG_MODELESS(pPlotMagIn);
	pPlotMagIn->AddData(pdTimeIn, pdMagIn, iOriginalDataBlockx2, _T("Measurement"));
	pPlotMagIn->SetDataProperties(0, RGB(255, 0, 0), 1, CNtgraphctrl1::Solid, CNtgraphctrl1::None);
	pPlotMagIn->AddData(pdTimeInLagrange, pdMagInLagrange, iLagrangeLength, _T("Lagrange interpolation"));
	pPlotMagIn->SetDataProperties(1, RGB(0, 0, 255), 1, CNtgraphctrl1::Solid, CNtgraphctrl1::None);
	pPlotMagIn->SetXLabel(_T("Time (us)"));
	pPlotMagIn->SetYLabel(_T("Magnitude"));
	pPlotMagIn->SetTitle(_T("Input Waveform"));
	pPlotMagIn->MyAutoFixAxis();
	pPlotMagIn->ShowLegend(0);

	CREATE_NTDIALOG_MODELESS(pPlotMagOut);
	pPlotMagOut->AddData(pdTimeIn, pdMagOut, iOriginalDataBlockx2, _T("Measurement"));
	pPlotMagOut->SetDataProperties(0, RGB(255, 0, 0), 1, CNtgraphctrl1::Solid, CNtgraphctrl1::None);
	pPlotMagOut->AddData(pdTimeOutLagrange, pdMagOutLagrange, iLagrangeLength, _T("Lagrange interpolation"));
	pPlotMagOut->SetDataProperties(1, RGB(0, 0, 255), 1, CNtgraphctrl1::Solid, CNtgraphctrl1::None);
	pPlotMagOut->SetXLabel(_T("Time (us)"));
	pPlotMagOut->SetYLabel(_T("Magnitude"));
	pPlotMagOut->SetTitle(_T("Output Waveform"));
	pPlotMagOut->MyAutoFixAxis();
	pPlotMagOut->ShowLegend(0);

	CREATE_NTDIALOG_MODELESS(pPlotLagrange);
	pPlotLagrange->AddData(pdTimeInLagrange, pdMagInLagrange, iLagrangeLength, _T("Input"));
	pPlotLagrange->SetDataProperties(0, RGB(255, 0, 0), 1, CNtgraphctrl1::Solid, CNtgraphctrl1::None);
	pPlotLagrange->AddData(pdTimeOutLagrange, pdMagOutLagrange, iLagrangeLength, _T("Output"));
	pPlotLagrange->SetDataProperties(1, RGB(0, 0, 255), 1, CNtgraphctrl1::Solid, CNtgraphctrl1::None);
	pPlotLagrange->SetXLabel(_T("Time (us)"));
	pPlotLagrange->SetYLabel(_T("Magnitude"));
	pPlotLagrange->SetTitle(_T("Measured Waveforms before Delay Adjustment"));
	pPlotLagrange->MyAutoFixAxis();
	pPlotLagrange->ShowLegend(0);

	// Perform finer tuning
	// did not need floor, just copy from MATLAB code.
	iMaxLags = int(floor(double(iLagrangeLength) / 2.0));
	ixcovLength = MatlabXcov(pdMagOutLagrange, pdMagInLagrange, iLagrangeLength, iMaxLags, 'c', pdCxy, piLags, 1);

	pdCxy = new double[ixcovLength];
	piLags = new int[ixcovLength];
	
	// Call xcov equivalent
	MatlabXcov(pdMagOutLagrange, pdMagInLagrange, iLagrangeLength, iMaxLags, 'c', pdCxy, piLags);

	// Find the maximum of the cross correlation
	iMaxCxyIndex = 0;
	for(iIndex = 1; iIndex < ixcovLength; iIndex++)
	{
		if(pdCxy[iIndex] > pdCxy[iMaxCxyIndex])
			iMaxCxyIndex = iIndex;
	}

	// Reference back to the lags vector
	iMaxCxyLag = piLags[iMaxCxyIndex];

	// Compute Delay
	dTimeDelay = (double(iMaxCxyLag) / ((m_iInterpRate + 1) * m_dSampleFrequency)) * 1000.0;

	// Plot the correlation results
	pdLags = new double[ixcovLength];
	for(iIndex = 0; iIndex < ixcovLength; iIndex++)
	{
		pdLags[iIndex] = (double) piLags[iIndex];
	}

	CREATE_NTDIALOG_MODELESS(pPlotCorrelation);
	pPlotCorrelation->AddData(pdLags, pdCxy, ixcovLength, _T("Cross-Correlation"));
	pPlotCorrelation->SetDataProperties(0, RGB(255, 0, 0), 1, CNtgraphctrl1::Solid, CNtgraphctrl1::None);
	pPlotCorrelation->SetXLabel(_T("Lags"));
	pPlotCorrelation->SetYLabel(_T("Rxy"));
	pPlotCorrelation->SetTitle(_T("Input and Output Waveform Cross Correlation"));
	pPlotCorrelation->MyAutoFixAxis();
	pPlotCorrelation->ShowLegend(0);


	// Cleanup
	SAFE_DELETE_ARRAY(pdMagIn);
	SAFE_DELETE_ARRAY(pdMagOut);

	SAFE_DELETE_ARRAY(pdTimeInLagrange);
	SAFE_DELETE_ARRAY(pdTimeOutLagrange);
	SAFE_DELETE_ARRAY(pdMagInLagrange);
	SAFE_DELETE_ARRAY(pdMagOutLagrange);
	SAFE_DELETE_ARRAY(pdTimeIn);
	SAFE_DELETE_ARRAY(pdCxy);
	SAFE_DELETE_ARRAY(piLags);
	SAFE_DELETE_ARRAY(pdLags);

	return dTimeDelay;
}


////////////////////////////////////////////////////////////////////////////////////
///
/// Sets the parameters for Time Delay Estimation/Adjustment
///
/// \param dSampleFrequency the sampling frequency
/// \param iBlockSize block size for lagrange interpolation
/// \param iInterpOrder interpolation order
/// \param iInterpRate interpolation rate
///
////////////////////////////////////////////////////////////////////////////////////
void TimeDelayAlgorithms::SetParameters(double dSampleFrequency, int iBlockSize, int iInterpOrder, int iInterpRate)
{
	m_dSampleFrequency = dSampleFrequency;
	m_iBlockSize = iBlockSize;
	m_iInterpOrder = iInterpOrder;
	m_iInterpRate = iInterpRate;

	return;
}


////////////////////////////////////////////////////////////////////////////////////
///
/// Sets the input vectors and length
///
/// \param pdInputI pointer to the array containing inphase input
/// \param pdInputQ pointer to the array containing quadraturephase input
/// \param pdOutputI pointer to the array containing inphase output
/// \param pdOutputQ pointer to the array containing quadraturephase output
/// \param iLength length of the arrays (this length applies to the previous 4 inputs)
///
////////////////////////////////////////////////////////////////////////////////////
void TimeDelayAlgorithms::SetInputs(double* pdInputI, double* pdInputQ, double* pdOutputI, double* pdOutputQ, int iLength)
{
	if(m_bOpenedFiles == false)
	{
		m_pdInputI = pdInputI;
		m_pdInputQ = pdInputQ;
		m_pdOutputI = pdOutputI;
		m_pdOutputQ = pdOutputQ;
		m_iLength = iLength;
	}
	else
	{
		AfxMessageBox(_T("This class already has opened files and allocated pointers, cannot overwrite."));
	}

	return;
}


////////////////////////////////////////////////////////////////////////////////////
///
/// Does a lagrange interpolation on the input data, pdInput
///
/// \param pdYin pointer to the input data
/// \param iLength length of the input data
/// \param iGranul granularity (larger than 1)
/// \param iOrder order
/// \param [out] pdXout linear array output (preallocated pointer)
/// \param [out] pdYout interpolated array output (preallocated pointer)
/// \param iComputeLength if != 0, returns the number of elements to allocate for Lagrange
///
/// \return the number of elements to allocate for Lagrange Interpolation if iComputeLength != 0
///
////////////////////////////////////////////////////////////////////////////////////
int TimeDelayAlgorithms::LagrangeInterpolation(double* pdYin, int iLength, int iGranul, int iOrder, double* pdXout, double* pdYout, int iComputeLength)
{
	int* piXin = NULL;
	double dXinStep;
	double dXoutStep;
	int iM;
	int iTotalLength;
	int* piDenominator = NULL;
	double* pdNominator = NULL;
	int iCounter;
	int iXTemp;
	int iIndex;

	iM = (iGranul + 1) * (iLength - 1) + 1;	
	
	// if iComputeLength != 0, return the length needed for lagrange
	if(iComputeLength != 0)
	{
		iTotalLength = iM - (3 * iGranul * (iOrder + 1)) - 1;

		// Need to round up to nearest granul * order + 1
		iTotalLength = ((iTotalLength / (iGranul + 1)) + 1) * (iGranul + 1);

		return iTotalLength;
	}


	piXin = new int[iLength];

	for(iIndex = 0; iIndex < iLength; iIndex++)
	{
		piXin[iIndex] = iIndex + 1;
	}

	dXinStep = piXin[1] - piXin[0];
	dXoutStep = dXinStep / (iGranul + 1);

	// Points after interpolation


	// Calculate denominators of lagrange polynomial
	piDenominator = new int[iOrder + 1];
	pdNominator = new double[iOrder + 1];
	for(iIndex = 0; iIndex < (iOrder + 1); iIndex++)
	{
		piDenominator[iIndex] = CrossDifferenceMultiplier(iIndex, piXin, (iOrder + 1));

	}
	
	// Make Lagrange Interpolation
	iCounter = 0;
	iTotalLength = iM - (3 * iGranul * (iOrder + 1)) - 1;

	for(iIndex = 0; iIndex < iTotalLength; iIndex += (iGranul + 1))
	{
		for(int ikIndex = 0; ikIndex < (iGranul + 1); ikIndex++)
		{
			pdXout[iIndex + ikIndex] = piXin[iCounter] + (ikIndex) * dXoutStep;
			pdYout[iIndex + ikIndex] = 0.0;

			for(int ijIndex = 0; ijIndex < (iOrder + 1); ijIndex++)
			{
				iXTemp = piXin[iCounter] + (ijIndex);
				pdNominator[ijIndex] = 0;
				
				// Calculate nominators of Lagrange polynomial
				pdNominator[ijIndex] = CoCrossDifferenceMultiplier(iXTemp, pdXout[iIndex + ikIndex], &piXin[iCounter], (iOrder+1));

				if(piDenominator[ijIndex] != 0)
				{
					pdYout[iIndex + ikIndex] = pdYout[iIndex + ikIndex] + (pdNominator[ijIndex] * pdYin[iCounter + ijIndex] / double(piDenominator[ijIndex]));
				}
			}
		}

		iCounter = iCounter + 1;
	}

	SAFE_DELETE_ARRAY(piXin);
	SAFE_DELETE_ARRAY(piDenominator);
	SAFE_DELETE_ARRAY(pdNominator);

	return 0;
}


////////////////////////////////////////////////////////////////////////////////////
///
/// Sets the time delay
///
/// \param dTimeDelay the time delay (specified by user or computed from
///					  TimeDelayEstimation() algorithm)
///
////////////////////////////////////////////////////////////////////////////////////
void TimeDelayAlgorithms::SetTimeDelay(double dTimeDelay)
{
	m_dTimeDelay = dTimeDelay;

	return;
}


////////////////////////////////////////////////////////////////////////////////////
///
/// Applies the time delay to the specified files.
/// This function will apply time delay adjustment to the input/output vectors, and
/// save them into a new file defined by \ref DELAY_ADJUSTMENT_FILE_SUFFIX.
/// \bug In the MATLAB file, I_in(1:NN) does nothing if NN < 0.
///
/// Reference MATLAB files:
/// - <a href = "../../NonlinearSystemModeler/Delay/timedelayadjustFunc.m"> timedelayadjustFunc.m </a>
///
////////////////////////////////////////////////////////////////////////////////////
void TimeDelayAlgorithms::TimeDelayAdjustment()
{
	int iCoarseIndex;
	int iFineIndex;
	double dFineDelay;
	double dTimeResolution;
	int iBlocksize;
	int iTotalLength;
	int iFinalLength;
	int iIndex;
	double* pdInputIShifted = NULL;
	double* pdInputQShifted = NULL;
	double* pdOutputIShifted = NULL;
	double* pdOutputQShifted = NULL;

	double* pdInputICoarse = NULL;
	double* pdInputQCoarse = NULL;
	double* pdOutputICoarse = NULL;
	double* pdOutputQCoarse = NULL;

	CStdioFile cfOutput;
	CStdioFile cfInput;
	CString csWriteFile;

	// Perform coarse adjustment
	iCoarseIndex = (int) floor(m_dTimeDelay * m_dSampleFrequency / 1000.0);
	iTotalLength = m_iLength - abs(iCoarseIndex);

	pdInputICoarse = new double[m_iLength];
	pdInputQCoarse = new double[m_iLength];
	pdOutputICoarse = new double[m_iLength];
	pdOutputQCoarse = new double[m_iLength];

	// Negative stands for the output delay iCoarseIndex samples
	if(iCoarseIndex < 0)
	{
		// Rotate Input
		for(iIndex = -iCoarseIndex; iIndex < m_iLength; iIndex++)
		{
			pdInputICoarse[iIndex + iCoarseIndex] = m_pdInputI[iIndex];
			pdInputQCoarse[iIndex + iCoarseIndex] = m_pdInputQ[iIndex];
		}
		for(iIndex = 0; iIndex < m_iLength; iIndex++)
		{
			pdOutputICoarse[iIndex] = m_pdOutputI[iIndex];
			pdOutputQCoarse[iIndex] = m_pdOutputQ[iIndex];
		}
	}
	// Positive stands for the output precedes iCoarseIndex samples
	else if(iCoarseIndex > 0)
	{
		// Rotate Output
		for(iIndex = iCoarseIndex; iIndex < m_iLength; iIndex++)
		{
			pdOutputICoarse[iIndex - iCoarseIndex] = m_pdOutputI[iIndex];
			pdOutputQCoarse[iIndex - iCoarseIndex] = m_pdOutputQ[iIndex];
		}
		for(iIndex = 0; iIndex < m_iLength; iIndex++)
		{
			pdInputICoarse[iIndex] = m_pdInputI[iIndex];
			pdInputQCoarse[iIndex] = m_pdInputQ[iIndex];
		}
	}
	// Else, do nothing because they are aligned for coarse adjustment.
	else
	{
		for(iIndex = 0; iIndex < m_iLength; iIndex++)
		{
			pdInputICoarse[iIndex] = m_pdInputI[iIndex];
			pdInputQCoarse[iIndex] = m_pdInputQ[iIndex];
			pdOutputICoarse[iIndex] = m_pdOutputI[iIndex];
			pdOutputQCoarse[iIndex] = m_pdOutputQ[iIndex];
		}
	}


	// Perform fine adjustment
	// FineDelay is smaller than the normal sample
	dFineDelay = m_dTimeDelay - (double(iCoarseIndex) * 1000.0 / m_dSampleFrequency);

	// Should this be +1???
	// Oct 5, 2007. It is in fact +1, otherwise, delay adjustments are off.
	dTimeResolution = 1000.0 / m_dSampleFrequency / (m_iInterpRate + 1);

	iFineIndex = round(dFineDelay / dTimeResolution);

	iBlocksize = DELAY_ADJUSTMENT_BLOCKSIZE;


	// Preinitialize iFinalLength
	pdInputIShifted = new double[iTotalLength];
	pdInputQShifted = new double[iTotalLength];
	pdOutputIShifted = new double[iTotalLength];
	pdOutputQShifted = new double[iTotalLength];
	iFinalLength = iTotalLength;

	// Output is delayed
	if(iFineIndex > 0)
	{
		datashift(pdOutputICoarse, pdOutputIShifted, iFineIndex, iBlocksize, iTotalLength, &iFinalLength);
		datashift(pdOutputQCoarse, pdOutputQShifted, iFineIndex, iBlocksize, iTotalLength, &iFinalLength);

		// Copy Input I/Q (just truncate the length basically)
		for(iIndex = 0; iIndex < iFinalLength; iIndex++)
		{
			pdInputIShifted[iIndex] = pdInputICoarse[iIndex];
			pdInputQShifted[iIndex] = pdInputQCoarse[iIndex];
		}
	}
	// Will NEVER fall into this block...due to floor() function for coarse index.
	else if(iFineIndex < 0)
	{
		datashift(pdInputICoarse, pdInputIShifted, -iFineIndex, iBlocksize, iTotalLength, &iFinalLength);
		datashift(pdInputQCoarse, pdInputQShifted, -iFineIndex, iBlocksize, iTotalLength, &iFinalLength);

		// Copy Input I/Q (just truncate the length basically)
		for(iIndex = 0; iIndex < iFinalLength; iIndex++)
		{
			pdOutputIShifted[iIndex] = pdOutputICoarse[iIndex];
			pdOutputQShifted[iIndex] = pdOutputQCoarse[iIndex];
		}		
	}
	// Do nothing...
	else
	{
		// Already Aligned, copy.
		for(iIndex = 0; iIndex < iFinalLength; iIndex++)
		{
			pdInputIShifted[iIndex] = pdInputICoarse[iIndex];
			pdInputQShifted[iIndex] = pdInputQCoarse[iIndex];
			pdOutputIShifted[iIndex] = pdOutputICoarse[iIndex];
			pdOutputQShifted[iIndex] = pdOutputQCoarse[iIndex];
		}		

	}

	// Write to text file, with modifications to filename
	m_csInputFilename.Replace(_T(".txt"), _T(DELAY_ADJUSTMENT_FILE_SUFFIX));
	m_csOutputFilename.Replace(_T(".txt"), _T(DELAY_ADJUSTMENT_FILE_SUFFIX));

	// Open files for writing
	cfInput.Open(m_csInputFilename, CFile::modeCreate | CFile::modeWrite);
	cfOutput.Open(m_csOutputFilename, CFile::modeCreate | CFile::modeWrite);

	// Iterate over iFinalLength times
	for(iIndex = 0; iIndex < iFinalLength; iIndex++)
	{
		/// \note should be save2col function?
		csWriteFile.Format(_T("%.16f \t %.16f"), pdInputIShifted[iIndex], pdInputQShifted[iIndex]);
		csWriteFile = csWriteFile + "\n";
		cfInput.WriteString(csWriteFile);

		csWriteFile.Format(_T("%.16f \t %.16f"), pdOutputIShifted[iIndex], pdOutputQShifted[iIndex]);
		csWriteFile = csWriteFile + "\n";
		cfOutput.WriteString(csWriteFile);
	}

	// Close files.
	cfInput.Close();
	cfOutput.Close();


	// Plot Delay Adjustment validation
	/// \note refactor this code later.
	{
		// Compare the waveform
		double* pdTimeVector = NULL;
		double* pdInput = NULL;
		double* pdOutput = NULL;
		double* pdInputAligned = NULL;
		double* pdOutputAligned = NULL;
		doublecomplex* pcdDataIn = NULL;
		doublecomplex* pcdDataDes = NULL;
		doublecomplex* pcdDataSim = NULL;

		// These are hardcoded as 80:120 in matlab.
		int iMinTime = 79;
		int iMaxTime = 120;
		int iTimeLength = iMaxTime - iMinTime;

		CNtGraphDialog* pInputvsAlignedOutput;
		CNtGraphDialog* pOutputvsAlignedOutput;
		CNtGraphDialog* pAlignedInputvsAlignedOutput;


		SPECTRUM_PARAMS sSpectrumParams;

		// Create complex double vectors
		pcdDataIn = new doublecomplex[iFinalLength];
		pcdDataDes = new doublecomplex[iFinalLength];
		pcdDataSim = new doublecomplex[iFinalLength];

		// Copy
		for(iIndex = 0; iIndex < iFinalLength; iIndex++)
		{
			pcdDataIn[iIndex].r = m_pdInputI[iIndex];
			pcdDataIn[iIndex].i = m_pdInputQ[iIndex];
			pcdDataDes[iIndex].r = m_pdOutputI[iIndex];
			pcdDataDes[iIndex].i = m_pdOutputQ[iIndex];
			pcdDataSim[iIndex].r = pdOutputIShifted[iIndex];
			pcdDataSim[iIndex].i = pdOutputQShifted[iIndex];
		}

		pdTimeVector = new double[iTimeLength];
		pdInput = new double[iTimeLength];
		pdOutput = new double[iTimeLength];
		pdInputAligned = new double[iTimeLength];
		pdOutputAligned = new double[iTimeLength];
		
		// Create Time Vector
		for(iIndex = iMinTime; iIndex < iMaxTime; iIndex++)
		{
			pdTimeVector[iIndex - iMinTime] = iIndex + 1;
			//pdInput[iIndex - iMinTime] = sqrt(pow(m_pdInputI[iIndex], 2.0) + pow(m_pdInputQ[iIndex], 2.0));
			//pdOutput[iIndex - iMinTime] = sqrt(pow(m_pdOutputI[iIndex], 2.0) + pow(m_pdOutputQ[iIndex], 2.0));
			// Input Shifted Absolute
			pdInputAligned[iIndex - iMinTime] = sqrt(pow(pdInputIShifted[iIndex], 2.0) + pow(pdInputQShifted[iIndex], 2.0));
		}
		
		// Absolute of original input vector
		MatlabAbs(&pcdDataIn[iMinTime], pdInput, iTimeLength);
		// Absolute of original output vector
		MatlabAbs(&pcdDataDes[iMinTime], pdOutput, iTimeLength);
		// Absolute of shifted output
		MatlabAbs(&pcdDataSim[iMinTime], pdOutputAligned, iTimeLength);

		// Plot Input vs Aligned Output
		CREATE_NTDIALOG_MODELESS(pInputvsAlignedOutput);
		pInputvsAlignedOutput->AddData(pdTimeVector, pdInput, iTimeLength, _T("Input"));
		pInputvsAlignedOutput->SetDataProperties(0, RGB(0, 255, 0), 1, CNtgraphctrl1::Solid, CNtgraphctrl1::DownTriangles);
		pInputvsAlignedOutput->AddData(pdTimeVector, pdInputAligned, iTimeLength, _T("Aligned Input"));
		pInputvsAlignedOutput->SetDataProperties(1, RGB(0, 0, 255), 1, CNtgraphctrl1::Solid, CNtgraphctrl1::UpTriangles);
		pInputvsAlignedOutput->SetTitle(_T("Input Waveforms"));
		pInputvsAlignedOutput->MyAutoFixAxis();
		pInputvsAlignedOutput->ShowLegend(0);

		// Plot Original Output vs Aligned Output
		CREATE_NTDIALOG_MODELESS(pOutputvsAlignedOutput);
		pOutputvsAlignedOutput->AddData(pdTimeVector, pdOutput, iTimeLength, _T("Output"));
		pOutputvsAlignedOutput->SetDataProperties(0, RGB(0, 255, 0), 1, CNtgraphctrl1::Solid, CNtgraphctrl1::DownTriangles);
		pOutputvsAlignedOutput->AddData(pdTimeVector, pdOutputAligned, iTimeLength, _T("Aligned Output"));
		pOutputvsAlignedOutput->SetDataProperties(1, RGB(0, 0, 255), 1, CNtgraphctrl1::Solid, CNtgraphctrl1::UpTriangles);
		pOutputvsAlignedOutput->SetTitle(_T("Output Waveforms"));
		pOutputvsAlignedOutput->MyAutoFixAxis();
		pOutputvsAlignedOutput->ShowLegend(0);

		// Plot Aligned Input vs Aligned Output
		CREATE_NTDIALOG_MODELESS(pAlignedInputvsAlignedOutput);
		pAlignedInputvsAlignedOutput->AddData(pdTimeVector, pdInputAligned, iTimeLength, _T("Aligned Input"));
		pAlignedInputvsAlignedOutput->SetDataProperties(0, RGB(0, 255, 0), 1, CNtgraphctrl1::Solid, CNtgraphctrl1::DownTriangles);
		pAlignedInputvsAlignedOutput->AddData(pdTimeVector, pdOutputAligned, iTimeLength, _T("Aligned Output"));
		pAlignedInputvsAlignedOutput->SetDataProperties(1, RGB(0, 0, 255), 1, CNtgraphctrl1::Solid, CNtgraphctrl1::UpTriangles);
		pAlignedInputvsAlignedOutput->SetTitle(_T("Aligned Waveforms"));
		pAlignedInputvsAlignedOutput->MyAutoFixAxis();
		pAlignedInputvsAlignedOutput->ShowLegend(0);

		// Compare Spectrum
		sSpectrumParams.m_dCenterFreq = 0;
		sSpectrumParams.m_dSpan = m_dSampleFrequency / 1.28;
		sSpectrumParams.m_dFstart = sSpectrumParams.m_dCenterFreq - sSpectrumParams.m_dSpan / 2;
		sSpectrumParams.m_dFstop = sSpectrumParams.m_dCenterFreq + sSpectrumParams.m_dSpan / 2;
		sSpectrumParams.m_dSampleRate = m_dSampleFrequency;
		sSpectrumParams.m_eWindow = WINDOW_HAMMING;
		
		// Mar 30/08
		// Comment out plotting of spectrums
		/*
		PlotUtility pPlot;
		pPlot.SetInputs(pcdDataIn, pcdDataDes, pcdDataSim, iFinalLength);
		pPlot.SetSpectrumParameters(sSpectrumParams);
		pPlot.PlotInSpectralDomain();
		*/

		SAFE_DELETE_ARRAY(pcdDataIn);
		SAFE_DELETE_ARRAY(pcdDataDes);
		SAFE_DELETE_ARRAY(pcdDataSim);
		SAFE_DELETE_ARRAY(pdTimeVector);
		SAFE_DELETE_ARRAY(pdInput);
		SAFE_DELETE_ARRAY(pdOutput);
		SAFE_DELETE_ARRAY(pdInputAligned);
		SAFE_DELETE_ARRAY(pdOutputAligned);
	}

	SAFE_DELETE_ARRAY(pdInputICoarse);
	SAFE_DELETE_ARRAY(pdInputQCoarse);
	SAFE_DELETE_ARRAY(pdOutputICoarse);
	SAFE_DELETE_ARRAY(pdOutputQCoarse);

	SAFE_DELETE_ARRAY(pdInputIShifted);
	SAFE_DELETE_ARRAY(pdInputQShifted);
	SAFE_DELETE_ARRAY(pdOutputIShifted);
	SAFE_DELETE_ARRAY(pdOutputQShifted);

	return;
}


////////////////////////////////////////////////////////////////////////////////////
///
/// Shifts the data using a spline interpolation method.
/// This function does not exactly replicate the piecewise cubic spline interpolation
/// as used by MATLAB's interp1 function. Instead, it uses the one found at 
/// http://people.scs.fsu.edu/~burkardt/cpp_src/spline/spline.html .
///
/// \param pdInputData data to be interpolated and shifted
/// \param [out] pdOutputData shifted data (preallocated)
/// \param iNumShift number of interpolated points to shift
/// \param iBlockSize size of block to use
/// \param iLength total length to write
/// \param [out] iWriteLength total written length
///
/// Reference MATLAB files:
/// - <a href = "../../NonlinearSystemModeler/Delay/datashift.m"> datashift.m </a>
/// - <a href = "../../NonlinearSystemModeler/Delay/blockdatashift.m"> blockdatashift.m </a>
///
////////////////////////////////////////////////////////////////////////////////////
void TimeDelayAlgorithms::datashift(double* pdInputData, double* pdOutputData, int iNumShift, int iBlockSize, int iLength, int* iWriteLength)
{
	int iIndex;//, iOutputIndex;//, iShiftedIndex;
	int iBlockIndex;
	int iInterpLength;

	double* pdInterpolatedBlock = NULL;
	double* pdInterpolatedBlockTime = NULL;
	double* pdTime = NULL;

	int ibcbeg;
	int ibcend;
	double ybcbeg;
	double ybcend;
	double *ypp;
	double yppval = 0;
	double ypval = 0;


	// Allocate Memory
	iInterpLength = (iBlockSize - 1) * (m_iInterpRate + 1);
	pdInterpolatedBlock = new double[iInterpLength];
	pdInterpolatedBlockTime = new double[iInterpLength];
	pdTime = new double[iBlockSize];

	// Initialize Time and Interpolated Block Time...
	for(iIndex = 0; iIndex < iInterpLength; iIndex++)
	{
		pdInterpolatedBlockTime[iIndex] = iIndex + 1;
	}

	for(iIndex = 0; iIndex < (iBlockSize); iIndex++)
	{
		pdTime[iIndex] = (iIndex * (m_iInterpRate)) + 1;
	}

	ibcbeg = 0;
	ybcbeg = 0.0;
	ibcend = 0;
	ybcend = 0.0;


	// Using a block size of BlockSize-1, as per MATLAB
	for(iIndex = 0; iIndex < iLength - iBlockSize; iIndex = (iIndex + iBlockSize - 1))
	{
		// Constant time index over all...

		// Setup Spline function
		ypp = spline_cubic_set(iBlockSize, pdTime, &pdInputData[iIndex], ibcbeg, ybcbeg, ibcend, ybcend);

		// Have to execute over ALL iInterpLength
//		for(iInterpIndex = 0; iInterpIndex < iInterpLength; iInterpIndex++)
//		{
//			pdInterpolatedBlock[iInterpIndex] = spline_cubic_val(iBlockSize-1, pdTime, pdInterpolatedBlockTime[iInterpIndex], &pdInputData[iIndex], ypp, &ypval, &yppval);
//		}

		// This is offset by 1, since our Initialized Time/Interpolated Time starts at 1.
		double dTvalIndex = (double) iNumShift + 1;
		
		// Modification, since we know the correct offsets, only call for the offsets and write directly into block
		// Note that in the MATLAB function, the blockdatashift function will return iBlockSize - 1 samples
		for(iBlockIndex = iIndex; iBlockIndex < (iIndex + iBlockSize); iBlockIndex++)
		{
			if(iNumShift > 0)
			{
				// Fixed Sept 22, condition was iBlockSize - 1.
				if(iBlockIndex == (iBlockSize + iIndex - 1))
				{
					break;
				}
				else
				{
					pdOutputData[iBlockIndex] = spline_cubic_val(iBlockSize, pdTime, dTvalIndex, &pdInputData[iIndex], ypp, &ypval, &yppval);
				}
			}
			else
			{
				if(iBlockIndex == 0)
				{
					// Do nothing
				}
				else
				{
					pdOutputData[iBlockIndex - 1] = spline_cubic_val(iBlockSize, pdTime, dTvalIndex, &pdInputData[iIndex], ypp, &ypval, &yppval);
				}
			}

			dTvalIndex = dTvalIndex + ((m_iInterpRate));
		}

		// Delete temporary ypp
		delete [] ypp;
	}

	*iWriteLength = iIndex;

	// Free
	SAFE_DELETE_ARRAY(pdTime);
	SAFE_DELETE_ARRAY(pdInterpolatedBlockTime);
	SAFE_DELETE_ARRAY(pdInterpolatedBlock);


// This was with lagrange, does not work because the lagrange function does not interpolate over all points.
#if 0
	iLagrangeLength = ComputeLagrangeLength(iBlockSize, m_iInterpRate, m_iInterpOrder);

	pdBlockData = (double*) malloc(sizeof(double) * iLagrangeLength);
	pdBlockDataTimeIndex = (double*) malloc(sizeof(double) * iLagrangeLength);

	if(iNumShift > 0)
	{
		for(iIndex = 0; iIndex < iLength - iBlockSize; iIndex = (iIndex + iBlockSize - 1))
		{
			// Use lagrange interpolation??
			LagrangeInterpolation(&pdInputData[iIndex], iBlockSize - 1, m_iInterpRate, m_iInterpOrder, pdBlockDataTimeIndex, pdBlockData);

			// Copy Back to pdOutputData (downsample)
			for(iOutputIndex = iIndex; iOutputIndex < iIndex + iBlockSize - 1; iOutputIndex++)
			{
				pdOutputData[iOutputIndex] = pdBlockData[iNumShift + ((iOutputIndex-iIndex) * (m_iInterpRate))];
			}
		}
	}
	else if(iNumShift < 0)
	{
		AfxMessageBox(_T("Shift number should be positive"));
	}

	free(pdBlockData);
	free(pdBlockDataTimeIndex);
#endif


	return;
}


////////////////////////////////////////////////////////////////////////////////////
///
/// Validates input files exist and assigns member pointers
///
/// \param csFilenameInI input filename string for input I vector (or input I/Q vector)
/// \param csFilenameInQ input filename string for input Q vector (optional)
/// \param bInputIQ flag to tell if input I/Q vectors are in the same file
/// \param csFilenameOutI output filename string for output I vector (or output I/Q vector)
/// \param csFilenameOutQ output filename string for output Q vector (optional)
/// \param bOutputIQ flag to tell if output I/Q vectors are in the same file
///
/// \return true if opening of files are successful and all pointers are assigned,
///			false otherwise.
///
////////////////////////////////////////////////////////////////////////////////////
bool TimeDelayAlgorithms::OpenFiles(CString csFilenameInI, CString csFilenameInQ, BOOL bInputIQ, CString csFilenameOutI, CString csFilenameOutQ, BOOL bOutputIQ)
{
	CFile cfTest;
	int iNumValuesInI, iNumValuesInQ;
	int iNumValuesOutI, iNumValuesOutQ;
	double* pdInputI = NULL;
	double* pdInputQ = NULL;
	double* pdOutputI = NULL;
	double* pdOutputQ = NULL;

	// Validate Files
	// Input I
	if(!cfTest.Open(csFilenameInI, CFile::readOnly))
	{
		AfxMessageBox(_T("Invalid Name for Input I file specified"));
		return false;
	}
	cfTest.Close();

	// If I/Q are in separate files, check Q
	if(bInputIQ == BST_UNCHECKED)
	{
		// Input Q
		if(!cfTest.Open(csFilenameInQ, CFile::readOnly))
		{
			AfxMessageBox(_T("Invalid Name for Input Q file specified"));
			return false;
		}
		cfTest.Close();
	}

	// Output I
	if(!cfTest.Open(csFilenameOutI, CFile::readOnly))
	{
		AfxMessageBox(_T("Invalid Name for Output I file specified"));
		return false;
	}
	cfTest.Close();

	// If I/Q are in separate files, check Q
	if(bOutputIQ == BST_UNCHECKED)
	{
		// Output Q
		if(!cfTest.Open(csFilenameOutQ, CFile::readOnly))
		{
			AfxMessageBox(_T("Invalid Name for Output Q file specified"));
			return false;
		}
		cfTest.Close();
	}	

	// Read in the data using the prebuilt function
	// Read Input Values
	if(bInputIQ  == BST_UNCHECKED)
	{
		iNumValuesInI = GetNumDataInTxtFile(FALSE, csFilenameInI);
		iNumValuesInQ = GetNumDataInTxtFile(FALSE, csFilenameInQ);

		if(iNumValuesInI != iNumValuesInQ)
		{
			AfxMessageBox(_T("Warning: Input I/Q Files do not have same lengths, taking smallest of two"));
			if(iNumValuesInI < iNumValuesInQ)
				iNumValuesInQ = iNumValuesInI;
			else
				iNumValuesInI = iNumValuesInQ;
		}

		// Cap values, running out of memory
		if(iNumValuesInI > DELAY_MAX_MEMORY_LENGTH)
		{
			iNumValuesInI = DELAY_MAX_MEMORY_LENGTH;
			iNumValuesInQ = DELAY_MAX_MEMORY_LENGTH;
		}


		pdInputI = new double[iNumValuesInI];
		pdInputQ = new double[iNumValuesInQ];

		ReadDataFromTxtFile(pdInputI, iNumValuesInI, csFilenameInI);
		ReadDataFromTxtFile(pdInputQ, iNumValuesInQ, csFilenameInQ);
	}
	else
	{
		iNumValuesInI = GetNumDataInTxtFile(TRUE, csFilenameInI);

		// Cap values, running out of memory
		if(iNumValuesInI > DELAY_MAX_MEMORY_LENGTH)
		{
			iNumValuesInI = DELAY_MAX_MEMORY_LENGTH;
			iNumValuesInQ = DELAY_MAX_MEMORY_LENGTH;
		}

		pdInputI = new double[iNumValuesInI];
		pdInputQ = new double[iNumValuesInI];

		ReadDataFromTxtFile(pdInputI, pdInputQ, iNumValuesInI, csFilenameInI);
	}
	
	// Read Output Values
	if(bOutputIQ  == BST_UNCHECKED)
	{
		iNumValuesOutI = GetNumDataInTxtFile(FALSE, csFilenameOutI);
		iNumValuesOutQ = GetNumDataInTxtFile(FALSE, csFilenameOutQ);

		if(iNumValuesOutI != iNumValuesOutQ)
		{
			AfxMessageBox(_T("Warning: Output I/Q Files do not have same lengths, taking smallest of two"));
			if(iNumValuesOutI < iNumValuesOutQ)
				iNumValuesOutQ = iNumValuesOutI;
			else
				iNumValuesOutI = iNumValuesOutQ;
		}

		// Cap values, running out of memory
		if(iNumValuesOutI > DELAY_MAX_MEMORY_LENGTH)
		{
			iNumValuesOutI = DELAY_MAX_MEMORY_LENGTH;
			iNumValuesOutQ = DELAY_MAX_MEMORY_LENGTH;
		}

		pdOutputI = new double[iNumValuesOutI];
		pdOutputQ = new double[iNumValuesOutQ];

		ReadDataFromTxtFile(pdOutputI, iNumValuesOutI, csFilenameOutI);
		ReadDataFromTxtFile(pdOutputQ, iNumValuesOutQ, csFilenameOutQ);
	}
	else
	{
		iNumValuesOutI = GetNumDataInTxtFile(TRUE, csFilenameOutI);

		// Cap values, running out of memory
		if(iNumValuesOutI > DELAY_MAX_MEMORY_LENGTH)
		{
			iNumValuesOutI = DELAY_MAX_MEMORY_LENGTH;
			iNumValuesOutQ = DELAY_MAX_MEMORY_LENGTH;
		}

		pdOutputI = new double[iNumValuesOutI];
		pdOutputQ = new double[iNumValuesOutI];

		ReadDataFromTxtFile(pdOutputI, pdOutputQ, iNumValuesOutI, csFilenameOutI);
	}
	
	// If Lengths of Input and Output do not match, take the lowest length and display a warning.
	if(iNumValuesInI != iNumValuesOutI)
	{
		AfxMessageBox(_T("Warning: Input and Output Data do not match length"));
		if(iNumValuesInI > iNumValuesOutI)
			iNumValuesInI = iNumValuesOutI;
	}


	// Assign to internal member variables
	m_pdInputI = pdInputI;
	m_pdInputQ = pdInputQ;
	m_pdOutputI = pdOutputI;
	m_pdOutputQ = pdOutputQ;
	m_iLength = iNumValuesInI;

	m_csInputFilename = csFilenameInI;
	m_csOutputFilename = csFilenameOutI;

	// Set flag
	m_bOpenedFiles = true;

	return true;
}