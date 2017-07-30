/// \file TimeDelayAlgorithms.h
/// Header file for TimeDelayAlgorithms class

#pragma once

/// \class TimeDelayAlgorithms
/// User class for estimation and compensation of the time delay between
/// input and output files
class TimeDelayAlgorithms
{
public:
	TimeDelayAlgorithms(void);
	~TimeDelayAlgorithms(void);
	
	double TimeDelayEstimation();
	void TimeDelayAdjustment();
	void SetParameters(double fSampleFrequency, int iBlockSize, int iInterpOrder, int iInterpRate);
	void SetInputs(double* pdInputI, double* pdInputQ, double* pdOutputI, double* pdOutputQ, int iLength);
	void SetTimeDelay(double dTimeDelay);
	bool OpenFiles(CString csFilenameInI, CString csFilenameInQ, BOOL bInputIQ, CString csFilenameOutI, CString csFilenameOutQ, BOOL bOutputIQ);

private:	
	double* m_pdInputI;			///< Pointer to the double floating point input (I)
	double* m_pdInputQ;			///< Pointer to the double floating point input (Q)
	double* m_pdOutputI;		///< Pointer to the double floating point output (I)
	double* m_pdOutputQ;		///< Pointer to the double floating point output (Q)
	int m_iLength;				///< Length of arrays (all arrays must be equal)
	double m_dSampleFrequency;	///< Sampling Frequency
	int m_iBlockSize;			///< Interpolation Block Size
	int m_iInterpOrder;			///< Interpolation Order
	int m_iInterpRate;			///< Interpolation Rate
	double m_dTimeDelay;		///< Computed Time Delay
	bool m_bOpenedFiles;		///< Flag to tell if this class opened files and allocated pointers
	CString m_csInputFilename;	///< Internal input filename for Delay Adjustment
	CString m_csOutputFilename;	///< Internal output filename for Delay Adjustment

	int LagrangeInterpolation(double* pdYin, int iLength, int iGranul, int iOrder, double* pdXout, double* pdYout, int iComputeLength = 0);
	void datashift(double* pdInputData, double* pdOutputData, int iNumShift, int iBlockSize, int iLength, int* iWriteLength);
};
