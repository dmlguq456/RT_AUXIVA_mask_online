#include <stdio.h>
#include "AUXIVA_mask_Online.h"
#include "header.h"
#include "sigproc.h"
#include <stdlib.h>

AUXIVA_MASK::AUXIVA_MASK()
{
	int i, j, k, freq, ch;
	int re, im;

	nfft = nWin;
	nshift = nWin / 4;
	nol = 3 * nWin / 4;
	nfreq = nfft / 2 + 1;
	//epsi = 0.000001;
	epsi = 2.220446049250313 * 1E-16;
	eps_floor = 1e-3;
	f_alpha = 0.98;
	gamma_t = 0.3;
	gamma_n = 0.9;
	gamma = new double[Nch];
	eta = new double *[Nch];
	for (i = 0; i < Nch; i++)
	{
		eta[i] = new double [nfreq];
	}
	for (i = 0; i < Nch; i++)
	{
		if (i == 0)
		{
			gamma[i] = gamma_t;
		}
		else
		{
			gamma[i] = gamma_n;
		}
	}
	Nrank = 2;
	double max = 32767;

	X = new double* [Nch]; // Nch X Nfreq(Complex)
	for (i = 0; i < Nch; i++)
	{
		X[i] = new double[nfreq * 2];
	}
	X_r = new double* [Nch]; // Nch X Nfreq(Complex)
	for (i = 0; i < Nch; i++)
	{
		X_r[i] = new double[nfreq * 2];
	}
	Y = new double* [Nch]; // Nch X Nfreq(Complex)
	for (i = 0; i < Nch; i++)
	{
		Y[i] = new double[nfreq * 2];
	}
	Pwr = new double* [Nch]; // Nch X Nfreq
	for (i = 0; i < Nch; i++)
	{
		Pwr[i] = new double[nfreq];
	}
	D = new double[nfreq]; // Nfreq

	V_nmf = new double* [Nch]; // Nch X Nrank
	for (i = 0; i < Nch; i++)
	{
		V_nmf[i] = new double[Nrank];
	}
	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nrank; j++)
		{
			V_nmf[i][j] = (rand() / max) + epsi;
		}
	}
	T_nmf = new double** [Nch]; // Nch X Nrank X Nfreq
	for (i = 0; i < Nch; i++)
	{
		T_nmf[i] = new double* [Nrank];
		for (j = 0; j < Nrank; j++)
		{
			T_nmf[i][j] = new double[nfreq];
		}
	}
	A_T_nmf = new double** [Nch]; // Nch X Nrank X Nfreq
	for (i = 0; i < Nch; i++)
	{
		A_T_nmf[i] = new double* [Nrank];
		for (j = 0; j < Nrank; j++)
		{
			A_T_nmf[i][j] = new double[nfreq];
		}
	}
	B_T_nmf = new double** [Nch]; // Nch X Nrank X Nfreq
	for (i = 0; i < Nch; i++)
	{
		B_T_nmf[i] = new double* [Nrank];
		for (j = 0; j < Nrank; j++)
		{
			B_T_nmf[i][j] = new double[nfreq];
		}
	}
	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nrank; j++)
		{
			for (k = 0; k < nfreq; k++)
			{
				T_nmf[i][j][k] = (rand() / max) + epsi;
				A_T_nmf[i][j][k] = 0.0;
				B_T_nmf[i][j][k] = 0.0;
			}
		}
	}

	lambda = new double* [Nch]; // Nch X Nfreq
	for (i = 0; i < Nch; i++)
	{
		lambda[i] = new double[nfreq];
	}
	lambda_n = new double[nfreq]; // Nch X Nfreq
	lambda_t = new double[nfreq]; // Nch X Nfreq
	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < nfreq; j++)
		{
			lambda[i][j] = 1.0;
		}
	}
	invWDE = new double* [Nch]; // Nch X Nfreq(Complex)
	for (i = 0; i < Nch; i++)
	{
		invWDE[i] = new double[nfreq * 2];
	}
	diag_WV = new double* [Nch]; // Nch X Nfreq(Complex)
	for (i = 0; i < Nch; i++)
	{
		diag_WV[i] = new double[nfreq * 2];
	}
	win_STFT = new double[nWin];
	for (i = 0; i < nWin; i++)
	{
		win_STFT[i] = sqrt((double)2 / 3) * 0.5 * (1.0 - cos(2.0 * (double)M_PI * (double)(i) / (nWin)));
	}
	W = new double** [Nch];
	for (i = 0; i < Nch; i++)
	{
		W[i] = new double* [Nch];
		for (j = 0; j < Nch; j++)
		{
			W[i][j] = new double[nfreq * 2];
		}
	}
	A = new double** [Nch];
	for (i = 0; i < Nch; i++)
	{
		A[i] = new double* [Nch];
		for (j = 0; j < Nch; j++)
		{
			A[i][j] = new double[nfreq * 2];
		}
	}
	V = new double*** [Nch];
	for (i = 0; i < Nch; i++)
	{
		V[i] = new double** [Nch];
		for (j = 0; j < Nch; j++)
		{
			V[i][j] = new double* [nfreq * 2];
			for (k = 0; k < nfreq * 2; k++)
			{
				V[i][j][k] = new double[Nch];
			}
		}
	}
	U = new double*** [Nch];
	for (i = 0; i < Nch; i++)
	{
		U[i] = new double** [Nch];
		for (j = 0; j < Nch; j++)
		{
			U[i][j] = new double* [nfreq * 2];
			for (k = 0; k < nfreq * 2; k++)
			{
				U[i][j][k] = new double[Nch];
			}
		}
	}

	//W Á¤ÀÇ
	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nch; j++)
		{
			for (freq = 0; freq < nfreq; freq++)
			{
				re = freq + freq;
				im = re + 1;
				if (i == j)
				{
					W[i][j][re] = 1.0;
					W[i][j][im] = 0.0;
					A[i][j][re] = 1.0;
					A[i][j][im] = 0.0;
				}
				else
				{
					W[i][j][re] = 0.0;
					W[i][j][im] = 0.0;
					A[i][j][re] = 0.0;
					A[i][j][im] = 0.0;
				}
			}
		}
	}

	//frameInd over 2
	p = new double* [Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		p[ch] = new double[nfreq];
	}
	Unumer = new double[nfreq * 2];
	Udenom = new double** [Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		Udenom[ch] = new double* [Nch];
		for (i = 0; i < Nch; i++)
		{
			Udenom[ch][i] = new double[nfreq * 2];
		}
	}
	p_U_X = new double* [Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		p_U_X[ch] = new double[nfreq * 2];
	}
	X_T_U = new double* [Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		X_T_U[ch] = new double[nfreq * 2];
	}
	p_U_X_X = new double*** [Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		p_U_X_X[ch] = new double** [Nch];
		for (i = 0; i < Nch; i++)
		{
			p_U_X_X[ch][i] = new double* [nfreq * 2];
			for (j = 0; j < nfreq * 2; j++)
			{
				p_U_X_X[ch][i][j] = new double[Nch];
			}
		}
	}

	//normalizing
	normCoef = new double[nfreq * 2];
	sqnorm = new double[nfreq * 2];

	WDE_V = new double* [Nch];
	for (i = 0; i < Nch; i++)
	{
		WDE_V[i] = new double[nfreq * 2];
	}
	w = new double* [Nch];
	for (i = 0; i < Nch; i++)
	{
		w[i] = new double[nfreq * 2];
	}
	dW = new double** [Nch];
	for (i = 0; i < Nch; i++)
	{
		dW[i] = new double* [Nch];
		for (j = 0; j < Nch; j++)
		{
			dW[i][j] = new double[nfreq * 2];
		}
	}

	//Calculate A
	Anumer = new double[nfreq * 2];
	AdW = new double** [Nch];
	for (i = 0; i < Nch; i++)
	{
		AdW[i] = new double* [Nch];
		for (j = 0; j < Nch; j++)
		{
			AdW[i][j] = new double[nfreq * 2];
		}
	}
	Adenom = new double** [Nch];
	for (i = 0; i < Nch; i++)
	{
		Adenom[i] = new double* [Nch];
		for (j = 0; j < Nch; j++)
		{
			Adenom[i][j] = new double[nfreq * 2];
		}
	}

	//Result
	Wbp = new double** [Nch];
	for (i = 0; i < Nch; i++)
	{
		Wbp[i] = new double* [Nch];
		for (j = 0; j < Nch; j++)
		{
			Wbp[i][j] = new double[nfreq * 2];
		}
	}
	Ytmp = new double* [Nch];
	for (i = 0; i < Nch; i++)
	{
		Ytmp[i] = new double[nfreq * 2];
	}
	Ybuff = new double* [Nch];
	for (i = 0; i < Nch; i++)
	{
		Ybuff[i] = new double[nWin];
	}
	int sample;
	for (i = 0; i < Nch; i++)
	{
		for (sample = 0; sample < nfreq * 2; sample++)
		{
			Ytmp[i][sample] = 0.0;
		}
		for (sample = 0; sample < nWin; sample++)
		{
			Ybuff[i][sample] = 0.0;
		}
	}

}

AUXIVA_MASK::~AUXIVA_MASK()
{
	int i, j, k;
	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nch; j++)
		{
			for (k = 0; k < nfreq * 2; k++)
			{
				delete[] V[i][j][k];
				delete[] U[i][j][k];
				delete[] p_U_X_X[i][j][k];
			}
			delete[] V[i][j];
			delete[] U[i][j];
			delete[] p_U_X_X[i][j];
		}
		delete[] V[i];
		delete[] U[i];
		delete[] p_U_X_X[i];
	}
	delete[] V;
	delete[] U;
	delete[] p_U_X_X;


	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nch; j++)
		{
			delete[] W[i][j];
		}
		for (j = 0; j < Nrank; j++)
		{
			delete[] T_nmf[i][j];
			delete[] A_T_nmf[i][j];
			delete[] B_T_nmf[i][j];
		}
		delete[] V_nmf[i];
		delete[] T_nmf[i];
		delete[] A_T_nmf[i];
		delete[] B_T_nmf[i];
		delete[] X[i];
		delete[] X_r[i];
		delete[] Y[i];
		delete[] Pwr[i];
		delete[] W[i];
		delete[] invWDE[i];
		delete[] diag_WV[i];
		delete[] lambda[i];

	}
	delete[] X;
	delete[] X_r;
	delete[] Y;
	delete[] Pwr;
	delete[] lambda;
	delete[] T_nmf;
	delete[] A_T_nmf;
	delete[] B_T_nmf;
	delete[] V_nmf;
	delete[] W;
	delete[] invWDE;
	delete[] diag_WV;
	delete[] win_STFT;
	delete[] lambda_n;
	delete[] lambda_t;
	//frameInd over 2
	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nch; j++)
		{
			delete[] Udenom[i][j];
		}
		delete[] Udenom[i];
		delete[] p_U_X[i];
		delete[] X_T_U[i];
		delete[] p[i];
	}
	delete[] Udenom;
	delete[] p_U_X;
	delete[] X_T_U;
	delete[] p;
	delete[] Unumer;

	//normalizing
	delete[] normCoef;
	delete[] sqnorm;
	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nch; j++)
		{
			delete[] A[i][j];
			delete[] dW[i][j];
		}
		delete[] A[i];
		delete[] dW[i];
		delete[] WDE_V[i];
		delete[] w[i];
	}
	delete[] A;
	delete[] dW;
	delete[] WDE_V;
	delete[] w;

	//Calculate A
	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nch; j++)
		{
			delete[] AdW[i][j];
			delete[] Adenom[i][j];
		}
		delete[] AdW[i];
		delete[] Adenom[i];
	}
	delete[] AdW;
	delete[] Adenom;
	delete[] Anumer;

	//result
	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nch; j++)
		{
			delete[] Wbp[i][j];
		}
		delete[] Wbp[i];
		delete[] Ytmp[i];
		delete[] Ybuff[i];
	}
	delete[] Wbp;
	delete[] Ytmp;
	delete[] Ybuff;
}

void AUXIVA_MASK::AUXIVA_MASK_lemma(double** input, int frameInd, double** output, double *Mask)
{
	
	double f_alpha_real = f_alpha;
	f_alpha = 0.94;
	if (frameInd == 50)
	{
		f_alpha = f_alpha_real;
	}
	int i, j, k, ch, freq, freqInd;
	int ch1, ch2;
	int re, im;
	for (ch = 0; ch < Nch; ch++)
	{
		for (i = 0; i < nWin; i++)
		{
			X[ch][i] = win_STFT[i] * input[ch][i];
		}
		hfft3(X[ch], nfft, 1);
	}

	for (freq = 0; freq < nfreq; freq++)
	{
		re = freq + freq;
		im = re + 1;
		for (ch1 = 0; ch1 < Nch; ch1++)
		{
			Y[ch1][re] = 0.0;
			Y[ch1][im] = 0.0;
			for (ch2 = 0; ch2 < Nch; ch2++)
			{
				Y[ch1][re] = Y[ch1][re] + W[ch1][ch2][re] * X[ch2][re] - W[ch1][ch2][im] * X[ch2][im];
				Y[ch1][im] = Y[ch1][im] + W[ch1][ch2][re] * X[ch2][im] + W[ch1][ch2][im] * X[ch2][re];
			}
		}
	}
	// Pwr
	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < nfreq; j++)
		{
			re = j + j;
			im = re + 1;
			Pwr[i][j] = (Y[i][re] * Y[i][re]) + (Y[i][im] * Y[i][im]);
			if (Pwr[i][j] < epsi)
			{
				Pwr[i][j] = epsi;
			}
		}
	}
	// optimize time activation at current frame according to bases
	// update bases frame by frame
#if ILRMA_OPT == 1
	for (i = 0; i < Nch; i++)
	{
		for (k = 0; k < nfreq; k++)
		{
			lambda[i][k] = 0.0;
			for (j = 0; j < Nrank; j++)
			{
				lambda[i][k] += V_nmf[i][j] * T_nmf[i][j][k];
			}
		}
		for (j = 0; j < Nrank; j++)
		{
			double Numer_V = 0.0;
			double Denom_V = 0.0;
			for (k = 0; k < nfreq; k++)
			{
				Numer_V += Pwr[i][k] * T_nmf[i][j][k] / (lambda[i][k] * lambda[i][k]);
				Denom_V += T_nmf[i][j][k] / (lambda[i][k]);
			}
			V_nmf[i][j] = V_nmf[i][j] * sqrt(Numer_V / Denom_V);
			if (V_nmf[i][j] < epsi)
			{
				V_nmf[i][j] = epsi;
			}
		}
		for (k = 0; k < nfreq; k++)
		{
			lambda[i][k] = 0.0;
			for (j = 0; j < Nrank; j++)
			{
				lambda[i][k] += V_nmf[i][j] * T_nmf[i][j][k];
			}
		}
		if (frameInd == 3)
		{
			for (k = 0; k < nfreq; k++)
			{
				for (j = 0; j < Nrank; j++)
				{
					A_T_nmf[i][j][k] = (T_nmf[i][j][k] * T_nmf[i][j][k]) * Pwr[i][k] * V_nmf[i][j] / (lambda[i][k] * lambda[i][k]);
					B_T_nmf[i][j][k] = V_nmf[i][j] / (lambda[i][k]);
					T_nmf[i][j][k] = sqrt(A_T_nmf[i][j][k] / B_T_nmf[i][j][k]);
					if (T_nmf[i][j][k] < epsi)
					{
						T_nmf[i][j][k] = epsi;
					}
				}
			}
		}
		else
		{
			for (k = 0; k < nfreq; k++)
			{
				for (j = 0; j < Nrank; j++)
				{
					A_T_nmf[i][j][k] = gamma[i] * A_T_nmf[i][j][k] + (1 - gamma[i]) * (T_nmf[i][j][k] * T_nmf[i][j][k]) * Pwr[i][k] * V_nmf[i][j] / (lambda[i][k] * lambda[i][k]);
					B_T_nmf[i][j][k] = gamma[i] * B_T_nmf[i][j][k] + (1 - gamma[i]) * V_nmf[i][j] / (lambda[i][k]);
					T_nmf[i][j][k] = sqrt(A_T_nmf[i][j][k] / B_T_nmf[i][j][k]);
					if (T_nmf[i][j][k] < epsi)
					{
						T_nmf[i][j][k] = epsi;
					}
				}
			}
		}
		for (k = 0; k < nfreq; k++)
		{
			lambda[i][k] = 0.0;
			for (j = 0; j < Nrank; j++)
			{
				lambda[i][k] += V_nmf[i][j] * T_nmf[i][j][k];
			}
		}
		for (j = 0; j < Nrank; j++)
		{
			double Numer_V = 0.0;
			double Denom_V = 0.0;
			for (k = 0; k < nfreq; k++)
			{
				Numer_V += Pwr[i][k] * T_nmf[i][j][k] / (lambda[i][k] * lambda[i][k]);
				Denom_V += T_nmf[i][j][k] / (lambda[i][k]);
			}
			V_nmf[i][j] = V_nmf[i][j] * sqrt(Numer_V / Denom_V);
			if (V_nmf[i][j] < epsi)
			{
				V_nmf[i][j] = epsi;
			}
		}
		for (k = 0; k < nfreq; k++)
		{
			lambda[i][k] = 0.0;
			for (j = 0; j < Nrank; j++)
			{
				lambda[i][k] += V_nmf[i][j] * T_nmf[i][j][k];
			}
		}
	}
#else
	// scaled variance by Target Mask
	for (k = 0; k < nfreq; k++)
	{
		if (Mask[k] < 1e-2)
		{
			Mask[k] = 1e-2;
		}
		if (Mask[k] > (1-1e-2))
		{
			Mask[k] = 1 - 1e-2;
		}
		//  p with Variances by Mask
		if (frameInd == 3)
		{
			eta[0][k] = Mask[k];
			eta[1][k] = 1 - Mask[k];
			lambda_t[k] = Pwr[0][k];
			double temp = 0;
			for (j = 1; j < Nch; j++)
			{
				temp += Pwr[j][k];
			}
			lambda_n[k] = temp / (Nch - 1);
		}
		else
		{
			eta[0][k] = f_alpha * eta[0][k] + (1 - f_alpha) * Mask[k];
			eta[1][k] = f_alpha * eta[1][k] + (1 - f_alpha) * (1 - Mask[k]);
			lambda_t[k] = gamma_t * lambda_t[k] + (1 - gamma_t) * Pwr[0][k];
			double temp = 0;
			for (j = 1; j < Nch; j++)
			{
				temp += Pwr[j][k];
			}
			lambda_n[k] = gamma_n * lambda_n[k] + (1 - gamma_n) * temp / (Nch - 1);

		}
		if (Mask[k] * lambda_t[k] < eps_floor)
		{
			p[0][k] = (1 - f_alpha) / eps_floor;
		}
		else
		{
			p[0][k] = (1 - f_alpha) / (Mask[k] * lambda_t[k] / eta[0][k]);
		}
		if ((1 - Mask[k]) * lambda_n[k] < eps_floor)
		{
			p[1][k] = (1 - f_alpha) / eps_floor;
		}
		else
		{
			p[1][k] = (1 - f_alpha) / ((1 - Mask[k]) * lambda_n[k] / eta[1][k]);
		}
	}
#endif

	for (ch = 0; ch < 2; ch++)
	{
		// Calculate V and U
		if (frameInd == 3)
		{
			for (freqInd = 0; freqInd < nfreq; freqInd++)
			{
				re = freqInd + freqInd;
				im = re + 1;
				for (ch1 = 0; ch1 < Nch; ch1++)
				{
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						V[ch1][ch2][re][ch] = (X[ch1][re] * X[ch2][re] + X[ch1][im] * X[ch2][im]) * p[ch][freqInd] / (1 - f_alpha);
						V[ch1][ch2][im][ch] = (X[ch1][im] * X[ch2][re] - X[ch1][re] * X[ch2][im]) * p[ch][freqInd] / (1 - f_alpha);
						U[ch1][ch2][re][ch] = 0.0;
						U[ch1][ch2][im][ch] = 0.0;
						if (ch1 == ch2)
						{
							U[ch1][ch2][re][ch] = V[ch1][ch2][re][ch] / ((V[ch1][ch2][re][ch] * V[ch1][ch2][re][ch]) + (V[ch1][ch2][im][ch] * V[ch1][ch2][im][ch]));
							U[ch1][ch2][im][ch] = -V[ch1][ch2][im][ch] / ((V[ch1][ch2][re][ch] * V[ch1][ch2][re][ch]) + (V[ch1][ch2][im][ch] * V[ch1][ch2][im][ch]));
						}
					}
				}
			}
		}
		else // over 3th frames
		{
			for (freqInd = 0; freqInd < nfreq; freqInd++)
			{
				re = freqInd + freqInd;
				im = re + 1;
				for (ch1 = 0; ch1 < Nch; ch1++)
				{
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						V[ch1][ch2][re][ch] = f_alpha * V[ch1][ch2][re][ch] + p[ch][freqInd] * (X[ch1][re] * X[ch2][re] + X[ch1][im] * X[ch2][im]);
						V[ch1][ch2][im][ch] = f_alpha * V[ch1][ch2][im][ch] + p[ch][freqInd] * (X[ch1][im] * X[ch2][re] - X[ch1][re] * X[ch2][im]);
					}
				}
				Unumer[re] = f_alpha * f_alpha;
				Unumer[im] = 0.0;
				for (ch1 = 0; ch1 < Nch; ch1++)
				{
					p_U_X[ch1][re] = 0.0;
					p_U_X[ch1][im] = 0.0;
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						p_U_X[ch1][re] += p[ch][freqInd] * (U[ch1][ch2][re][ch] * X[ch2][re] - U[ch1][ch2][im][ch] * X[ch2][im]);
						p_U_X[ch1][im] += p[ch][freqInd] * (U[ch1][ch2][im][ch] * X[ch2][re] + U[ch1][ch2][re][ch] * X[ch2][im]);
					}
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						p_U_X_X[ch1][ch2][re][ch] = p_U_X[ch1][re] * X[ch2][re] + p_U_X[ch1][im] * X[ch2][im];
						p_U_X_X[ch1][ch2][im][ch] = p_U_X[ch1][im] * X[ch2][re] - p_U_X[ch1][re] * X[ch2][im];
					}
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						int ch3;
						Udenom[ch1][ch2][re] = 0.0;
						Udenom[ch1][ch2][im] = 0.0;
						for (ch3 = 0; ch3 < Nch; ch3++)
						{
							Udenom[ch1][ch2][re] += p_U_X_X[ch1][ch3][re][ch] * U[ch2][ch3][re][ch] + p_U_X_X[ch1][ch3][im][ch] * U[ch2][ch3][im][ch];
							Udenom[ch1][ch2][im] += p_U_X_X[ch1][ch3][im][ch] * U[ch2][ch3][re][ch] - p_U_X_X[ch1][ch3][re][ch] * U[ch2][ch3][im][ch];
						}
					}
					X_T_U[ch1][re] = 0.0;
					X_T_U[ch1][im] = 0.0;
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						X_T_U[ch1][re] += X[ch2][re] * U[ch2][ch1][re][ch] + X[ch2][im] * U[ch2][ch1][im][ch];
						X_T_U[ch1][im] += X[ch2][re] * U[ch2][ch1][im][ch] - X[ch2][im] * U[ch2][ch1][re][ch];
					}
					Unumer[re] += (f_alpha * p[ch][freqInd]) * (X_T_U[ch1][re] * X[ch1][re] - X_T_U[ch1][im] * X[ch1][im]);
					Unumer[im] += (f_alpha * p[ch][freqInd]) * (X_T_U[ch1][re] * X[ch1][im] + X_T_U[ch1][im] * X[ch1][re]);
				}
				if (sqrt((Unumer[re] * Unumer[re]) + (Unumer[im] * Unumer[im])) < epsi)
				{
					Unumer[re] = 1E-6;
					Unumer[im] = 0.0;
				}
				for (ch1 = 0; ch1 < Nch; ch1++)
				{
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						U[ch1][ch2][re][ch] = (U[ch1][ch2][re][ch] / f_alpha) - (Udenom[ch1][ch2][re] * Unumer[re] + Udenom[ch1][ch2][im] * Unumer[im]) / ((Unumer[re] * Unumer[re]) + (Unumer[im] * Unumer[im]));
						U[ch1][ch2][im][ch] = (U[ch1][ch2][im][ch] / f_alpha) - (Udenom[ch1][ch2][im] * Unumer[re] - Udenom[ch1][ch2][re] * Unumer[im]) / ((Unumer[re] * Unumer[re]) + (Unumer[im] * Unumer[im]));
					}
				}
				for (ch1 = 0; ch1 < Nch; ch1++)
				{
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						U[ch1][ch2][re][ch] = (U[ch1][ch2][re][ch] + U[ch2][ch1][re][ch]) / 2;
						U[ch1][ch2][im][ch] = (U[ch1][ch2][im][ch] - U[ch2][ch1][im][ch]) / 2;
					}
				}
			}
		}
	}
	for (ch = 0; ch < Nch; ch++)
	{
		int ch4;
		if (ch == 0)
		{
			ch4 = 0;
		}
		else
		{
			ch4 = 1;
		}
		//Calculate invWV
		for (freqInd = 0; freqInd < nfreq; freqInd++)
		{
			re = freqInd + freqInd;
			im = re + 1;
			for (ch1 = 0; ch1 < Nch; ch1++)
			{
				invWDE[ch1][re] = 0.0;
				invWDE[ch1][im] = 0.0;
				for (ch2 = 0; ch2 < Nch; ch2++)
				{
					invWDE[ch1][re] += U[ch1][ch2][re][ch4] * A[ch2][ch][re] - U[ch1][ch2][im][ch4] * A[ch2][ch][im];
					invWDE[ch1][im] += U[ch1][ch2][im][ch4] * A[ch2][ch][re] + U[ch1][ch2][re][ch4] * A[ch2][ch][im];
				}
			}
		}
		// Normalizing
		for (freqInd = 0; freqInd < nfreq; freqInd++)
		{
			re = freqInd * 2;
			im = re + 1;
			normCoef[re] = 0.0;
			normCoef[im] = 0.0;
			for (ch1 = 0; ch1 < Nch; ch1++)
			{
				WDE_V[ch1][re] = 0.0;
				WDE_V[ch1][im] = 0.0;
				for (ch2 = 0; ch2 < Nch; ch2++)
				{
					WDE_V[ch1][re] += invWDE[ch2][re] * V[ch2][ch1][re][ch4] + invWDE[ch2][im] * V[ch2][ch1][im][ch4];
					WDE_V[ch1][im] += invWDE[ch2][re] * V[ch2][ch1][im][ch4] - invWDE[ch2][im] * V[ch2][ch1][re][ch4];
				}
				normCoef[re] += WDE_V[ch1][re] * invWDE[ch1][re] - WDE_V[ch1][im] * invWDE[ch1][im];
				normCoef[im] += WDE_V[ch1][re] * invWDE[ch1][im] + WDE_V[ch1][im] * invWDE[ch1][re];
			}
			sqnorm[re] = sqrt(normCoef[re]);
			sqnorm[im] = 0.0;
			if (sqnorm[re] < epsi)
			{
				sqnorm[re] = epsi;
			}

			for (ch1 = 0; ch1 < Nch; ch1++)
			{
				w[ch1][re] = invWDE[ch1][re] / sqnorm[re];
				w[ch1][im] = invWDE[ch1][im] / sqnorm[re];
			}
		}

		for (freq = 0; freq < nfreq; freq++)
		{
			re = freq * 2;
			im = re + 1;
			for (ch1 = 0; ch1 < Nch; ch1++)
			{
				dW[ch][ch1][re] = w[ch1][re] - W[ch][ch1][re];
				dW[ch][ch1][im] = -w[ch1][im] - W[ch][ch1][im];
				W[ch][ch1][re] = w[ch1][re];
				W[ch][ch1][im] = -w[ch1][im];
			}
		}

		// Calculate A
		for (freqInd = 0; freqInd < nfreq; freqInd++)
		{
			re = freqInd * 2;
			im = re + 1;
			Anumer[re] = 1;
			Anumer[im] = 0;
			for (ch1 = 0; ch1 < Nch; ch1++)
			{
				for (ch2 = 0; ch2 < Nch; ch2++)
				{
					AdW[ch1][ch2][re] = A[ch1][ch][re] * dW[ch][ch2][re] - A[ch1][ch][im] * dW[ch][ch2][im];
					AdW[ch1][ch2][im] = A[ch1][ch][re] * dW[ch][ch2][im] + A[ch1][ch][im] * dW[ch][ch2][re];
				}
				for (ch2 = 0; ch2 < Nch; ch2++)
				{
					int ch3;
					Adenom[ch1][ch2][re] = 0.0;
					Adenom[ch1][ch2][im] = 0.0;
					for (ch3 = 0; ch3 < Nch; ch3++)
					{
						Adenom[ch1][ch2][re] += AdW[ch1][ch3][re] * A[ch3][ch2][re] - AdW[ch1][ch3][im] * A[ch3][ch2][im];
						Adenom[ch1][ch2][im] += AdW[ch1][ch3][re] * A[ch3][ch2][im] + AdW[ch1][ch3][im] * A[ch3][ch2][re];
					}
				}
				Anumer[re] += dW[ch][ch1][re] * A[ch1][ch][re] - dW[ch][ch1][im] * A[ch1][ch][im];
				Anumer[im] += dW[ch][ch1][re] * A[ch1][ch][im] + dW[ch][ch1][im] * A[ch1][ch][re];
			}
			if (sqrt((Anumer[re] * Anumer[re]) + (Anumer[im] * Anumer[im])) < epsi)
			{
				Anumer[re] = epsi;
				Anumer[im] = 0.0;
			}
			for (ch1 = 0; ch1 < Nch; ch1++)
			{
				for (ch2 = 0; ch2 < Nch; ch2++)
				{
					A[ch1][ch2][re] = A[ch1][ch2][re] - (Adenom[ch1][ch2][re] * Anumer[re] + Adenom[ch1][ch2][im] * Anumer[im]) / ((Anumer[re] * Anumer[re]) + (Anumer[im] * Anumer[im]));
					A[ch1][ch2][im] = A[ch1][ch2][im] - (Adenom[ch1][ch2][im] * Anumer[re] - Adenom[ch1][ch2][re] * Anumer[im]) / ((Anumer[re] * Anumer[re]) + (Anumer[im] * Anumer[im]));
				}
			}
		}
	}

	// result - back projection using A
	for (freqInd = 0; freqInd < nfreq; freqInd++)
	{
		re = freqInd * 2;
		im = re + 1;
		for (ch1 = 0; ch1 < Nch; ch1++)
		{
			for (ch2 = 0; ch2 < Nch; ch2++)
			{
				Wbp[ch1][ch2][re] = A[ch1][ch1][re] * W[ch1][ch2][re] - A[ch1][ch1][im] * W[ch1][ch2][im];
				Wbp[ch1][ch2][im] = A[ch1][ch1][re] * W[ch1][ch2][im] + A[ch1][ch1][im] * W[ch1][ch2][re];
			}
			Ytmp[ch1][re] = 0.0;
			Ytmp[ch1][im] = 0.0;
			for (ch2 = 0; ch2 < Nch; ch2++)
			{
				Ytmp[ch1][re] += Wbp[ch1][ch2][re] * X[ch2][re] - Wbp[ch1][ch2][im] * X[ch2][im];
				Ytmp[ch1][im] += Wbp[ch1][ch2][re] * X[ch2][im] + Wbp[ch1][ch2][im] * X[ch2][re];
			}
		}
	}

	for (ch1 = 0; ch1 < Nch; ch1++)
	{
		hfft3(Ytmp[ch1], nfft, -1);
		for (i = 0; i < nWin - BufferSize; i++)
		{
			Ybuff[ch1][i] = Ybuff[ch1][BufferSize + i];
			Ybuff[ch1][i] += Ytmp[ch1][i] * win_STFT[i];
		}
		for (; i < nWin; i++)
		{
			Ybuff[ch1][i] = Ytmp[ch1][i] * win_STFT[i];
		}
	}
	for (ch = 0; ch < Nch; ch++)
	{
		for (i = 0; i < BufferSize; i++)
		{
			output[ch][i] = Ybuff[ch][i];
		}
	}
}

