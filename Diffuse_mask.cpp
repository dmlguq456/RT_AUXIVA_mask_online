#include <stdio.h>
#include <algorithm>
#include "AUXIVA_mask_Online.h"
#include "header.h"
#include "sigproc.h"

using namespace std;

CDR::CDR()
{
	nfft = nWin;
	nshift = nWin / 4;
	nol = 3 * nWin / 4;
	nfreq = nfft / 2 + 1;
	//epsi = 0.000001;
	epsi = 2.220446049250313*1E-16;

	int i, j;
	int Npair = Nch * (Nch - 1) / 2;

	win_STFT = new double[nWin];
	for (i = 0; i < nWin; i++)
	{
		win_STFT[i] = sqrt((double)2 / 3) * 0.5 * (1.0 - cos(2.0 * (double)M_PI * (double)(i) / (nWin)));
	}
	X = new double* [Nch];
	for (i = 0; i < Nch; i++)
	{
		X[i] = new double[nfreq * 2];
	}

	R_freq = new double [nfreq];
	
	Cnn = new double* [Npair];
	for (i = 0; i < Npair; i++)
	{
		Cnn[i] = new double[nfreq];
	}
	Cxx = new double* [Npair];
	for (i = 0; i < Npair; i++)
	{
		Cxx[i] = new double[2*nfreq];
	}
	Auto = new double** [Npair];
	for (i = 0; i < Npair; i++)
	{
		Auto[i] = new double* [2];
		for (j = 0; j < 2; j++)
		{
			Auto[i][j] = new double[2 * nfreq];
		}
	}
	Cross = new double* [Npair];
	for (j = 0; j < Npair; j++)
	{
		Cross[j] = new double[2 * nfreq];
	}
	cdr = new double* [Npair];
	for (j = 0; j < Npair; j++)
	{
		cdr[j] = new double[nfreq];
	}
	DIFF = new double* [nfreq];
	for (j = 0; j < nfreq; j++)
	{
		DIFF[j] = new double[Npair];
	}
	D = new double[nfreq];
	
	
}

CDR::~CDR()
{
	int i, j;
	int Npair = Nch * (Nch - 1) / 2;

	delete[] R_freq;
	for (i = 0; i < Nch; i++)
	{
		delete[] X[i];
	}
	delete[] X;
	delete[] win_STFT;

	for (i = 0; i < Npair; i++)
	{
		delete[] Cnn[i];
		delete[] Cxx[i];
		for (j = 0; j < 2; j++)
		{
			delete[] Auto[i][j];
		}
		delete[] Auto[i];
		delete[] Cross[i];
		delete[] cdr[i];

	}
	delete[] Auto;
	delete[] Cross;
	delete[] Cnn;
	delete[] Cxx;
	for (j = 0; j < nfreq; j++)
	{
		delete[] DIFF[j];
	}
	delete[] DIFF;
	delete[] D;

}

void CDR::CDR_mask(double **input, int frame_idx, double *Mask, double **mic_array)
{
	int i, ch, freq_idx, p_idx;
	int p1, p2;
	int re, im;
	int Npair = Nch * (Nch - 1) / 2;
	double lambda = 0.5;
	int kappa = 10;
	double threshold = 0.2;
	double mic_dist;

	for (ch = 0; ch < Nch; ch++)
	{
		for (i = 0; i < nWin; i++)
		{
			X[ch][i] = win_STFT[i] * input[ch][i];
		}
		hfft3(X[ch], nfft, 1);
	}


	for (i = 0; i < nfreq; i++)
	{
		R_freq[i] = i * SamplingFreq / double(nfft);
	}

	// Mic distance at every Mic pair and Coherence of Ideal Diffuse Noise
	p_idx = 0;
	for (p1 = 0; p1 < (Nch - 1); p1++)
	{
		for (p2 = p1+1; p2 < Nch; p2++)
		{
			double temp = 0;
			for (i = 0; i < 3; i++)
			{
				temp += (mic_array[p1][i] - mic_array[p2][i]) * (mic_array[p1][i] - mic_array[p2][i]);
			}
			mic_dist = sqrt(temp);
			for (freq_idx = 0; freq_idx < nfreq; freq_idx++)
			{
				temp = 2 * (double)M_PI * R_freq[freq_idx] * mic_dist / 343;
				if (temp == 0)
				{
					Cnn[p_idx][freq_idx] = 1;
				}
				else
				{
					Cnn[p_idx][freq_idx] = sin(temp) / temp;
				}
			}
			p_idx++;
		}
	}

	for (freq_idx = 0; freq_idx < nfreq; freq_idx++)
	{
		p_idx = 0;
		for (p1 = 0; p1 < (Nch - 1); p1++)
		{
			for (p2 = p1+1; p2 < Nch; p2++)
			{
				// Auto  and Cross Correlation
				if (frame_idx == 3)
				{
					re = freq_idx + freq_idx;
					im = re + 1;
					Auto[p_idx][0][freq_idx] = X[p1][re] * X[p1][re] + X[p1][im] * X[p1][im];
					Auto[p_idx][1][freq_idx] = X[p2][re] * X[p2][re] + X[p2][im] * X[p2][im];
					if (Auto[p_idx][0][freq_idx] < Cnn[p_idx][freq_idx] * Cnn[p_idx][freq_idx])
					{
						Auto[p_idx][0][freq_idx] = Cnn[p_idx][freq_idx] * Cnn[p_idx][freq_idx];
					}
					if (Auto[p_idx][1][freq_idx] < Cnn[p_idx][freq_idx] * Cnn[p_idx][freq_idx])
					{
						Auto[p_idx][1][freq_idx] = Cnn[p_idx][freq_idx] * Cnn[p_idx][freq_idx];
					}
					Cross[p_idx][re] = X[p1][re] * X[p2][re] + X[p1][im] * X[p2][im];
					Cross[p_idx][im] = X[p1][im] * X[p2][re] - X[p1][re] * X[p2][im];
				}
				else
				{
					re = freq_idx + freq_idx;
					im = re + 1;
					Auto[p_idx][0][freq_idx] = lambda * Auto[p_idx][0][freq_idx] + (1 - lambda) * (X[p1][re] * X[p1][re] + X[p1][im] * X[p1][im]);
					Auto[p_idx][1][freq_idx] = lambda * Auto[p_idx][1][freq_idx] + (1 - lambda) * (X[p2][re] * X[p2][re] + X[p2][im] * X[p2][im]);
					Cross[p_idx][re] = lambda * Cross[p_idx][re] + (1 - lambda) * (X[p1][re] * X[p2][re] + X[p1][im] * X[p2][im]);
					Cross[p_idx][im] = lambda * Cross[p_idx][im] + (1 - lambda) * (X[p1][im] * X[p2][re] - X[p1][re] * X[p2][im]);
				}

				Cxx[p_idx][re] = Cross[p_idx][re] / sqrt(Auto[p_idx][0][freq_idx] * Auto[p_idx][1][freq_idx]);
				Cxx[p_idx][im] = Cross[p_idx][im] / sqrt(Auto[p_idx][0][freq_idx] * Auto[p_idx][1][freq_idx]);

				double Cnn_sq = Cnn[p_idx][freq_idx] * Cnn[p_idx][freq_idx];
				double Cxx_sq = Cxx[p_idx][re] * Cxx[p_idx][re] + Cxx[p_idx][im] * Cxx[p_idx][im];
				double temp = sqrt(Cxx_sq + Cnn_sq * Cxx[p_idx][re] * Cxx[p_idx][re] - Cnn_sq * Cxx_sq - 2 * Cnn[p_idx][freq_idx] * Cxx[p_idx][re] + Cnn_sq);

				cdr[p_idx][freq_idx] = (Cnn[p_idx][freq_idx] * Cxx[p_idx][re] - Cxx_sq - temp) / (Cxx_sq - 1);
				DIFF[freq_idx][p_idx] = 1. / (1 + cdr[p_idx][freq_idx]);
				p_idx++;
			}
		}
		sort(DIFF[freq_idx], DIFF[freq_idx] + Npair);
		if (Npair % 2 == 1)
		{
			int median_idx = (Npair + 1) / 2;
			D[freq_idx] = DIFF[freq_idx][median_idx];
		}
		else
		{
			int median_idx = Npair / 2;
			D[freq_idx] = ( DIFF[freq_idx][median_idx] + DIFF[freq_idx][median_idx] ) / 2;
		}
		Mask[freq_idx] = 1 / (1 + exp(-kappa * (D[freq_idx] - threshold)));
		Mask[freq_idx] = 1 - Mask[freq_idx];
	}
}