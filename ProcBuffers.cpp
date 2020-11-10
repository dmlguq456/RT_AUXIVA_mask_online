#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include "ProcBuffers.h"
#include "sigproc.h"
#include <iostream>
#include <time.h>


using namespace std;

AUXIVA_MASK* iip_AUX;
CDR* iip_CDR;

#if MAKE_FILE == 1
FILE **IVA;
FILE** IN;

#endif

//For downsampling (48k to 16k)
double **InitCond, *XX_LP, *XX, **xx_lp, **x;

double **out_buff;
short **IVA_out;
double **input_temp;
double **output;

double** in_buff;
short** origin_out;
double** input;

double* Mask;
double** mic_array;

int ch_save;

ProcBuffers::ProcBuffers()
{
	if (SAVE_OPT == 1) ch_save = 1;
	else if (SAVE_OPT == 2) ch_save = Nch;

	int i, ch;

	iip_AUX = new AUXIVA_MASK();
	iip_CDR = new CDR();

	mic_array = new double* [Nch];
	for (i = 0; i < Nch; i++)
	{
		mic_array[i] = new double[3];
	}
	// ---------- Mic Configuration Setting ---------- //
	// 마이크의 좌표를 입력 
	// 모든 mic pair간의 거리를 찾기 위해서 필요한 것으로 따라서
	// 1. reference mic에 해당하는 mic_array[0]은 모두 0으로 설정
	// 2. x.y.z 기준은 임의로 잡고 진행해도 전혀 무방
	// mic_array[ch][0] = mic_array[0]과의 상대적인 x좌표입력
	// mic_array[ch][1] = mic_array[0]과의 상대적인 y좌표입력
	// mic_array[ch][2] = mic_array[0]과의 상대적인 z좌표입력
	
	mic_array[0][0] = 0;
	mic_array[0][1] = 0;
	mic_array[0][2] = 0;
	mic_array[1][0] = 0;
	mic_array[1][1] = 0;
	mic_array[1][2] = 0.2;
	mic_array[2][0] = 0;
	mic_array[2][1] = 0;
	mic_array[2][2] = -0.2;
	// ----------------------------------------------- //


	input_temp = new double*[Nch];
	output = new double*[Nch];
	out_buff = new double*[Nch];
	IVA_out = new short*[Nch];

	in_buff = new double* [Nch];
	origin_out = new short*[Nch];
	input = new double* [Nch];

	for (ch = 0; ch < Nch; ch++)
	{
		input_temp[ch] = new double[nWin];
		output[ch] = new double[BufferSize];
		out_buff[ch] = new double[BufferSize];
		IVA_out[ch] = new short[BufferSize];

		in_buff[ch] = new double[BufferSize];
		origin_out[ch] = new short[BufferSize];
		input[ch] = new double[BufferSize];


		for (i = 0; i < nWin; i++)
		{
			input_temp[ch][i] = 0.0;
		}
	}

	Mask = new double[nWin / 2 + 1];

	//For downsampling
	x = new double *[Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		x[ch] = new double[BufferSize];
	}
	InitCond = new double *[Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		InitCond[ch] = new double[BufferSize];
	}
	xx_lp = new double *[Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		xx_lp[ch] = new double[3 * BufferSize];
	}
	XX = new double[BufferSize * 2 + 2];
	XX_LP = new double[BufferSize * 2 + 2];
	for (ch = 0; ch < Nch; ch++)
	{
		for (i = 0; i < BufferSize; i++)
		{
			x[ch][i] = 0.0;
			InitCond[ch][i] = 0.0;
		}
		for (i = 0; i < 3 * BufferSize; i++)
		{
			xx_lp[ch][i] = 0.0;
		}
	}

	for (i = 0; i < BufferSize * 2 + 2; i++)
	{
		XX[i] = 0;
		XX_LP[i] = 0;
	}
	
#if MAKE_FILE == 1

	char file_name1[2][500];
	IVA = new FILE*[Nch];


	for (ch = 0; ch < ch_save; ch++)
	{
		if (ch == 0)
		{
			sprintf(file_name1[0], ".\\output\\IVA_Target.pcm");
		}
		else
		{
			sprintf(file_name1[0], ".\\output\\IVA_Noise_ch%d.pcm", ch);
		}
		IVA[ch] = fopen(file_name1[0], "wb");
	}

	char file_name2[2][500];
	IN = new FILE * [Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		sprintf(file_name2[0], ".\\input\\IN_ch%d.pcm", ch + 1);
		IN[ch] = fopen(file_name2[0], "wb");
	}
#endif

}

ProcBuffers::~ProcBuffers()
{
	int ch;
	for (ch = 0; ch < Nch; ch++)
	{
		delete[] xx_lp[ch];
		delete[] x[ch];
		delete[] InitCond[ch];
	}
	delete[] xx_lp;
	delete[] x;
	delete[] InitCond;
	delete[] XX_LP;
	delete[] XX;
	delete[] Mask;


	for (ch = 0; ch < Nch; ch++)
	{
		delete[] input_temp[ch];
		delete[] output[ch];
		delete[] out_buff[ch];
		delete[] IVA_out[ch];

		delete[] in_buff[ch];
		delete[] origin_out[ch];
		delete[] input[ch];
	}
	delete[] input_temp;
	delete[] output;
	delete[] out_buff;
	delete[] IVA_out;

	delete[] in_buff;
	delete[] origin_out;
	delete[] input;


#if MAKE_FILE == 1
	char file_name1[2][500];
	time_t t = time(NULL);
	struct tm tm = *localtime(&t);
	char file_name_time[500];
	

	for (ch = 0; ch < ch_save; ch++)
	{
		fclose(IVA[ch]);
		if (ch == 0)
		{
			sprintf(file_name1[0], ".\\output\\IVA_Target.pcm");
			sprintf(file_name1[1], ".\\output\\IVA_Target_%d%d%d_%d%d%d.wav", tm.tm_year - 100, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
		}
		else
		{
			sprintf(file_name1[0], ".\\output\\IVA_Noise_ch%d.pcm", ch);
			sprintf(file_name1[1], ".\\output\\IVA_Noise_ch%d_%d%d%d_%d%d%d.wav", ch, tm.tm_year - 100, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
		}
		pcm2wav(file_name1[0], file_name1[1], (long)(SamplingFreq));
		remove(file_name1[0]);
	}
	delete[] IVA;

	char file_name2[2][500];
	for (ch = 0; ch < Nch; ch++)
	{
		fclose(IN[ch]);
		sprintf(file_name2[0], ".\\input\\IN_ch%d.pcm", ch + 1);
		sprintf(file_name2[1], ".\\input\\IN_ch%d_%d%d%d_%d%d%d.wav", ch + 1, tm.tm_year - 100, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
		pcm2wav(file_name2[0], file_name2[1], (long)(SamplingFreq));
		remove(file_name2[0]);
	}
	delete[] IN;


#endif

	delete iip_AUX;
}

int ProcBuffers::Process(double **input, int Nframe, double **output)
{
	int i, j, ch;
	static int BuffCnt = 0, isNew16k = 0;
	//OCTA-Capture로 들어오는 input 48k이므로 process를 진행하기 위해 16k로 Down Sampling을 해야한다.
	isNew16k = (BuffCnt == 2);
	for (ch = 0; ch < Nch; ch++)
	{
		for (i = 0; i < BufferSize; i++)
		{
			XX[i + BufferSize] = 0;
			XX[i] = input[ch][i];
		}

		for (i = 0; i < BufferSize; i++)
		{
			xx_lp[ch][BuffCnt*BufferSize + i] = InitCond[ch][i] + XX[i];
			InitCond[ch][i] = XX_LP[BufferSize + i];
		}
		if (isNew16k == 1)
		{
			for (i = 0; i < BufferSize; i++)
			{
				x[ch][i] = xx_lp[ch][i + i + i];
			}
		}
	}
	BuffCnt = (BuffCnt + 1) % 3;

	if (isNew16k == 1)
	{
		for (ch = 0; ch < Nch; ch++)
		{
			for (i = 0; i < 3 * BufferSize; i++)
			{
				input_temp[ch][i] = input_temp[ch][BufferSize + i];
			}
			for (i = 0; i < BufferSize; i++)
			{
				input_temp[ch][3 * BufferSize + i] = x[ch][i];
			}
		}

		iip_CDR->CDR_mask(input_temp, Nframe, Mask, mic_array);
		iip_AUX->AUXIVA_MASK_lemma(input_temp, Nframe, output, Mask);

#if MAKE_FILE == 1
		for (i = 0; i < ch_save; i++)
		{
			for (j = 0; j < BufferSize; j++)
			{
				out_buff[i][j] = output[i][j] * 32768.0;
				IVA_out[i][j] = (short)(out_buff[i][j]);
			}
			fwrite(IVA_out[i], sizeof(short), BufferSize, IVA[i]);
		}

		for (i = 0; i < Nch; i++)
		{
			for (j = 0; j < BufferSize; j++)
			{
				in_buff[i][j] = input_temp[i][j] * 32768.0;
				origin_out[i][j] = (short)(in_buff[i][j]);
			}
			fwrite(origin_out[i], sizeof(short), BufferSize, IN[i]);
		}

#elif MAKE_FILE == 2

		for (j = 0; j < BufferSize; j++)
		{
			for (i = 0; i < ch_save; i++)
			{
				out_buff[i][j] = output[i][j] * 32768.0;
				IVA_out[i][j] = (short)(out_buff[i][j]);
				cout << IVA_out[i][j];
				if (i != Nch - 1)
					cout << "	";
			}
			cout << "\n";
		}

#else MAKE_FILE == 3

		for (j = 0; j < BufferSize; j++)
		{
			for (i = 0; i < ch_save; i++)
			{
				out_buff[i][j] = output[i][j] * 32768.0;
				IVA_out[i][j] = (short)(out_buff[i][j]);
				cout << IVA_out[i][j] << "	";
			}
			for (i = 0; i < Nch; i++)
			{
				in_buff[i][j] = input_temp[i][j] * 32768.0;
				origin_out[i][j] = (short)(in_buff[i][j]);
				cout << origin_out[i][j];
				if (i != Nch - 1)
					cout << "	";
			}
			cout << "\n";
		}

#endif
	}
	return 0;
}