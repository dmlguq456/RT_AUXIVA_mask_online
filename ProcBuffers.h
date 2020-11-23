#include "AUXIVA_mask_Online.h"
#define MAKE_FILE		1			//option 1 : wav 저장 (IVA출력 + 입력원본)		2: strout 출력(IVA출력)		3: strout 출력 (IVA출력 + 입력 원본)
#define SAVE_OPT		2			//option 1 : Target Output만 저장				2: Target 과 Noise Output 모두 저장
class ProcBuffers
{
private:
	int ch_save;
	int BuffCnt;
	int isNew16k;
	double** InitCond, * XX_LP, * XX, ** xx_lp, ** x;

	double** out_buff;
	short** IVA_out;
	double** input_temp;
	double** output;

	double** in_buff;
	short** origin_out;
	double** input;

	double* Mask;
	double** mic_array;


public:
	ProcBuffers();
	~ProcBuffers();
	int Process(double **input, int Nframe, double **output);
};