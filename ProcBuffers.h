#include "AUXIVA_mask_Online.h"
#define MAKE_FILE		1			//option 1 : wav 저장 (IVA출력 + 입력원본)		2: strout 출력(IVA출력)		3: strout 출력 (IVA출력 + 입력 원본)
#define SAVE_OPT		2			//option 1 : Target Output만 저장				2: Target 과 Noise Output 모두 저장
class ProcBuffers
{
private:
	

public:
	ProcBuffers();
	~ProcBuffers();
	static int Process(double **input, int Nframe, double **output);
};