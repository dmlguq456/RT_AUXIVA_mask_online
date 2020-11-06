#define Nch			3
#define nWin		2048
#define BufferSize		512
#define SamplingFreq	16000
#define ILRMA_OPT		0

class AUXIVA_MASK {

private:
	int nfft;
	int nshift;
	int nol;
	int nfreq;
	double epsi;
	double eps_floor;
	double f_alpha;
	double gamma_t;
	double gamma_n;
	double *gamma;
	double **eta;
	double *win_STFT;
	double **X;
	double **X_r;
	double **Y;
	double ***W;
	double** Pwr;
	double ****V;
	double ****U;
	double **diag_WV;
	double **invWDE;

	int Nrank;
	double** lambda; // Nch X Nrank
	double* lambda_n;
	double* lambda_t;
	double* D; // Nfreq
	double** V_nmf; // Nch X Nrank
	double*** T_nmf; // Nch X Nrank X Nfreq
	double*** A_T_nmf; // Nch X Nrank X Nfreq
	double*** B_T_nmf; // Nch X Nrank X Nfreq
	double Numer_V;
	double Denom_V;


	//frameInd over 2
	double **p;
	double **p_U_X;
	double ****p_U_X_X;
	double ***Udenom;
	double *Unumer;
	double **X_T_U;

	//normalizing
	double *normCoef;
	double *sqnorm;
	double ***A;
	double **WDE_V;
	double unW;
	double **w;
	double ***dW;

	//Calculate A
	double ***AdW;
	double ***Adenom;
	double *Anumer;

	//result
	double ***Wbp;
	double **Ytmp;
	double **Ybuff;

public:
	AUXIVA_MASK();
	~AUXIVA_MASK();
	void AUXIVA_MASK_lemma(double **input, int frameInd, double **output, double *Mask);
};


class CDR {

private:
	int nfft;
	int nshift;
	int nol;
	int nfreq;
	double epsi;

	double* win_STFT;
	double* R_freq;
	double** Cnn;
	double** Cxx;
	double** X;
	double*** Auto;
	double** Cross;
	double** cdr;
	double** DIFF;
	double* D;
public:
	CDR();
	~CDR();
	void CDR_mask(double** input, int frameInd, double* Mask, double** mic_array);
};
