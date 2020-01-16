#include "defines.h"
#include <stdarg.h>


extern PRECISION **PUNTEROS_CALCULOS_COMPARTIDOS;
extern int POSW_PUNTERO_CALCULOS_COMPARTIDOS;
extern int POSR_PUNTERO_CALCULOS_COMPARTIDOS;


extern PRECISION *dtaux, *etai_gp3, *ext1, *ext2, *ext3, *ext4;
extern PRECISION *gp4_gp2_rhoq, *gp5_gp2_rhou, *gp6_gp2_rhov;
extern PRECISION *gp1, *gp2, *dt, *dti, *gp3, *gp4, *gp5, *gp6, *etai_2;
extern PRECISION *dgp1, *dgp2, *dgp3, *dgp4, *dgp5, *dgp6, *d_dt;
extern PRECISION *d_ei, *d_eq, *d_eu, *d_ev, *d_rq, *d_ru, *d_rv;
extern PRECISION *dfi, *dshi;
extern PRECISION *fi_p, *fi_b, *fi_r, *shi_p, *shi_b, *shi_r;
extern PRECISION *spectra, *d_spectra,*spectra_mac;
extern PRECISION *etain, *etaqn, *etaun, *etavn, *rhoqn, *rhoun, *rhovn;
extern PRECISION *etai, *etaq, *etau, *etav, *rhoq, *rhou, *rhov;
extern PRECISION *parcial1, *parcial2, *parcial3;
extern PRECISION *nubB, *nupB, *nurB;
PRECISION **uuGlobalInicial;
PRECISION **HGlobalInicial;
PRECISION **FGlobalInicial;
extern int FGlobal, HGlobal, uuGlobal;
extern PRECISION *GMAC;
extern Cuantic *cuantic;

extern PRECISION * opa;

extern _Complex PRECISION *z,* zden, * zdiv;

void InitializePointerShareCalculation()
{

	PUNTEROS_CALCULOS_COMPARTIDOS = calloc(LONG_PUNTERO_CALCULOS_COMPARTIDOS, sizeof(PRECISION *));

	ResetPointerShareCalculation();
}


void DeleteSpectraCalculation()
{

	int i;
	for (i = 0; i < POSW_PUNTERO_CALCULOS_COMPARTIDOS; i++)
	{
		free(PUNTEROS_CALCULOS_COMPARTIDOS[i]);
	}

	ResetPointerShareCalculation();
}

void FreePointerShareCalculation()
{

	free(PUNTEROS_CALCULOS_COMPARTIDOS);
}
void ResetPointerShareCalculation()
{

	POSR_PUNTERO_CALCULOS_COMPARTIDOS = 0;
	POSW_PUNTERO_CALCULOS_COMPARTIDOS = 0;
}

void AsignPointerShareCalculation(int Numero, PRECISION *a, ...)
{

	va_list Param;

	PUNTEROS_CALCULOS_COMPARTIDOS[POSW_PUNTERO_CALCULOS_COMPARTIDOS] = a;
	POSW_PUNTERO_CALCULOS_COMPARTIDOS += 1;
	Numero--;

	va_start(Param, a);
	while (Numero > 0)
	{

		PUNTEROS_CALCULOS_COMPARTIDOS[POSW_PUNTERO_CALCULOS_COMPARTIDOS] = (PRECISION *)va_arg(Param, PRECISION *);
		POSW_PUNTERO_CALCULOS_COMPARTIDOS += 1;
		Numero--;
	}

	va_end(Param);
}

void ReadPointerShareCalculation(int Numero, PRECISION **a, ...)
{

	va_list Param;
	PRECISION **aux;

	aux = NULL;
	*aux = NULL;
	

	*a = (PRECISION *)PUNTEROS_CALCULOS_COMPARTIDOS[POSR_PUNTERO_CALCULOS_COMPARTIDOS];
	POSR_PUNTERO_CALCULOS_COMPARTIDOS += 1;

	va_start(Param, a);
	while (--Numero > 0)
	{

		aux = va_arg(Param, PRECISION **);
		*aux = (PRECISION *)PUNTEROS_CALCULOS_COMPARTIDOS[POSR_PUNTERO_CALCULOS_COMPARTIDOS];
		PUNTEROS_CALCULOS_COMPARTIDOS[POSR_PUNTERO_CALCULOS_COMPARTIDOS] = NULL;
		POSR_PUNTERO_CALCULOS_COMPARTIDOS += 1;
	}

	va_end(Param);
}


void AllocateMemoryDerivedSynthesis(int numl)
{

	/************* ME DER *************************************/
	dtaux = calloc(numl,sizeof(PRECISION));
	etai_gp3 = calloc(numl,sizeof(PRECISION));
	ext1 = calloc(numl,sizeof(PRECISION));
	ext2 = calloc(numl,sizeof(PRECISION));
	ext3 = calloc(numl,sizeof(PRECISION));
	ext4 = calloc(numl,sizeof(PRECISION));
	/**********************************************************/

	//***** VARIABLES FOR FVOIGT ****************************//
	z = malloc (numl * sizeof(_Complex PRECISION));
	zden = malloc (numl * sizeof(_Complex PRECISION));
	zdiv = malloc (numl * sizeof(_Complex PRECISION));
	/********************************************************/


	GMAC = calloc(numl, sizeof(PRECISION));

	spectra = calloc(numl * NPARMS, sizeof(PRECISION));
	spectra_mac = calloc(numl * NPARMS, sizeof(PRECISION));
	d_spectra = calloc(numl * NTERMS * NPARMS, sizeof(PRECISION));
	
	
	opa = calloc(numl,sizeof(PRECISION));

	gp4_gp2_rhoq = calloc(numl, sizeof(PRECISION));
	gp5_gp2_rhou = calloc(numl, sizeof(PRECISION));
	gp6_gp2_rhov = calloc(numl, sizeof(PRECISION));

	gp1 = calloc(numl, sizeof(PRECISION));
	gp2 = calloc(numl, sizeof(PRECISION));
	gp3 = calloc(numl, sizeof(PRECISION));
	gp4 = calloc(numl, sizeof(PRECISION));
	gp5 = calloc(numl, sizeof(PRECISION));
	gp6 = calloc(numl, sizeof(PRECISION));
	dt = calloc(numl, sizeof(PRECISION));
	dti = calloc(numl, sizeof(PRECISION));

	etai_2 = calloc(numl, sizeof(PRECISION));

	dgp1 = calloc(numl, sizeof(PRECISION));
	dgp2 = calloc(numl, sizeof(PRECISION));
	dgp3 = calloc(numl, sizeof(PRECISION));
	dgp4 = calloc(numl, sizeof(PRECISION));
	dgp5 = calloc(numl, sizeof(PRECISION));
	dgp6 = calloc(numl, sizeof(PRECISION));
	d_dt = calloc(numl, sizeof(PRECISION));

	d_ei = calloc(numl * 7, sizeof(PRECISION));
	d_eq = calloc(numl * 7, sizeof(PRECISION));
	d_eu = calloc(numl * 7, sizeof(PRECISION));
	d_ev = calloc(numl * 7, sizeof(PRECISION));
	d_rq = calloc(numl * 7, sizeof(PRECISION));
	d_ru = calloc(numl * 7, sizeof(PRECISION));
	d_rv = calloc(numl * 7, sizeof(PRECISION));
	dfi = calloc(numl * 4 * 3, sizeof(PRECISION));  //DNULO
	dshi = calloc(numl * 4 * 3, sizeof(PRECISION)); //DNULO

	fi_p = calloc(numl * 2, sizeof(PRECISION));
	fi_b = calloc(numl * 2, sizeof(PRECISION));
	fi_r = calloc(numl * 2, sizeof(PRECISION));
	shi_p = calloc(numl * 2, sizeof(PRECISION));
	shi_b = calloc(numl * 2, sizeof(PRECISION));
	shi_r = calloc(numl * 2, sizeof(PRECISION));

	etain = calloc(numl * 2, sizeof(PRECISION));
	etaqn = calloc(numl * 2, sizeof(PRECISION));
	etaun = calloc(numl * 2, sizeof(PRECISION));
	etavn = calloc(numl * 2, sizeof(PRECISION));
	rhoqn = calloc(numl * 2, sizeof(PRECISION));
	rhoun = calloc(numl * 2, sizeof(PRECISION));
	rhovn = calloc(numl * 2, sizeof(PRECISION));

	etai = calloc(numl, sizeof(PRECISION));
	etaq = calloc(numl, sizeof(PRECISION));
	etau = calloc(numl, sizeof(PRECISION));
	etav = calloc(numl, sizeof(PRECISION));
	rhoq = calloc(numl, sizeof(PRECISION));
	rhou = calloc(numl, sizeof(PRECISION));
	rhov = calloc(numl, sizeof(PRECISION));

	parcial1 = calloc(numl, sizeof(PRECISION));
	parcial2 = calloc(numl, sizeof(PRECISION));
	parcial3 = calloc(numl, sizeof(PRECISION));

	nubB = calloc(cuantic[0].N_SIG, sizeof(PRECISION));
	nurB = calloc(cuantic[0].N_SIG, sizeof(PRECISION));
	nupB = calloc(cuantic[0].N_PI, sizeof(PRECISION));

	uuGlobalInicial = calloc((int)(cuantic[0].N_PI + cuantic[0].N_SIG * 2), sizeof(PRECISION *));
	uuGlobal = 0;
	int i = 0;
	for (i = 0; i < (int)(cuantic[0].N_PI + cuantic[0].N_SIG * 2); i++)
	{
		uuGlobalInicial[i] = calloc(numl, sizeof(PRECISION));
	}

	HGlobalInicial = calloc((int)(cuantic[0].N_PI + cuantic[0].N_SIG * 2), sizeof(PRECISION *));
	HGlobal = 0;
	for (i = 0; i < (int)(cuantic[0].N_PI + cuantic[0].N_SIG * 2); i++)
	{
		HGlobalInicial[i] = calloc(numl, sizeof(PRECISION));
	}

	FGlobalInicial = calloc((int)(cuantic[0].N_PI + cuantic[0].N_SIG * 2), sizeof(PRECISION *));
	for (i = 0; i < (int)(cuantic[0].N_PI + cuantic[0].N_SIG * 2); i++)
	{
		FGlobalInicial[i] = calloc(numl, sizeof(PRECISION));
	}
	FGlobal = 0;
}

void FreeMemoryDerivedSynthesis()
{
	
	free(dtaux);
	free(etai_gp3);
	free(ext1);
	free(ext2);
	free(ext3);
	free(ext4);	

	free(zden);
	free(zdiv);
	free(z);

	free(gp1);
	free(gp2);
	free(gp3);
	free(gp4);
	free(gp5);
	free(gp6);
	free(dt);
	free(dti);

	free(etai_2);

	free(dgp1);
	free(dgp2);
	free(dgp3);
	free(dgp4);
	free(dgp5);
	free(dgp6);
	free(d_dt);

	free(d_ei);
	free(d_eq);
	free(d_ev);
	free(d_eu);
	free(d_rq);
	free(d_ru);
	free(d_rv);

	free(dfi);
	free(dshi);

	free(GMAC);
	free(opa);
	free(spectra);
	free(spectra_mac);
	free(d_spectra);
	
	

	free(fi_p);
	free(fi_b);
	free(fi_r);
	free(shi_p);
	free(shi_b);
	free(shi_r);

	free(etain);
	free(etaqn);
	free(etaun);
	free(etavn);
	free(rhoqn);
	free(rhoun);
	free(rhovn);

	free(etai);
	free(etaq);
	free(etau);
	free(etav);

	free(rhoq);
	free(rhou);
	free(rhov);

	free(parcial1);
	free(parcial2);
	free(parcial3);

	free(nubB);
	free(nurB);
	free(nupB);

	free(gp4_gp2_rhoq);
	free(gp5_gp2_rhou);
	free(gp6_gp2_rhov);

	int i;
	for (i = 0; i < (int)(cuantic[0].N_PI + cuantic[0].N_SIG * 2); i++)
	{
		free(uuGlobalInicial[i]);
	}

	for (i = 0; i < (int)(cuantic[0].N_PI + cuantic[0].N_SIG * 2); i++)
	{
		free(HGlobalInicial[i]);
	}

	for (i = 0; i < (int)(cuantic[0].N_PI + cuantic[0].N_SIG * 2); i++)
	{
		free(FGlobalInicial[i]);
	}

	free(uuGlobalInicial);
	free(HGlobalInicial);
	free(FGlobalInicial);
}
