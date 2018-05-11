/***************************************************************************
 *
 * Authors:    Amaya Jimenez            ajimenez@cnb.csic.es (2018)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <data/args.h>
#include <data/morphology.h>
#include <data/filters.h>
#include "volume_gain.h"

// Read arguments ==========================================================
void ProgVolumeGain::readParams()
{

    fn_vol = getParam("-i");
    fn_mask = getParam("--mask");

}

// Show ====================================================================
void ProgVolumeGain::show() const
{
    std::cout
    << "Input file   : " << fn_vol        << std::endl
	<< "Input mask   : " << fn_mask        << std::endl
    ;
}

// usage ===================================================================
void ProgVolumeGain::defineParams()
{
    addParamsLine("   -i <volume>                            : Volume to calculate the gain");
    addParamsLine("   [--mask <vol_file=\"\">]               : Mask defining the macromolecule");
}

// Produce side information ================================================
void ProgVolumeGain::produce_side_info()
{
	std::cout << "Producing side info..." << std::endl;

	V.read(fn_vol);
    V().setXmippOrigin();
    mask.read(fn_mask);
    mask().setXmippOrigin();

	FourierTransformer transformer;
	MultidimArray<double> &inputVol = V();
	VRiesz.resizeNoCopy(inputVol);

	transformer.FourierTransform(inputVol, fftV);
	iu.initZeros(fftV);

	// Calculate u (vector norm of frequency directions)
	double uz, uy, ux, uz2, u2, uz2y2;
	long n=0;
	for(size_t k=0; k<ZSIZE(fftV); ++k)
	{
		FFT_IDX2DIGFREQ(k,ZSIZE(inputVol),uz);
		uz2=uz*uz;

		for(size_t i=0; i<YSIZE(fftV); ++i)
		{
			FFT_IDX2DIGFREQ(i,YSIZE(inputVol),uy);
			uz2y2=uz2+uy*uy;

			for(size_t j=0; j<XSIZE(fftV); ++j)
			{
				FFT_IDX2DIGFREQ(j,XSIZE(inputVol),ux);
				u2=uz2y2+ux*ux;
				if ((k != 0) || (i != 0) || (j != 0))
					DIRECT_MULTIDIM_ELEM(iu,n) = 1.0/sqrt(u2);
				else
					DIRECT_MULTIDIM_ELEM(iu,n) = 1e38;
				++n;
			}
		}
	}

	// Prepare low pass filter
	lowPassFilter.FilterShape = RAISED_COSINE;
	lowPassFilter.raised_w = 0.01;
	lowPassFilter.do_generate_3dmask = false;
	lowPassFilter.FilterBand = LOWPASS;

	// Prepare mask
	MultidimArray<int> &pMask=mask();
	size_t R = (V().xdim)/2.0;

	NVoxelsOriginalMask = 0;
	FOR_ALL_ELEMENTS_IN_ARRAY3D(pMask)
	{
		if (A3D_ELEM(pMask, k, i, j) == 1)
			NVoxelsOriginalMask++;
		if (i*i+j*j+k*k > R*R)
			A3D_ELEM(pMask, k, i, j) = 0;
	}

	double u;
	int size_fourier = ZSIZE(fftV);
	freq_fourier.initZeros(size_fourier);

	int size = ZSIZE(pMask);

	VEC_ELEM(freq_fourier,0) = 1e-38;
	for(size_t k=0; k<size_fourier; ++k)
	{
		FFT_IDX2DIGFREQ(k,size, u);
		VEC_ELEM(freq_fourier,k) = u;
	}
}



void ProgVolumeGain::amplitudeMonogenicSignal3D(MultidimArray< std::complex<double> > &myfftV,
		MultidimArray<double> &amplitude)
{
	std::cout << "Calculating monogenic amplitude..." << std::endl;

	fftVRiesz.initZeros(myfftV);
	amplitude.resizeNoCopy(VRiesz);
	fftVRiesz_aux.initZeros(myfftV);
	std::complex<double> J(0,1);

	// Filter the input volume and add it to amplitude
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(myfftV)
	{
		double iun=DIRECT_MULTIDIM_ELEM(iu,n);
		double un=1.0/iun;
		DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
		DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = -J;
		DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= DIRECT_MULTIDIM_ELEM(fftVRiesz, n);
		DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= iun;
	}

	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);


	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n)=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate first component of Riesz vector
	fftVRiesz.initZeros(myfftV);
	double uz, uy, ux;
	long n=0;

	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				ux = VEC_ELEM(freq_fourier,j);
				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = ux*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate second component of Riesz vector
	fftVRiesz.initZeros(myfftV);
	n=0;

	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			uy = VEC_ELEM(freq_fourier,i);
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = uy*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate third component of Riesz vector
	fftVRiesz.initZeros(myfftV);
	n=0;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		uz = VEC_ELEM(freq_fourier,k);
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = uz*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	{
		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
		DIRECT_MULTIDIM_ELEM(amplitude,n)=sqrt(DIRECT_MULTIDIM_ELEM(amplitude,n));
	}

	/*// Low pass filter the monogenic amplitude
	double freq = 1.0;
	double aux_frequency;
	int fourier_idx;
	DIGFREQ2FFT_IDX(freq, ZSIZE(VRiesz), fourier_idx);
	FFT_IDX2DIGFREQ(fourier_idx, ZSIZE(VRiesz), aux_frequency);
	freq = aux_frequency;

	lowPassFilter.w1 = freq;
	amplitude.setXmippOrigin();
	lowPassFilter.applyMaskSpace(amplitude);*/

}

void ProgVolumeGain::calculateGlobalHistogram(MultidimArray<double> amplitude, MultidimArray<double> &histogram,
		MultidimArray<double> &cdfGlobal, MultidimArray<int> *pMask, int Nbins, double &step)
{
	/*double max = amplitude.computeMax();
	double min = amplitude.computeMin();

	std::cout << "NO MASKED Max = " << max << " Min = " << min << std::endl;*/

	/*FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude){
		DIRECT_MULTIDIM_ELEM(amplitude, n) *= DIRECT_MULTIDIM_ELEM(*pMask, n);
	}*/

	double max = amplitude.computeMax();
	double min = amplitude.computeMin();
	step = (max-min)/Nbins;
	std::cout << "Max = " << max << " Min = " << min << std::endl;

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	{
		if (DIRECT_MULTIDIM_ELEM(*pMask, n)==1)
		{
			int position = (int)floor(DIRECT_MULTIDIM_ELEM(amplitude, n)/step);
			DIRECT_MULTIDIM_ELEM(histogram, position)+=1;
		}
	}
	int histSum=histogram.sum();
	if (histSum!=0){
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(histogram){
			DIRECT_MULTIDIM_ELEM(histogram, n) /= histSum;
			//CDF global
			if(n==0)
				DIRECT_MULTIDIM_ELEM(cdfGlobal, n) = DIRECT_MULTIDIM_ELEM(histogram, n);
			else
				DIRECT_MULTIDIM_ELEM(cdfGlobal, n) = DIRECT_MULTIDIM_ELEM(histogram, n) + DIRECT_MULTIDIM_ELEM(cdfGlobal, n-1);
		}
	}

	//FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(histogram)
	//	std::cout << "GLOBAL: In " << n << " hist(n)= " << DIRECT_MULTIDIM_ELEM(histogram,n) << " cdf(n)= " << DIRECT_MULTIDIM_ELEM(cdfGlobal,n) << std::endl;

}


void ProgVolumeGain::matchingLocalHistogram(MultidimArray<double> amplitude, MultidimArray<double> &gainOut,
		std::vector<MultidimArray<double> > &cdfsLocal, MultidimArray<double> cdfGlobal,
		MultidimArray<int> *pMask, int Nbins, double step, int boxSize)
{

	MultidimArray<int> mask_aux;
	mask_aux.resize((*pMask));
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY((*pMask)){
		DIRECT_MULTIDIM_ELEM(mask_aux, n) = DIRECT_MULTIDIM_ELEM(*pMask, n);
	}

	MultidimArray<double> histogramAux, cdfAux;
	histogramAux.initZeros(Nbins);
	cdfAux.initZeros(Nbins);

	long p=0;
	long pp;

	int xdim= XSIZE(amplitude);
	int ydim= YSIZE(amplitude);
	int zdim= ZSIZE(amplitude);
	int total= ZYXSIZE(amplitude);
	//std::cout << "xdim = " << xdim << " ydim =  " << ydim << " zdim = " << zdim << " total = " << total << std::endl;
	for(size_t k=0; k<zdim; ++k){
		for(size_t i=0; i<ydim; ++i){
			for(size_t j=0; j<xdim; ++j){
				p = j+(i*xdim)+(k*xdim*ydim);
				if (DIRECT_MULTIDIM_ELEM(mask_aux, p)==0)
					continue;

				for (int kk=-boxSize; kk<=boxSize; kk++){
					for (int ii=-boxSize; ii<=boxSize; ii++){
						for (int jj=-boxSize; jj<=boxSize; jj++){
							pp = p + jj + ii*xdim + kk*xdim*ydim;
							//std::cout << "nn = " << nn << " n = " << n << " kk = " << kk << " ii = " << ii << " jj = " << jj << " ii*xdim = " << ii*xdim << " kk*xdim*ydim = " << kk*xdim*ydim << std::endl;
							if(pp<0 || pp>total)
								continue;

							if (DIRECT_MULTIDIM_ELEM(mask_aux, pp)==1){
								int position = (int)floor(DIRECT_MULTIDIM_ELEM(amplitude, pp)/step);
								DIRECT_MULTIDIM_ELEM(histogramAux, position)+=1;
								DIRECT_MULTIDIM_ELEM(mask_aux, pp)=0;
							}

						}
					}
				}
				int histSum=histogramAux.sum();
				if (histSum!=0){
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(histogramAux){
						DIRECT_MULTIDIM_ELEM(histogramAux, n) /= histSum;
						if(n==0)
							DIRECT_MULTIDIM_ELEM(cdfAux, n) = DIRECT_MULTIDIM_ELEM(histogramAux, n);
						else
							DIRECT_MULTIDIM_ELEM(cdfAux, n) = DIRECT_MULTIDIM_ELEM(histogramAux, n) + DIRECT_MULTIDIM_ELEM(cdfAux, n-1);
					}

					//FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(histogramAux)
					//	std::cout << "LOCAL: In " << n << " hist(n)= " << DIRECT_MULTIDIM_ELEM(histogramAux,n) << " cdf(n)= " << DIRECT_MULTIDIM_ELEM(cdfAux,n) << std::endl;

					//Matching
					for (int kk=-boxSize; kk<=boxSize; kk++){
						for (int ii=-boxSize; ii<=boxSize; ii++){
							for (int jj=-boxSize; jj<=boxSize; jj++){
								pp = p + jj + ii*xdim + kk*xdim*ydim;
								if(pp<0 || pp>total)
									continue;

								if (DIRECT_MULTIDIM_ELEM(*pMask, pp)==1){
									int position = (int)floor(DIRECT_MULTIDIM_ELEM(amplitude, pp)/step);
									double probLocal = DIRECT_MULTIDIM_ELEM(cdfAux, position);
									double diffProb=100000, diffAmp=100000;
									double newAmplitude;
									//std::cout << "LOCAL: In " << pp << std::endl;
									//std::cout << "LOCAL: DIRECT_MULTIDIM_ELEM(amplitude, pp) " << DIRECT_MULTIDIM_ELEM(amplitude, pp) << std::endl;
									//std::cout << "LOCAL: position " << position << std::endl;
									FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(cdfGlobal){
										//std::cout << "LOCAL: DIRECT_MULTIDIM_ELEM(cdfGlobal,n) " << DIRECT_MULTIDIM_ELEM(cdfGlobal,n) << std::endl;
										//std::cout << "LOCAL: probLocal " << probLocal << std::endl;
										if( fabs(DIRECT_MULTIDIM_ELEM(cdfGlobal,n)-probLocal)<=diffProb){
											diffProb = fabs(DIRECT_MULTIDIM_ELEM(cdfGlobal,n)-probLocal);
											//std::cout << "LOCAL: diffProb " << diffProb << std::endl;
											if( fabs(DIRECT_MULTIDIM_ELEM(amplitude, pp)-(step*n))< diffAmp){
												diffAmp = fabs(DIRECT_MULTIDIM_ELEM(amplitude, pp)-(step*n));
												//std::cout << "LOCAL: diffAmp " << diffAmp << std::endl;
												newAmplitude = step*n;
												//std::cout << "LOCAL: newAmplitude " << newAmplitude << std::endl;
												//std::cout << "LOCAL: step " << step << std::endl;
												//std::cout << "LOCAL: n " << n << std::endl;
											}
										}
									}
									DIRECT_MULTIDIM_ELEM(gainOut, pp) = newAmplitude/DIRECT_MULTIDIM_ELEM(amplitude,pp);
									//std::cout << "LOCAL: In " << pp << " position " << position << " step " << step << " prev_amp =  " << DIRECT_MULTIDIM_ELEM(amplitude, pp) << " new_amp= " << newAmplitude << std::endl;
								}


							}
						}
					}
					//exit(1);
					//cdfsLocal.push_back(cdfAux);
				}
				histogramAux.initZeros(Nbins);
				cdfAux.initZeros(Nbins);
			}
		}
	}

}


void ProgVolumeGain::run()
{
	MultidimArray<double> amplitude;
	MultidimArray<double> gainOut;
	int Nbins = 20;
	int boxSize = 20;

    produce_side_info();
    amplitudeMonogenicSignal3D(fftV, amplitude);
    Image<double> saveImg;
	saveImg = amplitude;
	saveImg.write("./amplitudeMono.vol");
	saveImg.clear();

	Image<int> saveMask;
	saveMask = mask();
	saveMask.write("./mask.vol");
	saveMask.clear();

	MultidimArray<double> histogramGlobal, cdfGlobal;
	histogramGlobal.initZeros(Nbins);
	cdfGlobal.initZeros(Nbins);
	MultidimArray<int> *pMask = &(mask());
	double step;

	calculateGlobalHistogram(amplitude, histogramGlobal, cdfGlobal, pMask, Nbins, step);

	std::vector<MultidimArray<double> > cdfsLocal;
	gainOut.initZeros(amplitude);
	matchingLocalHistogram(amplitude, gainOut, cdfsLocal, cdfGlobal, pMask, Nbins, step, boxSize);

	saveImg = gainOut;
	saveImg.write("./gainOut.vol");
	saveImg.clear();

	//Multiply volume and gain
	V() *= gainOut;
	saveImg = V();
	saveImg.write("./outputVol.vol");
	saveImg.clear();


	/*/AJ to check if the histogramsLocal contains all the information of the volume
	std::cout << "histogramLocal size = " << howManyHist << std::endl;
	MultidimArray<double> *aux;
	MultidimArray<double> total;
	total.initZeros(Nbins);
	for (int i=0; i<howManyHist; i++){
		aux = &histogramsLocal[i];
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(*aux){
			DIRECT_MULTIDIM_ELEM(total, n)+=DIRECT_MULTIDIM_ELEM(*aux, n);
		}
	}
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(total)
		std::cout << "Histogram total in " << n << " = " << total(n) << std::endl;
	*/



}
