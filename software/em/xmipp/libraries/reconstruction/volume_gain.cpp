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
    mono = checkParam("-mono");
    if (mono)
    	fn_mono = getParam("-mono");
    fn_mask = getParam("--mask");
    sampling = getDoubleParam("--sampling");
    boxSize = getIntParam("--boxSize");
    bandpass = checkParam("--bandpass");
    Nbands = getIntParam("--bandpass");
    iter = getIntParam("--iter");
    superposed = checkParam("--superposed");
    if(checkParam("-o"))
    	fn_out = getParam("-o");
    else
    	fn_out="./outputVol.vol";

}

// Show ====================================================================
void ProgVolumeGain::show() const
{
    std::cout
    << "Input file            : " << fn_vol        << std::endl
	<< "Input mask            : " << fn_mask        << std::endl
	<< "Box size              : " << boxSize        << std::endl
	<< "Number of iterations  : " << iter        << std::endl
    ;
}

// usage ===================================================================
void ProgVolumeGain::defineParams()
{
    addParamsLine("   -i <volume>                        : Volume to match the amplitudes");
    addParamsLine("   [-mono <monoResVolume>]            : MonoRes volume with the resolution in every voxel");
    addParamsLine("   [-o <output=\"\">]                 : Output volume filename");
	addParamsLine("   [--sampling <s=1>]                 : Sampling rate (A/px)");
    addParamsLine("   [--mask <mask=\"\">]               : Mask defining the macromolecule");
    addParamsLine("   [--boxSize <boxSize=5>]            : Size of the box in pixels per coordinate to calculate the histogram");
    addParamsLine("   [--bandpass <Nbands=5>]            : Carry out the matching in the whole frequency range or by bands. The integer value wiil be the number of band pass filter to apply");
    addParamsLine("   [--iter <iter=1>]                  : Number of iterations");
    addParamsLine("   [--superposed]                     : To allow superposed voxels in the histogram matching calculation");
}

// Produce side information ================================================
void ProgVolumeGain::produce_side_info()
{

	V.read(fn_vol);
    V().setXmippOrigin();
    mask.read(fn_mask);
    mask().setXmippOrigin();
    if (mono){
    	monoRes.read(fn_mono);
    	monoRes().setXmippOrigin();
    }

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
}


void ProgVolumeGain::calculateFFT()
{
	std::cout << "Calculating FFT..." << std::endl;

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
		FFT_IDX2DIGFREQ(k,ZSIZE(fftV),uz);
		uz2=uz*uz;

		for(size_t i=0; i<YSIZE(fftV); ++i)
		{
			FFT_IDX2DIGFREQ(i,YSIZE(fftV),uy);
			uz2y2=uz2+uy*uy;

			for(size_t j=0; j<XSIZE(fftV); ++j)
			{
				FFT_IDX2DIGFREQ(j,XSIZE(fftV),ux);
				u2=uz2y2+ux*ux;
				if ((k != 0) || (i != 0) || (j != 0))
					DIRECT_MULTIDIM_ELEM(iu,n) = 1.0/sqrt(u2);
				else
					DIRECT_MULTIDIM_ELEM(iu,n) = 1e38;
				++n;
			}
		}
	}

	double u;
	int size_fourier = ZSIZE(fftV);
	freq_fourier.initZeros(size_fourier);

	int size = ZSIZE(mask());

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

//	// Low pass filter the monogenic amplitude
//	double freq = 1.5;
//	double aux_frequency;
//	int fourier_idx;
//	DIGFREQ2FFT_IDX(freq, ZSIZE(amplitude), fourier_idx);
//	FFT_IDX2DIGFREQ(fourier_idx, ZSIZE(amplitude), aux_frequency);
//	std::cout << "Low pass filtering at " << freq << " " << aux_frequency << std::endl;
//	freq = aux_frequency;
//	lowPassFilter.w1 = freq;
//	amplitude.setXmippOrigin();
//	lowPassFilter.applyMaskSpace(amplitude);

}

void ProgVolumeGain::calculateGlobalHistogram(MultidimArray<double> amplitude, MultidimArray<double> &histogram,
		MultidimArray<double> &cdfGlobal, MultidimArray<int> *pMask, int Nbins, double &step)
{

	std::cout << "Calculating global histogram... " << std::endl;
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
	//AJ CUIDADO: que salen aqui valores negativos en el min a veces

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	{
		if (DIRECT_MULTIDIM_ELEM(*pMask, n)==1)
		{
			int position = (int)floor(DIRECT_MULTIDIM_ELEM(amplitude, n)/step);
			//if (position>=Nbins)
			//	std::cout << "ERROOOOORRRR " << position << " " << DIRECT_MULTIDIM_ELEM(amplitude, n) << " " << step << std::endl;
			if (position>=Nbins)
				position=Nbins-1;
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

	//FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(histogram){
		//std::cout << n << " " << DIRECT_MULTIDIM_ELEM(histogram,n) << " " << DIRECT_MULTIDIM_ELEM(cdfGlobal,n) << std::endl;
		//std::cout << "GLOBAL: In " << n << " hist(n)= " << DIRECT_MULTIDIM_ELEM(histogram,n) << " cdf(n)= " << DIRECT_MULTIDIM_ELEM(cdfGlobal,n) << std::endl;
	//}

}


void ProgVolumeGain::matchingLocalHistogram(MultidimArray<double> amplitude, MultidimArray<double> &gainOut,
		MultidimArray<double> cdfGlobal, MultidimArray<int> *pMask, int Nbins, double step, int boxSize)
{

	std::cout << "Matching local histograms... " << std::endl;
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
				//std::cout << "START MATCH INI " << p << std::endl;
				if (DIRECT_MULTIDIM_ELEM(mask_aux, p)==0)
					continue;

				//std::cout << "START MATCH " << p << std::endl;
				for (int kk=-boxSize; kk<=boxSize; kk++){
					for (int ii=-boxSize; ii<=boxSize; ii++){
						for (int jj=-boxSize; jj<=boxSize; jj++){
							pp = p + jj + ii*xdim + kk*xdim*ydim;
							//std::cout << "pp = " << pp << " p = " << p << " kk = " << kk << " ii = " << ii << " jj = " << jj << " ii*xdim = " << ii*xdim << " kk*xdim*ydim = " << kk*xdim*ydim << std::endl;
							if(pp<0 || pp>total)
								continue;

							if (DIRECT_MULTIDIM_ELEM(mask_aux, pp)==1){
								int position = (int)floor(DIRECT_MULTIDIM_ELEM(amplitude, pp)/step);
								if (position>=Nbins)
									position=Nbins-1;
								DIRECT_MULTIDIM_ELEM(histogramAux, position)+=1;
								DIRECT_MULTIDIM_ELEM(mask_aux, pp)=2;
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

					//FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(histogramAux){
						//std::cout << "LOCAL: In " << n << " hist(n)= " << DIRECT_MULTIDIM_ELEM(histogramAux,n) << " cdf(n)= " << DIRECT_MULTIDIM_ELEM(cdfAux,n) << std::endl;
						//std::cout << n << " " << DIRECT_MULTIDIM_ELEM(histogramAux,n) << " " << DIRECT_MULTIDIM_ELEM(cdfAux,n) << " " << DIRECT_MULTIDIM_ELEM(histNew, n) << std::endl;
					//}

					//Matching
					for (int kk=-boxSize; kk<=boxSize; kk++){
						for (int ii=-boxSize; ii<=boxSize; ii++){
							for (int jj=-boxSize; jj<=boxSize; jj++){
								pp = p + jj + ii*xdim + kk*xdim*ydim;
								if(pp<0 || pp>total)
									continue;

								if (DIRECT_MULTIDIM_ELEM(mask_aux, pp)==2){
									DIRECT_MULTIDIM_ELEM(mask_aux, pp)=0;
									int position = (int)floor(DIRECT_MULTIDIM_ELEM(amplitude, pp)/step);
									if (position==Nbins)
										position-=1;
									double probLocal = DIRECT_MULTIDIM_ELEM(cdfAux, position);
									double diffProb=100000, diffAmp=100000;
									double newAmplitude;
									//std::cout << "MATCH: In " << pp << std::endl;
									//std::cout << "MATCH: DIRECT_MULTIDIM_ELEM(amplitude, pp) " << DIRECT_MULTIDIM_ELEM(amplitude, pp) << std::endl;
									//std::cout << "MATCH: position " << position << std::endl;
									FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(cdfGlobal){
										//std::cout << "MATCH: DIRECT_MULTIDIM_ELEM(cdfGlobal,n) " << DIRECT_MULTIDIM_ELEM(cdfGlobal,n) << std::endl;
										//std::cout << "MATCH: probLocal " << probLocal << std::endl;
										if( fabs(DIRECT_MULTIDIM_ELEM(cdfGlobal,n)-probLocal)<=diffProb){
											diffProb = fabs(DIRECT_MULTIDIM_ELEM(cdfGlobal,n)-probLocal);
											//std::cout << "MATCH: diffProb " << diffProb << std::endl;
											if( fabs(DIRECT_MULTIDIM_ELEM(amplitude, pp)-(step*n))< diffAmp){
												diffAmp = fabs(DIRECT_MULTIDIM_ELEM(amplitude, pp)-(step*n));
												//std::cout << "MATCH: diffAmp " << diffAmp << std::endl;
												newAmplitude = step*n;
												//std::cout << "MATCH: newAmplitude " << newAmplitude << std::endl;
												//std::cout << "MATCH: step " << step << std::endl;
												//std::cout << "MATCH: n " << n << std::endl;
											}
										}
									}

									DIRECT_MULTIDIM_ELEM(gainOut, pp) = newAmplitude/DIRECT_MULTIDIM_ELEM(amplitude,pp);
									//std::cout << "MATCH: In " << pp << " position " << position << " step " << step << " prev_amp =  " << DIRECT_MULTIDIM_ELEM(amplitude, pp) << " new_amp= " << newAmplitude << " gain= " << DIRECT_MULTIDIM_ELEM(gainOut, pp) << std::endl;

								}

							}
						}
					}

					//FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(histogramAux){
						//std::cout << "LOCAL: In " << n << " hist(n)= " << DIRECT_MULTIDIM_ELEM(histogramAux,n) << " cdf(n)= " << DIRECT_MULTIDIM_ELEM(cdfAux,n) << std::endl;
						//std::cout << n << " " << DIRECT_MULTIDIM_ELEM(histogramAux,n) << " " << DIRECT_MULTIDIM_ELEM(cdfAux,n) << " " << DIRECT_MULTIDIM_ELEM(histNew, n) << std::endl;
					//}
					//exit(1);
					//cdfsLocal.push_back(cdfAux);
				}
				//std::cout << "FIN MATCH " << p << std::endl;
				histogramAux.initZeros(Nbins);
				cdfAux.initZeros(Nbins);
				//std::cout << "FIN MATCH 2 " << p << std::endl;

			}
		}
	}

}



void ProgVolumeGain::matchingLocalHistogram_new(MultidimArray<double> amplitude, MultidimArray<double> &gainOut,
		std::vector< double > globalAmplitudes, MultidimArray<int> *pMask, int boxSize, double freq)
{

	std::cout << "Matching local histograms... " << std::endl;
	MultidimArray<int> mask_aux;
	mask_aux.resize((*pMask));
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY((*pMask)){
		DIRECT_MULTIDIM_ELEM(mask_aux, n) = DIRECT_MULTIDIM_ELEM(*pMask, n);
	}

	long p=0;
	long pp;

	int xdim= XSIZE(amplitude);
	int ydim= YSIZE(amplitude);
	int zdim= ZSIZE(amplitude);
	int total= ZYXSIZE(amplitude);
	int lenGlobal= globalAmplitudes.size();

	//Calculate sort vector of local amplitudes values to obtain the cdf
	std::vector< double > localAmplitudes;

	//To allow re-visiting voxels
	MultidimArray<double> voxelCount;
	if(superposed){
		voxelCount.initZeros(amplitude);
	}

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

							if(pp<0 || pp>total)
								continue;

							if (DIRECT_MULTIDIM_ELEM(mask_aux, pp)==1){
								if(!mono || (mono && (sampling/DIRECT_MULTIDIM_ELEM(monoRes(), pp)) > freq)){
									localAmplitudes.push_back(DIRECT_MULTIDIM_ELEM(amplitude, pp));
									DIRECT_MULTIDIM_ELEM(mask_aux, pp)=2;
								}
							}

						}
					}
				}

				if (localAmplitudes.size()!=0){

					std::sort(localAmplitudes.begin(), localAmplitudes.end());

					//Matching
					for (int kk=-boxSize; kk<=boxSize; kk++){
						for (int ii=-boxSize; ii<=boxSize; ii++){
							for (int jj=-boxSize; jj<=boxSize; jj++){
								pp = p + jj + ii*xdim + kk*xdim*ydim;
								if(pp<0 || pp>total)
									continue;

								if (DIRECT_MULTIDIM_ELEM(mask_aux, pp)==2){
									if(superposed)
										DIRECT_MULTIDIM_ELEM(mask_aux, pp)=1;
									else
										DIRECT_MULTIDIM_ELEM(mask_aux, pp)=0;
									double valueLocal = DIRECT_MULTIDIM_ELEM(amplitude, pp);
									int position;
									for (int a=0; a<localAmplitudes.size(); a++){
										//std::cout << "MATCH: valueLocal " << valueLocal << " " << localAmplitudes[a] << std::endl;
										if (localAmplitudes[a]>valueLocal){
											position=a-1;
											break;
										}
									}
									double probLocal = (double)position/(double)localAmplitudes.size();
									double newAmplitude;
									//std::cout << "MATCH: In " << pp << std::endl;
									//std::cout << "MATCH: DIRECT_MULTIDIM_ELEM(amplitude, pp) " << DIRECT_MULTIDIM_ELEM(amplitude, pp) << std::endl;
									//std::cout << "MATCH: probLocal " << probLocal << std::endl;


									//std::cout << "MATCH: lenGlobal " << lenGlobal << std::endl;

									int posGlobal = (int)(probLocal*lenGlobal);
									//std::cout << "MATCH: posGlobal " << posGlobal << std::endl;

									if (posGlobal>=lenGlobal)
										posGlobal=lenGlobal-1;
									newAmplitude = globalAmplitudes[posGlobal];

									//std::cout << "MATCH: probGlobal " << (double)posGlobal/(double)lenGlobal << " " << (double)(posGlobal-1)/(double)lenGlobal << " " << (double)(posGlobal+1)/(double)lenGlobal << std::endl;
									//std::cout << "MATCH: newAmplitude " << newAmplitude << " " << globalAmplitudes[posGlobal-1] << " " << globalAmplitudes[posGlobal+1] << std::endl;

									if(superposed){
										DIRECT_MULTIDIM_ELEM(gainOut, pp) += newAmplitude;
										DIRECT_MULTIDIM_ELEM(voxelCount, pp) += 1;
									}else{
										DIRECT_MULTIDIM_ELEM(gainOut, pp) = newAmplitude;
									}
								}


							}
						}
					}

				}
				localAmplitudes.clear();

			}
		}
	}

	if(superposed){
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(gainOut){
			if (DIRECT_MULTIDIM_ELEM(voxelCount, n)!=0)
				DIRECT_MULTIDIM_ELEM(gainOut, n) = DIRECT_MULTIDIM_ELEM(gainOut, n)/DIRECT_MULTIDIM_ELEM(voxelCount, n);
		}
	}


}

void ProgVolumeGain::run_before()
{

	MultidimArray<double> amplitude, gainOut, histogramGlobal, cdfGlobal;
	MultidimArray<int> *pMask;

	int nBins=100;
	int iter=5;

	std::cout << "nBins " << nBins << std::endl;
	std::cout << "boxSize " << boxSize << std::endl;

    produce_side_info();
    pMask = &(mask());

    for (int i=0; i<iter; i++){
    	std::cout << "Iter " << i << std::endl;


		calculateFFT();
		amplitudeMonogenicSignal3D(fftV, amplitude);
		Image<double> saveImg;
		FileName name;
		name = formatString("./amplitudeMono%i.vol", i);
		saveImg = amplitude;
		saveImg.write(name);
		//saveImg.clear();

		Image<int> saveMask;
		name = formatString("./mask%i.vol", i);
		saveMask = mask();
		saveMask.write(name);
		saveMask.clear();

		histogramGlobal.initZeros(nBins);
		cdfGlobal.initZeros(nBins);
		double step;

		calculateGlobalHistogram(amplitude, histogramGlobal, cdfGlobal, pMask, nBins, step);

		gainOut.initZeros(amplitude);
		matchingLocalHistogram(amplitude, gainOut, cdfGlobal, pMask, nBins, step, boxSize);

		name = formatString("./gainOut%i.vol", i);
		saveImg = gainOut;
		saveImg.write(name);
		//saveImg.clear();

		//IDEAS: filtrar gain, bien calculada??

		/*// Low pass filter the gain
		double freq = 1.5;
		double aux_frequency;
		int fourier_idx;
		DIGFREQ2FFT_IDX(freq, ZSIZE(gainOut), fourier_idx);
		FFT_IDX2DIGFREQ(fourier_idx, ZSIZE(gainOut), aux_frequency);
		std::cout << "Low pass filtering at " << freq << " " << aux_frequency << std::endl;
		freq = aux_frequency;
		lowPassFilter.w1 = freq;
		gainOut.setXmippOrigin();
		lowPassFilter.applyMaskSpace(gainOut);

		name = formatString("./gainOutFilter%i.vol", i);
		saveImg = gainOut;
		saveImg.write(name);*/

		//Multiply volume and gain
		V() *= gainOut;

		/*//Low pass filtering the amplitude values
		double alpha=0.8;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(gainOut){
			DIRECT_MULTIDIM_ELEM(V(),n) = alpha*DIRECT_MULTIDIM_ELEM(V(),n) + (1-alpha)*DIRECT_MULTIDIM_ELEM(gainOut,n);
		}*/

		name = formatString("./outputVol%i.vol", i);
		saveImg = V();
		saveImg.write(name);
		saveImg.clear();

    }


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

void ProgVolumeGain::processing (MultidimArray<double> &V, MultidimArray<int> *pMask, double freq)
{

	MultidimArray<double> gainOut;
	std::vector< double > globalAmplitudes;
	double max = V.computeMax();
	double min = V.computeMin();
	std::cout << "Max = " << max << " Min = " << min << std::endl;

	//Calculate sort vector of global amplitudes values to obtain the cdf
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(V){
		if (DIRECT_MULTIDIM_ELEM(mask(), n)==1)
		{
			if(!mono || (mono && (sampling/DIRECT_MULTIDIM_ELEM(monoRes(), n)) > freq))
				globalAmplitudes.push_back(DIRECT_MULTIDIM_ELEM(V, n));
		}
	}
	std::sort(globalAmplitudes.begin(), globalAmplitudes.end());

	//Matching histograms
	gainOut.initZeros(V);
	matchingLocalHistogram_new(V, gainOut, globalAmplitudes, pMask, boxSize, freq);

	V=gainOut;

	//End
	globalAmplitudes.clear();

}


void ProgVolumeGain::run()
{

	MultidimArray<int> *pMask;
	produce_side_info();
	pMask = &(mask());


	Image<double> saveImg;
	FileName name;
	FileName nameRoot = fn_out.getRoot();

	if(bandpass){

		std::cout << "Selected band pass filtering with " << Nbands << " bands" << std::endl;

		MultidimArray<double> Vout;
		Vout.initZeros(V());

		MultidimArray<std::complex<double> > fftV, fftVaux;
		FourierTransformer transformer, transformer_inv;
		transformer.FourierTransform(V(),fftV);

		int filter_num = Nbands;
		double step = 0.5/filter_num;

		for(int j=0; j<iter; j++){
			std::cout << "Iteration " << j << std::endl;

			for (int i=0;i<filter_num;i++)
			{
				//if (i>0 && j==0){ //Primera iteracion, bandas de freq mayores que la primera
				//	V.read(fn_vol);
				//	V().setXmippOrigin();
				//}else if(j>0){ //Iteraciones mayor que 1, todas las bandas de freq
				if(j>0 && i==0){ //Iteraciones mayor que 1, primera banda de freq, rehacer la fft con el volumen de salida de la iteracion anterior
					//name = formatString("./outputVol_iter%i_new.vol",j-1);
					name.compose(nameRoot, j, "vol");
					V.read(name);
					V().setXmippOrigin();
					transformer.FourierTransform(V(),fftV);
				}

				fftVaux.initZeros(fftV);
				double w1=step*i;
				double w2=w1+step;
				double w12=w1*w1;
				double w22=w2*w2;
				std::cout << "Band pass filtering " << w1 << " " << w2 << " " << step << std::endl;
				double uz, uy, ux, uz2, u2, uz2y2, freqI;
				long n=0;
				for(size_t kk=0; kk<ZSIZE(fftV); ++kk)
				{
					FFT_IDX2DIGFREQ(kk,ZSIZE(fftV),uz);
					uz2=uz*uz;

					for(size_t ii=0; ii<YSIZE(fftV); ++ii)
					{
						FFT_IDX2DIGFREQ(ii,YSIZE(fftV),uy);
						uz2y2=uz2+uy*uy;

						for(size_t jj=0; jj<XSIZE(fftV); ++jj)
						{
							FFT_IDX2DIGFREQ(jj,XSIZE(fftV),ux);
							u2=uz2y2+ux*ux;
							//if ((kk != 0) || (ii != 0) || (jj != 0))
							//	freqI = sqrt(u2);
							//else
							//	freqI = 1e38;

							if(u2>=w12 && u2<w22){
								DIRECT_MULTIDIM_ELEM(fftVaux,n) = DIRECT_MULTIDIM_ELEM(fftV,n);
							}

							n++;
						}
					}
				}
				transformer_inv.inverseFourierTransform(fftVaux, V());

//				name = formatString("./filterVol%i_iter%i_new.vol", i, j);
//				saveImg = V();
//				saveImg.write(name);
//				saveImg.clear();

				//Calling to processing function
				processing(V(), pMask, w1);

				if (i==0){
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(V()){
						if (DIRECT_MULTIDIM_ELEM(mask(), n)==1)
							DIRECT_MULTIDIM_ELEM(Vout, n) = DIRECT_MULTIDIM_ELEM(V(), n);
					}
				}else{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(V()){
						if (DIRECT_MULTIDIM_ELEM(mask(), n)==1)
							DIRECT_MULTIDIM_ELEM(Vout, n) += DIRECT_MULTIDIM_ELEM(V(), n);
					}
				}

//				name = formatString("./filterVolProcess%i_iter%i_mono.vol", i, j);
//				saveImg = V();
//				saveImg.write(name);
//				saveImg.clear();

			}

			name.compose(nameRoot, j+1, "vol");
			saveImg = Vout;
			saveImg.write(name);
			saveImg.clear();

		}

		saveImg = Vout;
		saveImg.write(fn_out);
		saveImg.clear();


	}else{

		std::cout << "Selected the whole frequency range" << std::endl;

		for (int j=0; j<iter; j++){
			std::cout << "Iteration " << j << std::endl;

			processing(V(), pMask, 1.0);

			//name = formatString("./outputVol_iter%i.vol", j);
			name.compose(nameRoot, j+1, "vol");
			saveImg = V();
			saveImg.write(name);
			saveImg.clear();

		}

		saveImg = V();
		saveImg.write(fn_out);
		saveImg.clear();

	}


}
