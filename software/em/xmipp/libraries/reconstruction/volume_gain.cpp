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
    	Vweight.initZeros(monoRes());
    }

	// Prepare mask
	MultidimArray<int> &pMask=mask();
	size_t R = (V().xdim)/2.0;

	FOR_ALL_ELEMENTS_IN_ARRAY3D(pMask)
	{
		if (i*i+j*j+k*k > R*R)
			A3D_ELEM(pMask, k, i, j) = 0;
	}
}


int myBinarySearch(std::vector< double > amplitudes, double value){

	int posi;
	int ini=0;
	int fin=amplitudes.size()-1;
	int aux=0;
	int pos=-1;

	if(amplitudes.size()>2){
		while(pos==-1){
			posi = ((fin-ini+1)/2)+aux;
			if (amplitudes[posi]<=value){
				if(amplitudes[posi+1]>value){
					pos=posi;
					break;
				}else if(((fin-1)==ini || fin==ini) && fin==(amplitudes.size()-1)){
					pos=posi+1;
					break;
				}
				ini=posi+1;
				aux=posi;
			}else{
				if(amplitudes[posi-1]<=value){
					pos=posi-1;
					break;
				}
				fin=posi-1;
			}
			//std::cout << "My binary search " << amplitudes.size() << " " << posi << " " << ini << " " << fin << " " << amplitudes[posi-1] << " " << amplitudes[posi] << " " << amplitudes[posi+1] << " " << value << std::endl;
		}
	}
	else if(amplitudes.size()==1){
		pos=0;
	}
	else if(amplitudes.size()==2){
		pos=0;
		if(amplitudes[1]<=value){
			pos=1;
		}
	}
	//std::cout << "My binary search " << amplitudes.size() << " " << posi << " " << ini << " " << fin << " " << amplitudes[posi-1] << " " << amplitudes[posi] << " " << amplitudes[posi+1] << " " << value << std::endl;
	return pos;

}


void ProgVolumeGain::matchingLocalHistogram(MultidimArray<double> amplitude, MultidimArray<double> &gainOut,
		std::vector< double > globalAmplitudes, MultidimArray<int> mask, int boxSize, double freq)
{

	std::cout << "Matching local histograms... " << std::endl;
	/*MultidimArray<int> mask_aux;
	mask_aux.resize((*pMask));
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY((*pMask)){
		DIRECT_MULTIDIM_ELEM(mask_aux, n) = DIRECT_MULTIDIM_ELEM(*pMask, n);
	}*/

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
				if (DIRECT_MULTIDIM_ELEM(mask, p)==0)
					continue;

				for (int kk=-boxSize; kk<=boxSize; kk++){
					for (int ii=-boxSize; ii<=boxSize; ii++){
						for (int jj=-boxSize; jj<=boxSize; jj++){
							pp = p + jj + ii*xdim + kk*xdim*ydim;

							if(pp<0 || pp>total)
								continue;

							if (DIRECT_MULTIDIM_ELEM(mask, pp)==1){
								//if(!mono || (mono && (sampling/DIRECT_MULTIDIM_ELEM(monoRes(), pp)) > freq)){
									localAmplitudes.push_back(DIRECT_MULTIDIM_ELEM(amplitude, pp));
									DIRECT_MULTIDIM_ELEM(mask, pp)=2;
								//}
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

								if (DIRECT_MULTIDIM_ELEM(mask, pp)==2){
									if(superposed)
										DIRECT_MULTIDIM_ELEM(mask, pp)=1;
									else
										DIRECT_MULTIDIM_ELEM(mask, pp)=0;
									double valueLocal = DIRECT_MULTIDIM_ELEM(amplitude, pp);
									/*int position=-1;
									for (int a=0; a<localAmplitudes.size(); a++){
										if (localAmplitudes[a]>valueLocal){
											position=a-1;
											break;
										}
									}
									if(position==-1)
										position=localAmplitudes.size()-1;
									std::cout << "Normal search: " << valueLocal << " " << localAmplitudes[position] << " " << localAmplitudes[position+1] << std::endl;
									int posi = myBinarySearch(localAmplitudes, valueLocal);
									if(posi!=position){
										std::cout << "ERRORRRRR: " << posi << " " << position << std::endl;
										exit(0);
									}else{
										std::cout << "MATCH: " << posi << " " << position << std::endl;
									}
									std::cout << "---------------------------------------------------------" << std::endl;*/
									int position = myBinarySearch(localAmplitudes, valueLocal);
									double probLocal = (double)position/(double)localAmplitudes.size();
									double newAmplitude;
									int posGlobal = (int)(probLocal*lenGlobal);

									if (posGlobal>=lenGlobal)
										posGlobal=lenGlobal-1;
									newAmplitude = globalAmplitudes[posGlobal];

									if(superposed){
										DIRECT_MULTIDIM_ELEM(gainOut, pp) += newAmplitude;
										if(!mono)
											DIRECT_MULTIDIM_ELEM(voxelCount, pp) += 1;
										else
											DIRECT_MULTIDIM_ELEM(voxelCount, pp) += DIRECT_MULTIDIM_ELEM(Vweight, pp);
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

	/*Image<double> saveImg;
	FileName name;
	name = formatString("./voxelCount%i.vol", (int)(freq*1000));
	saveImg = voxelCount;
	saveImg.write(name);
	saveImg.clear();*/


}


void ProgVolumeGain::processing (MultidimArray<double> &V, MultidimArray<int> *pMask, double freq)
{

	MultidimArray<double> gainOut;
	std::vector< double > globalAmplitudes;
	double max = V.computeMax();
	double min = V.computeMin();
	//std::cout << "Max = " << max << " Min = " << min << std::endl;

	//Calculate sort vector of global amplitudes values to obtain the cdf
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(V){
		if (DIRECT_MULTIDIM_ELEM(mask(), n)==1)
		{
			//if(!mono || (mono && (sampling/DIRECT_MULTIDIM_ELEM(monoRes(), n)) > freq))
				globalAmplitudes.push_back(DIRECT_MULTIDIM_ELEM(V, n));
		}
	}
	std::sort(globalAmplitudes.begin(), globalAmplitudes.end());

	//Matching histograms
	gainOut.initZeros(V);
	matchingLocalHistogram(V, gainOut, globalAmplitudes, *pMask, boxSize, freq);

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
				if(j>0 && i==0){ //Iteraciones mayor que 1, primera banda de freq, rehacer la fft con el volumen de salida de la iteracion anterior
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

							if(u2>=w12 && u2<w22){
								DIRECT_MULTIDIM_ELEM(fftVaux,n) = DIRECT_MULTIDIM_ELEM(fftV,n);
							}

							n++;
						}
					}
				}
				transformer_inv.inverseFourierTransform(fftVaux, V());

//				name = formatString("./beforeMono_filterVol%i_iter%i.vol", i, j);
//				saveImg = V();
//				saveImg.write(name);
//				saveImg.clear();

				//To apply monores pushing down the frequencies above the maximum resolution (minimum value)
				double eval, weightGauss;
				double mean = sampling/w1;
				double sigma=7.0;
				//double weightAux = (1/(sigma*sqrt(2*PI)));
				double aux = 1.0/(2*sigma*sigma);
				if(mono){
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(V()){
						DIRECT_MULTIDIM_ELEM(Vweight, n)=1.0;
						if (DIRECT_MULTIDIM_ELEM(mask(), n)==1)
						{
							if(sampling/DIRECT_MULTIDIM_ELEM(monoRes(), n) < w1){
								eval = DIRECT_MULTIDIM_ELEM(monoRes(), n);
								weightGauss = (exp(-(eval-mean)*(eval-mean)*aux));
								//std::cout << "AAAAA mean " << mean << " eval " << eval << " weightGauss " << weightGauss << std::endl;
								DIRECT_MULTIDIM_ELEM(Vweight, n)=weightGauss;
								//std::cout << "1 DIRECT_MULTIDIM_ELEM(V, n) " << DIRECT_MULTIDIM_ELEM(V(), n) << std::endl;
								DIRECT_MULTIDIM_ELEM(V(), n) = DIRECT_MULTIDIM_ELEM(V(), n) * weightGauss;
								//std::cout << "2 DIRECT_MULTIDIM_ELEM(V, n) " << DIRECT_MULTIDIM_ELEM(V(), n) << std::endl;
							}
						}
					}
				}

//				name = formatString("./filterVol%i_iter%i.vol", i, j);
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

//				name = formatString("./filterVolProcess%i_iter%i_50.vol", i, j);
//				saveImg = V();
//				saveImg.write(name);
//				saveImg.clear();

			}

			name.compose(nameRoot, j+1, "vol");
			saveImg = Vout;
			saveImg.write(name);
			saveImg.clear();

		}

		fftV.clear();
		fftVaux.clear();
		if(mono)
			Vweight.clear();
		Vini.read(fn_vol);
	    Vini().setXmippOrigin();

		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Vout){
			if (DIRECT_MULTIDIM_ELEM(mask(), n)==0)
				DIRECT_MULTIDIM_ELEM(Vout, n) += DIRECT_MULTIDIM_ELEM(Vini(), n);
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

		Vini.read(fn_vol);
		Vini().setXmippOrigin();

		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(V()){
			if (DIRECT_MULTIDIM_ELEM(mask(), n)==0)
				DIRECT_MULTIDIM_ELEM(V(), n) += DIRECT_MULTIDIM_ELEM(Vini(), n);
		}

		saveImg = V();
		saveImg.write(fn_out);
		saveImg.clear();

	}


}
