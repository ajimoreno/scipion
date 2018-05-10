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
    addParamsLine("   -i <volume>              : Volume to calculate the gain");
    addParamsLine("   [--mask >]               : Mask defining the macromolecule");
}

// Produce side information ================================================
void ProgVolumeGain::produce_side_info()
{
	std::cout << "Calculate monogenic signal..." << std::endl;

	V.read(fn_vol);
    V().setXmippOrigin();
    mask.read(fn_mask);
    mask().setXmippOrigin();

	FourierTransformer transformer;
	MultidimArray<double> &inputVol = V();
	VRiesz.resizeNoCopy(inputVol);

	transformer.FourierTransform(inputVol, fftV);
	iu.initZeros(fftV);

	// Calculate u and first component of Riesz vector
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

	// Prepare mask
	MultidimArray<int> &pMask=mask();
	size_t R = (V().xdim)/2.0;

	NVoxelsOriginalMask = 0;
	FOR_ALL_ELEMENTS_IN_ARRAY3D(pMask)
	{
		if (A3D_ELEM(pMask, k, i, j) == 1)
			NVoxelsOriginalMask++;
		if (i*i+j*j+k*k > R*R)
			A3D_ELEM(pMask, k, i, j) = -1;
	}

	mask.write("mask.vol");

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

void ProgVolumeGain::run()
{
    produce_side_info();
}
