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
#ifndef _PROG_VOL_GAIN
#define _PROG_VOL_GAIN

#include <data/xmipp_funcs.h>
#include <data/multidim_array.h>
#include <data/xmipp_image.h>
#include <data/xmipp_program.h>
#include "fourier_filter.h"

class ProgVolumeGain: public XmippProgram
{
public:
    /// Input, monores, and output volume
    FileName fn_vol, fn_mono, fn_out;
    // Input mask
    FileName fn_mask;
    //Size of the box to calculate every histogram
    int boxSize;
    //Histogram matching by frequency bands or with the complete frequencies
    bool bandpass;
    //Number of bandpass filters
    int Nbands;
    //Number of iterations
    int iter;
    //To allow superposed voxel in histogram calculation
    bool superposed;
    //Sampling rate
    double sampling;
    //To allow monores input
    bool mono;

public:
    // Input volume
    Image<double> V, monoRes, Vini;
    Image<int> mask;
    MultidimArray<double> iu, VRiesz, Vweight;
    MultidimArray< std::complex<double> > fftV, fftVRiesz, fftVRiesz_aux;
    int NVoxelsOriginalMask;
    Matrix1D<double> freq_fourier;
    FourierTransformer transformer_inv;



public:
    /// Read arguments
    void readParams();

    /// Show
    void show() const;

    /// Define parameters
    void defineParams();

    /** Produce side info**/
    void produce_side_info();

    /** Run */
    void run();

    void matchingLocalHistogram(MultidimArray<double> amplitude, MultidimArray<double> &gainOut,
    		std::vector< double > cdfGlobal, MultidimArray<int> *pMask, int boxSize, double freq);

    void processing (MultidimArray<double> &V, MultidimArray<int> *pMask, double freq);

};
//@}
#endif
