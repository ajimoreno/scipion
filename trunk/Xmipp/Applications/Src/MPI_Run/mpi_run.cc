/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'                                  
 ***************************************************************************/
/* INCLUDES ---------------------------------------------------------------- */
#include <ReconstructionMPI/Programs/Prog_MPI_Run.hh>
#include <mpi.h>

int main (int argc, char *argv[]) {
   if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
      fprintf(stderr, "MPI initialization error\n"); 
      exit(EXIT_FAILURE);
  }  
  
  Prog_MPI_Run_Parameters prm(argc,argv);
    try {
    prm.read(argc,argv);
   }

   catch (Xmipp_error XE) {cout << XE; prm.usage(); exit(1);}

   try {
    if (prm.rank==0) prm.show();
      prm.run();
   }
   catch (Xmipp_error XE) {cout << XE; exit(1);}
   exit(0);
}


