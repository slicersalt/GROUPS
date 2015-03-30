#include "GroupwiseRegistration.h"

int main(int argc, char** argv)
{
	char sphere[1024] = "/home/hmali/Example/sphere/stx_noscale_996312_V12_t1w_label_pp_surf_tMeanSPHARM_procalign.sphere.vtk";
	
	char *tmpDepth[1024] = {"/home/hmali/Example/property/stx_noscale_961753_V12_t1w_label_pp_surf_tMeanSPHARM_procalign.C.txt",
							"/home/hmali/Example/property/stx_noscale_961753_V12_t1w_label_pp_surf_tMeanSPHARM_procalign.K.txt"};
	
	char *subjDepth[1024] = {"/home/hmali/Example/property/stx_noscale_961850_V12_t1w_label_pp_surf_tMeanSPHARM_procalign.C.txt",
							 "/home/hmali/Example/property/stx_noscale_961850_V12_t1w_label_pp_surf_tMeanSPHARM_procalign.K.txt",
							 "/home/hmali/Example/property/stx_noscale_963368_V12_t1w_label_pp_surf_tMeanSPHARM_procalign.C.txt",
							 "/home/hmali/Example/property/stx_noscale_963368_V12_t1w_label_pp_surf_tMeanSPHARM_procalign.K.txt",
							 "/home/hmali/Example/property/stx_noscale_963992_V12_t1w_label_pp_surf_tMeanSPHARM_procalign.C.txt",
							 "/home/hmali/Example/property/stx_noscale_963992_V12_t1w_label_pp_surf_tMeanSPHARM_procalign.K.txt",
							 "/home/hmali/Example/property/stx_noscale_965421_V12_t1w_label_pp_surf_tMeanSPHARM_procalign.C.txt",
							 "/home/hmali/Example/property/stx_noscale_965421_V12_t1w_label_pp_surf_tMeanSPHARM_procalign.K.txt",
	};
	
	/*char *coeff[1024] = {"../../data/freesurfer3/coeff/coeff_s40s8.lh.mid.txt", 
						"../../data/freesurfer3/coeff/coeff_s40s9.lh.mid.txt"
	};*/
	
	/*char *correspondence[1024] = {"../../data/freesurfer3/corr_s40s1.lh.mid.txt",      
								"../../data/freesurfer3/disp_corr_s40s1.lh.mid.txt", 
	};*/

	//GroupwiseRegistration *r = new GroupwiseRegistration(sphere, tmpDepth, subjDepth, coeff, correspondence, 42, 9, "log.txt", 1);
	//GroupwiseRegistration *r = new GroupwiseRegistration(sphere, tmpDepth, subjDepth, 2, 9, 1, "log.txt");

	//GroupwiseRegistration *r = new GroupwiseRegistration(sphere, tmpDepth, subjDepth, 4, 10, 2);
	
	// char *sphere, char **tmpDepth, char **subjDepth, int nSubj, int deg, int nProperties, char *coeffLog, char **coeff
	GroupwiseRegistration *r = new GroupwiseRegistration(sph, tempProp, propFileList, surfSize, degree, extSize, log, coeffFileList);
	
	//r->saveLDeformation("group.txt");
	/*for (int i = 0; i < 42; i++)
	{
		char fn[1024];
		sprintf(fn, "%s", &coeff[i][29]);
		r->saveLCoeff(fn, i);
	}*/

	return 0;
}
