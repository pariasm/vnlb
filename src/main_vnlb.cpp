/* Modified work: Copyright (c) 2019, Pablo Arias <pariasm@gmail.com>
 * Original work: Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * 
 * This program is free software: you can use, modify and/or redistribute it
 * under the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version. You should have received a copy of this license
 * along this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#include <string>
#include <sstream>
#include <float.h>

#include "VNLBayes/VideoNLBayes.h"
#include "cmd_option.h"

using namespace std;

/* main_vnlb.cpp
 * Main executable file
 */

enum Mode { BSIC_DENO, BSIC_ONLY, DENO_ONLY };

int main(int argc, char **argv)
{
	clo_usage("Video NL-Bayes video denoising");
	clo_help(" NOTE: Input (<) and output (>) sequences are specified by their paths in printf format.\n");

	// Paths to input/output sequences
	using std::string;
	const string  input_path = clo_option("-i"    , ""              , "< input sequence");
	const string  oracl_path = clo_option("-c"    , ""              , "< oracle sequence");
	const string  inbsc_path = clo_option("-b"    , ""              , "< input basic sequence");
	const string  final_path = clo_option("-deno" , "deno_%03d.png" , "> denoised sequence");
	const string  basic_path = clo_option("-bsic" , "bsic_%03d.png" , "> basic denoised sequence");

	      unsigned first_frame = clo_option("-f", 0, "< first frame");
	const unsigned last_frame  = clo_option("-l", 0, "< last frame");
	const unsigned frame_step  = clo_option("-s", 1, "< frame step");

	// Paths to optical flow
	const string  fflow_path = clo_option("-fof", "", "< input forward  optical flow");
	const string  bflow_path = clo_option("-bof", "", "< input backward optical flow");

	// General parameters
	const float sigma  = clo_option("-sigma", 0.f, "< standard deviation of the noise");
	const float sigmab = clo_option("-sigma_basic", 0.f, "< standard dev. of remanent noise in the basic estimate");
	const bool verbose    = (bool) clo_option("-verbose"     , true , "< verbose output");
	const unsigned print_prms = (unsigned) clo_option("-print-prms", 0, "< prints parameters for given channels");

	// Video NLB parameters
	const int time_search1  = clo_option("-wt1", 0  , "< search window temporal radius, step 1");
	const int time_search2  = clo_option("-wt2", 0  , "< search window temporal radius, step 2");
	const int space_search1 = clo_option("-wx1",-1  , "< search window spatial radius, step 1");
	const int space_search2 = clo_option("-wx2",-1  , "< search window spatial radius, step 2");
	const int patch_sizex1  = clo_option("-px1",-1  , "< spatial patch size, step 1");
	const int patch_sizex2  = clo_option("-px2",-1  , "< spatial patch size, step 2");
	const int patch_sizet1  = clo_option("-pt1", 1  , "< temporal patch size, step 1");
	const int patch_sizet2  = clo_option("-pt2", 1  , "< temporal patch size, step 2");
	const int num_patches1  = clo_option("-np1",-1  , "< number of similar patches, step 1");
	const int num_patches2  = clo_option("-np2",-1  , "< number of similar patches, step 2");
	const int rank1         = clo_option("-r1" ,-1  , "< rank of covariance matrix, step 1");
	const int rank2         = clo_option("-r2" ,-1  , "< rank of covariance matrix, step 2");
#if (defined(CLIPPED_VARIANCE) || defined(PAUL_SIMPLE_VARIANCE) || defined(PAUL_VARIANCE))
	const float thres1      = clo_option("-th1",-1.f, "< variance threshold for cov matrix, step 1");
#else
	const float thres1      = clo_option("-th1", -FLT_MAX, "< variance threshold for cov matrix, step 1");
#endif
	const float thres2      = clo_option("-th2",-1.f, "< variance threshold for cov matrix, step 2");
	const float beta1       = clo_option("-b1" ,-1.f, "< noise variance correction factor, step 1");
	const float beta2       = clo_option("-b2" ,-1.f, "< noise variance correction factor, step 2");
	const int patch_step1   = clo_option("-sp1",-1  , "< patch skipping step, step 1");
	const int patch_step2   = clo_option("-sp2",-1  , "< patch skipping step, step 2");
	const bool flat_area1 = (bool) clo_option("-flat-area1", false , "< use flat area trick, step 1");
	const bool flat_area2 = (bool) clo_option("-flat-area2", true  , "< use flat area trick, step 2");
	const bool couple_ch1 = (bool) clo_option("-couple-chnls1", false, "< joint Gaussian for channels, step 1");
	const bool couple_ch2 = (bool) clo_option("-couple-chnls2", false, "< joint Gaussian for channels, step 2");
	const bool no_paste1  = (bool) clo_option("-no-paste1", false , "< disable paste trick, step 1");
	const bool no_paste2  = (bool) clo_option("-no-paste2", false , "< disable paste trick, step 2");
//	const int  only_frame = clo_option("-only",-1, "< process only given frame");

	// Check inputs
	{
		if (input_path == "")
			return fprintf(stderr, "%s: no input sequence.\nTry `%s --help' for more information.\n",
					argv[0], argv[0]), EXIT_FAILURE;

		if (patch_sizex1 == 0 && patch_sizex2 > 0 && inbsc_path == "")
			return fprintf(stderr, "%s: if px1 = 0 and px2 > 0, a basic sequence path must be given.\n"
					"Try `%s --help' for more information.\n", argv[0], argv[0]), EXIT_FAILURE;

		if ((patch_sizex1 < 0 && patch_sizex1 != -1) ||
		    (patch_sizex2 < 0 && patch_sizex2 != -1) ||
		    (patch_sizet1 < 0 || patch_sizet2 <   0) )
			return fprintf(stderr, "%s: px1, px2, pt1 and pt2 cannot be negative.\n"
					"Try `%s --help' for more information.\n", argv[0], argv[0]), EXIT_FAILURE;

		if ((num_patches1 < 0 && num_patches1 != -1) ||
		    (num_patches2 < 0 && num_patches2 != -1))
			return fprintf(stderr, "%s: np1, np2, r1, r2, th1 and th2 cannot be negative.\n"
					"Try `%s --help' for more information.\n", argv[0], argv[0]), EXIT_FAILURE;

		if ((space_search1 < 0 && space_search1 != -1) ||
		    (space_search2 < 0 && space_search2 != -1) ||
		    ( time_search1 < 0 ||  time_search2 <   0) )
			return fprintf(stderr, "%s: wx1, wx2, wt1 and wt2 cannot be negative.\n"
					"Try `%s --help' for more information.\n", argv[0], argv[0]), EXIT_FAILURE;

		if (patch_sizex1 > 0 && inbsc_path != "")
			fprintf(stderr, "\x1b[33;1mWarning:\x1b[0m basic sequence path ignored since px1 > 0.\n");

		if ((fflow_path != "" && bflow_path == "") || 
		    (fflow_path == "" && bflow_path != ""))
			return fprintf(stderr, "Only one oflow path provided.\n"
					"Try `%s --help' for more information.\n", argv[0]), EXIT_FAILURE;

		if ((patch_sizex1 == 0) && (patch_sizex2 == 0))
			return fprintf(stderr, "Given patch sizes for both steps are 0. Nothing to do.\n"),
				EXIT_FAILURE;
	}

	// Determine mode
	Mode mode;
	if ((patch_sizex1 != 0) && (patch_sizex2 != 0)) mode = BSIC_DENO;
	if ((patch_sizex1 != 0) && (patch_sizex2 == 0)) mode = BSIC_ONLY;
	if ((patch_sizex1 == 0) && (patch_sizex2 != 0)) mode = DENO_ONLY;

	bool use_oflow = (fflow_path != "");

	// Only print parameters
	if (print_prms)
	{
		VideoSize tmp;
		tmp.channels = print_prms;

		// Compute denoising default parameters
		VideoNLB::nlbParams prms1, prms2;
		VideoNLB::defaultParameters(prms1, patch_sizex1, patch_sizet1, 1, sigma, tmp, false);
		VideoNLB::defaultParameters(prms2, patch_sizex2, patch_sizet2, 2, sigma, tmp, false);

		// Override with command line parameters
		if (space_search1 >= 0) VideoNLB::setSizeSearchWindow(prms1, (unsigned)space_search1);
		if (space_search2 >= 0) VideoNLB::setSizeSearchWindow(prms2, (unsigned)space_search2);
		if (num_patches1  >= 0) VideoNLB::setNSimilarPatches(prms1, (unsigned)num_patches1);
		if (num_patches2  >= 0) VideoNLB::setNSimilarPatches(prms2, (unsigned)num_patches2);
		if (rank1         >= 0) prms1.rank = rank1;
		if (rank2         >= 0) prms2.rank = rank2;
		if (thres1        >= 0) prms1.variThres = thres1;
		if (thres2        >= 0) prms2.variThres = thres2;
		if (beta1         >= 0) prms1.beta = beta1;
		if (beta2         >= 0) prms2.beta = beta2;
		prms1.flatAreas = flat_area1;
		prms2.flatAreas = flat_area2;
		prms1.coupleChannels = couple_ch1;
		prms2.coupleChannels = couple_ch2;

//		prms1.onlyFrame = only_frame - first_frame;
//		prms2.onlyFrame = only_frame - first_frame;

		if (no_paste1) prms1.aggreBoost = false;
		if (no_paste2) prms2.aggreBoost = false;

		if (patch_step1 >= 0) prms1.procStep = patch_step1;
		if (patch_step2 >= 0) prms2.procStep = patch_step2;

		prms2.sigmaBasic = sigmab;
 
		VideoNLB::printNlbParameters(prms1);
		VideoNLB::printNlbParameters(prms2);

		return EXIT_SUCCESS;
	}

	// Declarations
	Video<float> oracle, noisy, basic, final, diff;
	Video<float> fflow, bflow;

	// Load input videos
	                       noisy.loadVideo(input_path, first_frame, last_frame, frame_step);
	if (mode == DENO_ONLY) basic.loadVideo(inbsc_path, first_frame, last_frame, frame_step);
	if (use_oflow)         fflow.loadVideo(fflow_path, first_frame, last_frame, frame_step);
	if (use_oflow)         bflow.loadVideo(bflow_path, first_frame, last_frame, frame_step);

	if (oracl_path != "")
	{
		oracle.loadVideo(oracl_path, first_frame, last_frame, frame_step);
		if (verbose) printf("Clean video provided as oracle.\n");
	}

	// Denoising
	if (verbose) printf("Running Video NL-Bayes for noise %f\n", sigma);

	// Compute denoising default parameters
	VideoNLB::nlbParams prms1, prms2;
	VideoNLB::defaultParameters(prms1, patch_sizex1, patch_sizet1, 1, sigma, noisy.sz, verbose);
	VideoNLB::defaultParameters(prms2, patch_sizex2, patch_sizet2, 2, sigma, noisy.sz, verbose);

	// Override with command line parameters
	if (space_search1 >= 0) VideoNLB::setSizeSearchWindow(prms1, (unsigned)space_search1);
	if (space_search2 >= 0) VideoNLB::setSizeSearchWindow(prms2, (unsigned)space_search2);
	if (num_patches1  >= 0) VideoNLB::setNSimilarPatches(prms1, (unsigned)num_patches1);
	if (num_patches2  >= 0) VideoNLB::setNSimilarPatches(prms2, (unsigned)num_patches2);
	if (rank1         >= 0) prms1.rank = rank1;
	if (rank2         >= 0) prms2.rank = rank2;
	if (thres1        >= 0) prms1.variThres = thres1;
	if (thres2        >= 0) prms2.variThres = thres2;
	if (beta1         >= 0) prms1.beta = beta1;
	if (beta2         >= 0) prms2.beta = beta2;
	prms1.flatAreas = flat_area1;
	prms2.flatAreas = flat_area2;
	prms1.coupleChannels = couple_ch1;
	prms2.coupleChannels = couple_ch2;

	if (no_paste1) prms1.aggreBoost = false;
	if (no_paste2) prms2.aggreBoost = false;

	if (patch_step1 >= 0) prms1.procStep = patch_step1;
	if (patch_step2 >= 0) prms2.procStep = patch_step2;

//	prms1.onlyFrame = only_frame - first_frame;
//	prms2.onlyFrame = only_frame - first_frame;

	prms2.sigmaBasic = sigmab;

	// Percentage or processed groups of patches over total number of pixels
	std::vector<float> groupsRatio;

	// Run denoising algorithm
	groupsRatio = VideoNLB::runNLBayesThreads(noisy, fflow, bflow, basic, final,
			                                    prms1, prms2, oracle);

//	if (only_frame >= 0) // FIXME
//	{
//		Video<float> tmp(noisy.sz.width, noisy.sz.height, 1, noisy.sz.channels);
//		if (prms2.sizePatch)
//		{
//			VideoUtils::crop(final, tmp, only_frame - first_frame);
//			final = tmp;
//		}
//
//		VideoUtils::crop(basic , tmp, only_frame - first_frame); basic = tmp;
//		VideoUtils::crop(noisy , tmp, only_frame - first_frame); noisy = tmp;
//		VideoUtils::crop(oracle, tmp, only_frame - first_frame); oracle = tmp;
//		first_frame = only_frame;
//	}

	if (verbose)
		printf("Done. Processed %5.2f%% of possible patch groups in 1st step, and\n"
		       "%5.2f%% in 2nd step.\n", groupsRatio[0], groupsRatio[1]);

	// Save output sequences
	if (verbose) printf("Saving output sequences\n");
	basic.saveVideo(basic_path, first_frame, frame_step);
	if (prms2.sizePatch) final.saveVideo(final_path, first_frame, frame_step);

	return EXIT_SUCCESS;
}
