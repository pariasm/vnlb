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
#ifndef VIDEO_NL_BAYES_H_INCLUDED
#define VIDEO_NL_BAYES_H_INCLUDED

#define VNLB_H_VERSION
//#define VNLB_S_VERSION
#if defined(VNLB_H_VERSION)

//	#pragma message ( " Compiling VNLB_H_VERSION " )
	#define CLIPPED_VARIANCE
	#undef PAUL_VARIANCE
	#undef PAUL_SIMPLE_VARIANCE
	#undef FAT_ORIGINAL

#elif defined(VNLB_S_VERSION)

//	#pragma message ( " Compiling VNLB_S_VERSION " )
	#undef CLIPPED_VARIANCE
	#define PAUL_VARIANCE
	#undef PAUL_SIMPLE_VARIANCE
	#undef FAT_ORIGINAL

#endif

/* Variance estimation method.
 *
 * DEFAULT
 * Just substract sigma^2 from the diagonal. This might create negative filter
 * weights for the variances of the sample noisy covariance matrix smaller than
 * sigma2.
 *
 * CLIPPED_VARIANCE (VNLB-H in preprint)
 * Avoids negative weights in the empirical Wiener filter. When the estimated
 * variance of a certain component is lower than the noise variance, the filter
 * coefficient is set to zero. This applies whenever the Gaussian model is
 * estimated from the noisy patches.
 *
 * PAUL_VARIANCE (VLNB-S in preprint)
 * Uses the full inverse of the asymptotic limit of the expected eigenvalue
 * derived by Debashis Paul.
 *
 * PAUL_SIMPLE_VARIANCE
 * Uses a simplified inverse of the asymptotic limit of the expected eigenvalue
 * derived by Debashis Paul. 
 */
//#define CLIPPED_VARIANCE
//#undef  CLIPPED_VARIANCE
//#define PAUL_VARIANCE
//#undef  PAUL_VARIANCE
//#define PAUL_SIMPLE_VARIANCE
//#undef  PAUL_SIMPLE_VARIANCE

/* Flat area tricks : enable the original flat area trick. Otherwise default to
 * the one which uses the basic baricenter in the second step (and does nothing
 * in the first step)*/
//#define FAT_ORIGINAL
//#undef  FAT_ORIGINAL



#include "../VidUtils/LibVideoT.hpp"

namespace VideoNLB
{

/* Structures of parameters dedicated to NL-Bayes process
 *
 * NOTES:
 * - Decoupling the color channels of the 3D patches in the 2nd step: Each
 *   channel is considered independently of the others as in the first step.
 *   This makes the 2nd step faster, with results slitghly slower (0.1~0.2 drop
 *   in PSNR)
 * - sigmaBasic is used for the variance estimators based on the theory of
 *   Debashis Paul (PAUL_VARIANCE, PAUL_SIMPLE_VARIANCE)
 * - procStep and aggreBoost are speed ups by reducing the number of processed
 *   patch groups
 */
struct nlbParams
{
	float sigma;                // noise standard deviation
	float sigmaBasic;           // std. dev. of remanent noise in the basic estimate
	unsigned sizePatch;         // spatial patch size
	unsigned sizePatchTime;     // temporal patch size
	unsigned nSimilarPatches;   // number of similar patches
	unsigned sizeSearchWindow;  // spatial size of search window (w x w)
	unsigned sizeSearchTimeFwd; // how many forward  frames in search window
	unsigned sizeSearchTimeBwd; // how many backward frames in search window
	bool flatAreas;             // use flat area trick
	float gamma;                // threshold parameter to detect flat areas
	bool coupleChannels;        // joint Gaussian model for all channels
	float variThres;            // variance threshold
	unsigned rank;              // rank of covariance matrix
	float beta;                 // noise multiplier
	float tau;                  // patch distance threshold
	bool isFirstStep;           // which step is it?
	unsigned procStep;          // step used to skip reference patches
	bool aggreBoost;            // if true, patches near denoised patches will be skipped
	int onlyFrame;              // denoise only onlyFrame (-1 means process all frames)
	bool verbose;               // verbose output
};

/* Structure containing memory buffers used in the Bayesian estimation.
 */
struct matWorkspace
{
	std::vector<float> group;  // buffer to store the group of similar patches
	std::vector<float> center; // buffer to store the center of the group
	std::vector<float> covMat; // buffer to store the cov matrix
	std::vector<float> covEigVecs; // buffer to store the eigenvecs
	std::vector<float> covEigVals; // buffer to store the eigenvals
};

/* Initialize Parameters of the NL-Bayes algorithm.
 *
 * params : parameter structure to be filled
 * psz_x  : spatial patch size
 * psz_t  : temporal patch size
 * step   : 1 for first step, 2 for second step
 * sigma  : standard deviation of the noise
 * size   : size of the video
 * verbose: if true, print some informations
 *
 * NOTE: The default parameters are set based on sigma. We have determined
 * default parameters for several patch sizes:
 *
 * size around 200:  7 x  7 x 4 grayscale
 *                   8 x  8 x 3 grayscale
 *                  10 x 10 x 2 grayscale
 *                  14 x 14 x 1 grayscale
 * size 128          8 x  8 x 2 grayscale
 * size aroudn 100:  5 x  5 x 4 grayscale RGB
 *                   6 x  6 x 3 grayscale
 *                   7 x  7 x 2 grayscale RGB
 *                  10 x 10 x 1 grayscale RGB
 *
 * For patch sizes that are not in the list, the function will use the default
 * parameters corresponding a patch with the same number of frames from the set
 * with size around 100 pixels.
 *
 */
void defaultParameters(
	nlbParams &params,
	int psz_x,
	int psz_t,
	const unsigned step,
	const float sigma,
	const VideoSize &size,
	const bool verbose);

/* Sets size of spatial search window. It sets the border width accordingly,
 * and also ensures that the number of similar patches is not larger that the 
 * total number of available patches.
 */
void setSizeSearchWindow(nlbParams& prms, unsigned sizeSearchWindow);

/* Sets number of similar patches, ensuring that the number of similar
 * patches is not larger that the total number of available patches.
 */
void setNSimilarPatches(nlbParams& prms, unsigned nSimilarPatches);

/* Display parameters of the NL-Bayes algorithm.
 */
void printNlbParameters(const nlbParams &params);

/* Run the NL-Bayes algorithm in several parallel CPU threads by splitting
 * the video in tiles.
 *
 * imNoisy: contains the noisy video;
 * fflow  : forward optical flow;
 * bflow  : backward optical flow;
 * imBasic: will contain the basic estimate image after the first step;
 * imFinal: will contain the final denoised image after the second step;
 * params1: parameters for first step
 * params1: parameters for second step
 * imClean: original video used as an oracle
 *
 * Returns: Percentage of processed groups over number of pixels.
 */
std::vector<float> runNLBayesThreads(
	Video<float> &imNoisy,
	Video<float> const &fflow,
	Video<float> const &bflow,
	Video<float> &imBasic,
	Video<float> &imFinal,
	const nlbParams params1,
	const nlbParams params2,
	Video<float> &imClean);

/* Run a step of the NL-Bayes denoising for a tile.
 *
 * imNoisy: contains the noisy video;
 * fflow  : optical flow
 * bflow  : optical flow
 * imBasic: will contain the denoised image after the first step (basic estimation);
 * imFinal: will contain the denoised image after the second step;
 * params : parameters of the method
 * crop   : coordinates of the tile
 * imClean: original video to be used as an oracle
 *
 * Returns: Number of processed groups of patches
 */
unsigned processNLBayes(
	Video<float> const& imNoisy,
	Video<float> const& fflow,
	Video<float> const& bflow,
	Video<float> &imBasic,
	Video<float> &imFinal,
	nlbParams const& params,
	VideoUtils::TilePosition crop = VideoUtils::TilePosition(),
	Video<float> const &imClean = Video<float>());

/* Search for patches similar to a reference patch. The input optical flow is
 * used to center the search windows along the trajectory of the reference
 * patch.
 *
 * imNoisy: contains the original noisy image;
 * imBasic: contains the basic estimation;
 * fflow  : forward optical flow (optional)
 * bflow  : backward optical flow (optional)
 * groupNoisy: similar patches from noisy video (output)
 * groupBasic: similar patches from basic estimate(output 2nd step only)
 * indices: positions of similar patches (as indices of a vectorized video)
 * pidx   : index of the reference patch
 * params : see processStep2 for more explanations
 * imClean: original video to be used as an oracle for the search
 *
 * Returns: number of similar patches kept
 */
unsigned estimateSimilarPatches(
	Video<float> const& imNoisy,
	Video<float> const& imBasic,
	Video<float> const& fflow,
	Video<float> const& bflow,
	std::vector<float> &groupNoisy,
	std::vector<float> &groupBasic,
	std::vector<unsigned> &indices,
	const unsigned pidx,
	const nlbParams &params,
	Video<float> const &imClean);

/* Detect flat areas with a statistical test. Flat areas are given a
 * different processing if the flat-area1/2 parameters are set to true.
 *
 * groupNoisy: group of noisy patches
 * groupBasic: group of patches from basic estimate
 * params    : params structure
 * nSimP     : number of similar patches
 * channels  : number of channels in the video
 *
 * Returns: 1 if an homogeneous area is detected, 0 otherwise.
 */
int computeFlatArea(
	std::vector<float> & groupNoisy,
	std::vector<float> const &groupBasic,
	const nlbParams &params,
	const unsigned nSimP,
	const unsigned channels);

/* Bayesian estimation of the patches in the group (inplace) using a
 * Gaussian model learned from the patches. If patches from the basic
 * estimate are given, these are used to learn the parameters of the
 * Gaussian model (2nd step). If not, these are learnt from the noisy
 * patches (1st step).
 *
 * groupNoisy: noisy patches (estimated patches will be written here)
 * groupBasic: patches from the basic estimate
 * mat       : matrix workspace
 * params    : see processStep2 for more explanations;
 * nSimP     : number of similar patches.
 * chnls     : number of video channels
 * flatPatch : true for patches in a flat region
 *
 * Returns: variance of Gaussian model
 */
float computeBayesEstimate(
	std::vector<float> &groupNoisy,
	std::vector<float> &groupBasic,
	matWorkspace &mat,
	nlbParams const& params,
	const unsigned nSimP,
	const unsigned chnls,
	const bool flatPatch = false);

/* Aggregate estimates of similar patches into output video.
 *
 * im     : output video
 * weight : aggregation weights
 * mask   : processing mask (denoised patch are labeled as processed)
 * group  : group of patches
 * indices: positions of the patches
 * params : params structure
 * nSimP  : number of similar patches.
 *
 * Returns: number of processable pixels that were flaged non-processable.
 *
 * NOTE: the processing mask labels patches that are eligible to be
 * used as reference patches as 0. Denoised patches are labeled as 1,
 * and will not be used as reference patches. If aggreBoost is set to 1,
 * neighbors of denoised patches are also labeled as 1.
 */
int computeAggregation(
	Video<float> &im,
	Video<float> &weight,
	Video<char>  &mask,
	std::vector<float> const& group,
	std::vector<unsigned> const& indices,
	const nlbParams& params,
	const unsigned nSimP);

/* Normalize by aggregation weigths.
 *
 * im    : reference, used only when the weight are null
 * outim : aggregation video to normalize
 * weight: aggregation weights
 */
void computeWeightedAggregation(
	Video<float> const& im,
	Video<float> &outIm,
	Video<float> const& weight);

} // namespace

#endif // VIDEO_NL_BAYES_H_INCLUDED
