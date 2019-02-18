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
#ifndef LIB_VIDEOT_HPP_INCLUDED
#define LIB_VIDEOT_HPP_INCLUDED

#include <vector>
#include <string>
#include <stdexcept>
#include <cassert>
#include <climits>
#include <cstdio>
#include <cmath>
#include <unistd.h>

#include "mt19937ar.h"

// Boundary condition type
enum VideoBC {BC_SYMM, BC_CLIP};

/* Structure containing size informations of a video.
 *
 * width    : width of the image;
 * height   : height of the image;
 * channels : number of channels in the image;
 * frames   : number of frames in the video;
 * wh       : equal to width * height. Provided for convenience;
 * whc      : equal to width * height * channels. Provided for convenience.
 * whcf     : equal to width * height * frames * channels. Provided for convenience.
 * whf      : equal to width * height * frames. Provided for convenience.
 */
struct VideoSize
{
	unsigned width;
	unsigned height;
	unsigned frames;
	unsigned channels;
	unsigned wh;
	unsigned whc;
	unsigned whcf;
	unsigned whf;

	// Constuctors
	VideoSize(void)
		: width(0), height(0), frames(0), channels(0)
	{
		update_fields();
	}

	VideoSize(unsigned w, unsigned h, unsigned f, unsigned c)
		: width(w), height(h), frames(f), channels(c)
	{
		update_fields();
	}

	// Comparison operators
	inline bool operator == (const VideoSize& sz) const
	{
		return (width    == sz.width     &&
		        height   == sz.height    &&
		        channels == sz.channels  &&
		        frames   == sz.frames    );
	}

	inline bool operator != (const VideoSize& sz) const
	{ 
		return !operator==(sz);
	}

	// Updates products of dimensions
	inline void update_fields(void)
	{
		wh = width * height;
		whc = wh * channels;
		whcf = whc * frames;
		whf  = wh  * frames;
	}

	// Returns index
	inline unsigned index(unsigned x, unsigned y, unsigned t, unsigned c) const
	{
		assert(x < width && y < height && t < frames && c < channels);
		return t*whc + c*wh + y*width + x;
	}

	inline unsigned index(VideoBC bc_type, unsigned x, unsigned y, unsigned t, unsigned c) const
	{
		assert(c >= 0 && c < sz.channels);

		switch (bc_type)
		{
			case BC_SYMM:
				// NOTE: assumes that -width+1 < x < 2*(width -1)
				assert( -(int)width  < x && x < 2*(int)width  - 1 &&
				        -(int)height < y && y < 2*(int)height - 1 &&
				        -(int)frames < t && t < 2*(int)frames - 1 );

				x = (x < 0) ? -x : (x >= (int)width ) ? 2*(int)width  - 2 - x : x ;
				y = (y < 0) ? -y : (y >= (int)height) ? 2*(int)height - 2 - y : y ;
				t = (t < 0) ? -t : (t >= (int)frames) ? 2*(int)frames - 2 - t : t ;
				break;

			case BC_CLIP:
				x = (x < 0) ? 0 : (x >= (int)width ) ? (int)width  - 1 : x ;
				y = (y < 0) ? 0 : (y >= (int)height) ? (int)height - 1 : y ;
				t = (t < 0) ? 0 : (t >= (int)frames) ? (int)frames - 1 : t ;
				break;

			default:
				throw std::runtime_error("VideoSize::index(bc_type,x,y,t,c): unsupported bc_type");
		}

		return t*whc + c*wh + y*width + x;
	}

	// Returns index assuming the video has one channel
	inline unsigned index(unsigned x, unsigned y, unsigned t) const
	{
		assert(x < width && y < height && t < frames);
		return t*wh + y*width + x;
	}

	inline unsigned index(VideoBC bc_type, unsigned x, unsigned y, unsigned t) const
	{
		return index(bc_type, x, y, t, 0);
	}

	// Compute coordinates from index
	inline
	void coords(unsigned idx, unsigned& x, unsigned& y, unsigned& t, unsigned& c) const
	{
		assert(idx < whcf);
		t = (idx      ) / whc;
		c = (idx % whc) / wh ;
		y = (idx % wh ) / width;
		x = (idx % width  );
	}

	// Compute coordinates from index assuming the video has one channel
	inline
	void coords(unsigned idx, unsigned& x, unsigned& y, unsigned& t) const
	{
		assert(idx < whf);
		t = (idx      ) / wh;
		y = (idx % wh ) / width;
		x = (idx % width  );
	}
};

/* A video class template with very basic functionalities.
 *
 * NOTE: should not be used with T = bool, since current implementation
 * relies on std::vector, and std::vector<bool> cannot return a non-
 * constant reference to an element of the array.
 *
 * sz       : VideoSize structure with size of the video;
 * data     : pointer to an std::vector<T> containing the data
 */
template <class T>
class Video
{
	public:

		// Size
		VideoSize sz;

		// Data
		std::vector<T> data;

		// Constructors
		Video(void); // empty
		Video(const Video& in); // copy
		Video(const std::string pathToFiles,
		      unsigned firstFrame, unsigned lastFrame, unsigned frameStep = 1); // from filename
		Video(unsigned width, unsigned height, unsigned frames, unsigned channels = 1);  // alloc
		Video(unsigned width, unsigned height, unsigned frames, unsigned channels, T val);  // init
		Video(const VideoSize& size);  // alloc
		Video(const VideoSize& size, T val);  // init

		// Destructor
		~Video(void) { };

		void clear(void);
		void resize(unsigned width, unsigned height, unsigned frames, unsigned channels = 1);
		void resize(const VideoSize& size);

		// Read/write pixel access ~ inline for efficiency
		T& operator () (unsigned idx); // from coordinates
		T& operator () (unsigned x, unsigned y, unsigned t, unsigned c = 0); // from coordinates
		T& operator () (VideoBC bc_type, int x, int y, int t, int c = 0); // with boundary conditions

		// Read only pixel access ~ inline for efficiency
		T operator () (unsigned idx) const; // from index
		T operator () (unsigned x, unsigned y, unsigned t, unsigned c = 0) const; // from coordinates
		T operator () (VideoBC bc_type, int x, int y, int t, int c = 0) const; // with boundary conditions

		// Pixel access with special boundary conditions
		T& getPixelSymmetric(int x, int y, int t, unsigned c = 0);
		T  getPixelSymmetric(int x, int y, int t, unsigned c = 0) const;
		
		// I/O
		void loadVideo(const std::string pathToFiles, 
		               unsigned firstFrame, unsigned lastFrame, unsigned frameStep = 1);
		void saveVideo(const std::string pathToFiles, 
		               unsigned firstFrame, unsigned frameStep = 1) const;
		void saveVideoAscii(const std::string prefix, 
		                    unsigned firstFrame, unsigned frameStep = 1) const;
};

// Implementations

template <class T>
Video<T>::Video(void) : sz(), data(0) { }
	
template <class T>
Video<T>::Video(const Video& in) : sz(in.sz), data(in.data) { }

template <class T>
Video<T>::Video(
	const std::string pathToFiles,
	unsigned firstFrame,
	unsigned lastFrame,
	unsigned frameStep) : sz(), data(0)
{
	loadVideo(pathToFiles, firstFrame, lastFrame, frameStep);
}
	
template <class T>
Video<T>::Video(
	unsigned width,
	unsigned height,
	unsigned frames,
	unsigned channels)
	: sz(width, height, frames, channels) , data(sz.whcf) { }

template <class T>
Video<T>::Video(
	unsigned width,
	unsigned height,
	unsigned frames,
	unsigned channels,
	T val)
	: sz(width, height, frames, channels) , data(sz.whcf, val) { }

template <class T>
Video<T>::Video(const VideoSize& size, T val) : sz(size), data(sz.whcf, val) { }

template <class T>
Video<T>::Video(const VideoSize& size) : sz(size) , data(sz.whcf) { }
	
template <class T>
void Video<T>::clear(void)
{
	sz.width = 0;
	sz.height = 0;
	sz.frames = 0;
	sz.channels = 0;
	sz.update_fields();
	data.clear();
}

template <class T>
void Video<T>::resize(const VideoSize& size)
{
	if (sz != size)
	{
		clear();
		sz = size;
		data.resize(sz.whcf);
	}
}

template <class T>
void Video<T>::resize(
	unsigned width,
	unsigned height,
	unsigned frames,
	unsigned channels)
{
	resize(VideoSize(width, height, frames, channels));
}

template <class T>
void Video<T>::loadVideo(
	const std::string pathToFiles,
	unsigned firstFrame,
	unsigned lastFrame,
	unsigned frameStep)
{
	throw std::runtime_error("Video<T>::loadVideo(...) is only implemented "
	                         "for T = float");
}

template <> void Video<float>::loadVideo(
	const std::string pathToFiles,
	unsigned firstFrame,
	unsigned lastFrame,
	unsigned frameStep);

template <class T>
void Video<T>::saveVideo(
	const std::string pathToFiles,
	unsigned firstFrame,
	unsigned frameStep) const
{
	throw std::runtime_error("Video<T>::saveVideo(...) is only implemented "
	                         "for T = float");
}

template <> void Video<float>::saveVideo(
	const std::string pathToFiles,
	unsigned firstFrame,
	unsigned frameStep) const;

template <class T>
void Video<T>::saveVideoAscii(
	const std::string prefix,
	unsigned firstFrame,
	unsigned frameStep) const
{
	throw std::runtime_error("Video<T>::saveVideoAscii(...) is only implemented "
	                         "for T = float");
}

template <> void Video<float>::saveVideoAscii(
	const std::string prefix,
	unsigned firstFrame,
	unsigned frameStep) const;

template <class T>
inline T& Video<T>::operator () (unsigned idx)
{
	assert(idx < sz.whcf);
	return data[idx];
}

template <class T>
inline T& Video<T>::operator() (unsigned x, unsigned y, unsigned t, unsigned c)
{
	return data[sz.index(x,y,t,c)];
}

template <class T>
inline T& Video<T>::operator() (VideoBC bc_type, int x, int y, int t, int c)
{
	return data[sz.index(bc_type,x,y,t,c)];
}

template <class T>
inline T Video<T>::operator () (unsigned idx) const
{
	assert(idx < sz.whcf);
	return data[idx];
}

template <class T>
inline T Video<T>::operator() (unsigned x, unsigned y, unsigned t, unsigned c) const
{
	return data[sz.index(x,y,t,c)];
}

template <class T>
inline T Video<T>::operator() (VideoBC bc_type, int x, int y, int t, int c) const
{
	return data[sz.index(bc_type,x,y,t,c)];
}

template <class T>
inline T& Video<T>::getPixelSymmetric(int x, int y, int t, unsigned c)
{
	// NOTE: assumes that -width+1 < x < 2*(width -1)
	assert(-(int)sz.width   < x && x < 2*(int)sz.width -1&&
	       -(int)sz.height  < y && y < 2*(int)sz.height-1&&
	       -(int)sz.frames  < t && t < 2*(int)sz.frames-1);
	// symmetrize
	x = (x < 0) ? -x : (x >= (int)sz.width  ) ? 2*(int)sz.width  - 2 - x : x ;
	y = (y < 0) ? -y : (y >= (int)sz.height ) ? 2*(int)sz.height - 2 - y : y ;
	t = (t < 0) ? -t : (t >= (int)sz.frames ) ? 2*(int)sz.frames - 2 - t : t ;

	return data[sz.index(x,y,t,c)];
}

template <class T>
inline T Video<T>::getPixelSymmetric(int x, int y, int t, unsigned c) const
{
	// NOTE: assumes that -width+1 < x < 2*(width -1)
	assert(-(int)sz.width   < x && x < 2*(int)sz.width  - 1 &&
	       -(int)sz.height  < y && y < 2*(int)sz.height - 1 &&
	       -(int)sz.frames  < t && t < 2*(int)sz.frames - 1 );
	// symmetrize
	x = (x < 0) ? -x : (x >= (int)sz.width  ) ? 2*(int)sz.width  - 2 - x : x ;
	y = (y < 0) ? -y : (y >= (int)sz.height ) ? 2*(int)sz.height - 2 - y : y ;
	t = (t < 0) ? -t : (t >= (int)sz.frames ) ? 2*(int)sz.frames - 2 - t : t ;

	return data[sz.index(x,y,t,c)];
}


// Utilities for video
namespace VideoUtils
{
	/* Structure to store the position of a rectangular tile. It also has data
	 * to describe the position of the tile when it corresponds to a rectangular
	 * tiling (potentially with an added border) of a video. The tiles in a
	 * tiling do not overlap. The tile is contained within a crop. A crop is the
	 * union of a tile and a border surrounding it. This border is necessary for
	 * any boundary conditions that the processing applied to the tile might
	 * require.
	 *
	 * In this structure we store
	 * 1) coordinates of the crop with respect to the video coordinate system
	 * 2) position of the tile in a tiling (i.e. 2nd tile from the left, 3rd from the top)
	 * 3) coordinates of the tile with respect to the video coordinate system
	 */
	struct TilePosition
	{
		int origin_x; // x coordinate of top-left-front corner of crop
		int origin_y; // y coordinate of top-left-front corner of crop
		int origin_t; // t coordinate of top-left-front corner of crop

		int ending_x; // x coordinate of bottom-right-back corner of crop
		int ending_y; // y coordinate of bottom-right-back corner of crop
		int ending_t; // t coordinate of bottom-right-back corner of crop

		VideoSize source_sz; // size of source video

		int tile_x;    // x index of tile (0 <= tile_x < ntiles_x)
		int tile_y;    // y index of tile (0 <= tile_y < ntiles_y)
		int tile_t;    // t index of tile (0 <= tile_t < ntiles_t)

		int ntiles_x;  // total number of tiles in x direction
		int ntiles_y;  // total number of tiles in y direction
		int ntiles_t;  // total number of tiles in t direction

		int tile_origin_x; // x coordinate of top-left-front corner of tile
		int tile_origin_y; // y coordinate of top-left-front corner of tile
		int tile_origin_t; // t coordinate of top-left-front corner of tile

		int tile_ending_x; // x coordinate of bottom-right-back corner of tile
		int tile_ending_y; // y coordinate of bottom-right-back corner of tile
		int tile_ending_t; // t coordinate of bottom-right-back corner of tile
	};

	/* Add additive white Gaussian noise to video.
	 * vid     : input video
	 * vidNoisy: output noisy video
	 * sigma   : standard deviation of the noise.
	 */
	template <class T>
	void addNoise(Video<T> const& vid, Video<T> &vidNoisy, const float sigma)
	{
		vidNoisy = vid;

		// Random seed
//		mt_init_genrand((unsigned long int) time (NULL) +
//		                (unsigned long int) getpid());
		mt_init_genrand(0); printf("\x1b[33;1mWarning:\x1b[0m random generator seed is 0\n");

		// Add noise
		for (unsigned k = 0; k < vid.sz.whcf; k++)
		{
			const double a = mt_genrand_res53();
			const double b = mt_genrand_res53();
			vidNoisy(k) += (T) (sigma *
				sqrtl(-2.0l * log(a)) * cos(2.0l * M_PI * b));
		}
	}
	
	/* Compute PSNR and RMSE between vid1 and vid2.
	 * vid1: video 1;
	 * vid2: video 2;
	 * psnr: PSNR;
	 * rmse: RMSE;
	 */
	template <class T>
	void computePSNR(Video<T> const& vid1, Video<T> const& vid2, float &psnr, float &rmse)
	{
		if (vid1.sz != vid2.sz)
			throw std::runtime_error("VideoUtils::computePSNR: videos have different sizes");

#if 1
		// using doubles
		double sum = 0.f;
		for (unsigned k = 0; k < vid1.sz.whcf; k++)
			sum += ((double)vid1(k) - (double)vid2(k)) *
			       ((double)vid1(k) - (double)vid2(k));

		rmse = sqrtf(sum / (double) vid1.sz.whcf);
		psnr = 20.f * log10f(255.f / rmse);
#else
		// using doubles and recursive formula
		double mse = 0.f;
		for (unsigned k = 0; k < vid1.sz.whcf; k++)
		{
			double mse_prev = mse;
			double term_k = vid1(k) - vid2(k);
			mse += (term_k*term_k - mse_prev)/((double)(k+1));
		}

		rmse = sqrtf(mse);
		psnr = 20.f * log10f(255.f / rmse);
#endif


		return;
	}

	/* Compute PSNR and RMSE between vid1 and vid2 for each frame
	 * vid1: video 1;
	 * vid2: video 2;
	 * psnr: vector with PSNR of each frame
	 * rmse: vector with RMSE of each frame
	 */
	template <class T>
	void computePSNR(Video<T> const& vid1, Video<T> const& vid2,
			std::vector<float> &psnr, std::vector<float> &rmse)
	{
		if (vid1.sz != vid2.sz)
			throw std::runtime_error("VideoUtils::computeFramesPSNR: videos have different sizes");

		psnr.resize(vid1.sz.frames);
		rmse.resize(vid1.sz.frames);

		for (unsigned f = 0, i = 0; f < vid1.sz.frames; ++f)
		{
			double sum = 0.f;
			for (unsigned k = 0; k < vid1.sz.whc; k++, i++)
				sum += ((double)vid1(i) - (double)vid2(i)) *
						 ((double)vid1(i) - (double)vid2(i));

			rmse[f] = sqrtf(sum / (double) vid1.sz.whc);
			psnr[f] = 20.f * log10f(255.f / rmse[f]);
		}

		return;
	}

	/* 'Generalized' croping of a video (cropped video may be larger than original).
	 *
	 * vid1  : original video;
	 * vid2  : output video, already allocated to desired size;
	 * origin: vid1 coordinates of vid2 origin. Origin coordinates larger
	 *         than corresponding vid1 dimension are redefined to center the crop
	 *         in that dimension.
	 */
	template <class T>
	void crop(Video<T> const &vid1, Video<T> &vid2,
		int origin_t = INT_MAX, int origin_x = INT_MAX, int origin_y = INT_MAX)
	{
		assert(vid2.sz.channels == vid1.sz.channels);

		// Redefine invalid origin coordinates to default (centered crop)
		if (origin_t > (int)vid1.sz.frames) origin_t = ((int)vid1.sz.frames - (int)vid2.sz.frames)/2;
		if (origin_x > (int)vid1.sz.width ) origin_x = ((int)vid1.sz.width  - (int)vid2.sz.width )/2;
		if (origin_y > (int)vid1.sz.height) origin_y = ((int)vid1.sz.height - (int)vid2.sz.height)/2;

		// TODO: more efficient implementation
		for (int      f = 0; f < vid2.sz.frames  ; f++)
		for (unsigned c = 0; c < vid2.sz.channels; c++)
		for (int      y = 0; y < vid2.sz.height  ; y++)
		for (int      x = 0; x < vid2.sz.width   ; x++)
			vid2(x,y,f,c) = 
				vid1.getPixelSymmetric(x + origin_x, y + origin_y, f + origin_t, c);
	}

	/* 'Generalized' croping of a video (cropped video may be larger than original).
	 * vid1  : original video;
	 * vid2  : output video, already allocated to desired size;
	 * origin: vid1 coordinates of vid2 origin. Origin coordinates larger than
	 *         corresponding vid1 dimension are redefined to center the crop in
	 *         that dimension.
	 */
	template <class T>
	void crop(Video<T> const &vid1, Video<T> &vid2, const int * const origin)
	{
		crop(vid1, vid2, origin[2], origin[0], origin[1]);
	}

	/* 'Generalized' croping of a video (cropped video may be larger than original).
	 * vid1  : original video;
	 * vid2  : output video, already allocated to desired size;
	 * origin: vid1 coordinates of vid2 origin. Origin coordinates larger than
	 *         corresponding vid1 dimension are redefined to center the crop in
	 *         that dimension.
	 */
	template <class T>
	void crop(Video<T> const &vid1, Video<T> &vid2, TilePosition const &tile)
	{
		// Resize output video
		VideoSize cropSize;
		cropSize.width    = tile.ending_x - tile.origin_x;
		cropSize.height   = tile.ending_y - tile.origin_y;
		cropSize.frames   = tile.ending_t - tile.origin_t;
		cropSize.channels = vid1.sz.channels;
		cropSize.update_fields();

		vid2.resize(cropSize);

		crop(vid1, vid2, tile.origin_t, tile.origin_x, tile.origin_y);
	}

	/* Add/remove border by symmetryzing image content
	 * vid     : original image;
	 * vidSym  : output image
	 * border  : border size
	 * forward : 1 to add border, 0 to remove
	 */
	template <class T>
	void addBorder(Video<T> const& vid1, Video<T> &vid2, const unsigned border, const bool forward)
	{
		// Sizes
		const unsigned w2 = vid1.sz.width  + (forward ? 2*border : - 2*border);
		const unsigned h2 = vid1.sz.height + (forward ? 2*border : - 2*border);
		const unsigned f2 = vid1.sz.frames + (forward ? 2*border : - 2*border);
		const unsigned ch = vid1.sz.channels;

		// Position of vid2 origin in vid1 coordinates
		const int tx = forward ? -border : border;
		const int ty = forward ? -border : border;
		const int tf = forward ? -border : border;

		// Resize output image, if necessary
		vid2.resize(w2, h2, f2, ch);

		// Call generalized crop
		crop(vid1, vid2, tx, ty, tf);
	}

	/* Transform color space of a video, from RGB to YUV, or vice-versa (inplace).
	 * vid    : input/output video
	 * forward: if true, go from RGB to YUV, otherwise go from YUV to RGB.
	 */
	template <class T>
	void transformColorSpace(Video<T> &vid, const bool forward)
	{
		// If the image as only one channel, do nothing
		if (vid.sz.channels == 1) return;

		// Initialization
		const unsigned width  = vid.sz.width;
		const unsigned height = vid.sz.height;
		const unsigned chnls  = vid.sz.channels;
		const unsigned wh     = vid.sz.wh;

		for (int f = 0; f < vid.sz.frames; f++)
			if (forward) // RGB to YUV
			{
				if (chnls == 3)
				{
					const unsigned red   = f * vid.sz.whc;
					const unsigned green = red + wh;
					const unsigned blue  = green + wh;
					const float a = 1.f / sqrtf(3.f);
					const float b = 1.f / sqrtf(2.f);
					const float c = 2.f * a * sqrtf(2.f);

					float yuv[3];
					for (unsigned k = 0; k < wh; k++)
					{
						yuv[0] = a * ((float)vid(k + red) + (float)vid(k + green) + (float)vid(k + blue));
						yuv[1] = b * ((float)vid(k + red) - (float)vid(k + blue));
						yuv[2] = c * (0.25f * (float)vid(k + red ) - 0.5f * (float)vid(k + green)
								      + 0.25f * (float)vid(k + blue));

						vid(k + red  ) = (T)yuv[0];
						vid(k + green) = (T)yuv[1];
						vid(k + blue ) = (T)yuv[2];
					}
				}
			}
			else // YUV to RGB
			{
				if (chnls == 3)
				{
					const unsigned red   = f * vid.sz.whc;
					const unsigned green = red + wh;
					const unsigned blue  = green + wh;
					const float a = 1.f / sqrtf(3.f);
					const float b = 1.f / sqrtf(2.f);
					const float c = a / b;

					float rgb[3];
					for (unsigned k = 0; k < wh; k++)
					{
						rgb[0] = a * (float)vid(k + red) + b * (float)vid(k + green)
						                             + c * 0.5f * (float)vid(k + blue );
						rgb[1] = a * (float)vid(k + red) - c * (float)vid(k + blue );
						rgb[2] = a * (float)vid(k + red) - b * (float)vid(k + green)
						                      + c * 0.5f * (float)vid(k + blue );
						vid(k + red  ) = (T)rgb[0];
						vid(k + green) = (T)rgb[1];
						vid(k + blue ) = (T)rgb[2];
					}
				}
			}
	}

	/* Determine a and b such that : n = a * b, with a and b as greatest as possible.
	 * n : number to decompose;
	 * a : will contain a;
	 * b : will contain b.
	 */
	void determineFactor(const unsigned n, unsigned &a, unsigned &b)
	{
		if (n == 1)
		{
			a = 1;
			b = 1;
			return;
		}

		b = 2;
		while (n % b > 0) b++;
		a = n / b;

		if (b > a)
		{
			a = b;
			b = n / a;
		}
	}

	/* Subdivide a video into small sub-videos. This version does not add a
	 * border to the video, resulting in parts which might have different sizes.
	 *
	 * video : image to subdivide;
	 * videoSub : will contain all sub-videos;
	 * tiles : will store position of the tiles;
	 * border : boundary around sub-videos;
	 * ntiles : number of sub-videos wanted. Need to be a power of 2.
	 */
	template <class T>
	void subDivideTight(
		Video<T> const& vid,
		std::vector<Video<T> > &vidSub,
		std::vector<TilePosition> &tiles,
		const int border,
		const int ntiles)
	{
		/* FIXME current version splits the video only spatially.
		 *       The reason is to mantain consistency with Marc's
		 *       code. For its proper extension to video, we need
		 *       to determine how to split the video in space and
		 *       time. */

		// Determine number of sub-images
		unsigned u_nW, u_nH; // FIXME problem with unsigned and int
		determineFactor((unsigned)ntiles, u_nW, u_nH);
		int nW = (int)u_nW;
		int nH = (int)u_nH;

		const int wTmp = ceil(float(vid.sz.width ) / float(nW)); // sizes w/out
		const int hTmp = ceil(float(vid.sz.height) / float(nH)); //     borders

		tiles.resize(ntiles);
		vidSub.resize(ntiles);
		for (int p = 0, n = 0; p < nH; p++)
		for (int q = 0;        q < nW; q++, n++)
		{
			// Set crop information 
			tiles[n].source_sz = vid.sz;

			// origin of the crop. the crop contains the tile plus a border
			tiles[n].origin_x = std::max(0, q * wTmp - border);
			tiles[n].origin_y = std::max(0, p * hTmp - border);
			tiles[n].origin_t = 0;

			// end of the crop
			tiles[n].ending_x = std::min((int)vid.sz.width , (q+1) * wTmp + border);
			tiles[n].ending_y = std::min((int)vid.sz.height, (p+1) * hTmp + border);
			tiles[n].ending_t = vid.sz.frames;

			// crop using symmetric boundary conditions
			VideoUtils::crop(vid, vidSub[n], tiles[n]);

			// tile position with respect to the other tiles
			tiles[n].tile_x = q;
			tiles[n].tile_y = p;
			tiles[n].tile_t = 0;

			// number of tiles
			tiles[n].ntiles_x = nW;
			tiles[n].ntiles_y = nH;
			tiles[n].ntiles_t =  1;

			// coordinates of tile origin
			tiles[n].tile_origin_x = q * wTmp;
			tiles[n].tile_origin_y = p * wTmp;
			tiles[n].tile_origin_t = 0;

			// coordinates of tile end
			tiles[n].tile_ending_x = std::min((int)vid.sz.width , (q+1) * wTmp);
			tiles[n].tile_ending_y = std::min((int)vid.sz.height, (p+1) * hTmp);
			tiles[n].tile_ending_t = vid.sz.frames;
		}

		return;
	}

	/* Subdivide a video into smaller tiles
	 * video   : image to subdivide
	 * videoSub: output tiles
	 * tiles   : position of tiles
	 * border  : boundary around sub-videos
	 * ntiles  : number of tiles wanted (need to be a power of 2)
	 */
	template <class T>
	void subDivide(
		Video<T> const& vid,
		std::vector<Video<T> > &vidSub,
		std::vector<TilePosition> &tiles,
		const unsigned border,
		const unsigned ntiles)
	{
		/* FIXME current version splits the video only spatially. 
		 *       The reason is to mantain consistency with Marc's
		 *       code. For its proper extension to video, we need
		 *       to determine how to split the video in space and
		 *       time. */
		
		// Determine number of sub-images
		unsigned u_nW, u_nH; // FIXME problem with unsigned and int
		determineFactor((unsigned)ntiles, u_nW, u_nH);
		int nW = (int)u_nW;
		int nH = (int)u_nH;

		const int wTmp = ceil(float(vid.sz.width ) / float(nW)); // sizes w/out 
		const int hTmp = ceil(float(vid.sz.height) / float(nH)); //     borders

		// Obtain sub-images
		VideoSize imSubSize;
		imSubSize.width    = wTmp + 2 * border; // each sub-image has border
		imSubSize.height   = hTmp + 2 * border;
		imSubSize.frames   = vid.sz.frames; // NOTE: same frames as original
		imSubSize.channels = vid.sz.channels;
		imSubSize.update_fields();

		tiles.resize(ntiles);
		vidSub.resize(ntiles);
		for (int p = 0, n = 0; p < nH; p++)
		for (int q = 0;        q < nW; q++, n++)
		{
			vidSub[n].resize(imSubSize);

			// The origin is shifted -p_N to account for the subimage border
			int origin[3] = {q * wTmp - border, p * hTmp - border, 0};

			// Crop using symmetric boundary conditions
			VideoUtils::crop(vid, vidSub[n], origin);

			// Set crop information
			tiles[n].source_sz = vid.sz;
			tiles[n].origin_x  = origin[0];
			tiles[n].origin_y  = origin[1];
			tiles[n].origin_t  = origin[2];

			tiles[n].ending_x  = origin[0] + vidSub[n].sz.width ;
			tiles[n].ending_y  = origin[1] + vidSub[n].sz.height;
			tiles[n].ending_t  = origin[2] + vidSub[n].sz.frames;

			// Add information about the tiling
			tiles[n].tile_x = q;
			tiles[n].tile_y = p;
			tiles[n].tile_t = 0;

			tiles[n].ntiles_x = nW;
			tiles[n].ntiles_y = nH;
			tiles[n].ntiles_t =  1;

			tiles[n].tile_origin_x = q * wTmp;
			tiles[n].tile_origin_y = p * wTmp;
			tiles[n].tile_origin_t = 0;

			tiles[n].tile_ending_x = std::min((int)vid.sz.width , (q+1) * wTmp);
			tiles[n].tile_ending_y = std::min((int)vid.sz.height, (p+1) * hTmp);
			tiles[n].tile_ending_t = vid.sz.frames;
		}

		return;
	}

	/* Subdivide a video into small sub-videos
	 *
	 * video : image to subdivide;
	 * videoSub : will contain all sub-videos;
	 * border : boundary around sub-videos;
	 * ntiles : number of sub-videos wanted. Need to be a power of 2.
	 *
	 * none.
	 */
	template <class T>
	void subDivide(
		Video<T> const& vid,
		std::vector<Video<T> > &vidSub,
		const unsigned border,
		const unsigned ntiles)
	{
		std::vector<TilePosition> tmpCrops;
		subDivide(vid, vidSub, tmpCrops, border, ntiles);
		return;
	}

	/* Reconstruct an video from its small sub-videos
	 *
	 * vid : image to reconstruct;
	 * vidSub : will contain all sub-images;
	 * border : boundary around sub-videos.
	 *
	 * none.
	 */
	template <class T>
	void subBuild(std::vector<Video<T> > const& vidSub, Video<T> &vid, const unsigned border)
	{
		/* FIXME current version builds a video that has been split
		 *       only spatially by subDivide.
		 *       The reason is to mantain consistency with Marc's
		 *       code. For its proper extension to video, we need
		 *       to determine how to split the video in space and
		 *       time. */

		assert(vidSub.size());
		assert(vidSub[0].sz.whcf);
		assert(vid.sz.whcf);
		assert(vid.sz.frames   == vidSub[0].sz.frames  );
		assert(vid.sz.channels == vidSub[0].sz.channels);

		// Determine width and height composition
		unsigned nW, nH;
		determineFactor(vidSub.size(), nW, nH);
		const unsigned hTmp = vidSub[0].sz.height - 2 * border;
		const unsigned wTmp = vidSub[0].sz.width  - 2 * border;

		// Obtain inner image (containing boundaries)
		// TODO pending decision for video
		for (unsigned py = 0, n = 0; py < nH*hTmp; py += hTmp)
		for (unsigned px = 0       ; px < nW*wTmp; px += wTmp, n++)
		{
			/* Diagram for a 1D image with W = 8, covered
			 * with 2 sub images of w = 5, with border 2.
			 * Symmetrized pixels are indicated with an s.
			 *
			 * ori         0  1  2  3  4  5  6  7
			 * sub1 s0 s1  2  3  4  5  6 s7 s8
			 * sub2                s0 s1  2  3  4 s5 s6 s7 s8
			 * 
			 * px         0*w            1*w       <-- don't exceed 1*w + W-1*w
			 *
			 * Notation: [px,py] coords on big image of sub-image top-left point 
			 *           [qx,qy] point on big image
			 *           [sx,sy] corresponding point on sub-image
			 */
			unsigned wmax = std::min(wTmp, vid.sz.width  - px) + border;
			unsigned hmax = std::min(hTmp, vid.sz.height - py) + border;

			for (unsigned f = 0; f < vid.sz.frames  ; f++)
			for (unsigned c = 0; c < vid.sz.channels; c++)
			for (unsigned sy = border, qy = py; sy < hmax; sy++, qy++)
			for (unsigned sx = border, qx = px; sx < wmax; sx++, qx++)
				vid(qx, qy, f, c) = vidSub[n](sx, sy, f, c);
		}

		return;
	}

	/* Reconstruct an video from its small sub-videos
	 *
	 * vid : image to reconstruct;
	 * vidSub : will contain all sub-images;
	 * border : boundary around sub-videos.
	 *
	 * none.
	 */
	template <class T>
	void subBuildTight(std::vector<Video<T> > const& vidSub, Video<T> &vid, const int border)
	{
		/* FIXME current version builds a video that has been split
		 *       only spatially by subDivide. 
		 *       The reason is to mantain consistency with Marc's
		 *       code. For its proper extension to video, we need
		 *       to determine how to split the video in space and
		 *       time. */

		assert(vidSub.size());
		assert(vidSub[0].sz.whcf);
		assert(vid.sz.whcf);
		assert(vid.sz.frames   == vidSub[0].sz.frames  );
		assert(vid.sz.channels == vidSub[0].sz.channels);

		// Determine width and height composition
		unsigned nW, nH;
		determineFactor(vidSub.size(), nW, nH);
		const int wTmp = ceil(float(vid.sz.width ) / float(nW)); // sizes w/out 
		const int hTmp = ceil(float(vid.sz.height) / float(nH)); //     borders

		// Obtain inner image (containing boundaries)
		// TODO pending decision for video
		for (int p = 0, n = 0; p < nH; p++)
		for (int q = 0;        q < nW; q++, n++)
		{
			// top-left-front corner of crop
			int crop_ori_x = std::max(0, q * wTmp - border);
			int crop_ori_y = std::max(0, p * hTmp - border);
			int crop_ori_t = 0;

			// start of inner crop, removing the boundary
			int ori_x = q * wTmp;
			int ori_y = p * hTmp;
			int ori_t = 0       ;

			// end of inner crop
			int end_x = std::min((q+1) * wTmp, (int)vid.sz.width );
			int end_y = std::min((p+1) * hTmp, (int)vid.sz.height);
			int end_t = 0       ;

			// start of inner crop, with inner coordinates
			int in_ori_x = ori_x - crop_ori_x;
			int in_ori_y = ori_y - crop_ori_y;
			int in_ori_t = ori_t - crop_ori_t;

//			int  in_end_x = (q+1) * wTmp - out_ori_x;
//			int  in_end_y = (p+1) * hTmp - out_ori_y;
//			int  in_end_t = vid.sz.frames;

			for (unsigned f = 0; f < vid.sz.frames  ; f++)
			for (unsigned c = 0; c < vid.sz.channels; c++)
			for (unsigned sy = in_ori_y, qy = ori_y; qy < end_y; sy++, qy++)
			for (unsigned sx = in_ori_x, qx = ori_x; qx < end_x; sx++, qx++)
				vid(qx, qy, f, c) = vidSub[n](sx, sy, f, c);
		}

		return;
	}
}

#endif // LIB_VIDEOT_HPP_INCLUDED
