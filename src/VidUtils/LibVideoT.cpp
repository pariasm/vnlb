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

/* LibVideoT.cpp
 * Spetializations of template functions defined in LibVideoT.hpp
 *
 *
 * @author Pablo Arias <pariasm@gmail.com>
 */

#include "LibVideoT.hpp"

#include <cstdio>
#include <cstdlib> // EXIT_FAILURE
#include <iostream>


extern "C" {
#include "iio.h"
}

template <>
void Video<float>::loadVideo(
	const std::string pathToFiles,
	unsigned firstFrame,
	unsigned lastFrame,
	unsigned frameStep)
{
	clear();

	// open first frame to retrieve frame size and allocate mem
	std::vector<float>::iterator p_data;
	{
		char filename[1024];
		std::sprintf(filename, pathToFiles.c_str(), firstFrame);

		float *imTmp = NULL;
		int w, h, c;
		imTmp = iio_read_image_float_split(filename, &w, &h, &c);
		if (!imTmp)
			throw std::runtime_error("Video<float>::loadVideo: loading of " +
			                         std::string(filename) + " failed");

		// set video size
		sz.width     = w;
		sz.height    = h;
		sz.channels  = c;
		sz.frames    = (lastFrame - firstFrame + frameStep)/frameStep;
		sz.update_fields();

		// allocate memory for video
		data.resize(sz.whcf);

		// copy first frame
		p_data = data.begin();
		for (unsigned k = 0; k < w * h * c; ++k, ++p_data)
			*p_data = imTmp[k];

		free(imTmp);
	}

	// load rest of frames
	for (unsigned f = firstFrame + frameStep; f <= lastFrame; f += frameStep)
	{
		char filename[1024];
		std::sprintf(filename, pathToFiles.c_str(), f);

		float *imTmp = NULL;
		int w, h, c;
		imTmp = iio_read_image_float_split(filename, &w, &h, &c);
		if (!imTmp)
			throw std::runtime_error("Video<float>::loadVideo: loading of " +
			                         std::string(filename) + " failed");

		if (w != sz.width || h != sz.height || c != sz.channels)
			throw std::runtime_error("Video<T>::loadVideo: size of " + 
			                         std::string(filename) + " doesn't match video size");

		// copy frame data
		for (unsigned k = 0; k < w * h * c; ++k, ++p_data)
			*p_data = imTmp[k];
	}

	return;
}

template <>
void Video<float>::saveVideo(
	const std::string pathToFiles,
	unsigned firstFrame,
	unsigned frameStep) const
{
	// allocate memory for frame buffer
	float* imTmp = new float[sz.whc];

	unsigned lastFrame = firstFrame + sz.frames * frameStep;
	std::vector<float> frame(sz.whc);
	std::vector<float>::const_iterator p_data = data.begin();
	for (unsigned f = firstFrame; f < lastFrame; f += frameStep)
	{
		// copy frame data to buffer
		for (unsigned k = 0; k < sz.whc; ++k, ++p_data) imTmp[k] = *p_data;

		// file name
		char filename[1024];
		std::sprintf(filename, pathToFiles.c_str(), f);
		iio_save_image_float_split(filename, imTmp, sz.width, sz.height, sz.channels);
	}

	// free memory
	delete[] imTmp;

	return;
}

template <>
void Video<float>::saveVideoAscii(
	const std::string prefix,
	unsigned firstFrame,
	unsigned frameStep) const
{

	char channel[32], frame[32];

	int C = sz.channels;
	int F = sz.frames;
	int W = sz.width;
	int H = sz.height;

	for (size_t c = 0; c < C; ++c)
	for (size_t f = 0; f < F; ++f)
	{
		if (C > 1) std::sprintf(channel,"_ch%d",(int)c); else channel[0] = 0;
		if (F > 1) std::sprintf(frame  ,"_%03d",(int)f); else frame[0] = 0;

		std::string filename = prefix + channel + frame + ".asc";
		std::FILE *const nfile = std::fopen(filename.c_str(),"w");

		unsigned idx = sz.index(0,0,f,c);
		for (size_t y = 0; y < H; ++y)
		{
			for (size_t x = 0; x < W; ++x, ++idx)
				std::fprintf(nfile,"%.16g " , (double)data[idx]);
			std::fprintf(nfile, "\n");
		}

		std::fclose(nfile);
	}
	return;
}

