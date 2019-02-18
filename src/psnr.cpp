/* Copyright (c) 2019, Pablo Arias <pariasm@gmail.com>
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
#include <cstdlib>
#include <cstdio>
#include <string>
#include "VidUtils/LibVideoT.hpp"

int main(int argc, char **argv)
{
	// help message
	if (argc != 3 && argc != 5 && argc != 6)
		return printf("Compute PSNR and RMSE between two image sequences\n\n"
		       "USAGE: psnr input1-path input2-path [first-frame last-frame [output-file]]\n\n"
		       "- Paths for video frames have to be given in printf format\n"
		       "- First and last frame if not provided are set to zero (for an image)\n"
		       "- It prints to standard output the global PSNR and RMSE\n"
		       "- If an output filename is given, it appends to it the global and per-frame PSNR and RMSE\n"),
				 EXIT_FAILURE;

	// parse command line
	const std::string i1_path(argv[1]);
	const std::string i2_path(argv[2]);
	int first_frame = (argc > 3) ? atoi(argv[3]) : 0;
	int last_frame  = (argc > 4) ? atoi(argv[4]) : 0;
	std::string out_path = (argc > 5) ? argv[5] : "";

	// check inputs
	if (i1_path == "" || i2_path == "")
		return fprintf(stderr, "%s: input sequence missing.\nTry `%s' for more "
				"information.\n", argv[0], argv[0]), EXIT_FAILURE;

	// load video, add noise and save
	Video<float> i1, i2;
	i1.loadVideo(i1_path, first_frame, last_frame);
	i2.loadVideo(i2_path, first_frame, last_frame);

	// global psnr
	float psnr = -1, rmse = -1;
	VideoUtils::computePSNR(i1, i2, psnr, rmse);

	// per-frame psnr
	std::vector<float> frames_psnr, frames_rmse;
	VideoUtils::computePSNR(i1, i2, frames_psnr, frames_rmse);

	// output to std
	printf("PSNR %f\nRMSE %f\n", psnr, rmse);

	if (out_path != "")
	{
		if (out_path == "-") out_path = "stdout";
		FILE *f = fopen(out_path.c_str(), "a");
		fprintf(f, "Total PSNR %f\nTotal RMSE %f", psnr, rmse);

		if (last_frame > first_frame)
		{
			fprintf(f, "\nFrame PSNR");
			for (int i = 0; i < frames_psnr.size(); ++i) fprintf(f, " %f", frames_psnr[i]);
			fprintf(f, "\nFrame RMSE ");
			for (int i = 0; i < frames_rmse.size(); ++i) fprintf(f, " %f", frames_rmse[i]);
		}

		fprintf(f, "\n");
		fclose(f);
	}

	return EXIT_SUCCESS;
}
