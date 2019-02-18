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
#include <cstdio>
#include <cstdlib>
#include <string>
#include "VidUtils/LibVideoT.hpp"

int main(int argc, char **argv)
{
	// help message
	if (argc != 4 && argc != 6)
		return printf("Add additive white Gaussian noise to an image sequence\n\n"
		       "USAGE: awgn sigma input-path output-path [first-frame last-frame]\n\n"
		       "- Paths for video frames have to be given in printf format.\n"
		       "- First and last frame if not provided are set to zero (for an image)\n"),
				 EXIT_FAILURE;

	// parse command line
	float sigma = atof(argv[1]);
	const std::string input_path(argv[2]);
	const std::string noisy_path(argv[3]);
	int first_frame = (argc > 4) ? atoi(argv[4]) : 0;
	int last_frame  = (argc > 5) ? atoi(argv[5]) : 0;

	// check inputs
	if (input_path == "")
		return fprintf(stderr, "%s: no input sequence.\nTry `%s' for more "
				"information.\n", argv[0], argv[0]), EXIT_FAILURE;

	if (noisy_path == "")
		return fprintf(stderr, "%s: no output sequence.\nTry `%s' for more "
				"information.\n", argv[0], argv[0]), EXIT_FAILURE;

	// load video, add noise and save
	Video<float> input;
	input.loadVideo(input_path, first_frame, last_frame);
	VideoUtils::addNoise(input, input, sigma);
	input.saveVideo(noisy_path, first_frame);
	return EXIT_SUCCESS;
}
