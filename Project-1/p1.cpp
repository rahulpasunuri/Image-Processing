/*
 B490/B659 Project 1 Skeleton Code (1/2015)
 Be sure to read over the project document and this code (should you choose to use it) before
 starting to program.

 Compiling:
 A simple console command 'make' will excute the compilation instructions
 found in the Makefile bundled with this code. It will result in an executable
 named p1.

 Running:
 The executable p1 takes commands of the form:
 ./p1 problem_ID input_File ouput_File additional_Arguments
 Some examples:

 ./p1 2.1 input.png out.png
 This runs the 'averageGrayscale' function defined in this file and described
 in the project doc Part 2, problem 1 on the input. The output is saved as out.png.
 ----

 ./p1 4.2b input.jpg output.png 0.5

 This runs the Gaussian filter function (gaussianFilter) with sigma = 0.5 on input.jpg.
 This is problem 2b from Part 4 in the project documentation.

 */

//Link to the header file
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <random>
#include<vector>
#include<algorithm>
#define _USE_MATH_DEFINES
#include <math.h>

//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;

CImg<double> fastGaussian(CImg<double> input, int numCon, double sigma)
{
	
	CImg<double> out(input); //init the output image with input..
	int width = pow(((12*sigma*sigma)/(double)numCon)+1, 0.5);		
	cout<<width<<endl;
	//if(width%2==0)
	//{
		//width+=1; //make width odd..
	//}	
	//width=1; //TODO
	for(int i=0; i<numCon; i++)
	{
		CImg<double> temp(out);
		///*
		for(int p=0; p<input.height(); p++)
		{
			//apply row filter here..
			//compute the oth value..
			double val=out(0, p, 0, 0);
			for(int k=0; k<width-1; k++)
			{
				val+= out(k, p, 0, 0); //handle reflection..
			}	
			
			temp(0, p, 0, 0) = val;	 
			//temp(0, p, 0, 0) = out(0, p, 0, 0);	 
			//column other values in the row..
			for(int k=1; k<input.width(); k++)
			{
				int lIndex = k-width;
				if(lIndex<0)
				{
					lIndex = lIndex*-1 - 1;
				}
				val = temp(k-1, p, 0, 0) + out(k, p, 0, 0) - out(lIndex, p, 0, 0);
				
				temp(k, p, 0, 0) = val;							
			}	
			
			//copy the temp image to output image..
			for(int k=0; k<input.width(); k++)
			{
				out(k, p, 0, 0) = temp(k, p, 0, 0)/width;
			}							
		}
		//*/
		///*
		for(int p=0; p<input.width(); p++)
		{
			//apply column filter here..
			//apply row filter here..
			//compute the oth value..
			double val=out(p, 0, 0, 0);
			for(int k=0; k<width-1; k++)
			{
				val+= out(p, k, 0, 0); //handle reflection..
			}	

			temp(p, 0, 0, 0) = val;	 
			//temp(p, 0, 0, 0)=out(0, p, 0, 0);
			
			//column other values in the column..
			for(int k=1; k<input.height(); k++)
			{
				int lIndex = k-width;
				if(lIndex<0)
				{
					lIndex = lIndex*-1 - 1;
				}
				val = temp( p, k-1, 0, 0) + out(p, k, 0, 0) - out(p, lIndex, 0, 0);
				temp(p, k, 0, 0) = val;				
			}	
			
			//copy the temp image to output image..
			for(int k=0; k<input.height(); k++)
			{
				out(p, k, 0, 0) = temp(p, k, 0, 0)/width;
			}	
		}
		//*/					
	}		
	
	for(int m=0; m<input.width(); m++)
	{
		for(int n=0; n<input.height(); n++)
		{
			double val = out(m, n, 0, 0);
			if(val>255)
			{
				val=255;
			}
			else if(val<0)
			{
				val=0;
			}
			out(m, n, 0, 0)=val;
		}
	}
			
	return out;
}

//Part 2 - Basic image operations
CImg<double> averageGrayscale(CImg<double> input) {
	int width = input.width();
	int height = input.height();
	CImg<double> out(input.width(), input.height(), 1, 1);

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			double r = input(i, j, 0, 0); //value for red..
			double g = input(i, j, 0, 1); //va;ue for green
			double b = input(i, j, 0, 2); //value for blue..

			out(i, j, 0, 0) = 0.3 * r + 0.6 * g + 0.1 * b;
		}
	}

	return out; // out is the final image
}

CImg<double> simpleBW(CImg<double> input);
CImg<double> advancedBW(CImg<double> input);

//Part 3 - Adding noise
CImg<double> uniformNoise(CImg<double> input);
CImg<double> gaussianNoise(CImg<double> input, double sigma);
CImg<double> saltAndPepperNoise(CImg<double> input);

//Part 4 - Filtering
CImg<double> filter(CImg<double> input, CImg<double> filter);
CImg<double> meanFilter(CImg<double> input, int filterSize);
CImg<double> gaussianFilter(CImg<double> input, double sigma);
CImg<double> medianFilter(CImg<double> input, int size);

double meanSquareError(CImg<double> origInput, CImg<double> noiseNFilter);

CImg<double> separableMeanKernel(CImg<double> input, int filterSize);
CImg<double> separableGaussianKernel(CImg<double> input, double filterSize);
CImg<double> separableMeanAdaptiveKernel(CImg<double> input, int filterSize);
CImg<double> adaptiveFilter(CImg<double> input, CImg<double> h);

int main(int argc, char **argv) {

	if (argc < 4) {
		cout << "Insufficent number of arguments. Please see documentation"
				<< endl;
		cout << "p1 problemID inputfile outputfile" << endl;
		return -1;
	}

	char* inputFile = argv[2];
	char* outputFile = argv[3];
	cout << "In: " << inputFile << " Out: " << outputFile << endl;

	CImg<double> input(inputFile);
	CImg<double> output;

	if (!strcmp(argv[1], "2.1")) {
		cout << "# Problem 2.1 - Average Grayscale" << endl;
		if (input.spectrum() != 3) {
			cout << "INPUT ERROR: Input image is not a color image!" << endl;
			return -1;
		}
		output = averageGrayscale(input);
		output.save(outputFile);
	} else if (!strcmp(argv[1], "2.2a")) {
		cout << "# Problem 2.1a - Simple Threshold Black and White" << endl;
		if (input.spectrum() != 1) {
			cout << "INPUT ERROR: Input image is not a grayscale image!"
					<< endl;
			return -1;
		}
		output = simpleBW(input);
		output.save(outputFile);
	} else if (!strcmp(argv[1], "2.2b")) {
		cout << "# Problem 2.2b - Advanced Threshold Black and White" << endl;
		if (input.spectrum() != 1) {
			cout << "INPUT ERROR: Input image is not a grayscale image!"
					<< endl;
			return -1;
		}
		output = advancedBW(input);
		output.save(outputFile);
	} else if (!strcmp(argv[1], "3.1")) {
		cout << "# Problem 3.1 - Uniform Noise" << endl;
		if (input.spectrum() != 1) {
			cout << "INPUT ERROR: Input image is not a grayscale image!"
					<< endl;
			return -1;
		}
		output = uniformNoise(input);
		output.save(outputFile);
	} else if (!strcmp(argv[1], "3.2")) {
		cout << "# Problem 3.2 - Gaussian Noise" << endl;
		if (input.spectrum() != 1) {
			cout << "INPUT ERROR: Input image is not a grayscale image!"
					<< endl;
			return -1;
		}
		if (argc != 5) {
			cout << "INPUT ERROR: Provide sigma as additional argument!"
					<< endl;
			return -1;
		}
		output = gaussianNoise(input, atof(argv[4]));
		output.save(outputFile);
	} else if (!strcmp(argv[1], "3.3")) {
		cout << "# Problem 3.3 - Salt & Pepper Noise" << endl;
		if (input.spectrum() != 1) {
			cout << "INPUT ERROR: Input image is not a grayscale image!"
					<< endl;
			return -1;
		}
		output = saltAndPepperNoise(input);
		output.save(outputFile);
	} else if (!strcmp(argv[1], "4.2a")) {
		cout << "# Problem 4.2a - Mean Filter Noise" << endl;
		std::clock_t start = std::clock();
		if (input.spectrum() != 1) {
			cout << "INPUT ERROR: Input image is not a grayscale image!"
					<< endl;
			return -1;
		}
		if (argc != 5) {
			cout << "INPUT ERROR: Provide filter size as additional argument!"
					<< endl;
			return -1;
		}
		output = meanFilter(input, atoi(argv[4]));
		output.save(outputFile);
		std::cout << "Time Elapsed: "
				<< (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000)
				<< " ms" << std::endl;
	} else if (!strcmp(argv[1], "4.2b")) {
    	std::clock_t start = std::clock();
		cout << "# Problem 4.2b - Gaussian Filter " << endl;
		if (input.spectrum() != 1) {
			cout << "INPUT ERROR: Input image is not a grayscale image!"
					<< endl;
			return -1;
		}
		if (argc != 5) {
			cout << "INPUT ERROR: Provide sigma as additional argument!"
					<< endl;
			return -1;
		}
		output = gaussianFilter(input, atof(argv[4]));
		output.save(outputFile);
		std::cout << "Time Elapsed: "
		<< (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000)
		<< " ms" << std::endl;
		

	} else if (!strcmp(argv[1], "4.3")) {
		cout << "# Problem 4.3 - Median Noise" << endl;
		if (input.spectrum() != 1) {
			cout << "INPUT ERROR: Input image is not a grayscale image!"
					<< endl;
			return -1;
		}
		if (argc != 5) {
			cout << "INPUT ERROR: Provide filter size as additional argument!"
					<< endl;
			return -1;
		}
		output = medianFilter(input, atoi(argv[4]));
		output.save(outputFile);
	} else if (!strcmp(argv[1], "5")) {
		cout << "# Problem 5 - Noise Removal Analysis" << endl;
		//FIXME You will need to implement this section yourself
		//Basic approach::pass original and Image with noise convoluted with a filter to get MSE
		CImg<double> secondInput(outputFile);
		double mse = meanSquareError(input, secondInput);
		//double mse = meanSquareError(input, output);
		if (mse < 0) {
			cout << "Please enter original image and transformed image" << endl;
		} else {
			cout << "The mean square error is :: " << mse << endl;
		}

	} else if (!strcmp(argv[1], "6.1")) {
		cout << "# Problem 6.1 - Separable Kernel Convolutions" << endl;

		std::clock_t start = std::clock();
		//FIXME You will need to implement this section yourself

		if (input.spectrum() != 1) {
			cout << "INPUT ERROR: Input image is not a grayscale image!"
					<< endl;
			return -1;
		}
		if (argc != 6) {
			cout
					<< "INPUT ERROR: Provide filter size  and Filter Type  as additional argument!  :: FilterType::Sequence: (0-MeanFilter) (1-GaussianFilter) "
					<< endl;
			return -1;
		}
		int choice = atoi(argv[5]);
		//cout << "choice is " << choice << endl;

		if (choice != 0 && choice != 1) {
			cout
					<< "Error With the Filter Type :: FilterType::Sequence: (0-MeanFilter) (1-GaussianFilter) "
					<< endl;
			exit(1);
		}

		if (choice == 1)
			output = separableGaussianKernel(input, atof(argv[4]));
		else
			output = separableMeanKernel(input, atoi(argv[4]));
		output.save(outputFile);

		std::cout << "Time Elapsed: "
				<< (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000)
				<< " ms" << std::endl;

	} else if (!strcmp(argv[1], "6.2")) {
		cout << "# Problem 6.2 - Dynamic Box Filter" << endl;

		std::clock_t start = std::clock();
		//FIXME You will need to implement this section yourself

		if (input.spectrum() != 1) {
			cout << "INPUT ERROR: Input image is not a grayscale image!"
					<< endl;
			return -1;
		}
		if (argc != 5) {
			cout << "INPUT ERROR: Provide Filtersize as additional argument!"
					<< endl;
			return -1;
		}
		output = separableMeanAdaptiveKernel(input, atof(argv[4]));		
		output.save(outputFile);
		
		std::cout << "Time Elapsed: "
				<< (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000)
				<< " ms" << std::endl;

	} else if (!strcmp(argv[1], "6.3")) {
		cout << "# Problem 6.3 - Fast Gaussian Smoothing" << endl;
		if (argc != 5) 
		{
			cout << "INPUT ERROR: Provide sigma as additional argument"
					<< endl;
			return -1;
		}
		std::clock_t start = std::clock();
		output = fastGaussian(input, 4, atof(argv[4]));
		output.save(outputFile);
		std::cout << "Time Elapsed: "
				<< (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000)
				<< " ms" << std::endl;

	} 
	else {
		cout << "Unknown input command" << endl;
	}

	return 0;
}



CImg<double> simpleBW(CImg<double> input) {
//Creates a new grayscale image with same size as the input initialized to all 0s (black)
	CImg<double> output(input.width(), input.height(), 1, 1, 0);

	int width = input.width();
	int height = input.height();

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			double val = input(i, j, 0, 0);
			double out;

			if (val > 127) {
				out = 255;
			} else {
				out = 0;
			}
			output(i, j, 0, 0) = out;
		}
	}

	return output;
}

CImg<double> advancedBW(CImg<double> input) {

//Creates a new grayscale image with same size as the input initialized to all 0s (black)
	CImg<double> output(input.width(), input.height(), 1, 1, 0);

	int width = input.width();
	int height = input.height();

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {

			double val = input(i, j, 0, 0);
			double out;

			if (val > 127) {
				out = 255;
			} else {
				out = 0;
			}
			double e = val - out;
			if (j != height - 1) {
				if (i != 0) {
					input(i - 1, j + 1, 0, 0) += 0.2 * e;
				}
				if (i != width - 1) {
					input(i + 1, j + 1, 0, 0) += 0.1 * e;
				}
				input(i, j + 1, 0, 0) += 0.3 * e;
			}

			if (i != width - 1) {
				input(i + 1, j, 0, 0) += 0.4 * e;
			}

			output(i, j, 0, 0) = out;

		}
	}
	return output;
}

//Part 3 - Adding noise
CImg<double> uniformNoise(CImg<double> input) {

//Creates a new image with same size as the input initialized to all 0s (black)
	CImg<double> output(input.width(), input.height(), 1, 1, 0);

	int width = input.width();
	int height = input.height();

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			double val = input(i, j, 0, 0);

			int random = rand() % 21; // random number between 0 and 20

			int sign = rand() % 2; //sign 0 = > negative number, 1 => postive number

			if (sign == 0) {
				random *= -1;
			}

			val += random;

			if (val > 255) {
				val = 255;
			} else if (val < 0) {
				val = 0;
			}
			output(i, j, 0, 0) = val;
		}
	}

	return output;
}

CImg<double> gaussianNoise(CImg<double> input, double sigma) {

	//Creates a new image with same size as the input initialized to all 0s (black)
	CImg<double> output(input, "xyzc", 0);
	std::normal_distribution<int> distribution(0, sigma); // 0 mean and sigma standard deviation

	std::random_device rd;

	std::mt19937 e2(rd());

	int width = input.width();
	int height = input.height();

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			//cout<<"what??"<<endl;
			int random = distribution(e2);

			double val = input(i, j, 0, 0);

			val += random;

			if (val > 255) {
				val = 255;
			} else if (val < 0) {
				val = 0;
			}
			output(i, j, 0, 0) = val;

		}

	}

	return output;
}

CImg<double> saltAndPepperNoise(CImg<double> input) {

//Creates a new image with same size as the input initialized to all 0s (black)
	CImg<double> output(input, "xyzc", 0);

	int width = input.width();
	int height = input.height();

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			//TODO. .. should we take this probability from command line ???
			int random = rand() % 10; //generates numbers between 0 and 9
			output(i, j, 0, 0) = input(i, j, 0, 0);
			if (random < 2) //adds noise with 20% probability
					{
				//add noise..
				random = rand() % 2;
				if (random == 0) {
					output(i, j, 0, 0) = 0;
				} else {
					output(i, j, 0, 0) = 255;

				}
			}

		}

	}

	return output;
}

//Part 4 - Filtering
CImg<double> filter(CImg<double> input, CImg<double> h) {
	CImg<double> output(input.width(), input.height(), 1, 1, 0);
	int width = input.width();
	int height = input.height();
	int k = h.width() / 2;
	int l = h.height() / 2;

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			double val = 0;

			for (int p = -k; p <= k; p++) {
				for (int q = -l; q <= l; q++) {
					double orig;
					int modx = i - p;
					int mody = j - q;

					if (modx < 0 || mody < 0 || modx >= width
							|| mody >= height) {
						//cover out of bounds by reflection...
						if (modx < 0) {
							//reflect along y axis..
							modx = -1 * modx - 1;
						} else if (modx >= width) {
							int diff = modx - width + 1;
							modx = modx - diff;
						}
						if (mody < 0) {
							//reflect along y axis..
							mody = -1 * mody - 1;
						} else if (mody <= height) {
							int diff = mody - height + 1;
							mody = mody - diff;
						}
					}

					orig = input(modx, mody, 0, 0);
					val += (orig * h(p + k, q + l, 0, 0));
				}
			}

			//boundary checks..
			if (val > 255) {
				val = 255;
			} else if (val < 0) {
				val = 0;
			}
			output(i, j, 0, 0) = val;
		}
	}
	return output;
}

CImg<double> getMeanFilter(int filterSize, CImg<double> H) {
	//init the filter..
	for (int i = 0; i < filterSize; i++) {
		for (int j = 0; j < filterSize; j++) {
			H(i, j, 0, 0) = 1.0 / (filterSize * filterSize);
		}
	}
	return H;
}

CImg<double> meanFilter(CImg<double> input, int filterSize) {

	//Creates a new grayscale image (just a matrix) to be our filter
	CImg<double> H(filterSize, filterSize, 1, 1);
	H = getMeanFilter(filterSize, H);

	//init the filter..
	//getMeanFilter(filterSize, H);

	//Convole with filter and return
	return filter(input, H);
}

CImg<double> getGaussianFilter(int filterSize, double sigma, CImg<double>& H) {
	//init the gaussian filter..
	int k = filterSize / 2;
	double sigmaSquared = sigma * sigma;
	double factor = 1.0 / (2 * M_PI * sigma * sigma);
	for (int i = -k; i <= k; i++) {
		double xVal = exp((-1.0 * i * i) / (2 * sigmaSquared));
		for (int j = -k; j <= k; j++) {
			double yVal = exp((-1.0 * j * j) / (2 * sigmaSquared));
			H[i + k, j + k] = xVal * yVal * factor;
		}
	}
	return H;
}

CImg<double> gaussianFilter(CImg<double> input, double sigma) {

	//FIXME Determine filter size, see Part 4 2b for how
	int filterSize = 3 * sigma;

	//we only handle odd filtersizes..
	if (filterSize % 2 == 0) //if the filter size is even, then increment it
			{
		filterSize++;
	}

	//Creates a new grayscale image (just a matrix) to be our filter
	CImg<double> H(filterSize, filterSize, 1, 1);

	//init the gaussian filter..
	H = getGaussianFilter(filterSize, sigma, H);
	//Convole with filter and return
	return filter(input, H);

}

CImg<double> medianFilter(CImg<double> input, int size) {

	//Creates a new image with same size as the input initialized to all 0s (black)
	CImg<double> output(input.width(), input.height(), 1, 1, 0);

	int width = output.width();
	int height = output.height();
	int k = size / 2;
	//TODO even size is odd here??

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			vector<double> arr;

			for (int p = -k; p <= k; p++) {
				for (int q = -k; q <= k; q++) {
					int modx = i + p;
					int mody = j + q;

					if (modx < 0 || mody < 0 || modx >= width
							|| mody >= height) {
						//cover out of bounds by reflection...
						if (modx < 0) {
							//reflect along y axis..
							modx = -1 * modx - 1;
						} else if (modx >= width) {
							int diff = modx - width + 1;
							modx = modx - diff;
						}
						if (mody < 0) {
							//reflect along y axis..
							mody = -1 * mody - 1;
						} else if (mody <= height) {
							int diff = mody - height + 1;
							mody = mody - diff;
						}
					}
					arr.push_back(input(modx, mody, 0, 0));
				}
			}

			sort(arr.begin(), arr.end());
			int len = arr.size();
			double val;
			if (len & 1 == 0) //if len is even
			{
				//median is the average of the middle two values..
				val = (arr[len / 2] + arr[(len / 2) + 1]) / 2.0;
			} else //len is odd
			{
				val = arr[len / 2];
			}

			output(i, j, 0, 0) = val;
		}
	}

	return output;
}
double meanSquareError(CImg<double> origInput, CImg<double> noiseNFilter)

{
	int oWidth = origInput.width();
	int oHeight = origInput.height();
	int nFWidth = noiseNFilter.width();
	int nFHeight = noiseNFilter.height();
	//cout << oWidth << "," << oHeight << "," << nFWidth << "," << nFHeight << ","
	//		<< endl;
	//Assumption given:: original image size equals transformed image size::
	if (oWidth != nFWidth || oHeight != nFHeight) {
		cout << " Size of input and output image varies::exiting!" << endl;
		return -1.0;
	}
	if (oWidth < 1 || oHeight < 1) {

		cout << " Invalid image size " << endl;

		return -1.0;
	}

	double sme = 0.0;
	long totalPixels = oWidth * oHeight;

	for (int i = 0; i < oWidth; i++) {
		for (int j = 0; j < oHeight; j++) {
			//cout<<origInput(i,j)<<endl;
			sme += (double) pow((origInput(i, j) - noiseNFilter(i, j)), 2);
		}
	}

	return (sme / totalPixels); // mean square error
}

CImg<double> separableMeanKernel(CImg<double> input, int filterSize) {

	CImg<double> H(filterSize, filterSize, 1, 1);

	//for mean filter
//	H = getMeanFilter(filterSize, H);
	//cout << H.print() << endl;
	CImg<double> H_ROW(1, filterSize, 1, 1);
	CImg<double> H_COL(filterSize, 1, 1, 1);
	//init the Row,column filter..
	for (int i = 0; i < filterSize; i++) {

		H_ROW(0, i, 0, 0) = 1.0 / filterSize;
		H_COL(i, 0, 0, 0) = 1.0 / filterSize;
	}

	//cout << H_ROW.print() << endl << endl << H_COL.print() << endl;
	//return filter(filter(input, H_ROW), H_COL);
	return filter(input, H_COL);
}

CImg<double> separableGaussianKernel(CImg<double> input, double sigma) {

	int filterSize = sigma * 3;

	int k = filterSize / 2;

	if (filterSize % 2 == 0) {
		filterSize++;
	}

	CImg<double> H_ROW(1, filterSize, 1, 1);
	CImg<double> H_COL(filterSize, 1, 1, 1);

	double sigmaSquared = sigma * sigma;

	double factor = (1.0 / sqrt((2 * M_PI))) / (double) sigma;

	for (int i = -k; i <= k; i++) {

		double xVal = exp((-1.0 * i * i) / (2 * sigmaSquared));

		H_ROW(0, i + k, 0, 0) = (factor * xVal);
		H_COL(i + k, 0, 0, 0) = (factor * xVal);

	}

	//cout << H_ROW.print() << endl << endl << H_COL.print() << endl;
	return filter(filter(input, H_ROW), H_COL);

}

CImg<double> separableMeanAdaptiveKernel(CImg<double> input, int filterSize) {

	CImg<double> H_ROW(filterSize, 1, 1, 1);
	CImg<double> H_COL(1, filterSize, 1, 1);
	//init the Row,column filter..
	for (int i = 0; i < filterSize; i++) {

		H_ROW(i, 0, 0, 0) = 1.0 / filterSize;
		H_COL(0, i, 0, 0) = 1.0 / filterSize;
	}

	//return adaptiveFilter(adaptiveFilter(input, H_COL), H_ROW);
	//return adaptiveFilter(adaptiveFilter(input, H_ROW), H_COL);
	return adaptiveFilter(input, H_ROW);
	//return adaptiveFilter(input, H_COL);
}

CImg<double> adaptiveFilter(CImg<double> input, CImg<double> h) 
{
	CImg<double> output(input.width(), input.height(), 1, 1, 0);
	int width = input.width();
	int height = input.height();
	int k = h.width() / 2;
	int l = h.height() / 2;
	
    for(int j=0; j<height; j++)
    {
        double val=input(0, j, 0, 0);
        for(int p=1; p<=k; p++)
        {
            val+=2*input(p, j, 0, 0);
        }
        val = val*h(0,0,0,0);

        output(0, j, 0, 0) = val;
        for(int i=1; i<width; i++)
        {
			//check boundary condition for i...
			int last_i = i - k - 1;
			int next_i = i + k;

			if (last_i < 0) 
			{
				last_i = -1 * last_i - 1;
			} 
			else if (next_i >= width) 
			{
				next_i = next_i - (next_i - width + 1);
			}
			
			int orig_val_last = input(last_i, j, 0, 0);
	        int orig_val_next = input(next_i, j, 0, 0);

	        val = output(i-1, j, 0, 0) - orig_val_last * h(0, 0, 0, 0) + orig_val_next * h(0, 0, 0, 0);   

            output(i, j, 0, 0) = val;          
        }
    }
    
    for(int i=0; i<width; i++)
    {
        double val=input(i, 0, 0, 0);
        for(int p=1; p<=k; p++)
        {
            val+=2*input(i, p, 0, 0);
        }
        val *= h(0,0,0,0); 

        output(i, 0, 0, 0) = val;
        for(int j=1; j<height; j++)
        {
			int last_j = j - k-1;
			int next_j = j + k;

			if (last_j < 0) 
			{
				last_j = -1 * last_j - 1;
			} 
			else if (next_j >= height) 
			{
				next_j = next_j - (next_j - height + 1);
			}
		    int orig_val_last = input(i, last_j, 0, 0);
		    int orig_val_next = input(i, next_j, 0, 0);

		    val = output(i, j-1, 0, 0) - (orig_val_last * h(0, 0, 0, 0))+ (orig_val_next * h(0, 0, 0, 0));

		    output(i, j, 0, 0) = val;          
        }
    }
    

/*
	if (h.height() == 1) 
	{
		//adaptive for this row
		for (int j = 0; j < height; j++) 
		{
			double row_val = 0.0;
			bool isFirst = true;
			for (int i = 0; i < width; i++) 
			{
				//compute the row_sum first time by convolution kx1 on image
				if (isFirst) 
				{
				    for (int p = -k; p <= k; p++) 
				    {
					    double orig;
					    int modx = i - p;

					    //cover out of bounds by reflection...
					    if (modx < 0) 
					    {
						    //reflect along y axis..
						    modx = -1 * modx - 1;
					    } 
					    else if (modx >= width) 
					    {
						    int diff = modx - width + 1;
						    modx = modx - diff;
					    }					    
					    orig = input(modx, j, 0, 0);
					    row_val += (orig * h(0, 0, 0, 0));
			        }			        
					isFirst = false;
				} 
				else
				{

					//check boundary condition for i...
					int last_i = i - k - 1;
					int next_i = i + k;

					if (last_i < 0) 
					{
						last_i = -1 * last_i - 1;
					} 
					else if (next_i >= width) 
					{
						next_i = next_i - (next_i - width + 1);
					}
					
					int orig_val_last = input(last_i, j, 0, 0);
			        int orig_val_next = input(next_i, j, 0, 0);

			        row_val = output(i-1, j, 0, 0) - orig_val_last * h(0, 0, 0, 0) + orig_val_next * h(0, 0, 0, 0);
				}
		        if (row_val > 255) 
		        {
			        row_val = 255;
		        }
		        else if (row_val < 0)
		        {
			        row_val = 0;
		        }
		        output(i, j, 0, 0) = row_val;
			}
		}
	} 
	else if (h.width() == 1) // column(Y) filter:: vertical
	{	//adaptive for this row
		for (int i = 0; i < width; i++) 
		{
			double row_sum = 0.0;
			bool isFirst = true;
			for (int j = 0; j < height; j++) 
			{
				//compute the row_sum first time by convolution 1xk on image
				if (isFirst) 
				{

					for (int q = -l; q <= l; q++) 
					{
						double orig;
						int mody = j - q;
						if (mody < 0) 
						{
							//reflect along y axis..
							mody = -1 * mody - 1;
						}
						else if (mody <= height) 
						{
							int diff = mody - height + 1;
							mody = mody - diff;
						}
					
						orig = input(i, mody, 0, 0);
						row_sum += (orig * h(0, 0, 0, 0));

					}
					isFirst = false;
				} 
				else 
				{

					//check boundary condition for i...

					int last_j = j - k-1;
					int next_j = j + k;

					if (last_j < 0) {
						last_j = -1 * last_j - 1;
					} else if (next_j >= height) {
						next_j = next_j - (next_j - height + 1);
					}
				    int orig_val_last = input(i, last_j, 0, 0);
				    int orig_val_next = input(i, next_j, 0, 0);

				    row_sum = output(i, j-1, 0, 0) - (orig_val_last * h(0, 0, 0, 0))
						    + (orig_val_next * h(0, 0, 0, 0));
				}
			    if (row_sum > 255) 
			    {
				    row_sum = 255;
			    } 
			    else if (row_sum < 0) 
			    {

				    row_sum = 0;
			    }
			    output(i, j, 0, 0) = row_sum;

			}

		}
    }
    */
	return output;

}

