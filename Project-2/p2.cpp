/*
	B490/B659 Project 2 Skeleton Code    (2/2015)
	
	Be sure to read over the project document and this code (should you choose to use it) before 
	starting to program. 
	

	Compiling:
		A simple console command 'make' will excute the compilation instructions
		found in the Makefile bundled with this code. It will result in an executable
		named p2.

	Running:
		The executable p2 should take commands of the form:
			./p2 problem_ID input_File ouput_File additional_Arguments TODO
	
*/


//Link to the header file
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <fft.h>

/************* Parameters for Watermarking *************/

#define PI 3.14159265359
#define alpha 10
//#define radius 128
#define vlength 10
#define ccThreshold 0.3

/************* Parameters for Watermarking *************/

//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;


class ComplexClass
{
    public:
        double real;
        double imag;
        ComplexClass(double real, double imag)
        {
            this->real=real;
            this->imag=imag;
        }
        
        void add(ComplexClass d)
        {
        	this->real+=d.real;
        	this->imag+=d.imag;
        }
};

void resizeToSquare(CImg<double> &inp)
{	
	int n = inp.width();
	if(n<inp.height())
	{
		n=inp.height();
	}
	
	//find the next highest number which is a power of 2..
	int p =1;
	int m=pow(2, p);

	while(m<n)
	{
		p++;
		m=pow(2, p);
	}
	
	inp.resize(m, m, 1, 1, 0);
	return;
}

ComplexClass multiplyComplex(const ComplexClass c1, const ComplexClass c2)
{
    ComplexClass c(c1.real*c2.real-c1.imag*c2.imag, c1.real*c2.imag+c2.real*c1.imag);
    return c;    
}

// This code requires that input be a *square* image, and that each dimension
//  is a power of 2; i.e. that input.width() == input.height() == 2^k, where k
//  is an integer. You'll need to pad out your image (with 0's) first if it's
//  not a square image to begin with. (Padding with 0's has no effect on the FT!)
//
// Forward FFT transform: take input image, and return real and imaginary parts.
//
void fft(const CImg<double> &input, CImg<double> &fft_real, CImg<double> &fft_imag)
{
  fft_real = input;
  fft_imag = input;
  fft_imag = 0.0;

  FFT_2D(1, fft_real, fft_imag);
}

// Inverse FFT transform: take real and imaginary parts of fourier transform, and return
//  real-valued image.
//
void ifft(const CImg<double> &input_real, const CImg<double> &input_imag, CImg<double> &output_real)
{
  output_real = input_real;
  CImg<double> output_imag = input_imag;

  FFT_2D(0, output_real, output_imag);
}

CImg<double> getGaussianFilter(int k, double sigma) 
{	
	CImg<double> H(2*k+1, 2*k+1, 1, 1, 0);
	//init the gaussian filter..	
	double sigmaSquared = sigma * sigma;	
	double factor = 1.0 / (2 * M_PI * sigma * sigma);
	double sum=0;
	for (int i = -k; i <= k; i++) 
	{
		double xVal = exp((-1.0 * i * i) / (2 * sigmaSquared));
		for (int j = -k; j <= k; j++) 
		{
			double yVal = exp((-1.0 * j * j) / (2 * sigmaSquared));
			H(i + k, j + k, 0, 0) = xVal * yVal * factor;
			sum+=H(i + k, j + k, 0, 0);
		}
	}
	
	///*
	//normalize the values..
	for (int i = -k; i <= k; i++) 
	{
		for (int j = -k; j <= k; j++) 
		{
			H(i + k, j + k, 0, 0) = H(i + k, j + k, 0, 0)/sum;
		}
	}	
	//*/	
	
	return H;
}

//Problem: 2
CImg<double> morphing(CImg<double> &inp1, CImg<double> &inp2, double sigma, double sigma2)
{
	
	CImg<double> out(inp1.width(), inp1.height(), 1, 1, 0); 
	
	int f1 = 3*sigma;
	int f2 = 3*sigma2;
	
	CImg<double> g1=getGaussianFilter( f1, sigma);
	CImg<double> g2=getGaussianFilter( f2, sigma2);
	
	//make back up;
	CImg<double> inp2_bk(inp2);
	CImg<double> m1 = inp1.convolve(g1); 
	CImg<double> m2 = inp2.convolve(g2); 
	
	for(int i=0; i <inp1.width(); i++)
	{
		for(int j=0; j<inp1.height(); j++)
		{
			double val = inp2_bk(i, j, 0, 0) - inp2(i, j, 0, 0);		
			m2(i, j, 0, 0) = val;
			val+=inp1(i, j, 0, 0);
			if(val>255)
			{
				val=255;
			}
			else if(val<0)
			{
				val=0;
			}	
			
			out(i, j, 0, 0) = val;
		}
	}
	return out;
}

// Problem: 3.1
CImg<double> fft_magnitude(const CImg<double> &fft_real, const CImg<double> &fft_imag)
{
	CImg<double> out(fft_real.width(), fft_real.height(), 1, 1, 0); 
		
	for(int i=0; i<fft_real.width(); i++)
	{
		for(int j=0; j<fft_real.width(); j++)
		{
			double val = log(sqrt(fft_real(i, j, 0, 0)*fft_real(i, j, 0, 0) + fft_imag(i, j, 0, 0)*fft_imag(i, j, 0, 0)));
			out(i, j, 0, 0) = val;
		}	
	}
	
	return out;
}


void saveSpectrum(CImg<double> inp, const char* name)
{
	CImg<double> r(inp);
	CImg<double> i(inp);
	
	fft(inp, r, i);
	CImg<double> out = fft_magnitude(r, i);
	out.save(name);
}


//Problem: 3.2
CImg<double> remove_interference(const CImg<double> &inp)
{
    CImg<double> out(inp);
    CImg<double> real(inp), imag(inp);
    fft(inp, real, imag);
    
    //remove the left top error..
    //150-165 in x..
    //150-165 in y..    
    for(int i=150; i<165; i++)
    {
        for(int j=150; j<165; j++)
        {
        	real(i, j, 0, 0) =0;
        }
    }   

    //remove the right bottom error..
    //350-365 in x..
    //350-365 in y..    
    for(int i=350; i<365; i++)
    {
        for(int j=350; j<365; j++)
        {
        	real(i, j, 0, 0) =0;
        }
    }
    
    //do the inverse fft..
    ifft(real, imag, out);
    return out;
}

// Write this in Part 4 -- check if watermark N is in image
bool check_image(const CImg<double> &input, int N);


//convolve using fft and ifft
CImg<double> fftConvolve(CImg<double> inp1, CImg<double> inp2)
{
	CImg<double> real1(inp1);
	CImg<double> imag1(inp1); 
	CImg<double> real2(inp2);
	CImg<double> imag2(inp2);

	fft(inp1, real1, imag1);
	fft(inp2, real2, imag2);

    CImg<double> out_real(real1);
    CImg<double> out_imag(real1);    
    CImg<double> out(real1);    
    
    CImg<double> t1=fft_magnitude(real1, imag1);
    CImg<double> t2=fft_magnitude(real2, imag2);
    
    for(int i=0; i<out_imag.width(); i++)
    {
    	for(int j=0; j<out_imag.height(); j++)
    	{    		    		
			ComplexClass c1(real1(i, j, 0, 0), imag1(i, j, 0, 0));
			ComplexClass c2(real2(i, j, 0, 0), imag2(i, j, 0, 0));
    		    		
    		ComplexClass c=multiplyComplex(c1, c2);
			out_real(i, j, 0, 0) = c1.real*c2.real - c1.imag*c2.imag;
			out_imag(i, j, 0, 0) = c2.imag*c1.real+c1.imag*c2.real;
    	}
    }
    CImg<double> test = fft_magnitude(out_real, out_imag);
    //do the inverse transform..
    ifft(out_real, out_imag, out);
    return out;
}



//Problem: 3.3
CImg<double> fftConvolveGaussian(CImg<double> inp, double sigma, bool getOriginal = false)
{
	int origWidth = inp.width();
	int origHeight = inp.height();

	//construct the gaussian filter..	
	int filterSize = 3*sigma;
	int size = 2*filterSize+1;
	CImg<double> gauss=getGaussianFilter(filterSize, sigma);					;
	
	inp.resize(inp.width()+size, inp.height()+size, 1, 1, 0);
	resizeToSquare(inp);	
		
	//this is the standard resize logic..	
	gauss.resize(inp.width(), inp.width(),1, 1, 0);		//resize the gaussian image..	
	
	CImg<double> temp = fftConvolve(inp, gauss);
	
	if(getOriginal == true)
	{
		//translate the original image to the right location..
		CImg<double> out(origWidth, origHeight, 1, 1, 0);			
		for(int i=0; i<origWidth; i++)
		{
			for(int j=0; j<origHeight; j++)
			{
				out(i, j, 0, 0) = temp(i+filterSize, j+filterSize, 0, 0);
			}
		}
		temp = out;
		temp.normalize(0, 255);
	}			
	
	return temp;
}



CImg<double> getMeanFilter(int filterSize)
{
	CImg<double> t(filterSize, filterSize, 1, 1, 0);
	double val = 1/((double)filterSize*filterSize);
	for(int i=0; i<filterSize; i++)
	{
		for(int j=0; j<filterSize; j++)
		{
			t(i, j, 0, 0) = val;
		}
	}	
	return t;
}


CImg<double> inverseConvolve(CImg<double> inp1, CImg<double> inp2)
{
    CImg<double> real1(inp1);
    CImg<double> imag1(inp1);
    
    CImg<double> real2(inp2);
    CImg<double> imag2(inp2);    

    fft(inp1, real1, imag1);
    fft(inp2, real2, imag2);
    
    CImg<double> real(inp1);
    CImg<double> imag(inp1);
    CImg<double> out(inp1);
    bool isError=false;
    //do the complex number division..
    for(int i=0; i<inp1.width(); i++)
    {
        for(int j=0; j<inp1.height(); j++)
        {
            ComplexClass c1(real1(i, j, 0, 0,0), imag1(i, j, 0, 0));
            ComplexClass c2(real2(i, j, 0, 0,0), imag2(i, j, 0, 0));
            
            double d = pow(c2.real, 2)+pow(c2.imag,2);   
            if(d!=0)
            {
            	
		    	real(i, j,  0,  0) = (c1.real*c2.real + c1.imag*c2.imag)/d;
		        imag(i, j, 0, 0) = (c1.imag*c2.real-c1.real*c2.imag)/d;	
        	}
        	else
        	{
        		//TODO: what to do here??
        		cout<<"Kernel is not invertible !!! "<<endl;
        		
            	cout<<i<<"\t"<<j<<endl;
            	real(i, j,  0,  0) = c1.real;
            	imag(i, j, 0, 0) = c1.imag;
        		isError = true;
        		break;
        	}
        }
        
        if(isError)
        {
        	break;
        }
    }
    ifft(real, imag, out); //do the inverse fft....
    return out;
}


//Problem 3.4
CImg<double> inverseGaussian(CImg<double> inp, double sigma)
{
	//construct the gaussian filter..	
	int filterSize = 3*sigma;
	CImg<double> gauss = getGaussianFilter(filterSize, sigma);					
	
	gauss.resize(inp.width(), inp.width(),1, 1, 0);		//resize the gaussian image..	
	
	return inverseConvolve(inp, gauss);
}

//Problem 3.4a mean filter in question 3.4
CImg<double> inverseMean(CImg<double> inp, double size)
{
	//construct the gaussian filter..	
	CImg<double> mean = getMeanFilter(size);						
	mean.resize(inp.width(), inp.width(),1, 1, 0);		//resize the gaussian image..	
	
	return inverseConvolve(inp, mean);
}


//Problem: 5
CImg<double> fftConvolveMean(CImg<double> inp, int filterSize)
{
	
	//construct the gaussian filter..		
	CImg<double> l=getMeanFilter(filterSize);
	l.resize(inp.width(), inp.width(),1, 1, 0);		//resize the gaussian image..	
		
	return fftConvolve(inp, l);
}

CImg<double> mark_image(const CImg<double> &input, int N);


void transposeImage(CImg<double> &check)
{

	for(int i=0; i<check.width()/2; i++)
	{
		for(int j=0; j<check.height(); j++)
		{
			//swap to cause reflection...
			double temp = check(i, j, 0, 0);					
			check(i, j, 0, 0) = check(check.width()-i-1, j, 0, 0);
			check(check.width()-i-1, j, 0, 0) = temp;
		}
	}
	
	//reflect in the x direction..
	for(int i=0; i<check.width(); i++)
	{
		for(int j=0; j<check.height()/2; j++)
		{
			//swap to cause reflection...
			double temp = check(i, j, 0, 0);					
			check(i, j, 0, 0) = check(i, check.height()-j-1, 0, 0);
			check(i, check.height()-j-1, 0, 0) = temp;
		}
	}
}


void resizeKernel(CImg<double> &check)
{
	double left=0; double right=check.width()-1;
	double top=0; double bottom = check.height()-1;
	
	//get the nonzero left..
	for(int i=0; i<check.width(); i++)
	{
		bool isZero=true;
		for(int j=0; j<check.height(); j++)
		{
			if(check(i, j, 0, 0)!=0)
			{
				isZero=false;
				break;
			}
		}
		if(isZero)
		{
			
			left++;
		}
		else
		{
			break;
		}				
	}
	
	//get the non-zero right..
	for(int i=check.width()-1; i>=0; i--)
	{
		bool isZero=true;
		for(int j=0; j<check.height(); j++)
		{
			if(check(i, j, 0, 0)!=0)
			{
				isZero=false;
				break;
			}
		}
		
		if(isZero)
		{
			right--;
		}
		else
		{
			break;
		}				
	}
	
	//get the nonzero top..
	for(int j=0; j<check.height(); j++)
	{
		bool isZero=true;
		for(int i=0; i<check.width(); i++)
		{
			if(check(i, j, 0, 0)!=0)
			{
				isZero=false;
				break;
			}
		}
		
		if(isZero)
		{
			top++;
		}
		else
		{
			break;
		}				
	}
	
	//get the non-zero bottom..
	for(int j=check.height()-1; j>=0; j--)
	{
		bool isZero=true;
		for(int i=0; i<check.width(); i++)
		{
			if(check(i, j, 0, 0)!=0)
			{
				isZero=false;
				break;
			}
		}
		
		if(isZero)
		{
			bottom--;
		}
		else
		{
			break;
		}				
	}
	
	if(right-left+1<=0 || bottom-top+1<=0)
	{
	    cout<<"OMGGGGGGGGGGGGGGGGGGGGG!!"<<endl;
	    return;
	}
	
	CImg<double> modKernel(right-left+1, bottom-top+1, 1, 1, 0);
	
	for(int i=left; i<=right; i++)
	{
		for(int j=top; j<=bottom; j++)
		{
			modKernel(i-left, j-top, 0, 0) = check(i, j, 0, 0);
		}
	}
	
	check=modKernel;	
    cout<<"Pruning success"<<endl;
}


void normalize(CImg<double> &check)
{
	double sum=0;
	//normalize the kernel....
	for(int i=0; i<check.width(); i++)
	{
		for(int j=0; j<check.height(); j++)
		{
			sum+=check(i, j, 0, 0); 
		}
	}
	
	for(int i=0; i<check.width(); i++)
	{
		for(int j=0; j<check.height(); j++)
		{
			check(i, j, 0, 0)/=sum;
		}
	}	
}


string getStringRep(int n)
{
    if(n==0)
    {
        return "0";
    }

	string rep="";
	while(n !=0)
	{
		int p = n%10;
		
		rep = (char)(p+48)+rep;
		n=n/10;
	}
	return rep;
}
 
int main(int argc, char **argv)
{
    try 
    {

        if(argc < 4)
        {
            cout << "Insufficent number of arguments; correct usage:" << endl;
            cout << "    p2 problemID inputfile outputfile" << endl;
            return -1;
        }

        string part = argv[1];
        string inputFile = argv[2];
        //string outputFile = argv[3];
        //cout << "In: " << inputFile <<"  Out: " << outputFile << endl;

        CImg<double> input_image(inputFile.c_str());

        if(part == "2")
        {
            // do something here!
            
            if(argc < 7)
            {
                cout << "Insufficent number of arguments; correct usage:" << endl;
                cout << "p2 problemID inputfile1 inputfile2 outFile sigma1 sigma2" << endl;                
                return -1;
            }
        	CImg<double> inp1(argv[2]);
        	CImg<double> inp2(argv[3]);           
            double sigma =  atof(argv[5]);
			double sigma2 = atof(argv[6]);
            //can both input images be of different size ??? TODO
            if(inp1.width()!=inp2.width() || inp1.height()!=inp2.height())
            {
            	cout<<"INPUT ERROR: Input images are not of equal size"<<endl;
            	return -1;
            }
            
            //do error checking..
    		if (inp1.spectrum() != 1 || inp2.spectrum() != 1) 
    		{
				cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
				return -1;
			}		
            
			CImg<double> output = morphing(inp1, inp2, sigma, sigma2);
			output.save(argv[4]);            
        }
        else if(part == "3.1")
        {
            if(argc < 4)
            {
                cout << "Insufficent number of arguments; correct usage:" << endl;
                cout << "    p2 problemID inputfile1 outFile" << endl;                
                return -1;
            }
        	
        	CImg<double> inp(argv[2]);
        	resizeToSquare(inp);
        	
        	CImg<double> real(inp, "xyzc", 0);
        	CImg<double> imag(inp, "xyzc", 0);
        	        	
        	fft(inp, real, imag);
			CImg<double> output = fft_magnitude(real, imag);
			output.save(argv[3]);       
        }
        else if(part == "3.2")
        {
            if(argc < 4)
            {
                cout << "Insufficent number of arguments; correct usage:" << endl;
                cout << "    p2 problemID inputFileName outFileName" << endl;                
                return -1;
            }
        	CImg<double> inp(argv[2]);
        	resizeToSquare(inp);
			CImg<double> output = remove_interference(inp);
			output.save(argv[3]);      
        }
        else if(part == "3.3")
        {
            if(argc < 5)
            {
                cout << "Insufficent number of arguments; correct usage:" << endl;
                cout << "    p2 problemID inputFileName outFileName sigma" << endl;                
                return -1;
            }
                        
        	CImg<double> inp(argv[2]);             	
        	double sigma = atof(argv[4]);        	        	        	
			
			CImg<double> output = fftConvolveGaussian(inp, sigma, true);									
						
			output.save(argv[3]);  
        }
        else if (part=="3.4")
        {            
            if(argc < 5)
            {
                cout << "Insufficent number of arguments; correct usage:" << endl;
                cout << "    p2 problemID inputFileName outFileName sigma" << endl;                
                return -1;
            }
        	CImg<double> inp(argv[2]);
        	double sigma = atof(argv[4]);
        	
        	int w = inp.width();
        	int h = inp.height();
        	        	
        	//apply the gaussian filter first..
			CImg<double> temp = fftConvolveGaussian(inp, sigma);		
			
        	//apply the inverse convolution..
        	CImg<double> output = inverseGaussian(temp, sigma);
        	
        	//resize the image to its original dimensions..
			CImg<double> t(w, h, 1, 1, 0);
			
			for(int i=0; i<w; i++)
			{
			    for(int j=0; j<h; j++)
			    {
			        t(i, j, 0, 0) = output(i, j, 0, 0);
			    }
			}
			
			t.normalize(0, 255);       	        	
        	t.save(argv[3]);
        }
        else if (part=="3.4a")
        {            
            if(argc < 5)
            {
                cout << "Insufficent number of arguments; correct usage:" << endl;
                cout << "    p2 problemID inputFileName outFileName size" << endl;                
                return -1;
            }
        	CImg<double> inp(argv[2]);
        	double size = atof(argv[4]);
        	
        	int w = inp.width();
        	int h = inp.height();
        	resizeToSquare(inp);
        	        	
        	//apply the gaussian filter first..
			CImg<double> temp = fftConvolveMean(inp, size);		
        	//apply the inverse convolution..
        	CImg<double> output = inverseMean(temp, size);
        	
        	//resize the image to its original dimensions..
			CImg<double> t(w, h, 1, 1, 0);
			
			for(int i=0; i<w; i++)
			{
			    for(int j=0; j<h; j++)
			    {
			        t(i, j, 0, 0) = output(i, j, 0, 0);
			    }
			}
			
			t.normalize(0, 255);       	        	
        	t.save(argv[3]);
        }
        else if(part == "4")
        {
			if(argc < 5)
			{
				cout << "Insufficent number of arguments; correct usage:" << endl;
				cout << "    p2 problemID inputFileName outFileName Key" << endl;                
				return -1;
			}
		
			string todo = "";		
			cout << "Do you want to add watermark in \"" << argv[2] << "\" or check watermark in \"" << argv[3] << "\"? (mark / check) :" << endl;
			getline(cin, todo);
		
			if(todo == "mark") {
			    CImg<double> inp(argv[2]);		
			    if(inp.width() != inp.height()){
			        cout << "Please provide and square image with 2^n pixels" << endl;
			        return -1;
			    }
			    int x = 0;
			    int s = inp.width();
			    bool power = false;
			    while(x <= s) {
			        if(pow(2,x) == s) {
			            power = true;
			            break;
			        }
			        else
			            x++;
			    }
			    if(!power) {
			        cout << "Please provide and square image with 2^n pixels" << endl;
			        return -1;
			    }
				CImg<double> output = mark_image(inp, atoi(argv[4]));
				output.save(argv[3]);
			
				cout << "Watermark added" << endl;
			} else if(todo == "check") {
			    CImg<double> inp(argv[3]);		
			
				if(check_image(inp, atoi(argv[4]))) {
					cout << "Watermark is present" << endl;
				} else {
					cout << "No watermark found" << endl;
				}
			} else {
				cout << "Invalid input." << endl;
			}
        }
        else if(part == "5a")
        {
			//create 20 images with different values in the middle.. 
			int rounds = 20;
						
			//mark the (1, 1) position//
			int num = 20;
			for(int i=0; i< rounds; i++)
			{
				CImg<double> test(255, 255, 1, 3, 0);
				test(127, 127, 0, 0) =10*(i+1);
				
				string name = "noisfy_trials/inp";
				name += getStringRep(i);
				name+=".jpg";				
				test.save(name.c_str());
			}
											
			return 0;        
        }
        else if(part == "5b")
        {
			if(argc < 4)
			{
				cout << "Insufficent number of arguments; correct usage:" << endl;
				cout << "    p2 problemID inputFileName outFileName" << endl;                
				return -1;
			}
			cout<<"what??"<<endl;
			int rounds = 20;
			
			CImg<double> kernel(255, 255, 1, 1, 0);
						
			for(int i=0; i< rounds; i++)
			{												
				string name = "noisfy_trials/out";
				name += getStringRep(i);
				name+=".jpg";				
				CImg<double> test(name.c_str());
			    
			    for(int p=0; p<255; p++)
			    {
			        for(int q=0; q<255; q++)
			        {
			            kernel(p, q, 0, 0) += test(p, q, 0, 0);
			        }
			    }
			}
			
			normalize(kernel);
			kernel.normalize(0, 255);						
			
			//get the kernel..
			//check the spectrum of the image..
			CImg<double> check = kernel;
			cout<<"Before resize"<<endl;
			check.normalize(0, 255);
			
			for(int i=0; i<check.width(); i++)
			{
			    for(int j=0; j<check.height(); j++)
			    {
    			    if(i==0)
			        {
			            check(i, j, 0, 0) =0;			        
			        }
			        //removing noise..
			        else if(j<4)
			        {
			            check(i, j, 0, 0) =0;
			        }			  
			    }
			}
						
			//removing some noise in the image..
			for(int i=0; i<check.width(); i++)
			{
			    for(int j=0; j<check.height(); j++)
			    {
                    if(check(i, j, 0, 0 ) <20)
                    {
                        check(i, j, 0, 0 ) =0;
                    }
			    }
			}
			resizeKernel(check);
			
			//transposeImage(check);
			normalize(check);	
						
        	CImg<double> inp(argv[2]);
        	//inp.normalize(-1,1);
        	int w = inp.width();
        	int h = inp.height();
        	
        	inp.resize(w+check.width(), h+check.height(), 1, 1, 0);        	    		
        	resizeToSquare(inp);

    		check.resize(inp.width(), inp.height(), 1, 1, 0);    		

        	//apply the inverse convolution..
        	CImg<double> output = inverseConvolve(inp, check);
        	//CImg<double> output = inverseMean(inp, 3);		
        	
        	double min = 300;
        	for(int i=0; i<output.width(); i++)
        	{
        		for(int j=0; j<output.height();j++)
        		{
        			double val = output(i, j, 0, 0);
        			if(val < min)
        			{
        				min=val;
        			}
        		}
        	}
			output.normalize(0, 255); 						      	        	
        	output.save(argv[3]);						
			saveSpectrum(output, "f.jpg");
			return 0; 	
					 
        }                               
        else if(part == "5")
        {
			if(argc < 4)
			{
				cout << "Insufficent number of arguments; correct usage:" << endl;
				cout << "    p2 problemID inputFileName outFileName" << endl;                
				return -1;
			}

			CImg<double> blur("blur_standard.png");
			CImg<double> deblur("deblur_standard.png");

			
			CImg<double> kernel=inverseConvolve(blur, deblur);
			
			kernel.normalize(0, 255);
			kernel.save("omgkernel.jpg");
			
			normalize(kernel);
			
			
        	CImg<double> inp(argv[2]);
			CImg<double> out=inverseConvolve(inp, kernel); 
			out.normalize(0, 255);
			out.save(argv[3]);
			return 0;					 
        }
    } 
    catch(const string &err) 
    {
        cerr << "Error: " << err << endl;
    }
}

CImg<double> mark_image(const CImg<double> &input, int N) {
	CImg<double> real(input, "xyzc", 0);
	CImg<double> markedReal(input, "xyzc", 0);
	CImg<double> imag(input, "xyzc", 0);
	
	srand(N);
	
	int v[vlength];
	for(int i = 0; i < vlength; i++) {
	    v[i] = rand() % 2;
	}
	
	fft(input, real, imag);
	
	int cx = (int) real.width() / 2;
	int cy = (int) real.height() / 2;
	int radius = (int) real.width() / 4;
	double degree = 360 / (double) vlength;
	
	for(int i = 0; i < vlength; i++) {
	    	int x1 = (int) (cx + radius * cos(degree * PI / 180.0));
		int y1 = (int) (cy + radius * sin(degree * PI / 180.0));
		real(x1, y1, 0, 0) = real(x1, y1, 0, 0) + alpha * abs(real(x1, y1, 0, 0)) * v[i];
		
		int x2 = (int) (cx + radius * cos((degree * PI / 180.0) + 180.0));
		int y2 = (int) (cy + radius * sin((degree * PI / 180.0) + 180.0));
		real(x2, y2, 0, 0) = real(x2, y2, 0, 0) + alpha * abs(real(x2, y2, 0, 0)) * v[i];
		
		degree += 360 / (double) vlength;
	}
	
	ifft(real, imag, markedReal);
	
	return markedReal;
}

bool check_image(const CImg<double> &input, int N) {
	srand(N);
	
	int v[vlength];
	double c[vlength];
	for(int i = 0; i < vlength; i++) {
	    v[i] = rand() % 2;
	}
	
	CImg<double> real(input, "xyzc", 0);
	CImg<double> imag(input, "xyzc", 0);
	
	fft(input, real, imag);
	
	int cx = (int) real.width() / 2;
	int cy = (int) real.height() / 2;
	int radius = (int) real.width() / 4;
	double degree = 360 / (double) vlength;
	double vAvg = 0.0, cAvg = 0.0;
	for(int i = 0; i < vlength; i++) {
		int x = (int) (cx + radius * cos(degree * PI / 180.0));
		int y = (int) (cy + radius * sin(degree * PI / 180.0));
		c[i] = real(x, y, 0, 0);
		degree += 360 / (double) vlength;
		vAvg += v[i];
		cAvg += c[i];
	}
	
	vAvg /= vlength;
	cAvg /= vlength;
	double numerator = 0.0, denomV = 0.0, denomC = 0.0;
	for(int i = 0; i < vlength; i++) {
		numerator += (v[i] - vAvg) * (c[i] - cAvg);
		denomV += (v[i] - vAvg) * (v[i] - vAvg);
		denomC += (c[i] - cAvg) * (c[i] - cAvg);		
	}
	
	bool marked = false;
	if(numerator / (sqrt(denomV) * sqrt(denomC)) > ccThreshold) {
	    marked = true;
	}
	return marked;
}

