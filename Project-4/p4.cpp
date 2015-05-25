// B490/B659 Project 4 skeleton code
//
// Compile with: "make"
//
// See assignment handout for command line and project specifications.


//Link to the header file
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Sift.h>
#include <random>
#include <math.h>
#include<map>
#include <time.h>
#include <ctime>
#include<limits.h>
#include<algorithm>
//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;


//matrix code

int matrixInverse(float a[3][3], float (&out)[3][3])
{
    int i,j;
    float determinant=0;

    for(i=0;i<3;i++)
    {
        determinant = determinant + a[0][i]*(a[1][(i+1)%3]*a[2][(i+2)%3] - a[1][(i+2)%3]*a[2][(i+1)%3]);
    }

    if(determinant==0)
    {
        cout<<"Inverse does not exist (Determinant=0).\n";
        return -1;
    }
    
    for(i=0;i<3;i++)
    {
        for(j=0;j<3;j++)
        {
            out[j][i]=((a[(i+1)%3][(j+1)%3] * a[(i+2)%3][(j+2)%3]) - (a[(i+1)%3][(j+2)%3]*a[(i+2)%3][(j+1)%3]))/ determinant;            
        }      
    }
    
    return 1;
}

//end of matrix code



class Match
{   
	private:
        string fileName;
        int matchCount;		
             
    public:        
        int getMatchCount()
        {
            return matchCount;
        }

        string getFileName()
        {
            return fileName;
        }
        
        Match(string fileName, int matchCount)
        {
            this->fileName = fileName; 
            this->matchCount = matchCount;
        }
        
		//override the comparator operator
	    bool operator < (const Match& e1) const
    	{
	        return matchCount > e1.matchCount;
    	}
};

//#define MIN_SIFT_DISTANCE 10


#define MIN_SIFT_DISTANCE 100
#define displacementRatio 0.8

//part2
//This method maps directly from source coordinates to target coordinates
CImg<double> applyTransform(CImg<double> inp, float m[3][3])
{
    int width=inp.width();
    int height = inp.height();
    CImg<double> out(width, height, 1, inp.spectrum(), 0);
    
    for(int i=0; i<width; i++)
    {
        for(int j=0; j<height; j++)
        {
            //apply transformation of coordinates..
            
            float base = m[2][0]*i+m[2][1]*j+1;
            
            //how to handle non-int coordiantes??
            float newx = (m[0][0]*i+m[0][1]*j+m[0][2])/base;
            float newy = (m[1][0]*i+m[1][1]*j+m[1][2])/base;
            
            if(newx>=0 && newy>=0 && newx<width && newy<height)
            {
                //the below code is independent of number of channels,,,
                //so the input can be grayscale or color..
                for(int ch=0; ch<inp.spectrum(); ch++)
                {
                    out(newx, newy, 0, ch) = inp(i, j, 0, ch);                    
                }
            }
        }
    }

    return out;
}

//part2
//This method maps from target coordinates to source coordinates using matrix inverse
CImg<double> applyTransformInverse(CImg<double> inp, float n[3][3])
{
    float m[3][3] = {{1.1247, -0.3147, 222.94}, {0.1089, 0.685, -19.925}, {0.000265, -0.000597, 1.08278}};
    double status = matrixInverse(n, m);
    
    /*
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {
            cout<<m[i][j]<<"\t";
        }
        cout<<endl;
    }
    */
    
    if(status==-1)
    {
        //do a simple transform
        return applyTransform(inp, n);
    }
    
    int width=inp.width();
    int height = inp.height();
    CImg<double> out(width, height, 1, inp.spectrum(), 0);
    
    for(int i=0; i<width; i++)
    {
        for(int j=0; j<height; j++)
        {
            //apply transformation of coordinates..
            
            float base = m[2][0]*i+m[2][1]*j+1;
            
            //how to handle non-int coordiantes??
            float newx = (m[0][0]*i+m[0][1]*j+m[0][2])/base;
            float newy = (m[1][0]*i+m[1][1]*j+m[1][2])/base;
            
            if(newx>=0 && newy>=0 && newx<width && newy<height)
            {
                //the below code is independent of number of channels,,,
                //so the input can be grayscale or color..
                for(int ch=0; ch<inp.spectrum(); ch++)
                {
                    out(i, j, 0, ch) = inp(newx, newy, 0, ch);                    
                }
            }
        }
    }

    return out;
}




//part3.1
int matchSIFT(CImg<double> inp1, CImg<double> inp2, bool displayMatch)
{
    int match=0;
    
    CImg<double> gray1;
    // convert image to grayscale
    if(inp1.spectrum()==1)
    {
    	gray1=inp1;
    }
    else
    {
    	gray1 = inp1.get_RGBtoHSI().get_channel(2);
    }
    
    CImg<double> gray2;
    if(inp2.spectrum()==1)
    {
    	gray2=inp2;
    }
    else
    {
	    gray2 = inp2.get_RGBtoHSI().get_channel(2);
    }    
    
    vector<SiftDescriptor> d1 = Sift::compute_sift(gray1);
    vector<SiftDescriptor> d2 = Sift::compute_sift(gray2);
    
    
    //create a merged image..
    int maxHeight = inp1.height();
    if(inp2.height() > maxHeight)
    {
    	maxHeight = inp2.height();
    }
    
    //create a big image here..
    CImg<double> merge(inp1.width()+inp2.width(), maxHeight, 1, inp1.spectrum(), 0);
    if(displayMatch)
    {
		//place the first image in the left half..
		for(int i=0; i<inp1.width(); i++)
		{
			for(int j=0; j<inp1.height(); j++)
			{
		        for(int ch=0; ch<inp1.spectrum(); ch++)
				{
				    merge(i, j, 0, ch) = inp1(i, j, 0, ch);                    
				}
			}
		}

		//place the second image in the right half..
		for(int i=inp1.width(); i<merge.width(); i++)
		{
			for(int j=0; j<inp2.height(); j++)
			{
		        for(int ch=0; ch<inp2.spectrum(); ch++)
				{
				    merge(i, j, 0, ch) = inp2(i-inp1.width(), j, 0, ch);                    
				}
			}
		}
    }
    
    const unsigned char color[] = {0, 255, 0};
    for(int i=0; i<d1.size(); i++)
    {
        for(int j=0; j<d2.size(); j++)
        {
            //match both vectors here..
            double dist =0;
            for(int k=0; k<128; k++)
            {
                dist+= pow((d1[i].descriptor[k]-d2[j].descriptor[k]),2);
            }
            dist = sqrt(dist);
            if(dist<MIN_SIFT_DISTANCE)
            {
                //new coordinates for second image will be translated by inp1.width() amount..
   			    //merge.draw_line( d1[i].col, d1[i].row, d2[i].col+inp1.width(), d2[i].row, color, 1, 0U);	
                if(displayMatch)
                {
                	merge.draw_line( d1[i].col, d1[i].row, d2[j].col+inp1.width(), d2[j].row, color);
                }
                match++;
            }
        }
    }
    
    if(displayMatch)
    {
        merge.save("matchSift.png");   
    }
    
    return match; 
}

// part 4

CImg<double> stitch(CImg<double> input1, CImg<double> input2) {
	// convert image to grayscale
	CImg<double> gray1 = input1.get_RGBtoHSI().get_channel(2);
	vector<SiftDescriptor> d1 = Sift::compute_sift(gray1);

	// convert image to grayscale
	CImg<double> gray2 = input2.get_RGBtoHSI().get_channel(2);
	vector<SiftDescriptor> d2 = Sift::compute_sift(gray2);
	
	vector<int> inliers(d1.size(), 0), xArr(d1.size(), 0), yArr(d1.size(), 0);
	vector<bool> check(d1.size(), false);
	
	for(int i = 0; i < d1.size(); i++) {
		double min1 = MIN_SIFT_DISTANCE, min2 = MIN_SIFT_DISTANCE + 10;
		int j1 = 0, j2 = 0;
		for(int j = 0; j < d2.size(); j++) {
			double dist = 0;
			for(int k = 0; k < 128; k++) {
				dist += pow((d1[i].descriptor[k] - d2[j].descriptor[k]),2);
			}
			dist = sqrt(dist);
			if(dist < min2) {
				if(dist < min1) {
					min1 = dist;
					j1 = j;
				} else {
					min2 = dist;
					j2 = j;
				}
			}
		}
		double ratio = min1 / min2;
		if(ratio < displacementRatio) {
			xArr[i] = d2[j1].col - d1[i].col;
			yArr[i] = d2[j1].row - d1[i].row;
			check[i] = true;
		}
	}
	
	for(int i = 0; i < d1.size(); i++) {
		for(int j = 0; j < d1.size(); j++) {
			if(check[i] && check[j] && xArr[i] == xArr[j] && yArr[i] == yArr[j])
				inliers[i]++;
		}
	}	
	
	int x = 0, y = 0, maxCount = 0, x1, x2, y1, y2;
	
	for(int i = 0; i < d1.size(); i++) {
		if(inliers[i] > maxCount) {
			x = xArr[i];
			y = yArr[i];
			maxCount = inliers[i];
		}			   
	}
	
	if(x >= 0 && y >= 0) {
		x1 = x, y1 = y, x2 = 0, y2 = 0;
	} else if(x >= 0 && y < 0) {
		x1 = x, y1 = 0, x2 = 0, y2 = -y;
	} else if(x < 0 && y >= 0) {
		x1 = 0, y1 = y, x2 = -x, y2 = 0;
	} else {
		x1 = 0, y1 = 0, x2 = -x, y2 = -y;
	}
	
	CImg<double> output(input1.width() + input2.width() + abs(x), input1.height() + input2.height() + abs(y), 1, input1.spectrum(), 255);
	
	for(int i = x1; i < x1 + input1.width(); i++) {
		for(int j = y1; j < y1 + input1.height(); j++) {
			for(int k = 0; k < input1.spectrum(); k++) {					
				output(i, j, 0, k) = output(i, j, 0, k) == 255 ? input1(i - x1, j - y1, 0, k) : (output(i, j, 0, k) + input1(i - x1, j - y1, 0, k)) / 2;
			}
		}
	}
	
	for(int i = x2; i < x2 + input2.width(); i++) {
		for(int j = y2; j < y2 + input2.height(); j++) {
			for(int k = 0; k < input2.spectrum(); k++) {					
				output(i, j, 0, k) = output(i, j, 0, k) == 255 ? input2(i - x2, j - y2, 0, k) : (output(i, j, 0, k) + input2(i - x2, j - y2, 0, k)) / 2;
			}
		}
	}
	
	output.autocrop(255, "xy");
	
	return output;
}

//part 5

class HashMatch
{
	private:
		float fingerPrint;
		int index;
	
	public:		
		HashMatch(float fingerPrint, int index)
		{
			this->fingerPrint = fingerPrint;
			this->index = index;
		}
		
		bool isMatch(float fp)
		{
			return this->fingerPrint == fp;	
		}
		
		int getIndex()
		{
			return index; //this is the index of the sift in the first image's list..
		}
};

//solution for question 5.
int optMatchSIFT(CImg<double> inp1, CImg<double> inp2, bool displayMatch, bool verifyMatch)
{

	int correctMatches=0;
	
	clock_t begin = clock();

	time_t t;
	//seed the random number generator..
	srand((unsigned) time(&t));	
    int match=0;
    
    // convert image to grayscale
    CImg<double> gray1 = inp1.get_RGBtoHSI().get_channel(2);
    vector<SiftDescriptor> d1 = Sift::compute_sift(gray1);

    // convert image to grayscale
    CImg<double> gray2 = inp2.get_RGBtoHSI().get_channel(2);
    vector<SiftDescriptor> d2 = Sift::compute_sift(gray2);
    
    cout<<"Number of SIFT descriptors in Image-1: "<<d1.size()<<endl;
    cout<<"Number of SIFT descriptors in Image-2: "<<d2.size()<<endl;
    
    //pick only 10 from d1.
    //start of debug code
    bool isDebug=false;
    vector<SiftDescriptor> temp;
    
    if(isDebug)
    {
		for(int i=0; i<10; i++)
		{
			int index = rand()%d1.size();
			temp.push_back(d1[index]);    	
		}
		d2=temp;
    }
	//end of debug code..
    
    //create a merged image..
    int maxHeight = inp1.height();
    if(inp2.height() > maxHeight)
    {
    	maxHeight = inp2.height();
    }
    
    //create a big image here..
    CImg<double> merge(inp1.width()+inp2.width(), maxHeight, 1, inp1.spectrum(), 0);
    if(displayMatch)
    {
		//place the first image in the left half..
		for(int i=0; i<inp1.width(); i++)
		{
			for(int j=0; j<inp1.height(); j++)
			{
		        for(int ch=0; ch<inp1.spectrum(); ch++)
				{
				    merge(i, j, 0, ch) = inp1(i, j, 0, ch);                    
				}
			}
		}

		//place the second image in the right half..
		for(int i=inp1.width(); i<merge.width(); i++)
		{
			for(int j=0; j<inp2.height(); j++)
			{
		        for(int ch=0; ch<inp2.spectrum(); ch++)
				{
				    merge(i, j, 0, ch) = inp2(i-inp1.width(), j, 0, ch);                    
				}    
			}
		}
    }
    
    int k_size = 5; //the size of the reduced dimensional space,
    int bin_size = 100; //this is the bin size for each k.
    
    vector< vector<float> > transform;
	map< int, vector<HashMatch> > hash; //this is the final hash which we will use..
	vector<float> weights;
	vector<float> fp_weights; //weights for calculating the finger print..
	int fp_hash_size = 73; //hash table size for fingerprint calculation..
	
	int hash_size = 101; //reasonably big prime number ..?
	float offset; //should this be continous or discrete..

	//below are the params for normal distribution -  N(0,1)
	std::normal_distribution<int> distribution(0, 1);
	std::random_device rd;
	std::mt19937 e2(rd());
    
    //compute the transformation metrics..
    for(int i=0; i<k_size; i++)
    {
    	vector<float> curr;
    	for(int j=0; j<128; j++)
    	{
    		//pick an element from a normal distribution..
			curr.push_back(distribution(e2));
    	}
    	transform.push_back(curr);
    }
    
    //compute H values, which are the weights of individual transformations.
    double sum1=0, sum2=0;
    for(int i=0; i<k_size; i++)
    {
		int r = rand()%hash_size;
		sum1+=r;
    	//get weights in the range of 0-hash_size-1
    	weights.push_back(r);

		r = rand()%hash_size;
		sum2+=r;
    	fp_weights.push_back(r);
    }    
    
    //normalize the weights so that the hash value doesnt leave the int bounds..
    for(int i=0; i<k_size; i++)
    {
    	weights[i] /=sum1;
    	fp_weights[i] /=sum2; 
    }
    
    
    //pick offset...from uniform distribution of 0-bin_size
    offset = rand()%(bin_size+1);
    double test;
    const unsigned char color[] = {0, 255, 0};
    for(int i=0; i<d1.size(); i++)
    {
		//compute hash values for the first image..
		double hash_value=0;
		double fingerPrint = 0;
		for(int j=0; j<k_size; j++)
		{			
			//calculate k_i
			double k=offset; //init with offset..
 			for(int m=0; m<128; m++)
			{
				k+= (d1[i].descriptor[m]*transform[j][m]);
			}
			
			k = (long int)(k/bin_size); //discretize it to bins..						
			
			hash_value+= (k*weights[j]);
			fingerPrint+=(k*fp_weights[j]);	
		}
		
		if(hash_value<0)
		{
			hash_value*=-1;			
		}
		
		if(fingerPrint<0)
		{
			fingerPrint*=-1;			
		}
		//cout<<hash_value<<endl;
		hash_value = ((long int)hash_value)%hash_size;
		fingerPrint = ((long int)fingerPrint)%fp_hash_size;
		//cout<<hash_value<<endl;
		//insert into the map..
		if(hash.count((int)hash_value)==0)
		{
			vector<HashMatch> temp;
			HashMatch h(fingerPrint, i);
			temp.push_back(h);
			hash[(int)hash_value] = temp;
		}
		else
		{
			test=hash_value;
			HashMatch h(fingerPrint, i);
			hash[(int)hash_value].push_back(h);
		}
    }
    
	int single=0;
	int more=0;
    for(int i=0; i<hash_size; i++)
    {
    	if(hash.count(i)==1)
    	{
    		if(hash[i].size() == 1)
    		{
    			single++;
    		}
    		else
    		{
    			more++;	
    		}
    	}    	
    }
    
    cout<<"Number of Bins with Single Elements: "<<single<<endl;
    cout<<"Number of Bins with More than one Element: "<<more<<endl;
    int numMiss=0;    
    int numfpMiss=0;
    for(int i=0; i<d2.size(); i++)
    {
		//compute hash values for the first image..
		double hash_value=0;
		double fingerPrint = 0;
		for(int j=0; j<k_size; j++)
		{			
			//calculate k_i
			double k=offset; //init with offset..
 			for(int m=0; m<128; m++)
			{
				k+= (d2[i].descriptor[m]*transform[j][m]);
			}
		
			k = (long int)(k/bin_size); //discretize it to bins..						
		
			hash_value+= (k*weights[j]);
			fingerPrint+=(k*fp_weights[j]);	
		}
	
		if(hash_value<0)
		{
			hash_value*=-1;			
		}
	
		if(fingerPrint<0)
		{
			fingerPrint*=-1;			
		}
		//cout<<hash_value<<endl;
		hash_value = ((long int)hash_value)%hash_size;
		fingerPrint = ((long int)fingerPrint)%fp_hash_size;
		
		
		
		if(hash.count((int)hash_value)==0)
		{
			//miss case..where the bin has no other feature descriptors
			//solution: loop through all the features to get the match
			float minDistance=INT_MAX;
			int index=-1;
			for(int j=0; j<d1.size(); j++)
			{
				float dist=0;
				
				for(int k=0; k<128; k++)
				{
					dist += pow(d1[j].descriptor[k]-d2[i].descriptor[k],2);
				}			
				
				if(dist<minDistance)
				{
					minDistance = dist;
					index=j;
				}				
			}
			
			if (verifyMatch && sqrt(minDistance) < MIN_SIFT_DISTANCE)
			{
				correctMatches++;
			}
			
			//draw a line for the match..
	    	merge.draw_line( d1[index].col, d1[index].row, d2[i].col+inp1.width(), d2[i].row, color);
			match++;
			numMiss++;
		}
		else
		{
			//check for finger print..
			bool isMatchFound = false;
			for(int l=0; l<hash[(int)hash_value].size(); l++)
			{
				if(hash[(int)hash_value][l].isMatch(fingerPrint))
				{	
					int index = hash[(int)hash_value][l].getIndex();
					isMatchFound=true;				
					//match found..
		            if(displayMatch)
				    {
				    	merge.draw_line( d1[index].col, d1[index].row, d2[i].col+inp1.width(), d2[i].row, color);
				    }
				    
				    if(verifyMatch)
				    {
						//compute the distance..
						float minDistance=0;
						
						for(int p=0; p<128;p++)
						{
							minDistance += pow(d1[index].descriptor[p]-d2[i].descriptor[p],2);
						}					
							
						if(sqrt(minDistance) < MIN_SIFT_DISTANCE)
						{
							correctMatches++;
						}
				    }
				    
				    
				    match++;
					break;
				}
			}
			if(!isMatchFound)
			{
				//finger print miss..				
				
				float minDistance=INT_MAX;
				int index_match=-1;
				
				for(int l=0; l<hash[(int)hash_value].size(); l++)
				{
					int index = hash[(int)hash_value][l].getIndex();	
					
					float dist=0;
				
					for(int k=0; k<128; k++)
					{
						dist += pow(d1[index].descriptor[k]-d2[i].descriptor[k],2);
					}			
				
					if(dist<minDistance)
					{
						minDistance = dist;
						index_match=index;
					}						
				}				

				if (verifyMatch && sqrt(minDistance) < MIN_SIFT_DISTANCE)
				{
					correctMatches++;
				}

				//draw a line for the match..
				merge.draw_line( d1[index_match].col, d1[index_match].row, d2[i].col+inp1.width(), d2[i].row, color);
				match++;
				numfpMiss++;
			}
		}
		
    }
    
    cout<<"Number of Empty Bin Hash Misses: "<<numMiss<<endl;
    cout<<"Number of Finger Print Misses: "<<numfpMiss<<endl;
    
    if(verifyMatch)
    {
    	cout<<"Number of Accurate matches is"<<correctMatches<<endl;
    	cout<<"Percentage of correct matches is "<<(correctMatches*100/float(match))<<endl;
    }
    
    if(displayMatch)
    {
        merge.save("matchSift.png");  
    }
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout<<"Total Time taken (in seconds) for Matching is "<<elapsed_secs<<endl;
    return match; 
}



int main(int argc, char **argv)
{
    try 
    {

        if(argc < 2)
        {
            cout << "Insufficent number of arguments; correct usage:" << endl;
            cout << "    p4 part_id ..." << endl;
            return -1;
        }

        string part = argv[1];
        string inputFile = argv[2];

        if(part == "part2")
        {
            cout<<"Image Warping!!"<<endl<<endl;;
            if(argc < 4)
            {
                cout << "Insufficent number of arguments; correct usage:" << endl;
                cout << "    p4 part2 <input-File-Name> <output-File-Name>"  << endl;
                return -1;
            }
        
            //define the transformation matrix
            float m[3][3] = {{0.907,0.258,-182},{-0.153,1.44,58},{-0.000306,0.000731,1}};
            CImg<double> inp(inputFile.c_str());
            CImg<double> out = applyTransformInverse(inp, m);
            //CImg<double> out = applyTransform(inp, m);
            out.save(argv[3]);
        }
        else if(part == "part3.1")
        {
            if(argc < 4)
            {
                cout << "Insufficent number of arguments; correct usage:" << endl;
                cout << "    p4 part3.1 <input-File-Name1> <input-File-Name2>"  << endl;
                return -1;
            }

            clock_t begin = clock();
            CImg<double> input1(argv[2]);
            CImg<double> input2(argv[3]);

            double match=matchSIFT(input1, input2, true);
            cout<<"Number of matched SIFT descriptors is "<<match<<endl;            
            
            clock_t end = clock();
			double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			cout<<"Total Time taken (in seconds) for Matching is "<<elapsed_secs<<endl;
            
        }   
        else if(part == "part3.2")
        {
            if(argc <4)
            {                
                cout << "Insufficent number of arguments; correct usage:" << endl;
                cout << "    p4 part3.2 <Query-File-Name1> <input-File-Name(1)> ... <input-File-Name(N)>"  << endl;
                return -1;                
            }
            CImg<double> query(argv[2]);
            
            cout<<endl<<endl<<"Query Image is"<<argv[2]<<endl;
            
            vector<Match> matches; 
            for(int i=3; i<argc; i++)
            {
                CImg<double> temp(argv[i]);
                Match m(argv[i], matchSIFT(query, temp, false));
                matches.push_back(m);
            } 
                    
            //sort matches   
            std::sort(matches.begin(), matches.end());
            
            cout<<"Printing the top 10 Matches in the decreasing order of Match"<<endl;
            for(int i=0; i<matches.size() && i<10; i++)
            {
                cout<<matches[i].getFileName()<<": "<<matches[i].getMatchCount()<<endl;
            }
        }
        else if(part == "part4")
        {
            if(argc < 5)
            {                
                cout << "Insufficent number of arguments; correct usage:" << endl;
                cout << "    p4 part4 outputFile inputImg1 inputImg2"  << endl;
                return -1;                
            }
            
            CImg<double> input1(argv[3]);
            CImg<double> input2(argv[4]);
            
            CImg<double> output = stitch(input1, input2);
            output.save(argv[2]);
        }
        else if(part == "part5")
        {
            if(argc <4)
            {                
                cout << "Insufficent number of arguments; correct usage:" << endl;
                cout << "    p4 part3.1 <input-File-Name1> <input-File-Name2>"  << endl;
                return -1;           
            }
            
            CImg<double> input1(argv[2]);
            CImg<double> input2(argv[3]);

            double match=optMatchSIFT(input1, input2, true, false);
            cout<<"Number of matched SIFT descriptors is "<<match<<endl;
        }
        else if(part == "part5b")
        {
            if(argc <4)
            {                
                cout << "Insufficent number of arguments; correct usage:" << endl;
                cout << "    p4 part3.1 <input-File-Name1> <input-File-Name2>"  << endl;
                return -1;           
            }
            
            CImg<double> input1(argv[2]);
            CImg<double> input2(argv[3]);

            double match=optMatchSIFT(input1, input2, false, true);
            cout<<"Number of matched SIFT descriptors is "<<match<<endl;
        }
    } 
    catch(const string &err) 
    {
        cerr << "Error: " << err << endl;
    }
    return 0;
}








