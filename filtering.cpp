// filtering.cpp
// Assignment 3


#include "filtering.h"
#include <math.h>

using namespace std;

/**************************************************************
 //            IMAGE CONVOLUTION AND FILTERING               //
 *************************************************************/


// PS03 - 2.1.1 - convolve an image with a box filter of size k by k
Image boxBlur(const Image &im, const int &k, bool clamp) {
    
    return Image(0); // change this
}


// PS03 - 2.2.2 - reimeplement the box filter using the filter class.
// check that your results math those in the previous function "boxBlur"
Image boxBlur_filterClass(const Image &im, const int &k, bool clamp) {
  
    return Image(0); // change this
}


// PS03 - 2.3.1 - uses a Sobel kernel to compute the horizontal and vertical
// components of the gradient of an image and returns the gradient magnitude.
Image gradientMagnitude(const Image &im, bool clamp){
    
    return Image(0); // change this
    
}

// PS03 - 2.4.1 - create a vector containing the normalized values in a 1D Gaussian filter
vector<float> gauss1DFilterValues(float sigma, float truncate){
    
    return vector<float>(); // change this
}

// PS03 - 2.4.2 - blur across the rows of an image
Image gaussianBlur_horizontal(const Image &im, float sigma, float truncate, bool clamp){
    
    return Image(0); // change this
}

// PS03 - 2.4.3 - create a vector containing the normalized values in a 2D Gaussian filter
vector<float> gauss2DFilterValues(float sigma, float truncate){
    
    return vector<float>(); // change this
}


// PS03 - 2.4.4 - Blur an image with a full  full 2D rotationally symmetric Gaussian kernel
Image gaussianBlur_2D(const Image &im, float sigma, float truncate, bool clamp){
    
    return Image(0); // change this
}

// PS03 - 2.4.5 - Use principles of seperabiltity to blur an image using 2 1D Gaussian Filters
Image gaussianBlur_seperable(const Image &im, float sigma, float truncate, bool clamp){
    
    return Image(0); // change this
}


// PS03 - 2.5.1 - sharpen an image
Image unsharpMask(const Image &im, float sigma, float truncate, float strength, bool clamp){

    return Image(0); // change this
    
}


// PS03 - 3.0.1 -  Denoise an image using bilateral filtering
Image bilateral(const Image &im, float sigmaRange, float sigmaDomain, float truncateDomain, bool clamp){
    
    return Image(0); // change this
}


// PS03 - 3.1.1 - 6.865 only: Bilaterial Filter an image seperatly for
// the Y and UV components of an image
Image bilaYUV(const Image &im, float sigmaRange, float sigmaY, float sigmaUV, float truncateDomain, bool clamp){

    return Image(0); // change this
}


/**************************************************************
 //                 FILTER CLASS FUNCTIONS                  //
 *************************************************************/


// PS03 - 2.2.1 - write a convolution function for the filter class
Image Filter::Convolve(const Image &im, bool clamp){
    
    return Image(0); // change this
}


/**************************************************************
 //               DON'T EDIT BELOW THIS LINE                //
 *************************************************************/

// Create an image of 0's with a value of 1 in the middle. This function
// can be used to test that you have properly set the kernel values in your
// Filter object. Make sure to set k to be larger than the size of your kernel
Image impulseImg(const int &k){
    
    // initlize a kxkx1 image of all 0's
    Image impulse(k, k, 1);
    
    // set the center pixel to have intensity 1
    int center = floor(k/2);
    impulse(center,center,0) = 1;
    
    return impulse;
}


// Filter class constructor
Filter::Filter(const vector<float> &fData, const int &fWidth, const int &fHeight) {
    
    kernel = fData;
    width = fWidth;
    height = fHeight;
    
}

Filter::Filter(const int &fWidth, const int &fHeight) {
  width = fWidth;
  height = fHeight;
  kernel = std::vector<float>(width*height,0);
}

const float & Filter::operator()(int x, int y) const {
    if (x < 0 || x >= width)
        throw OutOfBoundsException();
    if ( y < 0 || y >= height)
        throw OutOfBoundsException();
    
    return kernel[x +y*width];
    
}


float & Filter::operator()(int x, int y) {
    if (x < 0 || x >= width)
        throw OutOfBoundsException();
    if ( y < 0 || y >= height)
        throw OutOfBoundsException();
    
    return kernel[x +y*width];
}


Filter::~Filter() {}
