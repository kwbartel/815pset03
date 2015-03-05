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
    Image output(im.width(), im.height(), im.channels());
    for (int x = 0; x < output.width(); x++) {
        for (int y = 0; y < output.height(); y++) {
            for (int z = 0; z < output.channels(); z++) {
                float accum = 0;
                for (int a = -((k-1)/2); a <= (k-1)/2; a++) {
                    for (int b = -((k-1)/2); b <= (k-1)/2; b++) {
                        accum += im.smartAccessor(x + a, y + b, z, clamp);
                    }
                }
                output(x, y, z) = accum/(k*k);
            }
        }
    }
    
    return output; // change this
}


// PS03 - 2.2.2 - reimeplement the box filter using the filter class.
// check that your results math those in the previous function "boxBlur"
Image boxBlur_filterClass(const Image &im, const int &k, bool clamp) {
    // Create the filter
    vector<float> boxData;
    float factor = 1.0/(k*k);

    for (int x = 0; x < k*k; x++)
        boxData.push_back(factor);
    
    Filter boxFilter(boxData, k, k);
    return boxFilter.Convolve(im, clamp);
}


// PS03 - 2.3.1 - uses a Sobel kernel to compute the horizontal and vertical
// components of the gradient of an image and returns the gradient magnitude.
Image gradientMagnitude(const Image &im, bool clamp){
    
    float fDataXArray[] = { -1.0, 0.0, 1.0, -2.0, 0.0, 2.0, -1.0, 0.0, 1.0 };
    vector<float> fDataX (fDataXArray, fDataXArray + sizeof(fDataXArray) / sizeof(float) );
    Filter sobelX(fDataX, 3, 3);
    
    float fDataYArray[] = { -1.0, -2.0, -1.0, 0.0, 0.0, 0.0, 1.0, 2.0, 1.0 };
    vector<float> fDataY (fDataYArray, fDataYArray + sizeof(fDataYArray) / sizeof(float) );
    Filter sobelY(fDataY, 3, 3);
    
    Image horizontalGradient = sobelX.Convolve(im);
    Image verticalGradient = sobelY.Convolve(im);
    Image output(im.width(), im.height(), im.channels());
    
    for (int x = 0; x < im.width(); x++)
        for (int y = 0; y < im.height(); y++)
            for (int z = 0; z < im.channels(); z++)
                output(x, y, z) = sqrt(pow(horizontalGradient(x, y, z), 2) + pow(verticalGradient(x, y, z), 2));
    return output;
}

// PS03 - 2.4.1 - create a vector containing the normalized values in a 1D Gaussian filter
vector<float> gauss1DFilterValues(float sigma, float truncate){
    
    // -1/(2 sigma**2)
    float factor = -1 / (2 * pow(sigma, 2));
    float length = 1 + 2 * (ceil(sigma * truncate));
    vector<float> values;
    
    // Calculate nonnormalized Gaussian values
    float accum = 0;
    for (int r = -ceil(sigma * truncate); r <= ceil(sigma * truncate); r++) {
        float v = exp(r*r * factor);
        values.push_back(v);
        accum += v;
    }
    
    // Normalize
    for (int i = 0; i < values.size(); i++)
        values[i] /= accum;
    
    return values;
}

// PS03 - 2.4.2 - blur across the rows of an image
Image gaussianBlur_horizontal(const Image &im, float sigma, float truncate, bool clamp){
    
    vector<float> gaussianValues = gauss1DFilterValues( sigma, truncate );
    Filter gaussianFilter( gaussianValues, gaussianValues.size(), 1 );
    return gaussianFilter.Convolve( im, clamp );
    
}

// PS03 - 2.4.3 - create a vector containing the normalized values in a 2D Gaussian filter
vector<float> gauss2DFilterValues(float sigma, float truncate){
    
    // push back pow(e,exponent) for 1 + 2 * ceil(sigma * truncate) times
    // -1/(2 sigma**2)
    float factor = -1 / (2 * pow(sigma, 2));
    float length = 1 + 2 * (ceil(sigma * truncate));
    vector<float> values;
    
    // Calculate nonnormalized Gaussian values
    float accum = 0;
    for (int x = -ceil(sigma * truncate); x <= ceil(sigma * truncate); x++) {
        for (int y = -ceil(sigma * truncate); y <= ceil(sigma * truncate); y++) {
            float r = sqrt(pow(x, 2) + pow(y, 2));
            float v = pow(2.718281828, r*r * factor);
            values.push_back(v);
            accum += v;
        }
    }
    
    // Normalize
    for (int i = 0; i < values.size(); i++)
        values[i] /= accum;
    
    return values;
}


// PS03 - 2.4.4 - Blur an image with a full  full 2D rotationally symmetric Gaussian kernel
Image gaussianBlur_2D(const Image &im, float sigma, float truncate, bool clamp){
    
    vector<float> gaussianValues = gauss2DFilterValues(sigma, truncate);
    float length = 1 + 2 * (ceil(sigma * truncate));
    Filter gaussianFilter2D( gaussianValues, length, length );
    return gaussianFilter2D.Convolve( im, clamp );
    
}

// PS03 - 2.4.5 - Use principles of seperabiltity to blur an image using 2 1D Gaussian Filters
Image gaussianBlur_seperable(const Image &im, float sigma, float truncate, bool clamp){
    
    vector<float> gaussValues = gauss1DFilterValues( sigma, truncate );
    Filter gaussHorizontal( gaussValues, gaussValues.size(), 1 );
    Filter gaussVertical( gaussValues, 1, gaussValues.size() );
    
    Image output = gaussHorizontal.Convolve( im, clamp );
    output = gaussVertical.Convolve( output, clamp );
    return output;
}


// PS03 - 2.5.1 - sharpen an image
Image unsharpMask(const Image &im, float sigma, float truncate, float strength, bool clamp){

    Image lowpass = gaussianBlur_seperable( im, sigma, truncate, clamp );
    Image highpass = im - lowpass;
    //Image output( im.width(), im.height(), im.channels() );
    //for (i = 0; i < strength; i++)
    return im + strength * highpass;
    
}


// PS03 - 3.0.1 -  Denoise an image using bilateral filtering
Image bilateral(const Image &im, float sigmaRange, float sigmaDomain, float truncateDomain, bool clamp){
    
    // Precompute spatial weight f (domain gaussian)
    vector<float> domainGaussianValues = gauss2DFilterValues( sigmaDomain, truncateDomain );
    float radius = ceil( sigmaDomain * truncateDomain );
    float length = 1 + 2 * radius;
    float range_weight_factor = -1.0/(2*sigmaRange * sigmaRange);
    
    
    Filter domainGaussianFilter( domainGaussianValues , length, length ) ;
    
    Image output( im.width(), im.height(), im.channels() );
    
    // For each output pixel
    for ( int x = 0; x < output.width(); x++ ) {
        for ( int y = 0; y < output.height(); y++ ) {
                
            // Set sums to 0
            float k_sum = 0;
            
            // For each neighborhood pixel
            for (int a = -radius; a <= radius; a++) {
                // Should b go from + to -?
                for (int b = -radius; b <= radius; b++) {
                    float range_weight = 0;
                    float gauss = 0;
                    
                    for (int c = 0; c < output.channels(); c++) {
                        // Dummy pixel
                        int x_prime = x + a;
                        int y_prime = y + b;
                        
                        // Compute color distance
                        float color_distance = 0;
                        for (int z = 0; z < im.channels(); z++)
                            color_distance += pow(im.smartAccessor(x, y, z, clamp) - im.smartAccessor(x_prime, y_prime, z, clamp), 2);
                        
                        // Compute range weight
                        //float range_weight = exp(range_weight_factor * color_distance);
                        range_weight = exp(-color_distance/(2*sigmaRange*sigmaRange));
                        gauss = exp(-(a*a+b*b) / (2 * sigmaDomain * sigmaDomain));
                        
                        output(x, y, c) += range_weight * gauss * im.smartAccessor(x_prime, y_prime, c, clamp);

                        //k_sum += range_weight * domainGaussianFilter(a + radius, b + radius);
                        //output(x, y, c) += range_weight * domainGaussianFilter(a + radius, b + radius) * im.smartAccessor(x_prime, y_prime, c, clamp);
                    }
                    k_sum += range_weight * gauss;
                }
                // After summing over all neighbors, set output pixel value
                /*cout << "------" << endl;
                cout << bilateral_sum << endl;
                cout << k_sum << endl;
                cout << bilateral_sum / k_sum - im(x, y, c) << endl;*/
            }
            for (int c = 0; c < output.channels(); c++) {
                output(x, y, c) = output(x, y, c) / k_sum;
            }
        }
    }
    
    return output;
    
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
    Image output(im.width(), im.height(), im.channels());
    for (int x = 0; x < output.width(); x++) {
        for (int y = 0; y < output.height(); y++) {
            for (int z = 0; z < output.channels(); z++) {
                float accum = 0;
                for (int a = ((width-1)/2); a >= -(width-1)/2; a--) {
                    for (int b = -((height-1)/2); b <= (height-1)/2; b++) {
                        int filter_x = a + (width-1)/2;
                        int filter_y = b + (height-1)/2;
                        accum += operator()(filter_x, filter_y) * im.smartAccessor(x + a, y + b, z, clamp);
                    }
                }
                output(x, y, z) = accum;
            }
        }
    }
    return output;
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
