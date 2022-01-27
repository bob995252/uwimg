#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>
#define PI 3.1415926
// Frees an array of descriptors.
// descriptor *d: the array.
// int n: number of elements in array.
void free_descriptors(descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        free(d[i].data);
    }
    free(d);
}

// Create a feature descriptor for an index in an image.
// image im: source image.
// int i: index in image for the pixel we want to describe.
// returns: descriptor for that index.
descriptor describe_index(image im, int i)
{
    int w = 5;
    descriptor d;
    d.p.x = i%im.w;
    d.p.y = i/im.w;
    d.data = calloc(w*w*im.c, sizeof(float));
    d.n = w*w*im.c;
    int c, dx, dy;
    int count = 0;
    // If you want you can experiment with other descriptors
    // This subtracts the central value from neighbors
    // to compensate some for exposure/lighting changes.
    for(c = 0; c < im.c; ++c){
        float cval = im.data[c*im.w*im.h + i];
        for(dx = -w/2; dx < (w+1)/2; ++dx){
            for(dy = -w/2; dy < (w+1)/2; ++dy){
                float val = get_pixel(im, i%im.w+dx, i/im.w+dy, c);
                d.data[count++] = cval - val;
            }
        }
    }
    return d;
}

// Marks the spot of a point in an image.
// image im: image to mark.
// ponit p: spot to mark in the image.
void mark_spot(image im, point p)
{
    int x = p.x;
    int y = p.y;
    int i;
    for(i = -9; i < 10; ++i){
        set_pixel(im, x+i, y, 0, 1);
        set_pixel(im, x, y+i, 0, 1);
        set_pixel(im, x+i, y, 1, 0);
        set_pixel(im, x, y+i, 1, 0);
        set_pixel(im, x+i, y, 2, 1);
        set_pixel(im, x, y+i, 2, 1);
    }
}

// Marks corners denoted by an array of descriptors.
// image im: image to mark.
// descriptor *d: corners in the image.
// int n: number of descriptors to mark.
void mark_corners(image im, descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        mark_spot(im, d[i].p);
    }
}

// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row image of the filter.
image make_1d_gaussian(float sigma)
{
    // TODO: optional, make separable 1d Gaussian.
    int filter_size, i;
    filter_size = ((int)floor(6*sigma) % 2 == 0) ? floor(6*sigma)+1:floor(6*sigma);
    image gaussian_1d_filter = make_image(filter_size, 1, 1);
    for (i = 0; i<filter_size; ++i) {
        set_pixel(gaussian_1d_filter, i, 0, 0, 
        1./(sqrt(2*PI) * sigma) * exp(-(pow(i-gaussian_1d_filter.w/2, 2.)) / (2 * pow(sigma, 2.))));
    }
    return gaussian_1d_filter;
}

// Smooths an image using separable Gaussian filter.
// image im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed image.
image smooth_image(image im, float sigma)
{
    if(0){
        image g = make_gaussian_filter(sigma);
        image s = convolve_image(im, g, 1);
        free_image(g);
        return s;
    } else {
        // TODO: optional, use two convolutions with 1d gaussian filter.
        // If you implement, disable the above if check.
        image g1 = make_1d_gaussian(sigma);
        image g2 = make_image(g1.h, g1.w, g1.c);
        for (int i = 0; i < g1.w; ++i) set_pixel(g2, 0, i, 0, get_pixel(g1, i, 0, 0));

        im = convolve_image(im, g1, 1);
        im = convolve_image(im, g2, 1);
        //feature_normalize(im);
        return im;
    }
}

// Calculate the structure matrix of an image.
// image im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.


image structure_matrix(image im, float sigma)
{
    int i, j;
    image S = make_image(im.w, im.h, 3);
    image S1 = make_image(im.w, im.h, 3);
    image S2 = make_image(im.w, im.h, 3);
    // TODO: calculate structure matrix for im.
    image Ix = make_gx_filter();
    image Iy = make_gy_filter();
    image gauss_filter = make_gaussian_filter(sigma);
    S1 = convolve_image(im, Ix, 0);
    S2 = convolve_image(im, Iy, 0);
    for (j = 0; j < im.h; ++j){
        for (i = 0; i < im.w; ++i){
            set_pixel(S, i, j, 0, get_pixel(S1, i, j, 0) * get_pixel(S1, i, j, 0));
            set_pixel(S, i, j, 1, get_pixel(S2, i, j, 0) * get_pixel(S2, i, j, 0));
            set_pixel(S, i, j, 2, get_pixel(S1, i, j, 0) * get_pixel(S2, i, j, 0));
        }
    }
    //S = convolve_image(S, gauss_filter, 1);
    //feature_normalize(S);
    S = smooth_image(S, sigma);
    return S;
}

// Estimate the cornerness of each pixel given a structure matrix S.
// image S: structure matrix for an image.
// returns: a response map of cornerness calculations.
// S = [Ix^2 Ixy; Ixy Iy^2]
// det(S) = Ix^2 * Iy^2 - Ixy * Ixy
// trace(S) = Ix^2 + Iy^2
image cornerness_response(image S)

{
    image R = make_image(S.w, S.h, 1);
    int i, j;
    float detS, traceS, alpha = 0.06;
    // TODO: fill in R, "cornerness" for each pixel using the structure matrix.
    // We'll use formulation det(S) - alpha * trace(S)^2, alpha = .06.
    for (j = 0; j < S.h; ++j){
        for (i = 0; i < S.w; ++i){
            detS = get_pixel(S, i, j, 0) * get_pixel(S, i, j, 1) - get_pixel(S, i, j, 2) * get_pixel(S, i, j, 2);
            traceS = get_pixel(S, i, j, 0) + get_pixel(S, i, j, 1);
            set_pixel(R, i, j, 0, detS - alpha * traceS * traceS);
        }
    }
    return R;
}
// Perform non-max supression on an image of feature responses.
// image im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: image with only local-maxima responses within w pixels.


void connectStrong(image r, int StrongRow, int StrongCol, int **flag)
{
    int x, y;
    for (int j = -1; j<2; ++j){
        for (int i = -1; i<2; ++i){
            x = (StrongRow + i < 0) ? 0: ( (StrongRow + i >= r.w) ? r.w-1: StrongRow + i);
            y = (StrongCol + j < 0) ? 0: ( (StrongCol + j >= r.h) ? r.h-1: StrongCol + j);
            //printf("x = %d, y = %d \n", x, y);
            if (get_pixel(r, StrongRow + i, StrongCol + j, 0) > 0 &&
                get_pixel(r, StrongRow + i, StrongCol + j, 0) < 1 && flag[x][y] == 0) {
                set_pixel(r, StrongRow + i, StrongCol + j, 0, 1);
                flag[x][y] = 1;
                connectStrong(r, StrongRow + i, StrongCol + j, flag);
            }
        }
    }
}

image nms_image(image im, int w)
{
    //image r = copy_image(im);
    
    image r = make_image(im.w, im.h, 1);
    //printf("r.w = %d, r.h = %d, r.c = %d\n", r.w, r.h, r.c);
    image *res, mag, theta;
    int i, j, neigh, yr[2], yl[2], xr[2], xl[2], strong_ind=0, weak_ind=0;
    int *StrongCol, *StrongRow, *WeakCol, *WeakRow, **flag;
    float theta_radian, right_pixel, left_pixel, target, HighThre = 0.2, LowThre = 0.05;
    res = sobel_image(im);
    mag = res[0];
    theta = res[1];
    feature_normalize(mag);
    //feature_normalize2(theta);

    // TODO: perform NMS on the response map.
    for (j = 0; j < im.h; ++j){
        for (i = 0; i < im.w; ++i){
            theta_radian = get_pixel(theta, i, j, 0);
            if ((theta_radian >=0 && theta_radian <= PI/4.) || (theta_radian <= -3./4.*PI && theta_radian > -PI)){
                //for (neigh = 0; neigh < w; ++neigh){
                xr[0] = i+1;
                yr[0] = j;
                xr[1] = i+1;
                yr[1] = j-1;
                xl[0] = i-1;
                yl[0] = j;
                xl[1] = i-1;
                yl[1] = j+1;
                target = get_pixel(mag, i, j, 0);
                right_pixel = get_pixel(mag, xr[0], yr[0], 0) 
                            + tan(theta_radian) * (get_pixel(mag, xr[1], yr[1], 0) - get_pixel(mag, xr[0], yr[0], 0));
                left_pixel = get_pixel(mag, xl[0], yl[0], 0) 
                            + tan(theta_radian) * (get_pixel(mag, xl[1], yl[1], 0) - get_pixel(mag, xl[0], yl[0], 0));
                if (target >= right_pixel && target >= left_pixel) set_pixel(r, i, j, 0, target);
                else set_pixel(r, i, j, 0, -999999.);
                //}
            }
            else if ((theta_radian > PI/4. && theta_radian < PI/2.) || (theta_radian < -PI/2 && theta_radian > - 3./4.*PI)){
                //for (neigh = 0; neigh < w; ++neigh){
                xr[0] = i;
                yr[0] = j-1;
                xr[1] = i+1;
                yr[1] = j-1;
                xl[0] = i;
                yl[0] = j+1;
                xl[1] = i-1;
                yl[1] = j+1;
                target = get_pixel(mag, i, j, 0);
                right_pixel = get_pixel(mag, xr[0], yr[0], 0) 
                            + tan(PI/2 - theta_radian) * (get_pixel(mag, xr[1], yr[1], 0) - get_pixel(mag, xr[0], yr[0], 0));
                left_pixel = get_pixel(mag, xl[0], yl[0], 0) 
                            + tan(PI/2 - theta_radian) * (get_pixel(mag, xl[1], yl[1], 0) - get_pixel(mag, xl[0], yl[0], 0));
                if (target >= right_pixel && target >= left_pixel) set_pixel(r, i, j, 0, target);
                else set_pixel(r, i, j, 0, -999999.);
                //}
            }
            else if ((theta_radian > PI/2. && theta_radian <= 3./4.*PI) || (theta_radian <= -PI/4. && theta_radian > -PI/2.)){
                //for (neigh = 0; neigh < w; ++neigh){
                xr[0] = i;
                yr[0] = j+1;
                xr[1] = i+1;
                yr[1] = j+1;
                xl[0] = i;
                yl[0] = j-1;
                xl[1] = i-1;
                yl[1] = j-1;
                target = get_pixel(mag, i, j, 0);
                right_pixel = get_pixel(mag, xr[0], yr[0], 0) 
                            + fabs(tan(PI/2 + theta_radian)) * (get_pixel(mag, xr[1], yr[1], 0) - get_pixel(mag, xr[0], yr[0], 0));
                left_pixel = get_pixel(mag, xl[0], yl[0], 0) 
                            + fabs(tan(PI/2 + theta_radian)) * (get_pixel(mag, xl[1], yl[1], 0) - get_pixel(mag, xl[0], yl[0], 0));
                if (target >= right_pixel && target >= left_pixel) set_pixel(r, i, j, 0, target);
                else set_pixel(r, i, j, 0, -999999.);
                //}
            }
            else if ((theta_radian > 3./4.*PI && theta_radian <= PI) || (theta_radian < 0 && theta_radian >= -PI/4.)){
                //for (neigh = 0; neigh < w; ++neigh){
                xr[0] = i+1;
                yr[0] = j;
                xr[1] = i+1;
                yr[1] = j+1;
                xl[0] = i-1;
                yl[0] = j;
                xl[1] = i-1;
                yl[1] = j-1;
                target = get_pixel(mag, i, j, 0);
                right_pixel = get_pixel(mag, xr[0], yr[0], 0) 
                            + fabs(tan(theta_radian)) * (get_pixel(mag, xr[1], yr[1], 0) - get_pixel(mag, xr[0], yr[0], 0));
                left_pixel = get_pixel(mag, xl[0], yl[0], 0) 
                            + fabs(tan(theta_radian)) * (get_pixel(mag, xl[1], yl[1], 0) - get_pixel(mag, xl[0], yl[0], 0));
                if (target >= right_pixel && target >= left_pixel) set_pixel(r, i, j, 0, target);
                else set_pixel(r, i, j, 0, -999999.);
                //}
            }
        }
    }
    
    for (j = 0; j < im.h; ++j){
        for (i = 0; i < im.w; ++i){
            if(get_pixel(r, i, j, 0) >= HighThre){
                set_pixel(r, i, j, 0, 1.);
                strong_ind = strong_ind + 1;
            }
            else if (get_pixel(r, i, j, 0) < LowThre) set_pixel(r, i, j, 0, 0.);
            else weak_ind = weak_ind + 1;
        }
    }
    StrongCol = malloc(strong_ind * sizeof(int));
    StrongRow = malloc(strong_ind * sizeof(int));
    WeakCol = malloc(weak_ind * sizeof(int));
    WeakRow = malloc(weak_ind * sizeof(int));
    strong_ind = 0;
    weak_ind = 0;
    for (j = 0; j < im.h; ++j){
        for (i = 0; i < im.w; ++i){
            if(get_pixel(r, i, j, 0) >= HighThre){
                StrongCol[strong_ind] = j;
                StrongRow[strong_ind] = i;
                strong_ind = strong_ind + 1;
            }
            else if (get_pixel(r, i, j, 0) >= LowThre && get_pixel(r, i, j, 0) < HighThre){
                WeakCol[weak_ind] = j;
                WeakRow[weak_ind] = i;
                weak_ind = weak_ind + 1;
            }
        }
    }
    
    flag = malloc(im.w * sizeof(*flag));
    for (i = 0; i < im.w; i++) flag[i] = malloc(im.h * sizeof(flag[0]));
    for (i = 0; i< im.w; ++i)
        for (j = 0; j<im.h; ++j) 
            flag[i][j] = 0;

    for (i = 0; i< strong_ind; ++i){
        connectStrong(r, StrongRow[i], StrongCol[i], flag);
    }
    for (i = 0; i< weak_ind; ++i){
        if (get_pixel(r, WeakRow[i], WeakCol[i], 0) != 1) set_pixel(r, WeakRow[i], WeakCol[i], 0, 0.);
    }
    // for every pixel in the image:
    //     for neighbors within w:
    //         if neighbor response greater than pixel response:
    //             set response to be very low (I use -999999 [why not 0??])
    //feature_normalize(r);
    for (i = 0; i < im.w; i++) free(flag[i]);
    free(flag);
    free(StrongRow);
    free(StrongCol);
    free(WeakRow);
    free(WeakCol);
    return r;
}

// Perform harris corner detection and extract features from the corners.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// int *n: pointer to number of corners detected, should fill in.
// returns: array of descriptors of the corners in the image.
descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
{
    // Calculate structure matrix
    image S = structure_matrix(im, sigma);

    // Estimate cornerness
    image R = cornerness_response(S);

    // Run NMS on the responses
    image Rnms = nms_image(R, nms);

    int i, j;
    //TODO: count number of responses over threshold
    int count = 0; // change this
    for (j = 0; j<im.h; ++j){
        for (i = 0; i<im.w; ++i){
            printf ("nms = %f", get_pixel(Rnms, i, j, 0));
            if (get_pixel(Rnms, i, j, 0) > thresh) count++;
        }
    }
    int *target_pixel = calloc(count, sizeof(int)), cc = 0;
    for (j = 0; j<im.h; ++j){
        for (i = 0; i<im.w; ++i){
            if (get_pixel(Rnms, i, j, 0) > thresh){
                target_pixel[cc] = i + j*im.w;
                cc++;
            }
        }
    }
    printf("cc = %d, count = %d", cc, count);

    *n = count; // <- set *n equal to number of corners in image.
    descriptor *d = calloc(count, sizeof(descriptor));
    //TODO: fill in array *d with descriptors of corners, use describe_index.
    for (i = 0; i < count; ++i){
        describe_index(im, target_pixel[i]);
    }
    free(target_pixel);
    free_image(S);
    free_image(R);
    free_image(Rnms);
    return d;
}

// Find and draw corners on an image.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
void detect_and_draw_corners(image im, float sigma, float thresh, int nms)
{
    int n = 0;
    descriptor *d = harris_corner_detector(im, sigma, thresh, nms, &n);
    mark_corners(im, d, n);
}