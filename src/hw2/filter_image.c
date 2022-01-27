#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853
#define MAX_INT 2147483647
#define MIN_INT -2147483647

void l1_normalize(image im)
{
    // TODO
    int i, j, k;
    float *sum;
    
    sum = (float*)malloc(im.c * sizeof(float));

    for (k = 0; k < im.c; ++k) {
        sum[k] = 0;
    }
    for(k = 0; k < im.c; ++k){
        for (j = 0; j < im.h; ++j){
            for (i = 0; i < im.w; ++i){
                sum[k] = sum[k] + get_pixel(im, i, j, k);
            }
        }
        for (j = 0; j < im.h; ++j){
            for (i = 0; i < im.w; ++i){
                set_pixel(im, i, j, k, get_pixel(im, i, j, k)/sum[k]);
            }
        }
    }
    /*
    float *c_max, *c_min;
    c_max = (float*)malloc(im.c * sizeof(float));
    c_min = (float*)malloc(im.c * sizeof(float));
    for (k = 0; k < im.c; ++k) {
        c_max[k] = MIN_INT;
        c_min[k] = MAX_INT;
    }

    for(k = 0; k < im.c; ++k){
        for(j = 0; j < im.h; ++j){
            for(i = 0; i < im.w; ++i){
                if (c_max[k] < get_pixel(im, i, j, k)) c_max[k] = get_pixel(im, i, j, k);
                if (c_min[k] > get_pixel(im, i, j, k)) c_min[k] = get_pixel(im, i, j, k);
            }
        }
    }
    */
    free(sum);
}

image make_box_filter(int w)
{
    // TODO
    int i, j;
    image box_filter = make_image(w, w, 1);
    for (j = 0; j < w; ++j){
        for (i = 0; i < w; ++i){
            set_pixel(box_filter, i, j, 0, 1./(w * w));
        }
    }
    return box_filter;
    
}

float do_conv_cal(image im, image filter, int x, int y, int z){
    int i, j, filter_c, count=0;
    float conv_sum=0;
    filter_c = (im.c > filter.c) ? 0 : z;
    for (j = 0; j < filter.h; ++j){
        for (i = 0; i < filter.w; ++i){
            conv_sum = conv_sum + get_pixel(im, x+i-filter.w/2, y+j-filter.h/2, z) * get_pixel(filter, i, j, filter_c);
            count = count + 1;
        }
    }
    return conv_sum;
}

image convolve_image(image im, image filter, int preserve)
{
    // TODO
    int i, j, k;
    float conv_sum, sum_pixel=0;
    assert (filter.c == 1 || filter.c == im.c);
    if (im.c == filter.c && !preserve){
        image pre_conv_im = make_image(im.w, im.h, im.c);
        image conv_im = make_image(im.w, im.h, 1);
        
        for (j = 0; j < im.h; ++j){
            for (i = 0; i < im.w; ++i){
                for(k = 0; k < im.c; ++k){
                    conv_sum = do_conv_cal(im, filter, i, j, k);
                    set_pixel(pre_conv_im, i, j, k, conv_sum);
                }
                for(k = 0; k < im.c; ++k) sum_pixel = sum_pixel + get_pixel(pre_conv_im, i, j, k);
                set_pixel(conv_im, i, j, 0, sum_pixel);
                sum_pixel = 0;
            }
        }
        return conv_im;
    }
    else if (im.c == filter.c && preserve){
        image conv_im = make_image(im.w, im.h, im.c);
        for(k = 0; k < im.c; ++k){
            for (j = 0; j < im.h; ++j){
                for (i = 0; i < im.w; ++i){
                    conv_sum = do_conv_cal(im, filter, i, j, k);
                    set_pixel(conv_im, i, j, k, conv_sum);
                }
            }
        }
        return conv_im;
    }
    else if (filter.c == 1 && preserve){
        image conv_im = make_image(im.w, im.h, im.c);
        for(k = 0; k < im.c; ++k){
            for (j = 0; j < im.h; ++j){
                for (i = 0; i < im.w; ++i){
                    conv_sum = do_conv_cal(im, filter, i, j, k);
                    set_pixel(conv_im, i, j, k, conv_sum);
                }
            }
        }
        return conv_im;
    }
    else {
        image pre_conv_im = make_image(im.w, im.h, im.c);
        image conv_im = make_image(im.w, im.h, 1);

        for (j = 0; j < im.h; ++j){
            for (i = 0; i < im.w; ++i){
                for(k = 0; k < im.c; ++k){
                    conv_sum = do_conv_cal(im, filter, i, j, k);
                    set_pixel(pre_conv_im, i, j, k, conv_sum);
                }
                for(k = 0; k < im.c; ++k) sum_pixel = sum_pixel + get_pixel(pre_conv_im, i, j, k);
                set_pixel(conv_im, i, j, 0, sum_pixel);
                sum_pixel = 0;
            }
        }
        return conv_im;
    }
}

image make_highpass_filter()
{
    // TODO
    image highpass_filter = make_image(3, 3, 1);
    int filter_element[9];
    for (int i = 0; i<9; ++i) filter_element[i] = 0;
    filter_element[1] = -1;
    filter_element[3] = -1;
    filter_element[4] = 4;
    filter_element[5] = -1;
    filter_element[7] = -1;
    for (int j = 0; j<3; ++j)
        for (int i = 0; i<3; ++i) set_pixel(highpass_filter, i, j, 0, filter_element[i+j*3]);
    return highpass_filter;
}

image make_sharpen_filter()
{
    // TODO
    image sharpen_filter = make_image(3, 3, 1);
    int filter_element[9];
    for (int i = 0; i<9; ++i) filter_element[i] = 0;
    filter_element[1] = -1;
    filter_element[3] = -1;
    filter_element[4] = 5;
    filter_element[5] = -1;
    filter_element[7] = -1;
    for (int j = 0; j<3; ++j)
        for (int i = 0; i<3; ++i) set_pixel(sharpen_filter, i, j, 0, filter_element[i+j*3]);
    return sharpen_filter;
}


image make_emboss_filter()
{
    // TODO
    image emboss_filter = make_image(3, 3, 1);
    int filter_element[9];
    for (int i = 0; i<9; ++i) filter_element[i] = 1;
    filter_element[0] = -2;
    filter_element[1] = -1;
    filter_element[2] = 0;
    filter_element[3] = -1;
    filter_element[6] = 0;
    filter_element[8] = 2;
    for (int j = 0; j<3; ++j)
        for (int i = 0; i<3; ++i) set_pixel(emboss_filter, i, j, 0, filter_element[i+j*3]);
    return emboss_filter;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: We should preserve the channels for emboss and sharpen filters because the purpose of these two filters is 
//         changing the style of pics. As a result, the color of pics should look like one before change. On the other hand, 
//         there is no need to preserve channels for the highpass filter, which is used to detect the edge of objects where 
//         possesses high frequencies. Accordingly, we only care about the skeleton of objects instead of the color of pics.


// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: Each of them should do the clamping after we convolute the pics. This is because when we employ the filter to convolute
//         the pics, it is possible the pixels become larger or less than the bounds ([0, 1] or [0, 255]). As a result, each of 
//         filters should do post-processing. 

image make_gaussian_filter(float sigma)
{
    // TODO
    int filter_size, i, j;
    filter_size = ((int)floor(6*sigma) % 2 == 0) ? floor(6*sigma)+1:floor(6*sigma);
    //filter_size = 7;
    image gaussian_filter = make_image(filter_size, filter_size, 1);
    for (j = 0; j < gaussian_filter.h; ++j){
        for (i = 0; i < gaussian_filter.w; ++i){
            set_pixel(gaussian_filter, i, j, 0, 
            1./(TWOPI * pow(sigma, 2.)) * exp(-(pow(i-gaussian_filter.w/2, 2.) + pow(j-gaussian_filter.h/2, 2.))/ (2 * pow(sigma, 2.))));
        }
    }

    l1_normalize(gaussian_filter);

    return gaussian_filter;
}

image add_image(image a, image b)
{
    // TODO
    assert(a.w == b.w && a.h == b.h && a.c == b.c);
    image add_im = make_image(a.w, a.h, a.c);
    int i, j, k;
    for(k = 0; k < a.c; ++k){
        for (j = 0; j < a.h; ++j){
            for (i = 0; i < a.w; ++i){
                set_pixel(add_im, i, j, k, get_pixel(a, i, j, k) + get_pixel(b, i, j, k));
            }
        }
    }
    return add_im;
}

image sub_image(image a, image b)
{
    // TODO
    assert(a.w == b.w && a.h == b.h && a.c == b.c);
    int i, j, k;
    image sub_im = make_image(a.w, a.h, a.c);
    for(k = 0; k < a.c; ++k){
        for (j = 0; j < a.h; ++j){
            for (i = 0; i < a.w; ++i){
                set_pixel(sub_im, i, j, k, get_pixel(a, i, j, k) - get_pixel(b, i, j, k));
            }
        }
    }
    return sub_im;
}

image make_gx_filter()
{
    // TODO
    image gx_filter = make_image(3, 3, 1);
    int filter_element[9];
    for (int i = 0; i<9; ++i) filter_element[i] = 0;
    filter_element[0] = -1;
    filter_element[2] = 1;
    filter_element[3] = -2;
    filter_element[5] = 2;
    filter_element[6] = -1;
    filter_element[8] = 1;
    for (int j = 0; j<3; ++j)
        for (int i = 0; i<3; ++i) set_pixel(gx_filter, i, j, 0, filter_element[i+j*3]);
    return gx_filter;
}

image make_gy_filter()
{
    // TODO
    image gy_filter = make_image(3, 3, 1);
    int filter_element[9];
    for (int i = 0; i<9; ++i) filter_element[i] = 0;
    filter_element[0] = -1;
    filter_element[1] = -2;
    filter_element[2] = -1;
    filter_element[6] = 1;
    filter_element[7] = 2;
    filter_element[8] = 1;
    for (int j = 0; j<3; ++j)
        for (int i = 0; i<3; ++i) set_pixel(gy_filter, i, j, 0, filter_element[i+j*3]);
    return gy_filter;
}

void feature_normalize(image im)
{
    // TODO
    int i, j, k;
    float c_max=im.data[0], c_min=im.data[0], range;

    for(k = 0; k < im.c; ++k){
        for(j = 0; j < im.h; ++j){
            for(i = 0; i < im.w; ++i){
                if (c_max < get_pixel(im, i, j, k)) c_max = get_pixel(im, i, j, k);
                if (c_min > get_pixel(im, i, j, k)) c_min = get_pixel(im, i, j, k);
            }
        }
    }
    range = c_max - c_min;
    for(k = 0; k < im.c; ++k){
        if (range == 0){
            for (j = 0; j < im.h; ++j){
                for (i = 0; i < im.w; ++i){
                    set_pixel(im, i, j, k, 0);
                }
            }
        }
        else{
            for (j = 0; j < im.h; ++j){
                for (i = 0; i < im.w; ++i){
                    set_pixel(im, i, j, k, ((get_pixel(im, i, j, k) - c_min) / range));
                }
            }
        }
    }
    /*
    float *c_max, *c_min, *range;
    c_max = (float*)malloc(im.c * sizeof(float));
    c_min = (float*)malloc(im.c * sizeof(float));
    range = (float*)malloc(im.c * sizeof(float));
    for (k = 0; k < im.c; ++k) {
        c_max[k] = MIN_INT;
        c_min[k] = MAX_INT;
        range[k] = 0;
    }
    for(k = 0; k < im.c; ++k){
        for(j = 0; j < im.h; ++j){
            for(i = 0; i < im.w; ++i){
                if (c_max[k] < get_pixel(im, i, j, k)) c_max[k] = get_pixel(im, i, j, k);
                if (c_min[k] > get_pixel(im, i, j, k)) c_min[k] = get_pixel(im, i, j, k);
            }
        }
    }
    for (k = 0; k < im.c; ++k) range[k] = c_max[k] - c_min[k];
    for(k = 0; k < im.c; ++k){
        if (range[k] == 0){
            for (j = 0; j < im.h; ++j){
                for (i = 0; i < im.w; ++i){
                    set_pixel(im, i, j, k, 0);
                }
            }
        }
        else{
            for (j = 0; j < im.h; ++j){
                for (i = 0; i < im.w; ++i){
                    set_pixel(im, i, j, k, (get_pixel(im, i, j, k) - c_min[k]) / (c_max[k] - c_min[k]));
                }
            }
        }
    }*/
}

image *sobel_image(image im)
{
    // TODO
    int i, j;
    image gx = make_gx_filter();
    image gy = make_gy_filter();
    image gx_im = convolve_image(im, gx, 0);
    image gy_im = convolve_image(im, gy, 0);
    image g = make_image(im.w, im.h, 1);
    image theta = make_image(im.w, im.h, 1);
    image *sobel;
    sobel = calloc(2, sizeof(image));
    for (j = 0; j < im.h; ++j){
        for (i = 0; i < im.w; ++i){
            set_pixel(g, i, j, 0, sqrt(pow(get_pixel(gx_im, i, j, 0), 2) + pow(get_pixel(gy_im, i, j, 0), 2)));
            set_pixel(theta, i, j, 0, atan2(get_pixel(gy_im, i, j, 0), get_pixel(gx_im, i, j, 0)));
        }
    }
    sobel[0] = g;
    sobel[1] = theta;
    return sobel;
}

image make_gaussian_filter_addsize(float sigma, int filter_size)
{
    // TODO
    int i, j;
    assert(filter_size % 2 != 0);
    image gaussian_filter = make_image(filter_size, filter_size, 1);
    for (j = 0; j < gaussian_filter.h; ++j){
        for (i = 0; i < gaussian_filter.w; ++i){
            set_pixel(gaussian_filter, i, j, 0, 
            1./(TWOPI * pow(sigma, 2.)) * exp(-(pow(i-gaussian_filter.w/2, 2.) + pow(j-gaussian_filter.h/2, 2.))/ (2 * pow(sigma, 2.))));
        }
    }
    l1_normalize(gaussian_filter);
    return gaussian_filter;
}

image colorize_sobel(image im)
{
    // TODO
    int i, j, k;
    image color_im = make_image(im.w, im.h, im.c);
    image mag, theta, g1, g2, f;
    image *res = sobel_image(im);
    mag = res[0];
    theta = res[1];
    feature_normalize(mag);
    feature_normalize(theta);
    //-(im.h/im.w) * x + im.h = y 
    for (k = 0; k < im.c; ++k){
        for (j = 0; j < im.h; ++j){
            for (i = 0; i < im.w; ++i){    
                if (-(float)im.h/(float)im.w * i + im.h > j){
                if (k==0) set_pixel(color_im, i, j, k, get_pixel(theta, i, j, 0));
                else set_pixel(color_im, i, j, k, get_pixel(mag, i, j, 0));}
            }
        }
    }
    //scale_image(color_im, 0, 2);
    //scale_image(color_im, 1, 2);
    //clamp_image(color_im);
    hsv_to_rgb(color_im);
    // make sure the size of gaussian filter is the same
    
    g1 = make_gaussian_filter_addsize(3, 7);
    g2 = make_gaussian_filter_addsize(2, 7);
    f = g1;
    for (j = 0; j < g1.h; ++j)
        for (i = 0; i < g1.w; ++i)
            set_pixel(f, i, j, 0, get_pixel(g1, i, j, 0) - get_pixel(g2, i, j, 0));
    color_im = convolve_image(color_im, f, 1);

    //clamp_image(color_im);
    
    for (k = 0; k < im.c; ++k){
        for (j = 0; j < im.h; ++j){
            for (i = 0; i < im.w; ++i){    
                if (-(float)im.h/(float)im.w * i + im.h < j) 
                set_pixel(color_im, i, j, k, get_pixel(im, i, j, k));
            }
        }
    }
    
    return color_im;
}
