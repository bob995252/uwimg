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
    /*
    image new = make_image(w, w, 3);
    for (int k = 0; k<3; ++k){
    for (j = 0; j < w; ++j){
        for (i = 0; i < w; ++i){
            set_pixel(new, i, j, k, 1./(w * w));
        }
    }}
    return new;
    */
    return box_filter;
    
}

float do_conv_cal(image im, image filter, int x, int y, int z){
    int i, j, filter_c;
    float conv_sum=0;
    filter_c = (im.c > filter.c) ? 0 : z;
    for (j = 0; j < filter.h; ++j){
        for (i = 0; i < filter.w; ++i){
            conv_sum = conv_sum + get_pixel(im, x+i-1, y+j-1, z) * get_pixel(filter, i, j, filter_c);
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
        for(k = 0; k < im.c; ++k){
            for (j = 0; j < im.h; ++j){
                for (i = 0; i < im.w; ++i){
                    if (k == 0 && j ==0) printf("pixel = %f\n", get_pixel(conv_im, i, j, k));
                }}}
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
    return make_image(1,1,1);
}

image make_sharpen_filter()
{
    // TODO
    return make_image(1,1,1);
}

image make_emboss_filter()
{
    // TODO
    return make_image(1,1,1);
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: TODO

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO

image make_gaussian_filter(float sigma)
{
    // TODO
    return make_image(1,1,1);
}

image add_image(image a, image b)
{
    // TODO
    return make_image(1,1,1);
}

image sub_image(image a, image b)
{
    // TODO
    return make_image(1,1,1);
}

image make_gx_filter()
{
    // TODO
    return make_image(1,1,1);
}

image make_gy_filter()
{
    // TODO
    return make_image(1,1,1);
}

void feature_normalize(image im)
{
    // TODO
}

image *sobel_image(image im)
{
    // TODO
    return calloc(2, sizeof(image));
}

image colorize_sobel(image im)
{
    // TODO
    return make_image(1,1,1);
}
