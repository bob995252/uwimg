#include <math.h>
#include "image.h"

float nn_interpolate(image im, float x, float y, int c)
{
    // TODO Fill in
    x = (x < 0) ? 0: ( (x >= im.w) ? im.w-1: x);
    y = (y < 0) ? 0: ( (y >= im.h) ? im.h-1: y);
    c = (c < 0) ? 0: ( (c > 2) ? 2: c);
    
    x = round(x);
    y = round(y);

    return get_pixel(im, x, y, c);
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix that first line)
    image resize_im = make_image(w, h, im.c);
    float i_scale, j_scale, m[2], b[2], del, del_m, del_b;
    int i, j, k;
    // Assume the linear transformation for changing desired size back to initial size is formulated as mx + b = y, 
    // where x -> inital and y -> desired.
    // Using Cramer's rule to solve m and b

    del = -0.5 - (im.w - 0.5);
    del_m = -0.5 - (w - 0.5);
    del_b = -0.5 * (w - 0.5) - (-0.5) * (im.w - 0.5);
    m[0] = del_m / del;
    b[0] = del_b / del;
    del = -0.5 - (im.h - 0.5);
    del_m = -0.5 - (h - 0.5);
    del_b = -0.5 * (h - 0.5) - (-0.5) * (im.h - 0.5);
    m[1] = del_m / del;
    b[1] = del_b / del;

    for(k = 0; k < im.c; ++k){
        for(j = 0; j < h; ++j){
            for(i = 0; i < w; ++i){
                i_scale = ((float)i - b[0]) / m[0];
                j_scale = ((float)j - b[1]) / m[1];
                set_pixel(resize_im, i, j, k, nn_interpolate(im, i_scale, j_scale, k));
            }
        }
    }
    return resize_im;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    float q1, q2, q;
    q1 = get_pixel(im, floor(x), floor(y), c) * (ceil(x) - x) + get_pixel(im, ceil(x), floor(y), c) * (x - floor(x));
    q2 = get_pixel(im, floor(x), ceil(y), c) * (ceil(x) - x) + get_pixel(im, ceil(x), ceil(y), c) * (x - floor(x));
    q = q1 * (ceil(y) - y) + q2 * (y - floor(y));
    return q;
}

image bilinear_resize(image im, int w, int h)
{
    // TODO
    image resize_im = make_image(w, h, im.c);
    float i_scale, j_scale, m[2], b[2], del, del_m, del_b;
    int i, j, k;
    // Assume the linear transformation for changing desired size back to initial size is formulated as mx + b = y, 
    // where x -> inital and y -> desired.
    // Using Cramer's rule to solve m and b
    
    del = -0.5 - (im.w - 0.5);
    del_m = -0.5 - (w - 0.5);
    del_b = -0.5 * (w - 0.5) - (-0.5) * (im.w - 0.5); 
    m[0] = del_m / del;
    b[0] = del_b / del;
    del = -0.5 - (im.h - 0.5);
    del_m = -0.5 - (h - 0.5);
    del_b = -0.5 * (h - 0.5) - (-0.5) * (im.h - 0.5);
    m[1] = del_m / del;
    b[1] = del_b / del;

    for(k = 0; k < im.c; ++k){
        for(j = 0; j < h; ++j){
            for(i = 0; i < w; ++i){
                i_scale = ((float)i - b[0]) / m[0];
                j_scale = ((float)j - b[1]) / m[1];
                set_pixel(resize_im, i, j, k, bilinear_interpolate(im, i_scale, j_scale, k));
            }
        }
    }
    return resize_im;
}

