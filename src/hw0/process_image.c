#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c)
{
    // TODO Fill this in
    int pixel_index;
    x = (x < 0) ? 0: ( (x >= im.w) ? im.w-1: x);
    y = (y < 0) ? 0: ( (y >= im.h) ? im.h-1: y);
    c = (c < 0) ? 0: ( (c > 2) ? 2: c);

    pixel_index = x + y * im.w + c * im.h * im.w;
    if (x < 0 || y< 0 || c < 0 || x >= im.w || y >= im.h || c > 2) return 0;
    else return im.data[pixel_index];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    // TODO Fill this in
    int im_h = im.h;
    int im_w = im.w;
    const int pixel_index = c * im_h * im_w + x + y * im_w;

    // Check if the coordinate is valid:
    if (x >= 0 && x <= im_w && y >= 0 && y <= im_h){
        //printf("h = %d\n", im.h);
        //printf("w = %d\n", im.w);
        /*for (int i = 0; i<im_h*im_w*im.c; i++){
            printf("%d value = %f\n", i, im.data[i]);
        }*/
        //printf("B4 value = %f\n", im.data[pixel_index]);
        //printf("af value = %f\n", im.data[pixel_index]);
        // set pixel value at a certain coordinate:
        im.data[pixel_index] = v;
    }

}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    // TODO Fill this in
    int i,j,k;
    for(k = 0; k < im.c; ++k){
        for(j = 0; j < im.h; ++j){
            for(i = 0; i < im.w; ++i){
                copy.data[i + j * im.w + k * im.h * im.w] = get_pixel(im, i, j, k);
            }
        }
    }
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    // TODO Fill this in
    int i,j;
    for(j = 0; j < im.h; ++j){
        for(i = 0; i < im.w; ++i){
            gray.data[i + j * im.w] = 0.299 * get_pixel(im, i, j, 0) + 0.587 * get_pixel(im, i, j, 1) + 0.114 * get_pixel(im, i, j, 2);
        }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    // TODO Fill this in
    int i,j;
    for(j = 0; j < im.h; ++j){
        for(i = 0; i < im.w; ++i){
            set_pixel(im, i, j, c, get_pixel(im, i, j, c) + v);
        }
    }
    
}

void clamp_image(image im)
{
    // TODO Fill this in
    int i,j,k;
    for (k = 0; k<im.c; ++k){
        for(j = 0; j < im.h; ++j){
            for(i = 0; i < im.w; ++i){
                im.data[i + j * im.w + k * im.w * im.h] = get_pixel(im, i, j, k) > 1 ? 1:get_pixel(im, i, j, k);
                im.data[i + j * im.w + k * im.w * im.h] = get_pixel(im, i, j, k) < 0 ? 0:get_pixel(im, i, j, k);
            }
        }
    }
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    // TODO Fill this in
    int i,j;
    //float min_Hp = 1000, max_H = -100, min_H = 100;
    float V, m, C, S, Hp, H, R, G, B;
    for(j = 0; j < im.h; ++j){
        for(i = 0; i < im.w; ++i){
            R = get_pixel(im, i, j, 0);
            G = get_pixel(im, i, j, 1);
            B = get_pixel(im, i, j, 2);
            V = three_way_max(R, G, B);
            m = three_way_min(R, G, B);
            C = V - m;
            if (R == 0 && G == 0 && B == 0) S = 0;
            else S = C / V;
            
            if (C == 0){
                Hp = 0;
            }    
            else {
                Hp = (R > G) ? ( (R > B) ? (G - B)/C : ((R - G) / C) + 4 ) : 
                    ( (G > B) ? (B - R)/C + 2 : ((R - G) / C) + 4);
            }
            H = (Hp < 0) ? Hp/6 + 1 : Hp/6;
            set_pixel(im, i, j, 0, H);
            set_pixel(im, i, j, 1, S);
            set_pixel(im, i, j, 2, V);
            
        }
    }
}

void hsv_to_rgb(image im)
{
    int i, j;
    float H, S, V, C, min, Hp, R, G, B;
    // TODO Fill this in

    for(j = 0; j < im.h; ++j){
        for(i = 0; i < im.w; ++i){
            H = get_pixel(im, i, j, 0);
            S = get_pixel(im, i, j, 1);
            V = get_pixel(im, i, j, 2);
            C = V * S;
            min = V - C;
            if (H >= 0 && H < 5./6.){
                Hp = 6 * H;
            }
            else {
                Hp = 6 * (H - 1);
            }
            //printf ("H = %f, Hp = %f\n", H, Hp);
            //scanf("123");
            if (Hp >= -1 && Hp < 0){
                R = V;
                G = min;
                B = G - Hp * C;
            }
            else if (Hp >= 0 && Hp < 1){
                R = V;
                B = min;
                G = B + Hp * C;
            }
            else if (Hp >= 1 && Hp < 2){
                G = V;
                B = min;
                R = B - (Hp - 2) * C;
            }
            else if (Hp >= 2 && Hp < 3){
                G = V;
                R = min;
                B = R + (Hp - 2) * C;
            }
            else if (Hp >= 3 && Hp < 4){
                B = V;
                R = min;
                G = R - (Hp - 4) * C;
            }
            else if (Hp >= 4 && Hp < 5){
                B = V;
                G = min;
                R = G + (Hp - 4) * C;
            }
            set_pixel(im, i, j, 0, R);
            set_pixel(im, i, j, 1, G);
            set_pixel(im, i, j, 2, B);
        }
    }
}

void scale_image(image im, int c, float v)
{
    int i, j;
    for(j = 0; j < im.h; ++j){
        for(i = 0; i < im.w; ++i){
            set_pixel(im, i, j, c, get_pixel(im, i, j, c) * v);
        }
    }
}
