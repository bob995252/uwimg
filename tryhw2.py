from uwimg import *
"""
im = load_image("data/dog.jpg")
f = make_box_filter(7)
blur = convolve_image(im, f, 1)
save_image(blur, "dog-box7")

im = load_image("data/dog.jpg")
f = make_box_filter(7)
blur = convolve_image(im, f, 1)
thumb = nn_resize(blur, blur.w//7, blur.h//7)
save_image(thumb, "dogthumb")

im = load_image("data/feng.jpg")
f = make_highpass_filter()
blur = convolve_image(im, f, 0)
clamp_image(blur)
save_image(blur, "fenghighpass")

im = load_image("data/riva.jpg")
f = make_emboss_filter()
blur = convolve_image(im, f, 1)
clamp_image(blur)
save_image(blur, "riva_emboss.jpg")

im = load_image("data/dog.jpg")
f = make_gaussian_filter(2)
blur = convolve_image(im, f, 1)
save_image(blur, "dog-gauss2")

im = load_image("data/dog.jpg")
f = make_gaussian_filter(2)
lfreq = convolve_image(im, f, 1)
hfreq = im - lfreq
reconstruct = lfreq + hfreq
save_image(lfreq, "low-frequency")
save_image(hfreq, "high-frequency")
save_image(reconstruct, "reconstruct")
"""
im = load_image("data/WEI.jpg")
res = sobel_image(im)
mag = res[0]
#theta = res[1]
feature_normalize(mag)
#save_image(theta, "theta")
save_image(mag, "dog_WEI.jpg")
"""


im = load_image("data/dog.jpg")
color_im = colorize_sobel(im)
save_image(color_im, "dog_color")


chiu = load_image("data/chiu_.png")
feng = load_image("data/feng.jpg")
chiu = rgb_to_grayscale(chiu)
feng = rgb_to_grayscale(feng)
chiu = bilinear_resize(chiu, feng.w, feng.h)
print(chiu.w, chiu.h, chiu.c)
print(feng.w, feng.h, feng.c)
f = make_gaussian_filter(2)
lfreq_c = convolve_image(chiu, f, 0)
lfreq_f = convolve_image(feng, f, 0)
hfreq_c = sub_image(chiu, lfreq_c)
reconstruct = add_image(lfreq_f , hfreq_c)
clamp_image(reconstruct)
save_image(reconstruct, "chiu2feng")
"""