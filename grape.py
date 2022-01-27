from uwimg import *

im = load_image("data/WEI.jpg")
nms = nms_image(im, 1)
#clamp_image(nms)
save_image(nms, "nms_test_WEI")