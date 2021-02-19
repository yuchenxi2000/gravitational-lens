import ctypes
import numpy as np
import cv2

# 自己改路径

background = cv2.imread('../src/milky-way-6000.png')

in_img = np.array(background)
in_width = ctypes.c_int(background.shape[1])
in_height = ctypes.c_int(background.shape[0])
channels = ctypes.c_int(background.shape[2])

lib = np.ctypeslib.load_library(libname='liblensing', loader_path='.')

img_scale = ctypes.c_double(1200.0)

pos = np.array([0.0, 0.0, 50.0], dtype=np.float64)
view = np.array([0.0, 0.0, -0.8], dtype=np.float64)
x_axis = np.array([1.0, 0.0, 0.0], dtype=np.float64)
sphere_r = ctypes.c_double(800.0)

lib.cal_image.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),
    ctypes.c_double,
    ctypes.c_double,
    np.ctypeslib.ndpointer(dtype=np.uint8, ndim=3, flags="C_CONTIGUOUS"),
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.uint8, ndim=3, flags="C_CONTIGUOUS"),
    ctypes.c_int,
    ctypes.c_int,
]

fps = 24  # 视频每秒24帧
size = (1920, 1080)  # 需要转为视频的图片的尺寸

video = cv2.VideoWriter('./Video.mp4', cv2.VideoWriter_fourcc('m', 'p', '4', 'v'), fps, size, True)

for i in range(400):
    out_img = np.zeros([size[0], size[1], 3], dtype=np.uint8)
    out_width = ctypes.c_int(out_img.shape[0])
    out_height = ctypes.c_int(out_img.shape[1])
    pos[0] = 40.0 - i * 0.2
    lib.cal_image(
        pos,
        view,
        x_axis,
        img_scale,
        sphere_r,
        in_img,
        in_width,
        in_height,
        channels,
        out_img,
        out_width,
        out_height
    )
    video.write(out_img.reshape([out_height.value, out_width.value, 3]))

video.release()
cv2.destroyAllWindows()
