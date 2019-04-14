#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_poly.h>
#include "SOIL.h"

struct OrbitInfo {
    double theta_m;
    double d;
    
    OrbitInfo(double theta_m, double d) : theta_m(theta_m), d(d) {}
};

OrbitInfo CalOrbit(double cc) {
    double e1, e2, e3;
    gsl_poly_solve_cubic(-0.5, 0.0, cc / 2.0, &e3, &e2, &e1);
    double w = sqrt((e1 - e2) / 2.0);
    double m = (e2 - e3) / (e1 - e2);
    double k = sqrt(m / (1 + m));
    double theta0 = acos(sqrt(e2 / (e2 - e3)));
    double theta_m = 1 / sqrt(1 + m) / w * (gsl_sf_ellint_F(M_PI_2, k, 0) - gsl_sf_ellint_F(theta0, k, 0));
    double d = (e2 - e3) / sqrt(-2.0 * e2 * e3 * (2 * e2 - e3) * ((e1 - e2) * (e2 - e3) + e2 * e2));
    return OrbitInfo(theta_m, d);
}

int main(int argc, const char * argv[]) {

    const char * read_path = "/Users/ycx/Desktop/star.png";
    int width, height, channels;
    unsigned char * data = SOIL_load_image(read_path, &width, &height, &channels, 0);
    
    if (channels != 3 && channels != 4) {
        return -1;
    }
    
    if (width != height) {
        // 注意，你的图片等长宽应当相等。如果使用天空盒等话，这是必须的。
        // 虽然我没有找到合适的天空盒，这导致黑洞周围出现黑色的环带，因为一些光线没有打到图片上。
        // 也许你找到了合适的天空盒，帮我改一下代码？
        return -1;
    }
    int scale_in = width;
    
    // 下面的参数的意义见README.md
    double scale = 50.0;
    double b = 50.0;
//    double gamma = atan(sqrt(2) * scale / b);
    double d = 20.0;
    const int scale_out = 800;
    double a = (b + d) / scale;
    
    unsigned char out_image[2 * scale_out * 2 * scale_out * 3];
    memset(out_image, 0, sizeof(out_image));
    
    // 逐点计算像素值
    for (int i = 0; i < 2 * scale_out; ++i) {
        for (int j = 0; j < 2 * scale_out; ++j) {
            double x, y;
            x = double(i) / double(scale_out) - 1.0;
            y = 1.0 - double(j) / double(scale_out);
//            x *= 1.2;
//            y *= 1.2;
            double r_2 = x * x + y * y;
            double r = sqrt(r_2);
            double alpha = atan(r / a);
            const double cc0 = 1.0 / (d * d) - 2.0 / (d * d * d);
            double cc = cc0 + a * a / (d * d * r_2);
            
            // 我们只求光子方程的外部解
            if (cc <= 0 || cc >= 1.0 / 27.0) {
                continue;
            }
            
            OrbitInfo orbit = CalOrbit(cc);
            double theta_m = orbit.theta_m;
            double D = orbit.d;
            
            double th = 2 * theta_m - alpha;
//            double remainder = fmod(2 * theta_m - alpha, 2 * M_PI);
            
            
//            if (remainder <= M_PI_2 - gamma || remainder >= M_PI_2 + gamma) {
//                continue;
//            }
            
            double rr = (-D - b * sin(th)) / cos(th);
            
            double xx = rr * x / r;
            double yy = rr * y / r;
            
            int xi = floor((xx / scale + 1.0) * 0.5 * scale_in);
            int yi = floor((-yy / scale + 1.0) * 0.5 * scale_in);
            
            // 由于没有天空盒。采用暴力的平铺处理。
            // 如果你有天空盒，应重写这些。
            if (xi < 0 || xi >= scale_in || yi < 0 || yi >= scale_in) {
                xi %= scale_in;
                yi %= scale_in;
                if (xi < 0) {
                    xi += scale_in;
                }
                if (yi < 0) {
                    yi += scale_in;
                }
            }
            
            out_image[3 * (i * 2 * scale_out + j)] = data[channels * (xi * scale_in + yi)];
            out_image[3 * (i * 2 * scale_out + j) + 1] = data[channels * (xi * scale_in + yi) + 1];
            out_image[3 * (i * 2 * scale_out + j) + 2] = data[channels * (xi * scale_in + yi) + 2];
        }
    }
    const char * save_path = "/Users/ycx/Desktop/photon.bmp";
    SOIL_save_image(save_path, SOIL_SAVE_TYPE_BMP, 2 * scale_out, 2 * scale_out, 3, out_image);
    
    return 0;
}
