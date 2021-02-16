/*
 * black hole lensing effect
 */
#include <iostream>

#include <glm/glm.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/transform.hpp>
#include <glm/gtx/vector_angle.hpp>

// gsl 不能并行。。
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <omp.h>

using std::cout;
using std::endl;
using std::ostream;

typedef glm::dvec3 Vec3;
typedef glm::dmat3 Mat3;

ostream & operator << (ostream & os, const Vec3 & v) {
    os << v.x << '\t' << v.y << '\t' << v.z;
    return os;
}

ostream & operator << (ostream & os, const Mat3 & m) {
    for (int i = 0; i < 3; ++i) {
        os << m[i][0] << '\t' << m[i][1] << '\t' << m[i][2] << '\n';
    }
    return os;
}

struct Line3D {
    Vec3 v; // direction vector
    Vec3 r; // a point on the line
};

struct ScatterInfo {
    double theta; // rotate angle
    double d; // distance between ray and center of hole
};

/*
 * solve cubic polynomial
 */
void poly_solve_cubic(double c2, double c1, double c0, double * e1, double * e2, double * e3) {
    Eigen::Matrix3d m;
    m(1, 0) = 1.0;
    m(2, 1) = 1.0;
    m(0, 2) = -c0;
    m(1, 2) = -c1;
    m(2, 2) = -c2;
    auto r = m.eigenvalues();
    *e3 = r(0).real();
    *e2 = r(1).real();
    *e1 = r(2).real();
}

/*
 * param cc: orbit constant in formula
 */
ScatterInfo scatter(double cc) {
    ScatterInfo res;
    
    res.d = 1 / sqrt(cc);
    
    double e1, e2, e3;
    poly_solve_cubic(-0.5, 0.0, cc / 2.0, &e1, &e2, &e3);
    double w = sqrt(8.0 / (e1 - e3));
    double k = sqrt((e2 - e3) / (e1 - e3));
    double theta0 = acos(sqrt(e2 / (e2 - e3)));
    res.theta = w * (std::ellint_1(k, M_PI_2) - std::ellint_1(k, theta0));
    
    return res;
}

/*
 * the ray after photon passes by the black hole
 * returns false if photon fell into the black hole
 */
bool transform(const Line3D & l1, Line3D & l2, const Vec3 & hole_pos) {
    // calculate cc
    double r = glm::distance(l1.r, hole_pos);
    Vec3 d = glm::normalize(l1.v);
    Vec3 vr = hole_pos - l1.r;
    double rd = glm::dot(vr, d);
    double cc =  1 / (r * r - rd * rd) - 2 / (r * r * r);
    
    // photon fell into the black hole
    if (cc <= 0 || cc >= 1.0 / 27.0) {
        return false;
    }
    
    Vec3 k = glm::normalize(glm::cross(l1.v, vr));
    ScatterInfo scatterinfo = scatter(cc);
    Mat3 m = glm::rotate(scatterinfo.theta-M_PI, k);
    
    l2.v = m * l1.v;
    l2.r = hole_pos + scatterinfo.d * glm::normalize(glm::cross(l2.v, k));
    
    return true;
};

class Image {
public:
    unsigned char * data;
    int width, height, channel;
    Image() {
        data = 0;
    }
    Image(unsigned char * data, int width, int height, int channel) {
        this->load(data, width, height, channel);
    }
    ~Image() {}

    void load(unsigned char * data, int width, int height, int channel) {
        this->data = data;
        this->width = width;
        this->height = height;
        this->channel = channel;
    }

    operator unsigned char * () {
        return data;
    }
    unsigned char & operator () (int w, int h, int ch) const {
        return data[channel * (width * h + w) + ch];
        // return data[channel * (height * w + h) + ch];
    }
    
    int size() {
        return width * height * channel;
    }
    
    // get a new blank image
    void blankImage(int w, int h, int ch) {
        this->width = w;
        this->height = h;
        this->channel = ch;
        if (data) delete [] data;
        data = new unsigned char[w * h * ch]();
    }
    
    bool isNull() {
        return data == 0;
    }
    
    // delete image buffer
    void clear() {
        if (data) delete [] data;
        this->data = 0;
        this->width = 0;
        this->height = 0;
        this->channel = 0;
    }
};

class ImageSampler {
    Image image;
    // distance between black hole and background image
    double background_d;
public:
    ImageSampler(unsigned char * data, int width, int height, int channel, int distance) : image(data, width, height, channel) {
        background_d = distance;
    }
    void sample_repeat(const Line3D & l, unsigned char outcolor[]) {
        double factor = background_d / l.v.z;
        Vec3 v = factor * l.v;
        int xi = (int)(v.x+image.width / 2) % image.width;
        int yi = (int)(image.height / 2-v.y) % image.height;
        if (xi < 0) xi += image.width;
        if (yi < 0) yi += image.height;
        for (int c = 0; c < image.channel; ++c) {
            outcolor[c] = image(xi, yi, c);
        }
    }
    void sample_mirrored_repeat(const Line3D & l, unsigned char outcolor[]) {
        double factor = background_d / l.v.z;
        Vec3 v = factor * l.v;
        int xi = (int)(v.x+image.width / 2);
        int yi = (int)(image.height / 2-v.y);
        xi %= 2 * image.width;
        yi %= 2 * image.height;
        if (xi < 0) {
            xi += 2 * image.width;
        }
        if (yi < 0) {
            yi += 2 * image.height;
        }
        if (xi >= image.width) {
            xi = 2 * image.width - 1 - xi;
        }
        if (yi >= image.height) {
            yi = 2 * image.height - 1 - yi;
        }
        for (int c = 0; c < image.channel; ++c) {
            outcolor[c] = image(xi, yi, c);
        }
    }
    void sample_normal(const Line3D & l, unsigned char outcolor[]) {
        if (l.v.z > 0) return;
        double t = (-background_d-l.r.z) / l.v.z;
        Vec3 v = l.r + t * l.v;
        int xi = v.x+image.width/2;
        int yi = image.height/2-v.y;
        if (xi >= 0 && xi < image.width && yi >= 0 && yi < image.height) {
            for (int c = 0; c < image.channel; ++c) {
                outcolor[c] = image(xi, yi, c);
            }
        }
    }
    void sample_sphere(const Line3D & l, unsigned char outcolor[]) {
        double rv = glm::dot(l.r, l.v);
        double vv = glm::length2(l.v);
        double rr = glm::length2(l.r);
        double dd = background_d * background_d;
        double t = (-rv+sqrt(rv*rv-rr*vv+dd))/vv;
        Vec3 R = l.r+t*l.v;
        
        // get xi, yi
        double phi = acos(R.z/sqrt(R.x*R.x+R.z*R.z));
        if (R.x < 0) {
            phi = 2*M_PI - phi;
        }
        double theta = M_PI_2-asin(R.y/background_d);
        int xi = phi/(2*M_PI)*(image.width+1);
        int yi = theta/(M_PI)*(image.height+1);
        for (int c = 0; c < image.channel; ++c) {
            outcolor[c] = image(xi, yi, c);
        }
    }
    int get_channel() {
        return image.channel;
    }
};

struct timer {
#ifdef _OPENMP
    double start;
#else
    clock_t start;
#endif
    
    const char * func_str;
    timer(const char * func_str) {
#ifdef _OPENMP
        start = omp_get_wtime();
#else
        start = clock();
#endif
        this->func_str = func_str;
    }
    ~timer() {
#ifdef _OPENMP
        double s = omp_get_wtime() - start;
#else
        double s = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
#endif
        printf("[timer] %s: %lf s\n", func_str, s);
    }
};

Vec3 d3_to_vec3(double d[3]) {
    return Vec3(d[0], d[1], d[2]);
}

extern "C" {

    void cal_image(double pos_d3[3], double view_d3[3], double x_axis_d3[3], double img_scale, double sphere_r, unsigned char * in_img, int in_width, int in_height, int in_channel, unsigned char * out_img, int out_width, int out_height) {
        timer timer0(__FUNCTION__);
        printf("[cal_image] openmp processors = %d\n", omp_get_num_procs());
        // black hole's position
        Vec3 hole_pos = {0, 0, 0};
        ImageSampler sampler(in_img, in_width, in_height, in_channel, sphere_r);
        Image outimage(out_img, out_width, out_height, in_channel);

        Vec3 view = d3_to_vec3(view_d3);
        Vec3 x_axis = glm::normalize(d3_to_vec3(x_axis_d3));
        Vec3 y_axis = glm::normalize(glm::cross(x_axis, view));

        // std::cout << d3_to_vec3(pos_d3) << std::endl;
        // std::cout << view << std::endl;
        // std::cout << x_axis << std::endl;
        // std::cout << img_scale << std::endl;
        // std::cout << sphere_r << std::endl;
        // std::cout << in_width << std::endl;
        // std::cout << in_height << std::endl;
        // std::cout << in_channel << std::endl;
        // std::cout << out_width << std::endl;
        // std::cout << out_height << std::endl;
        // std::cout << std::hex << (unsigned long long)in_img << std::endl;
        // std::cout << std::hex << (unsigned long long)out_img << std::endl;

    #pragma omp parallel for
        for (int i = 0; i < outimage.width; ++i) {
            for (int j = 0; j < outimage.height; ++j) {
                Line3D l1, l2;
                l1.r = d3_to_vec3(pos_d3);
                double x = ((double)i - 0.5 * (double)outimage.width) / img_scale;
                double y = ((double)j - 0.5 * (double)outimage.height) / img_scale;
                l1.v = view + x * x_axis + y * y_axis;
                if (transform(l1, l2, hole_pos)) {
                    sampler.sample_sphere(l2, &outimage(i, j, 0));
                }
            }
        }
    }

}
