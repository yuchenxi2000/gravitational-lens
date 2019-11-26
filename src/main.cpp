/*
 * black hole lensing effect
 */
#include <iostream>

#include "glm/glm.hpp"
#define GLM_ENABLE_EXPERIMENTAL
#include "glm/gtx/transform.hpp"
#include "glm/gtx/vector_angle.hpp"

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_poly.h>

#include "image.hpp"
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
 * param cc: orbit constant in formula
 */
ScatterInfo scatter(double cc) {
    ScatterInfo res;
    
    res.d = 1 / sqrt(cc);
    
    double e1, e2, e3;
    gsl_poly_solve_cubic(-0.5, 0.0, cc / 2.0, &e3, &e2, &e1);
    double w = sqrt(8.0 / (e1 - e3));
    double k = sqrt((e2 - e3) / (e1 - e3));
    double theta0 = acos(sqrt(e2 / (e2 - e3)));
    res.theta = w * (gsl_sf_ellint_Kcomp(k, 0) - gsl_sf_ellint_F(theta0, k, 0));
    
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

class ImageSampler {
    Image image;
    // distance between black hole and background image
    double background_d;
public:
    ImageSampler(const char * img_path, int distance) : image(img_path, 3) {
        if (image.isNull()) {
            cout << "cannot find image!" << endl;
            exit(-1);
        }
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
        if (R.x > 0) {
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

int main(int argc, const char * argv[]) {
    // black hole's position
    Vec3 hole_pos = {0, 0, 0};
    Line3D l1, l2;
    
    // sample from image
    ImageSampler sampler("milky-way-6000.png", 800);
    assert(sampler.get_channel() == 3);
    
    // the output image
    Image outimage;
    outimage.blankImage(2000, 2000, 3);
    
    // observer's position
    l1.r = {0, 0, 50};
    
    // adjust angle between view ray and z axis
    // (adjust the FOV)
    double factor = 0.8;
    l1.v.z = -factor;
    
    for (int i = 0; i < outimage.width; ++i) {
        for (int j = 0; j < outimage.height; ++j) {
            l1.v.x = (double)i/(double)outimage.width-0.5;
            l1.v.y = 0.5-(double)j/(double)outimage.height;
            if (transform(l1, l2, hole_pos)) {
                sampler.sample_sphere(l2, &outimage(i, j, 0));
            }
        }
    }
    
    // save image (bmp format only)
    outimage.save("example-out.bmp");
    return 0;
}
