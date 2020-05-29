#ifndef __Header_h__
#define __Header_h__

#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

template <size_t DIM, typename T> struct vec {
    vec() { for (size_t i = DIM; i--; data_[i] = T()); }
    T& operator[](const size_t i) { assert(i < DIM); return data_[i]; }
    const T& operator[](const size_t i) const { assert(i < DIM); return data_[i]; }
private:
    T data_[DIM];
};

typedef vec<3, float> Vec3f;

template <typename T> struct vec<3, T> {
    vec() : x(T()), y(T()), z(T()) {}
    vec(T X, T Y, T Z) : x(X), y(Y), z(Z) {}
    T& operator[](const size_t i) { assert(i < 3); return i <= 0 ? x : (1 == i ? y : z); }
    const T& operator[](const size_t i) const { assert(i < 3); return i <= 0 ? x : (1 == i ? y : z); }
    float norm() { return std::sqrt(x * x + y * y + z * z); }
    vec<3, T>& normalize(T l = 1) { *this = (*this) * (l / norm()); return *this; }
    T x, y, z;
};

template<size_t DIM, typename T> T operator*(const vec<DIM, T>& lhs, const vec<DIM, T>& rhs) {
    T ret = T();
    for (size_t i = DIM; i--; ret += lhs[i] * rhs[i]);
    return ret;
}

template<size_t DIM, typename T>vec<DIM, T> operator+(vec<DIM, T> lhs, const vec<DIM, T>& rhs) {
    for (size_t i = DIM; i--; lhs[i] += rhs[i]);
    return lhs;
}

template<size_t DIM, typename T>vec<DIM, T> operator-(vec<DIM, T> lhs, const vec<DIM, T>& rhs) {
    for (size_t i = DIM; i--; lhs[i] -= rhs[i]);
    return lhs;
}

template<size_t DIM, typename T, typename U> vec<DIM, T> operator*(const vec<DIM, T>& lhs, const U& rhs) {
    vec<DIM, T> ret;
    for (size_t i = DIM; i--; ret[i] = lhs[i] * rhs);
    return ret;
}

template<size_t DIM, typename T> vec<DIM, T> operator-(const vec<DIM, T>& lhs) {
    return lhs * T(-1);
}



class Colour {
public:
    Colour(const Vec3f& color) : dif_color(color) {}
    Colour() : dif_color() {}
    Vec3f dif_color;
};

class Sphere {
public:
    Vec3f center;
    float radius;
    Colour colour;

    Sphere(const Vec3f& c, const float r, const Colour col) : center(c), radius(r), colour(col) {}

    bool ray_intersect(const Vec3f& orig, const Vec3f& dir, float& t0) const {
        Vec3f L = center - orig;
        float tca = L * dir;
        float d2 = L * L - tca * tca;
        if (d2 > radius * radius) return false;
        float thc = sqrtf(radius * radius - d2);
        t0 = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    }
};

class Light {
public:
    Light(const Vec3f& p, const float& i) : position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};

class Picture {

public:
    Picture() {}

    std::vector<Sphere> spheres;
    void add(Sphere sphere) {
        spheres.push_back(sphere);
    }
private:
    bool pic_intersect(const Vec3f& orig, const Vec3f& dir, const std::vector<Sphere>& spheres, Vec3f& hit, Vec3f& N, Colour& colour) {
        float spheres_dist = std::numeric_limits<float>::max();
        for (size_t i = 0; i < spheres.size(); i++) {
            float dist_i;
            if (spheres[i].ray_intersect(orig, dir, dist_i) && dist_i < spheres_dist) {
                spheres_dist = dist_i;
                hit = orig + dir * dist_i;
                N = (hit - spheres[i].center).normalize();
                colour = spheres[i].colour;
            }
        }

        return spheres_dist < 1000;
    }

    Vec3f cast_ray(const Vec3f& orig, const Vec3f& dir, const std::vector<Sphere>& spheres, const Light lights) {
        Vec3f point, N;
        Colour colour;

        if (!pic_intersect(orig, dir, spheres, point, N, colour)) {
            return Vec3f(0, 0, 0); // fon
        }


        float dif_light_intensity = 0;
            Vec3f light_dir = (lights.position - point).normalize();
            dif_light_intensity += lights.intensity * std::max(0.f, light_dir * N);
        return colour.dif_color * dif_light_intensity;
    }
public:
    void render(const std::vector<Sphere>& spheres, const Light lights) {
        const int   width = 1024;
        const int   height = 768;
        const float fov = M_PI / 3.;
        std::vector<Vec3f> framebuffer(width * height);

#pragma omp parallel for
        for (size_t j = 0; j < height; j++) { // actual rendering loop
            for (size_t i = 0; i < width; i++) {
                float dir_x = (i + 0.5) - width / 2.;
                float dir_y = -(j + 0.5) + height / 2.;    // this flips the image at the same time
                float dir_z = -height / (2. * tan(fov / 2.));
                framebuffer[i + j * width] = cast_ray(Vec3f(0, 0, 0), Vec3f(dir_x, dir_y, dir_z).normalize(), spheres, lights);
            }
        }

        std::ofstream ofs; 
        ofs.open("Sun.ppm", std::ios::binary);
        ofs << "P6\n" << width << " " << height << "\n255\n";
        for (size_t i = 0; i < height * width; ++i) {
            Vec3f& c = framebuffer[i];
            float max = std::max(c[0], std::max(c[1], c[2]));
            if (max > 1) c = c * (1. / max);
            for (size_t j = 0; j < 3; j++) {
                ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
            }
        }
        ofs.close();
    }


};


#endif //__Header_h__
