
using namespace std;
class Methods {
protected:
 bool area_intersect(V& orig, V& dir, vector<Sphere>& spheres, V& hit, V& N, Skin& material) {
        double spheres_dist = numeric_limits<double>::max();
        for (size_t i = 0; i < spheres.size(); i++) {
            double dist_i;
            if (spheres[i].ray_intersect(orig, dir, dist_i) && dist_i < spheres_dist) {
                spheres_dist = dist_i;
                hit = orig + dir * dist_i;
                N = (hit - spheres[i].getCenter()).normalize();
                material = spheres[i].getSkin();
            }
        }
        return spheres_dist < 1000;
    }
V cast_ray(V& orig, V& dir, vector<Sphere>& spheres, vector<Light>& lights) {
        V point, N;
        Skin material;

        if (!area_intersect(orig, dir, spheres, point, N, material)) {
            return V(0.2, 0.4, 0.8);
        }

        double light_intensity = 0;
        for (size_t i = 0; i < lights.size(); i++) {
            V light_dir = (lights[i].getPosition() - point).normalize();
            double light_distance = (lights[i].getPosition() - point).norm();
            light_intensity += lights[i].getIntensity() * max(0.f, light_dir * N);
        }
        return material.getColor() * light_intensity;
    }
public:
 void render(vector<Sphere>& spheres, vector<Light>& lights) {
        const int   width = 1920;
        const int   height = 1080;
        const double fov = M_PI / 3.;
        vector<V> framebuffer(width * height);


        for (size_t j = 0; j < height; j++) {
            for (size_t i = 0; i < width; i++) {
                double x = (i + 0.5) - width / 2.;
                double y = -(j + 0.5) + height / 2.;
                double z = -height / (2. * tan(fov / 2.));
                V a = V(0, 0, 0);
                framebuffer[i + j * width] = cast_ray(a, V(x, y, z).normalize(), spheres, lights);
            }
        }

        ofstream image;
        image.open("./image.ppm", ios::binary);
        image << "P6\n" << width << " " << height << "\n255\n";
        for (size_t i = 0; i < height * width; ++i) {
            V& c = framebuffer[i];
            float maxi = max(c[0], max(c[1], c[2]));
            if (maxi > 1) c = c * (1. / maxi);
            for (size_t j = 0; j < 3; j++) {
                image << (char)(255 * max(0.f, min(1.f, framebuffer[i][j])));
            }
        }
        image.close();

    }
};
