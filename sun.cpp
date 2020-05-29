#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "Header.h"


int main()
{
    Colour   yellow(Vec3f(1, 1, 0)), green(Vec3f(0, 1, 0));

    Picture pic = Picture();

    Sphere m0 = Sphere(Vec3f(-7, 5, -15), 3, yellow), m1 = Sphere(Vec3f(9, -9, -19), 5, green);

    pic.add(m0);

    pic.add(m1);

    Light light = Light(Vec3f(0, 0, -10), 1);

    pic.render(pic.spheres, light);

    return 0;

}
