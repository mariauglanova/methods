#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "space.h"
#include "light.h"
#include "skin.h"
#include "sphere.h"
#include "methods.h"
 
using namespace std;


int main() {
    Skin one(V(0.7, 0.2, 0.3));
    Skin two(V(0.4, 0.5, 0.1));

    vector<Sphere> spheres;
    spheres.push_back(Sphere(V(-8, 0, -13), 2, one));
    spheres.push_back(Sphere(V(10, 5, -15), 3, two));

    vector<Light> lights;
    lights.push_back(Light(V(0, 25, 25), 2.5));

    Methods met;
    met.render(spheres, lights);
   

    return 0;
}


