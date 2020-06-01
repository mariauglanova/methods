
class Light {
public:
    Light(const V& p, const double i) : position(p), intensity(i) {}
    V getPosition() {
        return position;
    }
    double getIntensity() {
        return intensity;
    }
private:
    V position;
    double intensity;
};
