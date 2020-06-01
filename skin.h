
class Skin {
public:
    Skin(const V& color) : color(color) {}
    Skin() : color(){}
    V getColor() {
        return color;
    }
private:
    V color;
 
};