#include "nori/bitmap.h"
#include <iostream>

using namespace nori;

int main(int argc, char **argv) {

    Bitmap bitmap = Bitmap(Vector2i(4,2));
    bitmap(0) = Color3f(60, 60, 60);
    bitmap(1) = Color3f(20,20,20);
    bitmap(2) = Color3f(4,4,4);
    bitmap(3) = Color3f(16,16,16);
    bitmap(4) = Color3f(10, 10, 10);
    bitmap(5) = Color3f(30,30,30);
    bitmap(6) = Color3f(20,20,20);
    bitmap(7) = Color3f(2,2,2);
    // bitmap(8) = Color3f(60, 60, 60);
    // bitmap(9) = Color3f(20,20,20);
    // bitmap(10) = Color3f(4,4,4);
    // bitmap(11) = Color3f(16,16,16);

    bitmap.saveEXR("test");

}