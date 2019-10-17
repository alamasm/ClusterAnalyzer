#include "Interface.h"
using namespace std;

int main() {
    Control control;
    Interface interface(control);
    interface.start();
}