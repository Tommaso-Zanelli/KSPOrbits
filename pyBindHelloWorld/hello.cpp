#include <pybind11/pybind11.h>
#include <iostream>

// https://pybind11.readthedocs.io/en/stable/advanced/pycpp/object.html

namespace py = pybind11;

std::string hello(std::string name)
{
    return "Hello, " + name + "!";
}

void printHello(std::string name)
{
    std::cout << "Hello, " << name << "!" << '\n';
    return;
}

PYBIND11_MODULE(hello, m)
{
    m.doc() = "hello world module";
    m.def("hello", &hello, "A function that says hello");
    m.def("printHello", &printHello, "A function that prints hello");
}
