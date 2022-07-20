#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <vector>
#include <string>
#include "ChargeAna.h"
namespace py = pybind11;
using namespace std;


PYBIND11_MODULE(cChargeAnalysis, m) {
    py::options options;
    options.disable_function_signatures();

    //bindings and documentation for individual functions
    py::class_<ChargeAnalysis>(m, "ChargeAnalysis")
        //-----------------------------------------------------
        // Constructor, Destructor and Access functions
        //-----------------------------------------------------        
        .def(py::init< >())
        .def_readwrite("boxs", &ChargeAnalysis::boxs)
        .def_readwrite("natoms", &ChargeAnalysis::natoms)
        .def_readwrite("rho", &ChargeAnalysis::rho)
        .def_readwrite("nxyz", &ChargeAnalysis::nxyz)
        .def_readwrite("dx3s", &ChargeAnalysis::dx3s)
        .def_readwrite("dxvol", &ChargeAnalysis::dxvol)
        .def_readwrite("debug", &ChargeAnalysis::debug)
        .def_readwrite("gridtest", &ChargeAnalysis::gridtest)
        .def_readwrite("threadnum", &ChargeAnalysis::threadnum)
        .def("ChargeAnalysisInit", &ChargeAnalysis::ChargeAnalysisInit,"To Set the BoxSize AtomNumber and import ChargeDensity", py::arg("rhofilename"), py::arg("boxs"), py::arg("natoms"))
        .def("readrho", &ChargeAnalysis::readrho,"readrho",py::arg("rhofilename"))
        .def("SumMethod_voro", &ChargeAnalysis::SumMethod_voro,"Used a CenterPointMethod to calculate ChargeDensity integration in the voronoi boxes", py::arg("v3s"), py::arg("vertex_positions"))
        .def("SumMethod_cut", &ChargeAnalysis::SumMethod_cut, "Used a CenterPointMethod to calculate ChargeDensity integration in the cutoff sphere", py::arg("atompos"), py::arg("cutoff"))
        .def("getden", &ChargeAnalysis::getden,"getden",py::arg("i"), py::arg("j"), py::arg("k"))
        .def("CD_pdf", &ChargeAnalysis::CD_pdf, "CD_pdf", py::arg("atompos"), py::arg("cutoff"), py::arg("histlow")=0, py::arg("histbin")=1000)
        .def("calculate_q_sumYdotCharge", &ChargeAnalysis::calculate_q_sumYdotCharge, "calculate_q_sumYdotCharge",py::arg("qs"), py::arg("atompos"), py::arg("cutoff"), py::arg("histlow") = 0, py::arg("histbin") = 1000)

        ;


#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}