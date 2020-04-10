#include <vector>
#include <pybind11/stl.h>

using namespace std;
namespace py = pybind11;

template<typename T>
class Dataset{
    vector<T> data;
public:
    void populate(py::list){
        for(auto item : list){
            data.push_back()
        }
    }


}