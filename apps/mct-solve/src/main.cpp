#include <iostream>
#include <concepts>

template <typename T>
concept ReturnPolicy = requires{ typename T::Type; };

template <class T, template<typename> class RETURN_POLICY>
struct DataTransformer
{
    using ReturnType = RETURN_POLICY<T>::Type;

    ReturnType calculateStep(const T& previous) { return {};}

    
};

int main()
{
    std::cerr << "Hello mct-solve\n";
    return 0;
}