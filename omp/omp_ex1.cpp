#include <iostream>
#include <omp.h>

int main() {

#pragma omp parallel
    {
        int ID = omp_get_thread_num();
#pragma omp critical
        {
            std::cout << "hello(" << ID << ")" << std::endl;
            std::cout << " world(" << ID << ")" << std::endl;
        }
    }
    return 0;
}