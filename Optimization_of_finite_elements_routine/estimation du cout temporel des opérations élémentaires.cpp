#include <iostream>
#include <chrono>

using namespace std;
using namespace std::chrono;

// estimation du cout temporel des opérations élémentaires
int main() {
    const int iterations = 1000000; // Nombre d'itérations pour chaque opération
    double tempsTotalAddition = 0, tempsTotalMultiplication = 0, tempsTotalDivision = 0;

    // Mesurer le temps total pour les additions
    auto startAddition = high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        volatile int a = 10;
        volatile int b = 20;
        volatile int result = a + b;
    }
    auto endAddition = high_resolution_clock::now();
    tempsTotalAddition = duration_cast<nanoseconds>(endAddition - startAddition).count() / static_cast<double>(iterations);

    // Mesurer le temps total pour les multiplications
    auto startMultiplication = high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        volatile int a = 10;
        volatile int b = 20;
        volatile int result = a * b;
    }
    auto endMultiplication = high_resolution_clock::now();
    tempsTotalMultiplication = duration_cast<nanoseconds>(endMultiplication - startMultiplication).count() / static_cast<double>(iterations);

    // Mesurer le temps total pour les divisions
    auto startDivision = high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        volatile double a = 10.0;
        volatile double b = 2.0;
        volatile double result = a / b;
    }
    auto endDivision = high_resolution_clock::now();
    tempsTotalDivision = duration_cast<nanoseconds>(endDivision - startDivision).count() / static_cast<double>(iterations);

    cout << "Temps moyen pour l'addition : " << tempsTotalAddition << " nanosecondes" << endl;
    cout << "Temps moyen pour la multiplication : " << tempsTotalMultiplication << " nanosecondes" << endl;
    cout << "Temps moyen pour la division : " << tempsTotalDivision << " nanosecondes" << endl;

    return 0;
}
