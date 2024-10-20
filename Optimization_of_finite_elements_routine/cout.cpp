// bash : icpc cout.cpp ; ./a.out
// w10  : g++ cout.cpp ; .\cout.exe

#include <iostream>
#include <chrono>
#include <random>

using namespace std;
using namespace std::chrono;

int main() {
    const int iterations_outer = 1000; // Nombre d'itérations globales pour moyenner
    const int iterations_inner = 100000; // Nombre d'itérations pour chaque opération
    double tempsTotalAddition = 0, tempsTotalMultiplication = 0, tempsTotalDivision = 0;

    // Générateur de nombres aléatoires
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> disDouble(-10000000.0, 10000000.0);

    // Boucle globale pour moyenner le résultat
    for (int j = 0; j < iterations_outer; ++j) {
        // Définition des variables une seule fois
        double a = disDouble(gen);
        double b = disDouble(gen);
        double c = disDouble(gen);

        // Mesurer le temps total pour les additions
        auto startAddition = high_resolution_clock::now();
        for (int i = 0; i < iterations_inner; ++i) {
            volatile double result1 = a + b;
            volatile double result2 = a + c;
            volatile double result3 = c + b;
        }
        auto endAddition = high_resolution_clock::now();
        tempsTotalAddition += duration_cast<nanoseconds>(endAddition - startAddition).count() / static_cast<double>(iterations_inner) / 3;

        // Mesurer le temps total pour les multiplications
        auto startMultiplication = high_resolution_clock::now();
        for (int i = 0; i < iterations_inner; ++i) {
            volatile double result1 = a * b;
            volatile double result2 = a * c;
            volatile double result3 = c * b;
        }
        auto endMultiplication = high_resolution_clock::now();
        tempsTotalMultiplication += duration_cast<nanoseconds>(endMultiplication - startMultiplication).count() / static_cast<double>(iterations_inner) / 3;

        // Mesurer le temps total pour les divisions
        auto startDivision = high_resolution_clock::now();
        for (int i = 0; i < iterations_inner; ++i) {
            volatile double result1 = a / b;
            volatile double result2 = b / c;
            volatile double result3 = c / a;
        }
        auto endDivision = high_resolution_clock::now();
        tempsTotalDivision += duration_cast<nanoseconds>(endDivision - startDivision).count() / static_cast<double>(iterations_inner) / 3;
    }

    // Calculer les temps moyens sur les itérations globales
    tempsTotalAddition /= iterations_outer;
    tempsTotalMultiplication /= iterations_outer;
    tempsTotalDivision /= iterations_outer;

    cout << "Temps moyen pour l'addition : " << tempsTotalAddition << " nanosecondes" << endl;
    cout << "Performance arithmétique pour l'addition : " << 1 / (tempsTotalAddition * 1e-9) << " OPS" << endl;

    cout << "Temps moyen pour la multiplication : " << tempsTotalMultiplication << " nanosecondes" << endl;
    cout << "Performance arithmétique pour la multiplication : " << 1 / (tempsTotalMultiplication * 1e-9) << " OPS" << endl;

    cout << "Temps moyen pour la division : " << tempsTotalDivision << " nanosecondes" << endl;
    cout << "Performance arithmétique pour la division : " << 1 / (tempsTotalDivision * 1e-9) << " OPS" << endl;

    return 0;
}
