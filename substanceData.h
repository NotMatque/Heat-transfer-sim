#ifndef SUBSTANCEDATA_H
#define SUBSTANCEDATA_H

struct SubstanceData {
    std::string name;
    double conductivity;
    double alpha;
    double density;
    double specificHeat;

    friend std::ostream& operator<<(std::ostream& os, const SubstanceData& sd) {
        os << "Substance Name: " << sd.name << "\n";
        os << "Conductivity: " << sd.conductivity << " [W / (m * C)]\n";
        os << "Alpha: " << sd.alpha << " [W / (m^2 * C)]\n";
        os << "Density: " << sd.density << " [kg / m^3]\n";
        os << "SpecificHeat: " << sd.specificHeat << " [J / (kg * C)]\n";

        return os;
    }
};

#endif //SUBSTANCEDATA_H
