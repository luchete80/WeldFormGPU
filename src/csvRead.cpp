#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

void readCSV(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string cell;

        while (std::getline(ss, cell, ',')) {
            row.push_back(cell);
        }

        // Process the row (e.g., print)
        for (const auto& col : row) {
            std::cout << col << " ";
        }
        std::cout << std::endl;
    }

    file.close();
}

int main() {
    const std::string filename = "data.csv";
    readCSV(filename);
    return 0;
}
