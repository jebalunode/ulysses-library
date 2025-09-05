#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>

using namespace std;

std::vector<double> gfn2_xtb_opt(string & molecule, int charge, double Telec );

int main( int argc, char ** argv){
 string water=R"(3

O      0.0000000     0.0000000    -0.0559861    
H      0.7744122     0.0000000     0.5156460
H     -0.7744122     0.0000000     0.5156460)";



cout << water << endl;



std::vector<double> data =  gfn2_xtb_opt(water, 0, 300.0 );
cout << "Charges\n";

for (auto item : data)
	cout << item<< endl;


if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <filename>\n";
    return 1;
}

std::ifstream file(argv[1]);
if (!file) {
    std::cerr << "Error opening file: " << argv[1] << "\n";
    return 1;
}

std::stringstream buffer;
buffer << file.rdbuf(); // Read entire file into buffer

std::string fileContent = buffer.str(); // Convert buffer to string


data =  gfn2_xtb_opt(fileContent, 0, 300.0 );
cout << "Charges\n";

for (auto item : data)
        cout << item<< endl;




return 0;


}
