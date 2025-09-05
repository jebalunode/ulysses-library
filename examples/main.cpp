#include <iostream>
#include <vector>
#include <string>

using namespace std;

std::vector<double> gfn2_xtb_spe(string & molecule, int charge, double Telec );

int main( int argc, char ** argv){
 string water=R"(3

O      0.0000000     0.0000000    -0.0559861    
H      0.7744122     0.0000000     0.5156460
H     -0.7744122     0.0000000     0.5156460)";



cout << water << endl;



std::vector<double> data =  gfn2_xtb_spe(water, 0, 300.0 );
cout << "Charges\n";

for (auto item : data)
	cout << item<< endl;
 return 0;
}
