#include <algorithm>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <cmath>

using namespace std;

int main(int argc, char** argv)
{
    const double eps = 1e-18;
	if (argc != 3) {
		cerr << "Usage: compare <file1> <file2>" << endl;
        return 1;
	}

	ifstream file1( argv[1], ios::in | ios::binary );
	ifstream file2( argv[2], ios::in | ios::binary );

    assert( file1.good() );
    assert( file2.good() );

	int nErrs = 0;
    double norm = 0;
	while( file1.good() && file2.good() ) {
        double x1, x2;
        file1.read( reinterpret_cast<char*>( &x1 ), sizeof x1 );
        file2.read( reinterpret_cast<char*>( &x2 ), sizeof x2 );

        if( !file1.good() || !file1.good() ) {
            break;
        }
		if( abs( x1 - x2 ) / max( abs( x1 ), abs( x2 ) ) > eps ) {
			++nErrs;
			cout << x1 << " " << x2 << endl;
		}
        norm = max( norm, abs( x1 - x2 ) );
	}
    cerr << "norm = " << norm << endl;

    if( file1.good() || file2.good() ) {
        cerr << "Files are of different sizes" << endl;
        return 2;
    }

	if( nErrs == 0 ) {
        cerr << "Ok" << endl;
        return 0;
    } else {
        cerr << nErrs << " errors. failed" << endl;
        return 3;
    }
}
