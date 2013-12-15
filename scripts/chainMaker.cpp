/*
 * This is a small utility program to generate peptide chains. Usage is:
 * chainMaker -s <1-letter amino-acid sequence>
 *
 * Optional arguments are as follows (defaults in [] ):
 * 		-chain: 1 letter chain code [A]
 * 		-phi: phi angle (-57.8)
 * 		-psi: psi angle (-47)
 * 		-o: output pdb [chain.pdb]
 * 		--c-terminus: 0 for OXT,HO, 1 for OXT [0]
 * 		--n-terminus: 0 for NH2, 1 for NH3+ [0]
 * 		--r-stereo: use r-stereo [l-stereo]
 *
 * This code to make the chain was written by Geoff Hutchison for the Avogadro project
 * http://avogadro.openmolecules.net/
 *
 * That code was released under the GPL version 2, so this code is under the same.
 *
 * The code uses the openbabel library (http://openbabel.org/).
 *
 * To build it, the following CmakeLists.txt file was used:
 * cmake_minimum_required(VERSION 2.6)
 * add_executable(chainMaker chainMaker.cpp)
 * #target_link_libraries(chainMaker openbabel)
 * #target_link_libraries(chainMaker /home/jmht/Documents/avogadro/openbabel/build/lib/libopenbabel.so)
 * target_link_libraries(chainMaker /Users/jmht/Documents/avogadro/openbabel/OBINSTALL/lib/libopenbabel.dylib )
 * include_directories( /Users/jmht/Documents/avogadro/openbabel/OBINSTALL/include/openbabel-2.0 )
 *
*/

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/residue.h>
#include <openbabel/atom.h>
using namespace std;

class ChainMaker {

public:

	//Constructor
	ChainMaker() {};
	//Destructor
	~ChainMaker() {};

	void AddResidue(string residue, bool lStereo,
			OpenBabel::OBMol &mol, vector<OpenBabel::OBInternalCoord*> &vic,
			const char chain);

	void AddTerminus(int element, string atomID,
			int a, double distance,
			int b, double angle,
			int c, double dihedral,
			OpenBabel::OBMol &mol, vector<OpenBabel::OBInternalCoord*> &vic);

	void generateFragment( OpenBabel::OBMol &obfragment,
			std::string &sequence,
			int phi, int psi, bool lStereo,
			int nTerminus, int cTerminus, char chain
			);

	void setZmatDir( const char* directory ) { zmatDir=directory; };

private:

	// Path to the directory with the zmatrix fragments of the amino acids
	string zmatDir;

};


void ChainMaker::AddResidue(string residue, bool lStereo,
		OpenBabel::OBMol &mol, vector<OpenBabel::OBInternalCoord*> &vic,
		const char chain)
{
	string filename = zmatDir;

	if (residue != "gly") {
		if (lStereo)
			filename += "l-";
		else // D stereo
			filename += "d-";
	}
	filename += residue + ".zmat";

	ifstream ifs;
	ifs.open( filename.c_str() );

	if (!ifs) { // file doesn't exist
		cerr << " Cannot open residue file: " << filename;
		return;
	}

	// Offset:
	//  When we add the internal coordinates, we have to increment
	//  based on the size of the molecule so far
	unsigned int offset = mol.NumAtoms();

	// setup the parent residue
	int prevRes = mol.NumResidues() + 1;
	OpenBabel::OBResidue *res = mol.NewResidue();
	res->SetNum(prevRes);
	res->SetChain(chain);
	// needs to be in uppercase
	std::transform(residue.begin(), residue.end(),residue.begin(), ::toupper);
	res->SetName( residue );

	// Read in an amino z-matrix
	// similar to MOPAC zmat format
	char buffer[BUFF_SIZE];
	vector<string> vs;
	OpenBabel::OBAtom *atom;

	while (ifs.getline(buffer, BUFF_SIZE)) {
		OpenBabel::tokenize(vs, buffer);

		atom = mol.NewAtom();
		atom->SetAtomicNum(OpenBabel::etab.GetAtomicNum(vs[0].c_str()));
		atom->SetPartialCharge(atof(vs[7].c_str()));
		res->InsertAtom(atom);
		res->SetHetAtom(atom, false);
		res->SetSerialNum(atom, mol.NumAtoms());
		if (vs.size() == 9)
			res->SetAtomID(atom, vs[8]);

		OpenBabel::OBInternalCoord *coord = new OpenBabel::OBInternalCoord;
		coord->_dst = atof(vs[1].c_str());
		coord->_ang = atof(vs[2].c_str());
		coord->_tor = atof(vs[3].c_str());

		unsigned int index;
		// Set _a
		index = atoi(vs[4].c_str());
		if (index > 0 && index <= mol.NumAtoms())
			coord->_a = mol.GetAtom(index + offset);
		else
			coord->_a = NULL;
		// Set _b
		index = atoi(vs[5].c_str());
		if (index > 0 && index <= mol.NumAtoms())
			coord->_b = mol.GetAtom(index + offset);
		else
			coord->_b = NULL;
		// Set _c
		index = atoi(vs[6].c_str());
		if (index > 0 && index <= mol.NumAtoms())
			coord->_c = mol.GetAtom(index + offset);
		else
			coord->_c = NULL;

		vic.push_back(coord);
	}
}

void ChainMaker::AddTerminus(int element, string atomID,
		int a, double distance,
		int b, double angle,
		int c, double dihedral,
		OpenBabel::OBMol &mol, vector<OpenBabel::OBInternalCoord*> &vic)
{
	OpenBabel::OBResidue *res = mol.GetResidue(mol.NumResidues() - 1);
	if (!res || mol.NumResidues() == 0)
		return; // can't do anything -- we're in a weird state

	OpenBabel::OBAtom *atom;

	atom = mol.NewAtom();
	atom->SetAtomicNum(element);
	res->InsertAtom(atom);
	res->SetHetAtom(atom, false);
	res->SetSerialNum(atom, mol.NumAtoms());
	//res->SetAtomID(atom, atomID.toAscii().data());
	res->SetAtomID(atom, atomID.data());

	OpenBabel::OBInternalCoord *coord = new OpenBabel::OBInternalCoord;
	coord->_dst = distance;
	coord->_ang = angle;
	coord->_tor = dihedral;

	coord->_a = mol.GetAtom(a);
	coord->_b = mol.GetAtom(b);
	coord->_c = mol.GetAtom(c);

	// Add a bond between the recently created atom and our "a"
	mol.AddBond(mol.NumAtoms(), a, 1);

	vic.push_back(coord);
}

void ChainMaker::generateFragment( OpenBabel::OBMol &obfragment,
		std::string &sequence,
		int phi, int psi, bool lStereo,
		int nTerminus, int cTerminus, char chain
		)
{

	if (sequence.empty())
		return; // also nothing to do

	//OBMol obfragment;
	vector<OpenBabel::OBInternalCoord*> vic;
	vic.push_back((OpenBabel::OBInternalCoord*)NULL);
	OpenBabel::OBInternalCoord* ic;
	int lastN, lastCa, lastCac, lastO; // backbone atoms
	lastN = lastCa = lastCac = lastO = 0;
	int newN, newCa, newCac, newO;
	int lastAtom = 0; // last atom read from previous residue

	double omega=179.99;
	double amideLength = 1.34;
	double bondAngle = 120.0;

	// Now the work begins
	// Get the sequence (in lower case)
	obfragment.BeginModify();

	// Split the string using std library functions
	const char delim = '-';

	//http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
	std::vector<std::string> elems;
	std::stringstream ss(sequence);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}

	for(std::vector< std::string >::iterator it = elems.begin(); it != elems.end(); ++it) {
		AddResidue(*it, lStereo, obfragment, vic, chain);
		if (!obfragment.NumAtoms()) {
			// Residue was not added - bail
			cerr << "Problem adding new residues - file not read.";
			return;
		}

		newN = lastAtom + 1;
		newCa = lastAtom + 2;
		newCac = lastAtom + 3;
		newO = lastAtom + 4;

		if (lastAtom != 0) {
			// set the peptide bond to the previous residue
			// first the nitrogen
			ic = vic[newN];
			ic->_a = obfragment.GetAtom(lastCac);
			ic->_dst = amideLength;
			ic->_b = obfragment.GetAtom(lastCa);
			ic->_ang = bondAngle;
			ic->_c = obfragment.GetAtom(lastN);
			ic->_tor = psi;

			// fix the O=C from previous residue
			ic = vic[lastO];
			ic->_tor = 180.0 + psi;

			// now the Calpha
			ic = vic[newCa];
			ic->_b = obfragment.GetAtom(lastCac);
			ic->_ang = bondAngle;
			ic->_c = obfragment.GetAtom(lastCa);
			ic->_tor = omega;

			// now the new C=O
			ic = vic[newCac];
			ic->_c = obfragment.GetAtom(lastCac);
			ic->_tor = phi;

			// add the peptide bond
			obfragment.AddBond(lastCac, newN, 1);
		}
		else { // The first residue
			// Add the N-terminus modification
			switch ( nTerminus ) {
			case 0: // NH2
				AddTerminus(1, "HN", newN, 1.016, newCa, 120.0,
						newCac, 175.0, obfragment, vic);
				break;
			case 1: // NH3+
				AddTerminus(1, "HN", newN, 1.016, newCa, 109.5,
						newCac, 117.0, obfragment, vic);
				AddTerminus(1, "2HN", newN, 1.016, newCa, 109.5,
						newCac, -117.0, obfragment, vic);
				break;
			default:
				break;
			}
		}

		// add the known backbone bonds
		obfragment.AddBond(newN, newCa, 1);
		obfragment.AddBond(newCa, newCac, 1);
		obfragment.AddBond(newCac, newO, 2); // C=O

		lastN = newN;
		lastCa = newCa;
		lastCac = newCac;
		lastO = newO;
		lastAtom = obfragment.NumAtoms();
	}
	// Fix the final C=O if not straight-chain
	ic = vic[lastO];
	ic->_tor = 180.0 + psi;

	// Add the C-terminus end group
	switch ( cTerminus ) {
	case 0: // CO2H
		AddTerminus(8, "OXT", lastCac, 1.351, lastO, 120.0,
				lastCa, -180.0, obfragment, vic);
		obfragment.AddBond(obfragment.NumAtoms(), lastCac, 1);
		AddTerminus(1, "HO", obfragment.NumAtoms(), 1.064,
				lastCac, 120.0, lastO, 180.0, obfragment, vic);
		break;
	case 1: // CO2-
		AddTerminus(8, "OXT", lastCac, 1.351, lastO, 120.0,
				lastCa, -180.0, obfragment, vic);
		break;
	default:
		break;
	}

	obfragment.EndModify();
	if (obfragment.NumAtoms()) {
		// Don't do all this work, if there's nothing to do
		InternalToCartesian(vic,obfragment);
		OpenBabel::OBBitVec allAtoms;
		allAtoms.SetRangeOn(0, obfragment.NumAtoms());
		allAtoms.SetBitOff(obfragment.NumAtoms() - 1); // Don't add bonds for the terminus
		//jmht FIX!!!
		OpenBabel::resdat.AssignBonds(obfragment, allAtoms);

		// some of the fragments still miss bonds
		obfragment.ConnectTheDots();

		obfragment.SetPartialChargesPerceived();
	}
}

int main(int argc,char **argv)
{
	// Maps 1 AA code -> 3
	map <char, string> one2three;
	one2three['A'] = "ala";
	one2three['C'] = "cys";
	one2three['D'] = "asp";
	one2three['E'] = "glu";
	one2three['F'] = "phe";
	one2three['G'] = "gly";
	one2three['H'] = "his";
	one2three['I'] = "ile";
	one2three['K'] = "lys";
	one2three['L'] = "leu";
	one2three['M'] = "met";
	one2three['N'] = "asn";
	one2three['P'] = "pro";
	one2three['Q'] = "gln";
	one2three['R'] = "arg";
	one2three['S'] = "ser";
	one2three['T'] = "thr";
	one2three['V'] = "val";
	one2three['W'] = "trp";
	one2three['Y'] = "tyr";

	/* Default values */
	double phi=-57.8;
	double psi=-47;
	bool lStereo=true;
	int nTerminus=0; // 0=NH2, 1=NH3+
	int cTerminus=0; // 0=CO2H, 1=CO2-
	char chain='A';
	string pdbout = "chain.pdb";
	string sequence1 = "";

	/* Read in the command-line arguments */
	vector<string> args(argv + 1, argv + argc);
	for ( unsigned int i=0; i < args.size(); i++) {
		if (args[i] == "-h" || args[i] == "--help") {
			cout << "Usage: prog -s ASDSDS -phi -57.8 -psi -47 -o pdbout.pdb "<< endl;
			return 0;
		} else if (args[i] == "-chain") {
			chain = args[i+1].c_str()[0];
		} else if (args[i] == "-o") {
			pdbout = args[i+1];
		} else if (args[i] == "-phi") {
			phi = atoi(args[i+1].c_str());
		} else if (args[i] == "-psi") {
			psi = atoi(args[i+1].c_str());
		} else if (args[i] == "-s") {
			sequence1 = args[i+1];
		} else if (args[i] == "--c-terminus") {
			cTerminus = atoi(args[i+1].c_str());
		} else if (args[i] == "--n-terminus") {
			nTerminus = atoi(args[i+1].c_str());
		} else if (args[i] == "--r-stereo") {
			lStereo=false;
		}
	}

	assert( sequence1.size() > 0);
	cout << "Processing with phi: " << phi << " psi: " << psi << " sequence: " << sequence1 << endl;

	// Convert single AA sequence to triplets separted by -
	stringstream s3;
	char delim='-';
	for ( unsigned int i=0; i < sequence1.size(); i++ ) {
		if ( i == 0 ) {
			s3 << one2three.at( sequence1[i] );
		} else{
			s3 << delim << one2three.at( sequence1[i] );
		}
	}

	// Now create the chain
	string sequence3 = s3.str();
	OpenBabel::OBMol obfragment;

	ChainMaker chainMaker;
	//generateFragment( obfragment, sequence3, phi, psi, lStereo, nTerminus, cTerminus, chain );
	//string zmatDir="/Users/jmht/Documents/avogadro/avogadro/builder/amino/";
	chainMaker.setZmatDir( "/Users/jmht/Documents/avogadro/avogadro/builder/amino/" );
	chainMaker.generateFragment( obfragment, sequence3, phi, psi, lStereo, nTerminus, cTerminus, chain );

	// Write it out as a pdb
	OpenBabel::OBConversion conv;
	conv.SetOutFormat("PDB");
	conv.WriteFile(&obfragment,pdbout);

}
