#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "Polyhedron.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedronLibrary;
using namespace Gedim;

//g++ -std=c++17 main.cpp src/Utils.cpp src/UCDUtilities.cpp -I./src -I/usr/include/eigen3 -o geodesic

void EsportaUCD(const PolyhedronMesh& mesh, const std::string& outputDir)
{
    UCDUtilities utilities;

    const auto& coords = mesh.Cell0DsCoordinates;
    const auto& edges = mesh.Cell1DsExtrema;

    std::vector<double> markers0D(mesh.NumCell0Ds, 0.0);
    std::vector<double> markers1D(mesh.NumCell1Ds, 0.0);

    // === Proprietà per punti ===
    std::vector<UCDProperty<double>> cell0Ds_properties(1);
    cell0Ds_properties[0].Label = "Marker";
    cell0Ds_properties[0].UnitLabel = "-";
    cell0Ds_properties[0].NumComponents = 1;
    cell0Ds_properties[0].Data = markers0D.data();

    // === Proprietà per segmenti ===
    std::vector<UCDProperty<double>> cell1Ds_properties(1);
    cell1Ds_properties[0].Label = "Marker";
    cell1Ds_properties[0].UnitLabel = "-";
    cell1Ds_properties[0].NumComponents = 1;
    cell1Ds_properties[0].Data = markers1D.data();

    // === Esportazione ===
    utilities.ExportPoints(outputDir + "/Cell0Ds.inp", coords, cell0Ds_properties);
    utilities.ExportSegments(outputDir + "/Cell1Ds.inp", coords, edges, {}, cell1Ds_properties);

    std::cout << "Esportazione completata: Cell0Ds.inp, Cell1Ds.inp\n";
}

int main()
{
    int p, q, b, c;
    cout << "Inserisci il simbolo di Schläfli {p,q} (es. 3 4 per Ottaedro): ";
    cin >> p >> q;

    cout << "Inserisci i parametri di triangolazione b e c (classe I: b>0, c=0 oppure viceversa): ";
    cin >> b >> c;

    string name;
    if (p == 3 && q == 3) name = "Tetraedro";
    else if (p == 3 && q == 4) name = "Ottaedro";
    else if (p == 3 && q == 5) name = "Icosaedro";
    else {
        cerr << "Errore: solo i solidi platonici {3,3}, {3,4}, {3,5} sono supportati per ora.\n";
        return 1;
    }

    string path = "SolidiPlatonici/" + name;
    PolyhedronMesh base_mesh;
    if (!ImportPolyhedronMesh(base_mesh, path)) {
        cerr << "Errore: file non trovato in " << path << endl;
        return 1;
    }

    PolyhedronMesh geodesic_mesh;
    bool success = false;

    if ((b == 0 && c > 0) || (c == 0 && b > 0)) {
        int param = std::max(b, c);  // usare il parametro diverso da zero
        success = TriangolaClasseI(p, q, param, base_mesh, geodesic_mesh, 0);
    } else if (b == c && b >= 1) {
        success = TriangolaClasse2(p, q, b, base_mesh, geodesic_mesh, c);
    } else {
        cerr << "Errore: valori di b e c non validi. Permessi: (b > 0, c = 0), (b = 0, c > 0), (b = c >= 1).\n";
        return 1;
    }


    if (!success) {
        cerr << "Errore nella triangolazione classe I.\n";
        return 1;
    }


    string outputDir = "OutputGeodetico/" + name + "_b=" + to_string(b) + "_c=" + to_string(c);
    if (!EsportaMeshSuFile(geodesic_mesh, outputDir)) {
        cerr << "Errore nell'esportazione dei file di output.\n";
        return 1;
    }

    EsportaUCD(geodesic_mesh, outputDir);

    cout << "Triangolazione completata. File esportati in: " << outputDir << endl;
    
  /*  PolyhedronMesh dual_mesh;
    if (!CostruisciDualGeodetico(geodesic_mesh, dual_mesh)) {
        cerr << "Errore nella costruzione del duale.\n";
        return 1;
    }

    string dualOutputDir = outputDir + "_DUALE";
    if (!EsportaMeshSuFile(dual_mesh, dualOutputDir)) {
        cerr << "Errore nell'esportazione del duale.\n";
        return 1;
    }

    EsportaUCD(dual_mesh, dualOutputDir);
    cout << "Costruzione del duale completata. File esportati in: " << dualOutputDir << endl;
    return 0;*/
}
