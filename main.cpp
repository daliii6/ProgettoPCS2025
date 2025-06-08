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

    std::vector<UCDProperty<double>> no_properties; // Vuoto

    utilities.ExportPoints(outputDir + "/Cell0Ds.inp", coords, no_properties);
    utilities.ExportSegments(outputDir + "/Cell1Ds.inp", coords, edges, {}, no_properties);

    std::cout << "Esportazione completata.\n";
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

    string path = "../SolidiPlatonici/" + name;
    PolyhedronMesh base_mesh;
    if (!ImportPolyhedronMesh(base_mesh, path)) {
        cerr << "Errore: file non trovato in " << path << endl;
        return 1;
    }

    PolyhedronMesh geodesic_mesh;
    bool success = false;
    if ((b == 0 && c > 0) || (c == 0 && b > 0)) {
        int param = std::max(b, c);  // usare il parametro diverso da zero
        success = TriangolaClasseI(param, base_mesh, geodesic_mesh);
    } else if (b == c && b >= 1) {
        success = TriangolaClasseII(b, base_mesh, geodesic_mesh);
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

    unsigned int id_start, id_end;
    unsigned int n = geodesic_mesh.NumCell0Ds;
    cout << "id compreso tra 0 e " << n - 1 << endl;
    cout << "Inserisci id start e id fine: ";
    cin >> id_start >> id_end;

    if (id_start < n && id_end < n)
    {
        vector<list<pair<unsigned int, double>>> LA_geo = ListaAdiacenza(geodesic_mesh);
        vector<int> predecessori;
        double lunghezzaTotale;
        vector<unsigned int> cammino = DijkstraCamminoMinimo(LA_geo, id_start, id_end, lunghezzaTotale, predecessori);
        vector<unsigned int> lati = latiCamminoMinimo(geodesic_mesh, cammino, outputDir);

        cout << "Lunghezza totale cammino: " << lunghezzaTotale << "\n";
        cout << "Numero lati: " << lati.size() << "\n";
    }
    else
    {
        cout << "Errore: vertici fuori dal range. File UCD neutri verranno generati.\n";
        EsportaUCD(geodesic_mesh, outputDir);
    }

    // Calcolo e esportazione del duale
    PolyhedronMesh dual_mesh;
    if (!CostruisciDualMesh(geodesic_mesh, dual_mesh)) {
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
    return 0;

}