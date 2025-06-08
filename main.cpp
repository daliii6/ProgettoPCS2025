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

void EsportaUCD(const PolyhedronMesh& mesh, const std::string& outputDir)
{
    UCDUtilities utilities;

    const auto& coords = mesh.Cell0DsCoordinates;
    const auto& edges = mesh.Cell1DsExtrema;

    std::vector<UCDProperty<double>> no_properties;

    utilities.ExportPoints(outputDir + "/Cell0Ds.inp", coords, no_properties);
    utilities.ExportSegments(outputDir + "/Cell1Ds.inp", coords, edges, {}, no_properties);

    std::cout << "[INFO] File UCD per Paraview esportati: Cell0Ds.inp, Cell1Ds.inp\n";
}

int main()
{
    int p, q, b, c;
    int id_start = -1, id_end = -1;

    cout << "Inserisci i parametri (p, q, b, c, [id_start] [id_end]): ";
    string line;
    getline(cin >> ws, line);
    istringstream iss(line);

    iss >> p >> q >> b >> c;
    if (!(iss >> id_start >> id_end)) {
        id_start = -1;
        id_end = -1;
    }

    string name;
    if (p == 3 && q == 3) name = "Tetraedro";
    else if (p == 3 && q == 4) name = "Ottaedro";
    else if (p == 3 && q == 5) name = "Icosaedro";
    else {
        cerr << "[ERRORE] Supportati solo i solidi platonici {3,3}, {3,4}, {3,5}.\n";
        return 1;
    }

    string path = "../SolidiPlatonici/" + name;
    PolyhedronMesh base_mesh;
    if (!ImportPolyhedronMesh(base_mesh, path)) {
        cerr << "[ERRORE] File non trovato nella cartella: " << path << "\n";
        return 1;
    }

    PolyhedronMesh geodesic_mesh;
    bool success = false;
    if ((b == 0 && c > 0) || (c == 0 && b > 0)) {
        int param = std::max(b, c);
        success = TriangolaClasseI(param, base_mesh, geodesic_mesh);
    } else if (b == c && b >= 1) {
        success = TriangolaClasseII(b, base_mesh, geodesic_mesh);
    } else {
        cerr << "[ERRORE] Parametri non validi. Ammessi solo: (b > 0, c = 0), (b = 0, c > 0), (b = c >= 1).\n";
        return 1;
    }

    if (!success) {
        cerr << "[ERRORE] Triangolazione non riuscita: possibile presenza di facce non triangolari.\n";
        return 1;
    }

    string outputDir = "OutputGeodetico/" + name + "_b=" + to_string(b) + "_c=" + to_string(c);
    if (!EsportaMeshSuFile(geodesic_mesh, outputDir)) {
        cerr << "[ERRORE] Esportazione fallita dei file .txt delle celle.\n";
        return 1;
    }
    cout << "[INFO] File .txt generati in: " << outputDir << "\n";

    int n = geodesic_mesh.NumCell0Ds;
    cout << "[INFO] Vertici numerati da 0 a " << n - 1 << ".\n";

    if (id_start >= 0 && id_start < n && id_end >= 0 && id_end < n)
    {
        vector<list<pair<unsigned int, double>>> LA_geo = ListaAdiacenza(geodesic_mesh);
        vector<int> predecessori;
        double lunghezzaTotale;
        vector<unsigned int> cammino = DijkstraCamminoMinimo(LA_geo, id_start, id_end, lunghezzaTotale, predecessori);
        vector<unsigned int> lati = latiCamminoMinimo(geodesic_mesh, cammino, outputDir);

        cout << "\n[INFO] Cammino minimo tra i vertici " << id_start << " e " << id_end << ":\n";
        for (size_t i = 0; i < lati.size(); ++i) {
            unsigned int eid = lati[i];
            unsigned int v0 = geodesic_mesh.Cell1DsExtrema(0, eid);
            unsigned int v1 = geodesic_mesh.Cell1DsExtrema(1, eid);
            cout << "• Lato " << eid << ": collega vertici " << v0 << " e " << v1 << "\n";
        }
        cout << "Lunghezza totale: " << lunghezzaTotale << "\n";
        cout << "Numero lati nel cammino: " << lati.size() << "\n";
    }
    else
    {
        cout << "[AVVISO] Vertici non validi. Verranno generati solo file UCD neutri.\n";
        EsportaUCD(geodesic_mesh, outputDir);
    }

    PolyhedronMesh dual_mesh;
    if (!CostruisciDualMesh(geodesic_mesh, dual_mesh)) {
        cerr << "[ERRORE] Costruzione del duale fallita.\n";
        return 1;
    }

    string dualOutputDir = outputDir + "_DUALE";
    if (!EsportaMeshSuFile(dual_mesh, dualOutputDir)) {
        cerr << "[ERRORE] Esportazione fallita dei file del duale.\n";
        return 1;
    }

    EsportaUCD(dual_mesh, dualOutputDir);
    cout << "[INFO] Doppio geodetico esportato in: " << dualOutputDir << "\n";
    return 0;
}
