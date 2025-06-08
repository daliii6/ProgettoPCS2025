#pragma once

#include <iostream>
#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include "Utils.hpp"

#include <Eigen/Dense> // necessario se non incluso in Utils.hpp

using namespace PolyhedronLibrary;

const int b = 4;

// Calcola le proprietà teoriche attese di un solido geodetico di classe I
// p, q: simbolo di Schläfli del solido platonico (usa solo q)
// b: parametro di triangolazione (b ≥ 1)
std::array<unsigned int, 3> ProprietàSolidiI(int q, int b) {
    unsigned int T = b * b;

    unsigned int Exp_V = 0;
    unsigned int Exp_E = 0;
    unsigned int Exp_F = 0;

    if (q == 3) {
        Exp_V = 2 * T + 2;
        Exp_E = 6 * T;
        Exp_F = 4 * T;
    } else if (q == 4) {
        Exp_V = 4 * T + 2;
        Exp_E = 12 * T;
        Exp_F = 8 * T;
    } else if (q == 5) {
        Exp_V = 10 * T + 2;
        Exp_E = 30 * T;
        Exp_F = 20 * T;
    } else {
        std::cerr << "[ERRORE] q deve essere 3, 4 o 5 per i solidi platonici." << std::endl;
    }

    return {Exp_V, Exp_E, Exp_F};
}


// Calcola le proprietà teoriche attese di un solido geodetico di classe II
// p, q: simbolo di Schläfli del solido platonico di partenza
// b: parametro di triangolazione (b ≥ 1)
std::array<unsigned int, 3> ProprietàSolidiII(int q, int b) {
    unsigned int V = 0, E = 0, F = 0;

    if (q == 3) {
        V = 4; E = 6; F = 4;      // Tetraedro {3,3}
    } else if (q == 4) {
        V = 6; E = 12; F = 8;     // Ottaedro {3,4}
    } else if (q == 5) {
        V = 12; E = 30; F = 20;   // Icosaedro {3,5}
    } else {
        std::cerr << "[ERRORE] q deve essere 3, 4 o 5 per i solidi platonici." << std::endl;
        return {0, 0, 0};
    }

    unsigned int Exp_V = V + E * (2 * b - 1) + F * ((3 * b * b) / 2 - (3 * b / 2) + 1);
    unsigned int Exp_E = E * 2 * b + F * ((9 * b * b) / 2 + (3 * b / 2));
    unsigned int Exp_F = F * (3 * b * b + 3 * b);

    return {Exp_V, Exp_E, Exp_F};
}


// Calcola quanti vertici hanno un certo grado (numero di facce adiacenti)
// ExpectedDegree: grado atteso del vertice
// Vertices: ID dei vertici
// Faces: elenco delle facce come liste di ID vertici
int VertexDegree(int ExpectedDegree, const std::vector<unsigned int>& Vertices, const std::vector<std::vector<unsigned int>>& Faces){
    int NumVertexOfExpectedDegree = 0;
    for (const auto& idV : Vertices) {
        int CurrentDegree = 0;
        for (const auto& listF : Faces)
            if (std::find(listF.begin(), listF.end(), idV) != listF.end())
                CurrentDegree++;
        if (CurrentDegree == ExpectedDegree)
            NumVertexOfExpectedDegree++;
    }
    return NumVertexOfExpectedDegree;
}

TEST(TestPolyhedron, TestClasseITetraedro)
{
    PolyhedronMesh PlatonicPolyhedron;
    if (!ImportPolyhedronMesh(PlatonicPolyhedron, "../SolidiPlatonici/Tetraedro/"))
        FAIL() << "[ERRORE] File non trovato";

    PolyhedronMesh polyhedron;
    TriangolaClasseI(b, PlatonicPolyhedron, polyhedron);
    std::array<unsigned int, 3> expected = ProprietàSolidiI(3, b); // Tetraedro: {p=3, q=3}
    std::array<unsigned int, 3> actual = {
        polyhedron.NumCell0Ds,
        polyhedron.NumCell1Ds,
        polyhedron.NumCell2Ds
    };
    EXPECT_EQ(expected, actual);

    unsigned int T = b * b;
    EXPECT_EQ(VertexDegree(3, polyhedron.Cell0DsId, polyhedron.Cell2DsVertices), 4);
    EXPECT_EQ(VertexDegree(6, polyhedron.Cell0DsId, polyhedron.Cell2DsVertices), 2 * (T - 1));
}

TEST(TestPolyhedron, TestClasseIOttaedro)
{
    PolyhedronMesh PlatonicPolyhedron;
    if (!ImportPolyhedronMesh(PlatonicPolyhedron, "../SolidiPlatonici/Ottaedro/"))
        FAIL() << "[ERRORE] File non trovato";

    PolyhedronMesh polyhedron;
    TriangolaClasseI(b, PlatonicPolyhedron, polyhedron);
    std::array<unsigned int, 3> expected = ProprietàSolidiI(4, b); // Ottaedro: {p=3, q=4}
    std::array<unsigned int, 3> actual = {
        polyhedron.NumCell0Ds,
        polyhedron.NumCell1Ds,
        polyhedron.NumCell2Ds
    };
    EXPECT_EQ(expected, actual);

    unsigned int T = b * b;
    EXPECT_EQ(VertexDegree(4, polyhedron.Cell0DsId, polyhedron.Cell2DsVertices), 6);
    EXPECT_EQ(VertexDegree(6, polyhedron.Cell0DsId, polyhedron.Cell2DsVertices), 4 * (T - 1));
}

TEST(TestPolyhedron, TestClasseIITetraedro)
{
    PolyhedronMesh PlatonicPolyhedron;
    if (!ImportPolyhedronMesh(PlatonicPolyhedron, "../SolidiPlatonici/Tetraedro/"))
        FAIL() << "[ERRORE] File non trovato";

    PolyhedronMesh polyhedron;
    TriangolaClasseII(b, PlatonicPolyhedron, polyhedron);

    std::array<unsigned int, 3> expected = ProprietàSolidiII(3, b); // Tetraedro: {p=3, q=3}
    std::array<unsigned int, 3> actual = {
        polyhedron.NumCell0Ds,
        polyhedron.NumCell1Ds,
        polyhedron.NumCell2Ds
    };
    EXPECT_EQ(expected, actual);
}

TEST(TestPolyhedron, TestClasseIIOttaedro)
{
    PolyhedronMesh PlatonicPolyhedron;
    if (!ImportPolyhedronMesh(PlatonicPolyhedron, "../SolidiPlatonici/Ottaedro/"))
        FAIL() << "[ERRORE] File non trovato";

    PolyhedronMesh polyhedron;
    TriangolaClasseII(b, PlatonicPolyhedron, polyhedron);

    std::array<unsigned int, 3> expected = ProprietàSolidiII(4, b); // Ottaedro: {p=3, q=4}
    std::array<unsigned int, 3> actual = {
        polyhedron.NumCell0Ds,
        polyhedron.NumCell1Ds,
        polyhedron.NumCell2Ds
    };
    EXPECT_EQ(expected, actual);
}

TEST(TestPolyhedron, TestClasseIIIcosaedro)
{
    PolyhedronMesh PlatonicPolyhedron;
    if (!ImportPolyhedronMesh(PlatonicPolyhedron, "../SolidiPlatonici/Icosaedro/"))
        FAIL() << "[ERRORE] File non trovato";

    PolyhedronMesh polyhedron;
    TriangolaClasseII(b, PlatonicPolyhedron, polyhedron);

    std::array<unsigned int, 3> expected = ProprietàSolidiII(5, b); // Icosaedro: {p=3, q=5}
    std::array<unsigned int, 3> actual = {
        polyhedron.NumCell0Ds,
        polyhedron.NumCell1Ds,
        polyhedron.NumCell2Ds
    };
    EXPECT_EQ(expected, actual);
}


TEST(TestDualPolyhedron, TestClasseIDualeTetraedro)
{
    PolyhedronMesh PlatonicPolyhedron;
    if (!ImportPolyhedronMesh(PlatonicPolyhedron, "../SolidiPlatonici/Tetraedro/"))
        FAIL() << "[ERRORE] File non trovato";

    PolyhedronMesh polyhedron;
    TriangolaClasseI(b, PlatonicPolyhedron, polyhedron);
    PolyhedronMesh DualPolyhedron;
    CostruisciDualMesh(polyhedron, DualPolyhedron);

    int T = b * b;

    int ExpectedVertices = 4 * T;      // Una faccia della mesh originale diventa un vertice nel duale
    int ExpectedFaces = 2 * T + 2;     // Ogni vertice della mesh originale diventa una faccia nel duale

    EXPECT_EQ(DualPolyhedron.NumCell0Ds, ExpectedVertices);
    EXPECT_EQ(DualPolyhedron.NumCell2Ds, ExpectedFaces);
}


TEST(TestDualPolyhedron, TestClasseIIDualeTetraedro)
{
    PolyhedronMesh PlatonicPolyhedron;
    if (!ImportPolyhedronMesh(PlatonicPolyhedron, "../SolidiPlatonici/Tetraedro/"))
        FAIL() << "[ERRORE] File non trovato";

    PolyhedronMesh polyhedron;
    TriangolaClasseII(b, PlatonicPolyhedron, polyhedron);

    PolyhedronMesh DualPolyhedron;
    CostruisciDualMesh(polyhedron, DualPolyhedron);

    // Calcola proprietà della mesh di classe II (non del duale)
    std::array<unsigned int, 3> expected = ProprietàSolidiII(3, b); // Tetraedro: {p=3, q=3}

    int ExpectedVertices = expected[2]; // Le facce originali diventano vertici del duale
    int ExpectedFaces = expected[0];    // I vertici originali diventano facce del duale

    EXPECT_EQ(DualPolyhedron.NumCell0Ds, ExpectedVertices);
    EXPECT_EQ(DualPolyhedron.NumCell2Ds, ExpectedFaces);
}

TEST(TestOrdinaFacce, OrdinaFacceAttornoVertice_Ottaedro)
{
    using namespace PolyhedronLibrary;

    PolyhedronMesh PlatonicPolyhedron;
    ASSERT_TRUE(ImportPolyhedronMesh(PlatonicPolyhedron, "../SolidiPlatonici/Ottaedro/"));

    PolyhedronMesh polyhedron;
    ASSERT_TRUE(TriangolaClasseI(1, PlatonicPolyhedron, polyhedron));

    int vertex_id = 0; // scegli un ID valido, es: vertice al polo nord

    // Ottiene le facce che contengono il vertice e le ordina ciclicamente
    std::vector<int> facce_ordinate = OrdinaFacceAttornoAlVertice(vertex_id, polyhedron);

    // Verifica che le facce ordinate contengano effettivamente il vertice
    for (int fid : facce_ordinate) {
        const auto& face = polyhedron.Cell2DsVertices[fid];
        EXPECT_NE(std::find(face.begin(), face.end(), vertex_id), face.end())
            << "La faccia " << fid << " non contiene il vertice " << vertex_id;
    }

    // Verifica che ogni faccia sia adiacente alla successiva (2 vertici in comune)
    for (size_t i = 0; i + 1 < facce_ordinate.size(); ++i) {
        const auto& f1 = polyhedron.Cell2DsVertices[facce_ordinate[i]];
        const auto& f2 = polyhedron.Cell2DsVertices[facce_ordinate[i+1]];
        int count_comuni = 0;
        for (int v : f1)
            if (std::find(f2.begin(), f2.end(), v) != f2.end())
                count_comuni++;
        EXPECT_EQ(count_comuni, 2) << "Le facce " << facce_ordinate[i]
                                   << " e " << facce_ordinate[i+1]
                                   << " non sono adiacenti attorno al vertice " << vertex_id;
    }
}

TEST(TestShortestPath, CamminoMinimoClasseI)
{
    using namespace PolyhedronLibrary;

    PolyhedronMesh PlatonicPolyhedron;
    if (!ImportPolyhedronMesh(PlatonicPolyhedron, "../SolidiPlatonici/Ottaedro/"))
        FAIL() << "[ERRORE] File non trovato";

    PolyhedronMesh polyhedron;
    ASSERT_TRUE(TriangolaClasseI(2, PlatonicPolyhedron, polyhedron));

    auto LA = ListaAdiacenza(polyhedron);

    double lunghezzaTotale;
    std::vector<int> predecessori;
    std::vector<unsigned int> cammino = DijkstraCamminoMinimo(LA, 2, 8, lunghezzaTotale, predecessori);

    std::vector<unsigned int> expected_path = {2, 5, 8};

    EXPECT_EQ(cammino, expected_path);
    EXPECT_EQ(cammino.size() - 1, 2); // due spigoli
    EXPECT_NEAR(lunghezzaTotale, 1.531, 1e-3);
}

TEST(TestShortestPath, CamminoMinimoClasseII)
{
    using namespace PolyhedronLibrary;

    PolyhedronMesh PlatonicPolyhedron;
    if (!ImportPolyhedronMesh(PlatonicPolyhedron, "../SolidiPlatonici/Tetraedro/"))
        FAIL() << "[ERRORE] File non trovato";

    PolyhedronMesh polyhedron;
    ASSERT_TRUE(TriangolaClasseII(2, PlatonicPolyhedron, polyhedron));

    auto LA = ListaAdiacenza(polyhedron);

    double lunghezzaTotale;
    std::vector<int> predecessori;
    std::vector<unsigned int> cammino = DijkstraCamminoMinimo(LA, 4, 7, lunghezzaTotale, predecessori);

    // Qui metti un path valido secondo il grafo
    std::vector<unsigned int> expected_path = {4, 6, 7};

    EXPECT_EQ(cammino, expected_path);
    EXPECT_EQ(cammino.size() - 1, 2); // due spigoli
    EXPECT_NEAR(lunghezzaTotale, 1.45, 1e-2);
}
