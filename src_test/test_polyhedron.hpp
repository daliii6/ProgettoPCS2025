#pragma once

#include <iostream>
#include <gtest/gtest.h>
#include <algorithm>
#include "utils.hpp" // Assumiamo che qui ci sia PolyhedronMesh, TriangolaClasseI, ecc.

using namespace PolyhedronLibrary;

// Impostiamo il parametro di triangolazione b da usare in tutti i test
const int TRIANGULATION_PARAMETER = 3;  // Si può cambiare in 4, 5, ecc.

// Funzione di utilità per contare quanti vertici hanno una valenza (grado) specifica
// - ExpectedDegree è il grado atteso
// - Vertices è la lista di ID dei vertici
// - Faces è la lista delle facce, ciascuna è una tripletta di ID vertici
int VertexDegree(int ExpectedDegree, const std::vector<int>& Vertices, const std::vector<std::vector<int>>& Faces) {
	int NumVertexOfExpectedDegree = 0;

	for (const auto& idV : Vertices) {
		int CurrentDegree = 0;

		// Conta quante volte il vertice idV compare nelle facce
		for (const auto& listF : Faces)
			if (std::find(listF.begin(), listF.end(), idV) != listF.end())
				CurrentDegree++;

		// Se la valenza del vertice coincide con quella attesa, lo contiamo
		if (CurrentDegree == ExpectedDegree)
			NumVertexOfExpectedDegree++;
	}
	return NumVertexOfExpectedDegree;
}

// Funzione che restituisce il numero atteso di vertici, spigoli e facce
// a seconda del valore di q e del parametro b
std::array<int, 3> ExpectedCounts(int q, int b) {
	int T = b * b; // T = b^2, perché siamo in classe I (c = 0)
	std::array<int, 3> result;

	if (q == 3) result = {2*T + 2, 6*T, 4*T};    // Tetraedro
	if (q == 4) result = {4*T + 2, 12*T, 8*T};   // Ottaedro
	if (q == 5) result = {10*T + 2, 30*T, 20*T}; // Icosaedro

	return result;
}

// Primo test: costruiamo un solido geodetico a simmetria tetraedrica (q = 3)
TEST(TestGeodeticPolyhedronClasseI, Tetrahedron) {
	PolyhedronMesh PlatonicPolyhedron;

	// Importa la mesh del tetraedro da file
	ASSERT_TRUE(FileManagement::ImportPolyhedronMesh(PlatonicPolyhedron, "../SolidiPlatonici/Tetraedro/"));

	PolyhedronMesh GeodeticPolyhedron;

	// Costruisce la triangolazione geodetica di classe I
	ASSERT_TRUE(Generation::TriangolaClasseI(3, 3, TRIANGULATION_PARAMETER, PlatonicPolyhedron, GeodeticPolyhedron, 0));

	// Calcola i valori attesi di vertici, spigoli e facce
	auto [ExpectedV, ExpectedE, ExpectedF] = ExpectedCounts(3, TRIANGULATION_PARAMETER);

	// Verifica che il numero di vertici/spigoli/facce sia corretto
	EXPECT_EQ(GeodeticPolyhedron.NumCell0Ds, ExpectedV);
	EXPECT_EQ(GeodeticPolyhedron.NumCell1Ds, ExpectedE);
	EXPECT_EQ(GeodeticPolyhedron.NumCell2Ds, ExpectedF);

	int T = TRIANGULATION_PARAMETER * TRIANGULATION_PARAMETER;

	// Verifica il numero di vertici con valenza 3 (ai vertici del tetraedro)
	EXPECT_EQ(VertexDegree(3, GeodeticPolyhedron.Cell0DsId, GeodeticPolyhedron.Cell2DsVertices), 4);

	// Verifica il numero di vertici con valenza 6 (quelli interni)
	EXPECT_EQ(VertexDegree(6, GeodeticPolyhedron.Cell0DsId, GeodeticPolyhedron.Cell2DsVertices), 2*(T - 1));
}

// Secondo test: simmetria ottaedrica (q = 4)
TEST(TestGeodeticPolyhedronClasseI, Octahedron) {
	PolyhedronMesh PlatonicPolyhedron;
	ASSERT_TRUE(FileManagement::ImportPolyhedronMesh(PlatonicPolyhedron, "../SolidiPlatonici/Ottaedro/"));

	PolyhedronMesh GeodeticPolyhedron;
	ASSERT_TRUE(Generation::TriangolaClasseI(3, 4, TRIANGULATION_PARAMETER, PlatonicPolyhedron, GeodeticPolyhedron, 0));

	auto [ExpectedV, ExpectedE, ExpectedF] = ExpectedCounts(4, TRIANGULATION_PARAMETER);

	EXPECT_EQ(GeodeticPolyhedron.NumCell0Ds, ExpectedV);
	EXPECT_EQ(GeodeticPolyhedron.NumCell1Ds, ExpectedE);
	EXPECT_EQ(GeodeticPolyhedron.NumCell2Ds, ExpectedF);

	int T = TRIANGULATION_PARAMETER * TRIANGULATION_PARAMETER;

	EXPECT_EQ(VertexDegree(4, GeodeticPolyhedron.Cell0DsId, GeodeticPolyhedron.Cell2DsVertices), 6);
	EXPECT_EQ(VertexDegree(6, GeodeticPolyhedron.Cell0DsId, GeodeticPolyhedron.Cell2DsVertices), 4*(T - 1));
}

// Terzo test: simmetria icosaedrica (q = 5)
TEST(TestGeodeticPolyhedronClasseI, Icosahedron) {
	PolyhedronMesh PlatonicPolyhedron;
	ASSERT_TRUE(FileManagement::ImportPolyhedronMesh(PlatonicPolyhedron, "../SolidiPlatonici/Icosaedro/"));

	PolyhedronMesh GeodeticPolyhedron;
	ASSERT_TRUE(Generation::TriangolaClasseI(3, 5, TRIANGULATION_PARAMETER, PlatonicPolyhedron, GeodeticPolyhedron, 0));

	auto [ExpectedV, ExpectedE, ExpectedF] = ExpectedCounts(5, TRIANGULATION_PARAMETER);

	EXPECT_EQ(GeodeticPolyhedron.NumCell0Ds, ExpectedV);
	EXPECT_EQ(GeodeticPolyhedron.NumCell1Ds, ExpectedE);
	EXPECT_EQ(GeodeticPolyhedron.NumCell2Ds, ExpectedF);

	int T = TRIANGULATION_PARAMETER * TRIANGULATION_PARAMETER;

	EXPECT_EQ(VertexDegree(5, GeodeticPolyhedron.Cell0DsId, GeodeticPolyhedron.Cell2DsVertices), 12);
	EXPECT_EQ(VertexDegree(6, GeodeticPolyhedron.Cell0DsId, GeodeticPolyhedron.Cell2DsVertices), 10*(T - 1));
}
