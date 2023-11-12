//
// Created by manuel on 10/11/23.
//

#include <vector>
#include "gtest/gtest.h"
#include "sim/grid.cpp"


TEST(MiFuncionTest, Obtener_Indice) {
  std::vector<double> const vector = {0.00011779, 0.00830882};
  std::vector<particula> const particulas = {};
  grid const malla = grid(vector, particulas);

  EXPECT_EQ(malla.obtener_indice(0, 0, 0), 0);
  EXPECT_EQ(malla.obtener_indice(0, 0, 1), 1);
  EXPECT_EQ(malla.obtener_indice(0, 1, 1), 16);
  EXPECT_EQ(malla.obtener_indice(1, 1, 1), 331);
}

TEST(MiFuncionTest, Obtener_Coordenadas) {
  std::vector<double> const vector = {0.00011779, 0.00830882};
  std::vector<particula> const particulas = {};
  grid const malla = grid(vector, particulas);

  std::vector<int> esperado = {0, 0, 0};
  EXPECT_EQ(malla.obtener_coordenadas(0), esperado);

  esperado = {3, 2, 1};
  EXPECT_EQ(malla.obtener_coordenadas(976), esperado);
}

TEST(MiFuncionTest, Obtener_Contiguos) {
  std::vector<double> const vector  = {0.00011779, 0.00830882};
  std::vector<particula> const particulas = {};
  grid const malla = grid(vector, particulas);

  std::vector<int> const esperado = {1, 15, 16, 315, 316, 330, 331};
  EXPECT_EQ(malla.obtener_contiguos(0), esperado);
}


TEST(MiFuncionTest, RecolocarParticula) {
  std::vector<double> const vector  = {0.00011779, 0.00830882};
  std::vector<particula> const particulas = {};
  grid malla = grid(vector, particulas);

  malla.recolocar_particula(particula(0, -0.064, -0.079, -0.064, 0, 0, 0, 0, 0, 0));
  EXPECT_EQ(malla.acceder_bloque_part(0, 0), 0);

  malla.recolocar_particula(particula(1, -0.050, -0.079, -0.064, 0, 0, 0, 0, 0, 0));
  EXPECT_EQ(malla.acceder_bloque_part(315, 0), 1);

  malla.recolocar_particula(particula(2, -0.050, -0.070, -0.050, 0, 0, 0, 0, 0, 0));
  EXPECT_EQ(malla.acceder_bloque_part(331, 0), 2);

  malla.recolocar_particula(particula(3, -0.1, -0.1, -0.1, 0, 0, 0, 0, 0, 0));
  EXPECT_EQ(malla.acceder_bloque_part(0, 1), 3);
}
