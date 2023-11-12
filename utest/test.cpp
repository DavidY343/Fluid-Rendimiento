//
// Created by manuel on 10/11/23.
//

#include <vector>
#include "gtest/gtest.h"
#include "sim/grid.cpp"


TEST(MiFuncionTest, ObtenerIndice) {
  std::vector<double> const vector = {0.00011779, 0.00830882};
  std::vector<particula> const particulas = {};
  grid const malla = grid(vector, particulas);

  EXPECT_EQ(malla.obtener_indice(0, 0, 0), 0);
  EXPECT_EQ(malla.obtener_indice(0, 0, 1), 1);
  EXPECT_EQ(malla.obtener_indice(0, 1, 1), 16);
  EXPECT_EQ(malla.obtener_indice(1, 1, 1), 331);
}

TEST(MiFuncionTest, ObtenerCoordenadas) {
  std::vector<double> const vector = {0.00011779, 0.00830882};
  std::vector<particula> const particulas = {};
  grid const malla = grid(vector, particulas);

  std::vector<int> esperado = {0, 0, 0};
  EXPECT_EQ(malla.obtener_coordenadas(0), esperado);

  esperado = {3, 2, 1};
  EXPECT_EQ(malla.obtener_coordenadas(976), esperado);
}

TEST(MiFuncionTest, ObtenerContiguos) {
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
  EXPECT_EQ(malla.acceder_bloque_part(0, 0).getid(), 0);

  malla.recolocar_particula(particula(1, -0.050, -0.079, -0.064, 0, 0, 0, 0, 0, 0));
  EXPECT_EQ(malla.acceder_bloque_part(315, 0).getid(), 1);

  malla.recolocar_particula(particula(2, -0.050, -0.070, -0.050, 0, 0, 0, 0, 0, 0));
  EXPECT_EQ(malla.acceder_bloque_part(331, 0).getid(), 2);

  malla.recolocar_particula(particula(3, -0.1, -0.1, -0.1, 0, 0, 0, 0, 0, 0));
  EXPECT_EQ(malla.acceder_bloque_part(0, 1).getid(), 3);
}

TEST(MiFuncionTest, CalculoDensidad) {
  std::vector<double> const vector  = {0.00011779, 0.00830882};
  std::vector<particula> const particulas = {};
  grid malla = grid(vector, particulas);

  malla.recolocar_particula(particula(0, -0.064, -0.079, -0.064, 0, 0, 0, 0, 0, 0));
  malla.recolocar_particula(particula(1, -0.063, -0.079, -0.064, 0, 0, 0, 0, 0, 0));

  malla.inicializar_densidades();
  malla.calcular_densidades();

  EXPECT_EQ(malla.acceder_bloque_part(0, 0).getp(), 3.1493845807576875e-13);
  EXPECT_EQ(malla.acceder_bloque_part(0, 1).getp(), 3.1493845807576875e-13);

  malla.recolocar_particula(particula(2, -0.063, -0.077, -0.063, 0, 0, 0, 0, 0, 0));


  malla.inicializar_densidades();
  malla.calcular_densidades();

  EXPECT_EQ(malla.acceder_bloque_part(0, 0).getp(), 5.6542019373689712e-13);
  EXPECT_EQ(malla.acceder_bloque_part(0, 1).getp(), 5.7753110034231635e-13);
  EXPECT_EQ(malla.acceder_bloque_part(0, 2).getp(), 5.1307437792767598e-13);

  malla.recolocar_particula(particula(3, -0.048, -0.077, -0.063, 0, 0, 0, 0, 0, 0));
  malla.recolocar_particula(particula(4, -0.047, -0.077, -0.063, 0, 0, 0, 0, 0, 0));

  malla.inicializar_densidades();
  malla.calcular_densidades();

  EXPECT_EQ(malla.acceder_bloque_part(315, 0).getp(), 3.1493845807576875e-13);
  EXPECT_EQ(malla.acceder_bloque_part(630, 0).getp(), 3.1493845807576875e-13);

}
