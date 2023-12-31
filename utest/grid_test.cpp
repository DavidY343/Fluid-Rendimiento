#include "sim/block.cpp"
#include "sim/grid.cpp"

#include "gtest/gtest.h"
#include <vector>

TEST(ObtenerIndiceTest, ComprobarIndiceBase) {
  std::vector<double> const vector        = {0.00011779, 0.00830882};
  std::vector<particula> const particulas = {};
  grid const malla                        = grid(vector, particulas);

  EXPECT_EQ(malla.obtener_indice(0, 0, 0), 0);
}

TEST(ObtenerIndiceTest, ComprobarIndiceZ) {
  std::vector<double> const vector        = {0.00011779, 0.00830882};
  std::vector<particula> const particulas = {};
  grid const malla                        = grid(vector, particulas);

  EXPECT_EQ(malla.obtener_indice(0, 0, 1), 1);
}

TEST(ObtenerIndiceTest, ComprobarIndiceY) {
  std::vector<double> const vector        = {0.00011779, 0.00830882};
  std::vector<particula> const particulas = {};
  grid const malla                        = grid(vector, particulas);

  EXPECT_EQ(malla.obtener_indice(0, 1, 1), 16);
}

TEST(ObtenerIndiceTest, ComprobarIndiceX) {
  std::vector<double> const vector        = {0.00011779, 0.00830882};
  std::vector<particula> const particulas = {};
  grid const malla                        = grid(vector, particulas);

  EXPECT_EQ(malla.obtener_indice(1, 1, 1), 331);
}

TEST(ObtenerCoordenadasTest, ComprobarCoordenadasBase) {
  std::vector<double> const vector        = {0.00011779, 0.00830882};
  std::vector<particula> const particulas = {};
  grid const malla                        = grid(vector, particulas);

  std::vector<int> const esperado = {0, 0, 0};
  EXPECT_EQ(malla.obtener_coordenadas(0), esperado);
}

TEST(ObtenerCoordenadasTest, ComprobarCoordenadasXYZ) {
  std::vector<double> const vector        = {0.00011779, 0.00830882};
  std::vector<particula> const particulas = {};
  grid const malla                        = grid(vector, particulas);

  std::vector<int> const esperado = {1, 1, 1};
  EXPECT_EQ(malla.obtener_coordenadas(331), esperado);
}

TEST(ObtenerContiguosTest, ComprobarContiguosBase) {
  std::vector<double> const vector        = {0.00011779, 0.00830882};
  std::vector<particula> const particulas = {};
  grid const malla                        = grid(vector, particulas);

  std::vector<int> const esperado = {1, 15, 16, 315, 316, 330, 331};
  EXPECT_EQ(malla.obtener_contiguos(0), esperado);
}

TEST(RecolocarParticulaTest, ComprobarParticulaBase) {
  std::vector<double> const vector        = {0.00011779, 0.00830882};
  std::vector<particula> const particulas = {};
  grid malla                              = grid(vector, particulas);

  malla.recolocar_particula(particula(0, -0.064, -0.079, -0.064, 0, 0, 0, 0, 0, 0));
  EXPECT_EQ(malla.acceder_bloque_part(0, 0).getid(), 0);
}

TEST(RecolocarParticulaTest, ComprobarParticulaMalla) {
  std::vector<double> const vector        = {0.00011779, 0.00830882};
  std::vector<particula> const particulas = {};
  grid malla                              = grid(vector, particulas);

  malla.recolocar_particula(particula(1, -0.050, -0.079, -0.064, 0, 0, 0, 0, 0, 0));
  EXPECT_EQ(malla.acceder_bloque_part(315, 0).getid(), 1);
}

TEST(RecolocarParticulaTest, ComprobarParticulaFuera) {
  std::vector<double> const vector        = {0.00011779, 0.00830882};
  std::vector<particula> const particulas = {};
  grid malla                              = grid(vector, particulas);

  malla.recolocar_particula(particula(3, -0.1, -0.1, -0.1, 0, 0, 0, 0, 0, 0));
  EXPECT_EQ(malla.acceder_bloque_part(0, 1).getid(), 3);
}

TEST(CalculoDensidadTest, ComprobarDensidadBase) {
  std::vector<double> const vector        = {0.00011779, 0.00830882};
  std::vector<particula> const particulas = {};
  grid malla                              = grid(vector, particulas);

  malla.recolocar_particula(particula(0, -0.064, -0.079, -0.064, 0, 0, 0, 0, 0, 0));
  malla.recolocar_particula(particula(1, -0.063, -0.079, -0.064, 0, 0, 0, 0, 0, 0));

  malla.inicializar_densidades();
  malla.calcular_densidades();

  EXPECT_EQ(malla.acceder_bloque_part(0, 0).getp(), 3.1493845807576875e-13);
  EXPECT_EQ(malla.acceder_bloque_part(0, 1).getp(), 3.1493845807576875e-13);
}

TEST(CalculoDensidadTest, ComprobarDensidad3Particulas) {
  std::vector<double> const vector        = {0.00011779, 0.00830882};
  std::vector<particula> const particulas = {};
  grid malla                              = grid(vector, particulas);

  malla.recolocar_particula(particula(0, -0.064, -0.079, -0.064, 0, 0, 0, 0, 0, 0));
  malla.recolocar_particula(particula(1, -0.063, -0.079, -0.064, 0, 0, 0, 0, 0, 0));
  malla.recolocar_particula(particula(2, -0.063, -0.077, -0.063, 0, 0, 0, 0, 0, 0));

  malla.inicializar_densidades();
  malla.calcular_densidades();

  EXPECT_EQ(malla.acceder_bloque_part(0, 0).getp(), 5.6542019373689712e-13);
  EXPECT_EQ(malla.acceder_bloque_part(0, 1).getp(), 5.7753110034231635e-13);
  EXPECT_EQ(malla.acceder_bloque_part(0, 2).getp(), 5.1307437792767598e-13);
}

TEST(CalculoDensidadTest, ComprobarDensidad5Particulas) {
  std::vector<double> const vector        = {0.00011779, 0.00830882};
  std::vector<particula> const particulas = {};
  grid malla                              = grid(vector, particulas);

  malla.recolocar_particula(particula(0, -0.064, -0.079, -0.064, 0, 0, 0, 0, 0, 0));
  malla.recolocar_particula(particula(1, -0.063, -0.079, -0.064, 0, 0, 0, 0, 0, 0));
  malla.recolocar_particula(particula(2, -0.063, -0.077, -0.063, 0, 0, 0, 0, 0, 0));
  malla.recolocar_particula(particula(3, -0.048, -0.077, -0.063, 0, 0, 0, 0, 0, 0));
  malla.recolocar_particula(particula(4, -0.047, -0.077, -0.063, 0, 0, 0, 0, 0, 0));

  malla.inicializar_densidades();
  malla.calcular_densidades();

  EXPECT_EQ(malla.acceder_bloque_part(0, 0).getp(), 5.6542019373689712e-13);
  EXPECT_EQ(malla.acceder_bloque_part(0, 1).getp(), 5.7753110034231635e-13);
  EXPECT_EQ(malla.acceder_bloque_part(0, 2).getp(), 5.1307437792767598e-13);
  EXPECT_EQ(malla.acceder_bloque_part(315, 0).getp(), 3.1493845807576875e-13);
  EXPECT_EQ(malla.acceder_bloque_part(630, 0).getp(), 3.1493845807576875e-13);
}

TEST(TransfDensidadTest, ComprobarTransfDensidadBase) {
  double const m                    = 0.00011779;
  double const h                    = 0.00830882;
  std::vector<double> const vector  = {m, h};
  std::vector<particula> particulas = {};

  particula const part = particula(0, -0.064, -0.079, -0.064, 0, 0, 0, 0, 0, 0);

  particulas.push_back(part);
  grid malla = grid(vector, particulas);
  malla.transformar_densidades();

  EXPECT_EQ(malla.acceder_bloque_part(0, 0).getp(),
            (315 * pow(h, 6) * m) / (64 * std::numbers::pi * pow(h, 9)));
}

TEST(TransfDensidadTest, ComprobarTransfDensidadConPSet) {
  double const m                    = 0.00011779;
  double const h                    = 0.00830882;
  std::vector<double> const vector  = {m, h};
  std::vector<particula> particulas = {};

  particula part = particula(1, -0.063, -0.079, -0.064, 0, 0, 0, 0, 0, 0);
  part.setp(4);
  particulas.push_back(part);
  grid malla = grid(vector, particulas);

  malla.transformar_densidades();
  EXPECT_EQ(malla.acceder_bloque_part(0, 1).getp(),
            (315 * (pow(h, 6) + 4) * m) / (64 * std::numbers::pi * pow(h, 9)));
}

TEST(TransfDensidadTest, ComprobarTransfDensidadConPSetAlto) {
  double const m                    = 0.00011779;
  double const h                    = 0.00830882;
  std::vector<double> const vector  = {m, h};
  std::vector<particula> particulas = {};

  particula part = particula(2, -0.063, -0.078, -0.064, 0, 0, 0, 0, 0, 0);
  part.setp(40000);
  particulas.push_back(part);
  grid malla = grid(vector, particulas);

  malla.transformar_densidades();

  EXPECT_EQ(malla.acceder_bloque_part(0, 2).getp(),
            (315 * (pow(h, 6) + 40000) * m) / (64 * std::numbers::pi * pow(h, 9)));
}

TEST(TransfAccTest, ComprobarTransfAccBase) {
  double const m                    = 0.00011779;
  double const h                    = 0.00830882;
  std::vector<double> const vector  = {m, h};
  std::vector<particula> particulas = {};

  particula part1 = particula(0, -0.064, -0.079, -0.064, 0, 0, 0, 0, 0, 0);
  particula part2 = particula(1, -0.063, -0.079, -0.064, 0, 0, 0, 0, 0, 0);
  particula part3 = particula(2, -0.063, -0.078, -0.064, 0, 0, 0, 0, 0, 0);
  particula part4 = particula(3, -0.063, -0.078, -0.047, 0, 0, 0, 0, 0, 0);

  particulas.push_back(part1);
  particulas.push_back(part2);
  particulas.push_back(part3);
  particulas.push_back(part4);

  grid malla = grid(vector, particulas);

  malla.inicializar_densidades();
  malla.transformar_densidades();
  malla.transferir_aceleraciones();

  EXPECT_EQ(malla.acceder_bloque_part(0, 0).getax(), 8774.1149213445824);
  EXPECT_EQ(malla.acceder_bloque_part(0, 0).getay(), 3378.878339769858);
  EXPECT_EQ(malla.acceder_bloque_part(0, 0).getaz(), 0);

  EXPECT_EQ(malla.acceder_bloque_part(2, 0).getax(), 0);
  EXPECT_EQ(malla.acceder_bloque_part(2, 0).getay(), -9.8);
  EXPECT_EQ(malla.acceder_bloque_part(2, 0).getaz(), 0);
}
